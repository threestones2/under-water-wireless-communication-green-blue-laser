%% 1. WCI-MC vs Distance (Wavefront-Coupled + Adaptive HG-PIS, With Turbulence)
clc; clear; close all;

% ================= 参数初始化 =================
dist_array = 10:10:50;     % 传输距离数组 (m)
N_packets = 1e5;           % 仿真光子数
n_max = 10;                % 最大散射阶数
lambda = 514e-9; k_wave = 2 * pi / lambda;
w0 = 0.1; 
div_angle = 3 * pi / 180;

Rx_Aperture = 0.2; Rx_FOV = 10 * pi / 180; Rx_Area = pi * (Rx_Aperture / 2)^2;
c_water = 2.25e8;

param.coef_a = 0.114; param.coef_b = 0.0374; param.coef_c = 0.1514;      
param.albedo = param.coef_b / param.coef_c;
param = calc_haltrin_params(param);

p_peak = pdf_Empirical(cos(1e-6), param);
try param.theta_water = fzero(@(th) pdf_Empirical(cos(th), param) - p_peak*exp(-2), [1e-6, 0.1], optimset('Display','off'));
catch, param.theta_water = 0.02; end
C_val = 4 * pi * p_peak; g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);

% 弱湍流参数
T_avg = 20; S_avg = 35; H_ratio = -20; 
epsilon = 1e-9; chi_T = 1e-7; eta = 1e-3;
[Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
N_screens = 20; D_screen = 1; N_grid = 2^8; dx = D_screen / N_grid; x_axis = (-N_grid/2 : N_grid/2-1) * dx;

fprintf('=== 启动 WCI-MC 距离演化仿真 ===\n');

for d_idx = 1:length(dist_array)
    Link_Dist = dist_array(d_idx);
    Tx_Pos = [0, 0, 0]; Rx_Pos = [0, Link_Dist, 0];
    Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
    z_R = (pi * w0^2) / lambda; 
    
    delta_z_screen = Link_Dist / N_screens;
    [rho0_Link, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, Link_Dist);
    
    t_min = Link_Dist / c_water; dt = 1e-10; t_max = t_min + 1e-6;
    T_bins = t_min : dt : t_max; N_bins = length(T_bins); h_time = zeros(1, N_bins);
    
    u_vec = [1, 0, 0]; v_vec = [0, 0, 1];
    Screen_Chain = repmat(struct('Center',[], 'Normal',[0,1,0], 'u_vec',u_vec, 'v_vec',v_vec, 'grad_x',[], 'grad_y',[]), 1, N_screens);
    for i = 1:N_screens
        Screen_Chain(i).Center = Tx_Pos + i * delta_z_screen * Link_Dir;
        phi_screen = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
        [gx, gy] = gradient(phi_screen, dx);
        Screen_Chain(i).grad_x = gx; Screen_Chain(i).grad_y = gy;
    end
    
    tic;
    for p = 1:N_packets
        U_init = div_angle * sqrt(-0.5 * log(rand())); psi_init = 2 * pi * rand();
        dir_init = rotate_direction(Link_Dir, U_init, psi_init);
        pos = Tx_Pos; dir = dir_init; weight = 1.0; dist_travel = 0;
        
        % --- [核心修正]: 0阶直射路径恢复强制穿透湍流相位屏的射线追踪 ---
        Huge_Aperture = 1e5; % 使用虚拟大平面确保命中
        [pos_end, dir_end, plane_hit, path_len_ballistic] = ray_march_generic(pos, dir, 1e9, ...
            Rx_Pos, Huge_Aperture, pi, Rx_Normal, true, ...
            Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
        
        if plane_hit
            r_wander = norm(pos_end - Rx_Pos); % 真实的波束漂移距离
            cos_rx_tilt = abs(dot(dir_end, Rx_Normal));
            
            % 波束展宽综合计算 (衍射 + 湍流短期扩展)
            W_diff = w0 * sqrt(1 + (path_len_ballistic / z_R)^2);
            rho0_ballistic = rho0_Link * (Link_Dist / path_len_ballistic)^(3/5);
            W_turb_LT = 2 * path_len_ballistic / (k_wave * rho0_ballistic);
            cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*w0))^(1/3));
            W_turb_ST = W_turb_LT * cf;
            
            % 总有效光斑半径
            W_spot = sqrt(W_diff^2 + W_turb_ST^2);
            
            geom_loss = (2 * Rx_Area * cos_rx_tilt) / (pi * W_spot^2); 
            if acos(cos_rx_tilt) <= Rx_FOV / 2
                pointing_loss = exp(-2 * r_wander^2 / W_spot^2);
                w_ballistic = exp(-param.coef_c * path_len_ballistic) * min(1, geom_loss * pointing_loss);
                b_idx = floor((path_len_ballistic/c_water - t_min)/dt) + 1;
                if b_idx >= 1 && b_idx <= N_bins, h_time(b_idx) = h_time(b_idx) + w_ballistic; end
            end
        end
        
        % --- 1~n阶散射路径 ---
        d_step = -log(rand()) / param.coef_c;
        [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
        dist_travel = dist_travel + step_len;
        
        for order = 1 : n_max
            if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist, break; end
            
            vec2rx = Rx_Pos - pos; dist2rx = norm(vec2rx); dir2rx = vec2rx / dist2rx;
            if acos(dot(-dir2rx, Rx_Normal)) <= Rx_FOV/2
                base_w = weight * param.albedo * min(1, pdf_Empirical(dot(dir, dir2rx), param) * (Rx_Area / dist2rx^2 * abs(dot(dir2rx, Rx_Normal)))) * exp(-param.coef_c * dist2rx);
                if base_w > 1e-15
                    [pos_v, ~, v_hit, v_len] = ray_march_generic(pos, dir2rx, dist2rx+1e-1, Rx_Pos, Huge_Aperture, pi, Rx_Normal, true, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
                    if v_hit
                        r_wander_scat = norm(pos_v - Rx_Pos);
                        rho_0 = rho0_Link * (Link_Dist / v_len)^(3/5);
                        th_turb_ST = 2 / (k_wave * rho_0) * max(0, 1 - 0.37 * (rho_0 / Rx_Aperture)^(1/3));
                        penalty = (param.theta_water^2) / (param.theta_water^2 + th_turb_ST^2);
                        point_loss_scat = exp(-2 * r_wander_scat^2 / ((v_len^2) * (param.theta_water^2 + th_turb_ST^2)));
                        b_idx = floor(((dist_travel + v_len)/c_water - t_min)/dt) + 1;
                        if b_idx >= 1 && b_idx <= N_bins, h_time(b_idx) = h_time(b_idx) + base_w * penalty * point_loss_scat; end
                    end
                end
            end
            
            weight = weight * param.albedo; if weight < 1e-9, break; end
            xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi);
            cos_th_i = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1); th_i = acos(cos_th_i);
            weight = weight * min(1e3, pdf_Empirical(cos_th_i, param) / ((1 - g_prop^2) / (4 * pi * (1 + g_prop^2 - 2 * g_prop * cos_th_i)^1.5)));
            dir = rotate_direction(dir, th_i, 2*pi*rand());
            
            d_step = -log(rand()) / param.coef_c;
            [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
            dist_travel = dist_travel + step_len;
        end
    end
    t_run = toc;
    
    h_norm = h_time / N_packets; P_rx = sum(h_norm); tau_rms = 0;
    if P_rx > 0, t_mean = sum(T_bins .* h_norm) / P_rx; tau_rms = sqrt( sum(((T_bins - t_mean).^2) .* h_norm) / P_rx ); end
    fprintf('L = %2d m | Time: %5.2fs | Path Loss: %7.2f dB | RMS Delay: %.4f ns\n', Link_Dist, t_run, 10*log10(max(P_rx, 1e-300)), tau_rms*1e9);
end

% ---------- 辅助函数 ----------
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000, t=max(th(i)*180/pi,1e-6); val(i)=exp(p.q_e*(1-p.k1*t^0.5+p.k2*t^1.0-p.k3*t^1.5+p.k4*t^2.0-p.k5*t^2.5)); end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(max(min(cos_theta, 1), -1)) * 180 / pi, 1e-6); term = 1 - param.k1*t_deg^0.5 + param.k2*t_deg^1.0 - param.k3*t_deg^1.5 + param.k4*t_deg^2.0 - param.k5*t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end
function nd = rotate_direction(d, t, p)
    denom = sqrt(1 - d(3)^2);
    if denom < 1e-10, nd = [sin(t)*cos(p), sin(t)*sin(p), sign(d(3))*cos(t)]; else, nd = [sin(t)/denom*(d(1)*d(3)*cos(p)-d(2)*sin(p))+d(1)*cos(t), sin(t)/denom*(d(2)*d(3)*cos(p)+d(1)*sin(p))+d(2)*cos(t), -sin(t)*cos(p)*denom+d(3)*cos(t)]; end
    nd = nd / norm(nd);
end
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen)
    hit_flag = false; total_len = 0; remaining_dist = dist_limit; N_screens = length(Screen_Chain);
    while remaining_dist > 1e-6
        min_dist = remaining_dist; event_type = 'none'; target_idx = -1;
        current_layer = floor(dot(pos - Tx_Pos, Link_Dir) / delta_z_screen);
        if dot(dir, Link_Dir) > 1e-6, target_scr = current_layer + 1; elseif dot(dir, Link_Dir) < -1e-6, target_scr = current_layer; else, target_scr = -1; end
        if target_scr >= 1 && target_scr <= N_screens
            scr = Screen_Chain(target_scr); denom = dot(dir, scr.Normal);
            if abs(denom) > 1e-6, t = dot(scr.Center - pos, scr.Normal) / denom; if t > 1e-6 && t < min_dist, min_dist = t; event_type = 'screen'; target_idx = target_scr; end; end
        end
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx) > 1e-6, t_rx = dot(Rx_Pos - pos, Rx_Normal) / denom_rx; if t_rx > 1e-6 && t_rx <= min_dist, min_dist = t_rx; event_type = 'rx'; end; end
        end
        pos = pos + dir * min_dist; remaining_dist = remaining_dist - min_dist; total_len = total_len + min_dist; 
        if strcmp(event_type, 'rx')
            if norm(pos - Rx_Pos) <= Rx_Aperture/2 && acos(dot(-dir, Rx_Normal)) <= Rx_FOV/2, hit_flag = true; return; end
        elseif strcmp(event_type, 'screen')
            scr = Screen_Chain(target_idx); vec_p = pos - scr.Center;
            idx_x = mod(round((dot(vec_p, scr.u_vec) - x_axis(1))/dx), N_grid) + 1; idx_y = mod(round((dot(vec_p, scr.v_vec) - x_axis(1))/dx), N_grid) + 1;
            dir = dir + (scr.grad_x(idx_y, idx_x)*scr.u_vec + scr.grad_y(idx_y, idx_x)*scr.v_vec)/k_wave; dir = dir / norm(dir);
        end
    end
end
function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    k_wave = 2 * pi / lambda; A = -1.05e-6*S + 2*1.6e-8*T*S + 2*-2.02e-6*T -4.23e-3/(lambda*1e9); B = 1.779e-4 + -1.05e-6*T + 1.6e-8*T^2 + 1.155e-2/(lambda*1e9); 
    s_f = S*1e-3; Pr = 7; Sc = 700; c_T = 0.072^(4/3)/Pr; c_S = 0.072^(4/3)/Sc; c_TS = 0.072^(4/3)*(Pr+Sc)/(2*Pr*Sc);
    d_r = 0.15 * (2.6e-4*abs(H_ratio)/7.6e-4); chi_S = chi_T * d_r / (H_ratio^2); chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    Phi_Hill = @(K, chi_M, c_M) 0.0573 * chi_M * epsilon^(-1/3) .* (K.^2 + 1e-25).^(-11/6) .* exp(-176.9 * (K*eta).^2 .* c_M^(0.96)) .* (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + B^2 * Phi_Hill(K, chi_S, c_S) + 2*A*B * Phi_Hill(K, chi_TS, c_TS));
end
function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2*pi/D; dx = D/N; kx = (-N/2 : N/2-1)*dk; [KX, KY] = meshgrid(kx, kx); K_grid = sqrt(KX.^2 + KY.^2); K_grid(N/2+1, N/2+1) = 1e-10;
    F_phi = 2*pi*k_wave^2*delta_z*Phi_n_func(K_grid); F_phi(N/2+1, N/2+1) = 0;
    phase_high = real(ifftshift(ifft2(ifftshift((randn(N)+1i*randn(N)).*sqrt(F_phi)*dk)))) * N^2;
    phase_low = zeros(N, N); [xx, yy] = meshgrid((-N/2 : N/2-1)*dx);
    for p = 1:3, dk_p = dk/(3^p);
        for m = -1:1, for n = -1:1, if m==0&&n==0, continue; end
            kx_p = m*dk_p; ky_p = n*dk_p; k_p = sqrt(kx_p^2+ky_p^2);
            phase_low = phase_low + real((randn(1)+1i*randn(1))*sqrt(2*pi^2*k_wave^2*delta_z*Phi_n_func(k_p))*dk_p * exp(1i*(kx_p*xx+ky_p*yy)));
        end; end
    end
    phase_screen = phase_high + phase_low;
end
function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L)
    kappa_eval = 100; Cn2_eq = Phi_n_func(kappa_eval) / (0.033 * kappa_eval^(-11/3)); rho0_exact = (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5);
end