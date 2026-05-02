%% 1. 水下光通信仿真: MC-MPS (仅有湍流版本)
clc; clear; close all;

%% ================= 参数初始化 =================
param = struct();
param.phase_func = 'Empirical';  
param.n_max = 10;          

% --- 波长与光束 ---
lambda = 514 * 1e-9;
k_wave = 2*pi/lambda;
w0 = 0.1;                  
div_angle = 10*pi/180;     

% --- 空间布局 ---
Tx_Pos = [0, 0, 0];         
Rx_Pos = [0, 100, 0];       
Link_Vec = Rx_Pos - Tx_Pos;
Link_Dist = norm(Link_Vec);
Link_Dir = Link_Vec / Link_Dist;

mu_T = Link_Dir;            
Rx_Normal = -Link_Dir;

if abs(mu_T(3)) < 0.9, up_temp_Tx = [0, 0, 1]; else, up_temp_Tx = [1, 0, 0]; end
u_vec_Tx = cross(up_temp_Tx, mu_T); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(mu_T, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

Rx_Aperture = 0.2;          
Rx_FOV = 20 * pi/180;        
Rx_Area = pi*(Rx_Aperture/2)^2;

% --- 水体介质参数 ---
param.c_water = 2.25e8;     
param.coef_c = 0.1514;      
param.coef_a = 0.114;       
param.coef_b = 0.0374; 
param.albedo = param.coef_b / param.coef_c;

% --- Petzold 相函数参数 ---
param.q_e = 2.598 + 17.748*sqrt(param.coef_b) - 16.722*param.coef_b + 5.932*param.coef_b*sqrt(param.coef_b);
param.k1 = 1.188 - 0.688*param.albedo; param.k2 = 0.1 * (3.07 - 1.90*param.albedo);
param.k3 = 0.01 * (4.58 - 3.02*param.albedo); param.k4 = 0.001 * (3.24 - 2.25*param.albedo);
param.k5 = 0.0001 * (0.84 - 0.61*param.albedo);

th_test = linspace(0, pi, 2000); val_test = zeros(size(th_test));
for i=1:length(th_test)
    t_deg = max(th_test(i) * 180 / pi, 1e-6);
    term = 1 + (-1)^1*param.k1*t_deg^0.5 + (-1)^2*param.k2*t_deg^1.0 + ...
           (-1)^3*param.k3*t_deg^1.5 + (-1)^4*param.k4*t_deg^2.0 + (-1)^5*param.k5*t_deg^2.5;
    val_test(i) = exp(param.q_e * term);
end
param.b_emp_norm = 2 * pi * trapz(th_test, val_test .* sin(th_test));

% 计算水体发散角
obj_func = @(th) pdf_Empirical(cos(th), param) - pdf_Empirical(cos(1e-6), param) * exp(-2);
try opts = optimset('Display', 'off'); param.theta_water = fzero(obj_func, [1e-6, 0.1], opts); 
catch, param.theta_water = 0.02; end

% --- 湍流参数 ---
T_avg = 20; S_avg = 35; H_ratio = -20; epsilon = 1e-9; chi_T = 1e-7; eta = 1e-3;
N_screens = 20; D_screen = 0.2; N_grid = 2^10; delta_z_screen = Link_Dist / N_screens; 
dt = 1e-10; t_min = dt; t_max= 1e-6; param.T_bins = t_min : dt : t_max;
N_bins = length(param.T_bins);
N_packets = 1e5;

%% ================= 仿真核心 =================
[Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
[rho0_Link, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, Link_Dist);
dx = D_screen / N_grid; x_axis = (-N_grid/2 : N_grid/2-1) * dx;

if abs(Link_Dir(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
u_vec = cross(up_temp, Link_Dir); u_vec = u_vec / norm(u_vec);
v_vec = cross(Link_Dir, u_vec);   v_vec = v_vec / norm(v_vec);

Screen_Chain = repmat(struct('Center', [], 'Normal', Link_Dir, 'u_vec', u_vec, 'v_vec', v_vec, 'grad_x', [], 'grad_y', []), 1, N_screens);
for i = 1:N_screens
    Screen_Chain(i).Center = Tx_Pos + i * delta_z_screen * Link_Dir;
    phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
    [gx, gy] = gradient(phi, dx);
    Screen_Chain(i).grad_x = gx; Screen_Chain(i).grad_y = gy;
end

h_time = 1e-12*ones(1, N_bins); 
tic;
for p = 1:N_packets
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos_init = Tx_Pos + r0*cos(phi0)*u_vec_Tx + r0*sin(phi0)*v_vec_Tx;
    dir_init = rotate_direction(mu_T, div_angle * sqrt(-0.5 * log(rand())), 2 * pi * rand);
    
    % [直射光]
    [pos_end, dir_end, plane_hit, path_len] = ray_march_generic(pos_init, dir_init, 1e9, Rx_Pos, 1e5, pi, Rx_Normal, true, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
    if plane_hit
        W_spot = calc_beam_spot_size(w0, lambda, path_len, Cn2_eq);
        cos_tilt = abs(dot(dir_end, Rx_Normal));
        geom_loss = (2 * Rx_Area * cos_tilt) / (pi * W_spot^2); 
        if acos(cos_tilt) > Rx_FOV / 2, geom_loss = 0; end
        
        weight = exp(-param.coef_c * path_len) * min(1, geom_loss * exp(-2 * norm(pos_end - Rx_Pos)^2 / W_spot^2));
        idx = floor((path_len/param.c_water - t_min) / dt) + 1;
        if idx >= 1 && idx <= N_bins, h_time(idx) = h_time(idx) + weight; end
    end

    % [散射光]
    pos = pos_init; dir = dir_init; weight = 1.0; dist_travel = 0;
    d_step = -log(rand()) / param.coef_c;
    [pos, dir, ~, s_len] = ray_march_generic(pos, dir, d_step, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
    dist_travel = dist_travel + s_len;
    
    for order = 1 : param.n_max
        if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist, break; end
        v_rx = Rx_Pos - pos; d_rx = norm(v_rx); d_dir = v_rx / d_rx;
        
        if acos(dot(-d_dir, Rx_Normal)) <= Rx_FOV/2
            base_eng = weight * param.albedo * min(1, pdf_Empirical(dot(dir, d_dir), param) * (Rx_Area/d_rx^2*abs(dot(d_dir, Rx_Normal)))) * exp(-param.coef_c * d_rx);
            if base_eng > 1e-15
                [p_end, ~, v_hit, v_len] = ray_march_generic(pos, d_dir, d_rx + 1e-1, Rx_Pos, 1e5, pi, Rx_Normal, true, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
                if v_hit
                    r_0 = rho0_Link * (Link_Dist / v_len)^(3/5);
                    theta_ST = (2 / (k_wave * r_0)) * max(0, 1 - 0.37 * (r_0 / Rx_Aperture)^(1/3));
                    W_sq = (v_len^2) * (param.theta_water^2 + theta_ST^2);
                    
                    eng = base_eng * (param.theta_water^2 / (param.theta_water^2 + theta_ST^2)) * exp(-2 * norm(p_end - Rx_Pos)^2 / W_sq);
                    idx = floor(((dist_travel + v_len)/param.c_water - t_min) / dt) + 1;
                    if idx >= 1 && idx <= N_bins, h_time(idx) = h_time(idx) + eng; end
                end
            end
        end
        
        weight = weight * param.albedo; if weight < 1e-9, break; end
        th_i = pi * rand; p_scat = 2 * pi * rand;
        weight = weight * (2 * pi^2 * pdf_Empirical(cos(th_i), param) * sin(th_i));
        dir = rotate_direction(dir, th_i, p_scat);
        
        [pos, dir, ~, s_len] = ray_march_generic(pos, dir, -log(rand())/param.coef_c, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
        dist_travel = dist_travel + s_len;
    end
end
fprintf('Completed in %.2f s\n', toc);

% 绘图
figure('Color', 'w'); plot(param.T_bins*1e9, 10*log10(h_time/N_packets), 'r-'); grid on;
xlabel('Time (ns)'); ylabel('Received Power (dB)'); title('CIR - With Turbulence');

%% ================= 辅助函数 =================
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen)
    hit_flag = false; total_len = 0; remaining = dist_limit; N_screens = length(Screen_Chain);
    while remaining > 1e-6
        min_dist = remaining; event = 'none'; t_idx = -1;
        c_layer = floor(dot(pos - Tx_Pos, Link_Dir) / delta_z_screen);
        cos_t = dot(dir, Link_Dir);
        if cos_t > 1e-6, t_scr = c_layer + 1; elseif cos_t < -1e-6, t_scr = c_layer; else, t_scr = -1; end
        
        if t_scr >= 1 && t_scr <= N_screens
            denom = dot(dir, Screen_Chain(t_scr).Normal);
            if abs(denom) > 1e-6
                t = dot(Screen_Chain(t_scr).Center - pos, Screen_Chain(t_scr).Normal) / denom;
                if t > 1e-6 && t < min_dist, min_dist = t; event = 'screen'; t_idx = t_scr; end
            end
        end
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx) > 1e-6
                t_rx = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
                if t_rx > 1e-6 && t_rx <= min_dist, min_dist = t_rx; event = 'rx'; end
            end
        end
        pos = pos + dir * min_dist; remaining = remaining - min_dist; total_len = total_len + min_dist;
        
        if strcmp(event, 'rx')
            if norm(pos - Rx_Pos) <= Rx_Aperture/2 && acos(dot(-dir, Rx_Normal)) <= Rx_FOV/2
                hit_flag = true; return;
            end
        elseif strcmp(event, 'screen')
            scr = Screen_Chain(t_idx); vec = pos - scr.Center;
            idx_x = mod(round((dot(vec, scr.u_vec) - x_axis(1))/dx), N_grid) + 1;
            idx_y = mod(round((dot(vec, scr.v_vec) - x_axis(1))/dx), N_grid) + 1;
            delta = (scr.grad_x(idx_y, idx_x)*scr.u_vec + scr.grad_y(idx_y, idx_x)*scr.v_vec) / k_wave;
            dir = (dir + delta) / norm(dir + delta);
        end
    end
end
function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    k_wave = 2 * pi / lambda; A = -2.02e-6*2*T - 4.23e-3/(lambda*1e9); B = 1.779e-4 + 1.155e-2/(lambda*1e9);
    Phi_Hill = @(K, c) 0.0573 * chi_T * epsilon^(-1/3) .* (K.^2 + 1e-25).^(-11/6) .* exp(-176.9*(K*eta).^2.*c^0.96);
    Phi_n_func = @(K) A^2 * Phi_Hill(K, 1);
end
function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2 * pi / D; [KX, KY] = meshgrid((-N/2 : N/2-1) * dk);
    K_grid = sqrt(KX.^2 + KY.^2); K_grid(N/2+1, N/2+1) = 1e-10;
    C_nm = (randn(N) + 1i * randn(N)) .* sqrt(2*pi*k_wave^2*delta_z*Phi_n_func(K_grid)) * dk;
    phase_screen = real(ifftshift(ifft2(ifftshift(C_nm)))) * N^2;
end
function [rho0, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, L)
    Cn2_eq = Phi_func(100) / (0.033 * 100^(-11/3));
    rho0 = (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5);
end
function W_L = calc_beam_spot_size(w0, lambda, L, Cn2)
    k = 2 * pi / lambda; z_R = (pi * w0^2) / lambda;
    W_diff = w0 * sqrt(1 + (L / z_R)^2);
    W_turb = 2 * L / (k * (0.545 * k^2 * Cn2 * L)^(-3/5));
    W_L = sqrt(W_diff^2 + W_turb^2);
end
function new_dir = rotate_direction(dir, theta, psi)
    mz = dir(3); denom = sqrt(1 - mz^2);
    if denom < 1e-10, new_dir = [sin(theta)*cos(psi), sin(theta)*sin(psi), sign(mz)*cos(theta)];
    else, new_dir = [sin(theta)/denom*(dir(1)*mz*cos(psi) - dir(2)*sin(psi))+dir(1)*cos(theta), ...
                     sin(theta)/denom*(dir(2)*mz*cos(psi) + dir(1)*sin(psi))+dir(2)*cos(theta), -sin(theta)*cos(psi)*denom+mz*cos(theta)];
    end
    new_dir = new_dir / norm(new_dir);
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(max(min(cos_theta, 1), -1)) * 180 / pi, 1e-6);
    term = 1 + (-1)^1*param.k1*t_deg^0.5 + (-1)^2*param.k2*t_deg^1.0 + ...
           (-1)^3*param.k3*t_deg^1.5 + (-1)^4*param.k4*t_deg^2.0 + (-1)^5*param.k5*t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end