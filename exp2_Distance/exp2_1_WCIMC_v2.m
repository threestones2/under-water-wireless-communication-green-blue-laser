%% Exp2: WCI-MC (Wavefront-Coupled + Adaptive HG-PIS, Optimized Version)
% 优化架构: 
% 1. 物理降维: 仅0阶弹道光子执行级联相位屏穿透, 散射光子采用纯解析展宽评估
% 2. 内存扁平化: 废弃 Struct 结构体寻址, 采用 3D 连续内存矩阵提取波前梯度
clc; clear; close all;

% ================= 物理参数初始化 =================
dist_cell = {10:10:60, 10:10:60, 5:5:20}; 
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr = [0.1514, 0.398, 2.190];
num_W = length(water_types);

N_packets = 1e5; n_max = 200; 
lambda = 514e-9; k_wave = 2 * pi / lambda;
w0 = 0.002; 
div_angle = 20 * pi / 180; theta_half_div = div_angle / 2; % 宽波束
Rx_Aperture = 0.2; Rx_FOV = 90 * pi / 180; Rx_Area = pi * (Rx_Aperture / 2)^2; % 宽视场

% --- OTOPS 强光学湍流参数 ---
T_avg = 20; S_avg = 35; H_ratio = -1; 
epsilon = 1e-10; chi_T = 1e-5; 

N_screens = 20; D_screen = 1; N_grid = 2^8; dx = D_screen / N_grid; x_axis = (-N_grid/2 : N_grid/2-1) * dx;

PL_Cell_WCI = cell(1, num_W);

fprintf('--- 运行 Exp2: WCI-MC (深度优化版, 宽波束配置) ---\n');

% 预计算湍流谱及物理内尺度
[Phi_func, ~, eta_physical] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, H_ratio);

for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx); param.coef_b = coef_b_arr(w_idx); param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c;
    param = calc_haltrin_params(param);
    
    P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max;
    g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);
    
    dist_arr = dist_cell{w_idx};
    num_D = length(dist_arr);
    PL_arr_temp = zeros(1, num_D);
    
    for d_idx = 1:num_D
        L = dist_arr(d_idx);
        Tx_Pos = [0, 0, 0]; Rx_Pos = [0, L, 0]; Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
        u_vec_Tx = [1, 0, 0]; v_vec_Tx = [0, 0, 1];
        
        % 动态计算当前距离的相干参数
        [rho0_Link, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, L, eta_physical);
        delta_z_screen = L / N_screens;
        
        % ========================================================
        % 【优化方案二】：数据结构扁平化 (Data Flattening)
        % 废弃 struct，构建连续内存的 3D 矩阵以消除寻址瓶颈
        % ========================================================
        Grad_X_3D = zeros(N_grid, N_grid, N_screens);
        Grad_Y_3D = zeros(N_grid, N_grid, N_screens);
        Screen_Z_1D = zeros(1, N_screens);
        
        for i = 1:N_screens
            phi_screen = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
            [gx, gy] = gradient(phi_screen, dx);
            Grad_X_3D(:, :, i) = gx;
            Grad_Y_3D(:, :, i) = gy;
            Screen_Z_1D(i) = i * delta_z_screen; % 相对 Tx_Pos 的一维 Z 坐标
        end
        
        P_rx_accum = 0; tic;
        
        rng(123456, 'twister'); % 注入控制变量种子
        
        for p = 1:N_packets
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
            pos = Tx_Pos + r0*cos(phi0)*u_vec_Tx + r0*sin(phi0)*v_vec_Tx;
            U_init = theta_half_div * sqrt(-0.5 * log(rand()));
            dir = rotate_direction(Link_Dir, U_init, 2*pi*rand());
            
            weight = 1.0; P_packet = 0;
            
            % ========================================================
            % 【优化方案一】：物理降维 
            % 0 阶直射路径: 维持严密的级联相位屏穿透，捕获相干波束漂移
            % ========================================================
            [pos_end, dir_end, plane_hit, path_len_ballistic] = ray_march_flat(pos, dir, 1e9, Rx_Pos, 1e5, pi, Rx_Normal, true, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, u_vec_Tx, v_vec_Tx, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
            
            if plane_hit
                cos_rx_tilt = abs(dot(dir_end, Rx_Normal));
                if acos(cos_rx_tilt) <= Rx_FOV / 2
                    r_wander_perp = norm(cross(pos_end - Rx_Pos, dir_end));
                    r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);
                    
                    W_geo = w0 + path_len_ballistic * tan(theta_half_div);
                    rho0_ballistic = rho0_Link * (L / path_len_ballistic)^(3/5);
                    W_turb_LT = 2 * path_len_ballistic / (k_wave * rho0_ballistic);
                    cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*w0))^(1/3));
                    W_spot = sqrt(W_geo^2 + (W_turb_LT * cf)^2);
                    
                    if W_spot > 1e-6, recv_frac = 1 - marcumq(2*r_wander_perp/W_spot, 2*r_eff/W_spot); else, recv_frac = double(r_wander_perp <= r_eff); end
                    P_packet = P_packet + exp(-param.coef_c * path_len_ballistic) * recv_frac;
                end
            end
            
            % --- 多重散射演化 (非相干漫射场) ---
            for ord = 1:n_max
                % 1. 物理空间步进 (关闭非相干光子的相位屏扰动)
                d_step = -log(rand()) / param.coef_b; 
                pos = pos + dir * d_step; 
                weight = weight * exp(-param.coef_a * d_step);
                
                if dot(pos - Tx_Pos, Link_Dir) >= L || weight < 1e-15, break; end 
                
                % 2. 虚拟强制接收 (Local Estimation)
                vec2rx = Rx_Pos - pos; dist2rx = norm(vec2rx); dir2rx = vec2rx / dist2rx;
                cos_inc = dot(-dir2rx, Rx_Normal);
                
                if cos_inc >= cos(Rx_FOV/2)
                    omega = Rx_Area / (dist2rx^2) * cos_inc;
                    base_w = weight * min(1, pdf_Empirical(dot(dir, dir2rx), param) * omega) * exp(-param.coef_c * dist2rx);
                    
                    if base_w > 1e-15
                        % 纯解析计算次级散射点源的湍流短曝光展宽 (摒弃极耗时的光线步进)
                        r_eff_v = (Rx_Aperture/2) * sqrt(cos_inc);
                        rho_0_sph = rho0_Link * (L / dist2rx)^(3/5);
                        W_ST = (2 * dist2rx / (k_wave * rho_0_sph)) * max(0, 1 - 0.37*(rho_0_sph/(2*r_eff_v))^(1/3));
                        
                        % 由于虚拟光子精准指向接收机，离轴漂移 r_wp = 0，Marcum-Q 积分物理退化为纯指数项
                        if W_ST > 1e-6
                            point_loss = 1 - exp(-2 * r_eff_v^2 / W_ST^2);
                        else
                            point_loss = 1.0; 
                        end
                        P_packet = P_packet + base_w * point_loss;
                    end
                end
                
                % 3. 自适应 HG-PIS 采样下一跳极角
                xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); cos_th_i = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1); th_i = acos(cos_th_i);
                weight = weight * min(1e3, pdf_Empirical(cos_th_i, param) / ((1 - g_prop^2) / (4 * pi * (1 + g_prop^2 - 2 * g_prop * cos_th_i)^1.5)));
                dir = rotate_direction(dir, th_i, 2*pi*rand());
            end
            P_rx_accum = P_rx_accum + P_packet;
        end
        PL_arr_temp(d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
        fprintf('  %s | L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', water_types{w_idx}, L, toc, -PL_arr_temp(d_idx));
    end
    PL_Cell_WCI{w_idx} = PL_arr_temp;
end
save('data_exp2_WCIMC.mat', 'dist_cell', 'PL_Cell_WCI', 'water_types');

% ================= 核心辅助函数区域 =================

% 【重构】扁平化矩阵版光线步进函数
function [pos, dir, hit_flag, total_len] = ray_march_flat(pos, dir, dist_limit, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, u_vec_Tx, v_vec_Tx, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen)
    hit_flag = false; 
    step_actual = dist_limit; 
    hit_rx_plane = false; 
    N_screens = length(Screen_Z_1D);
    
    if enable_hit_check
        denom_rx = dot(dir, Rx_Normal);
        if abs(denom_rx) > 1e-6
            t_rx = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
            if t_rx > 1e-6 && t_rx <= dist_limit
                step_actual = t_rx; 
                hit_rx_plane = true; 
            end
        end
    end
    
    total_len = step_actual;
    dir_z = dot(dir, Link_Dir); 
    z_start = dot(pos - Tx_Pos, Link_Dir); 
    z_end = z_start + dir_z * step_actual;
    
    if dir_z > 1e-6
        crossed_screens = (floor(z_start / delta_z_screen) + 1) : floor(z_end / delta_z_screen);
    elseif dir_z < -1e-6
        crossed_screens = (ceil(z_start / delta_z_screen) - 1) : -1 : ceil(z_end / delta_z_screen);
    else
        crossed_screens = []; 
    end
    
    crossed_screens = crossed_screens(crossed_screens >= 1 & crossed_screens <= N_screens);
    dir_old = dir; 
    
    if ~isempty(crossed_screens)
        delta_dir = [0, 0, 0]; 
        for i = crossed_screens
            % 运用 3D 连续矩阵寻址，突破计算瓶颈
            t_i = (Screen_Z_1D(i) - z_start) / dir_z;
            pos_i = pos + dir_old * t_i;
            
            loc_u = dot(pos_i - Tx_Pos, u_vec_Tx); 
            loc_v = dot(pos_i - Tx_Pos, v_vec_Tx);
            
            idx_x = mod(round((loc_u - x_axis(1)) / dx), N_grid) + 1; 
            idx_y = mod(round((loc_v - x_axis(1)) / dx), N_grid) + 1;
            
            delta_dir = delta_dir + (Grad_X_3D(idx_y, idx_x, i) * u_vec_Tx + Grad_Y_3D(idx_y, idx_x, i) * v_vec_Tx);
        end
        dir = dir_old + delta_dir / k_wave; 
        dir = dir / norm(dir);
    end
    
    pos = pos + dir_old * step_actual; 
    
    if hit_rx_plane
        if norm(pos - Rx_Pos) <= Rx_Aperture/2 && acos(dot(-dir, Rx_Normal)) <= Rx_FOV/2
            hit_flag = true; 
        end
    end
end

function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000
        t=max(th(i)*180/pi,1e-6); 
        val(i)=exp(p.q_e*(1-p.k1*t^0.5+p.k2*t^1.0-p.k3*t^1.5+p.k4*t^2.0-p.k5*t^2.5)); 
    end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); 
    param = p;
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

function [Phi_n_func, k_wave, eta_physical] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, H_ratio)
    lambda_nm = lambda * 1e9; k_wave = 2 * pi / lambda; 
    A = -1.05e-6*S + 2*1.6e-8*T*S - 2*2.02e-6*T - 4.23e-3/lambda_nm; B = 1.779e-4 - 1.05e-6*T + 1.6e-8*T^2 + 1.155e-2/lambda_nm; 
    s_f = S*1e-3; T_k = T+273.15;
    cp = 1000 * ((5.328 - 9.76e-2 * S + 4.04e-4 * S^2) + (-6.913e-3 + 7.351e-4 * S - 3.15e-6 * S^2) * T + (9.6e-6 - 1.927e-6 * S + 8.23e-9 * S^2) * T^2 + (2.5e-9 + 1.666e-9 * S - 7.125e-12 * S^2) * T^3); 
    rho = (9.9992293295e2 + 2.0341179217e-2 * T - 6.1624591598e-3 * T^2 + 2.2614664708e-5 * T^3 - 4.6570659168e-8 * T^4) + s_f * (8.0200240891e2 - 2.0005183488 * T + 1.6771024982e-2 * T^2 - 3.0600536746e-5 * T^3 - 1.6132224742e-5 * T * S);
    mu = ((0.15700386464 * (T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5) * (1 + (1.5409136040 + 1.9981117208e-2 * T - 9.5203865864e-5 * T^2) * s_f + (7.9739318223 - 7.561456881e-2 * T + 4.7237011074e-4 * T^2) * s_f^2);
    sigma_T = 10^(log10(240 + 0.0002 * (S / 1.00472)) - 3 + 0.434 * (2.3 - (343.5 + 0.037 * (S / 1.00472)) / (1.00024 * T + 273.15)) * (1 - (1.00024 * T + 273.15) / (647.3 + 0.03 * (S / 1.00472)))^(1/3)); 
    Pr = mu * cp / sigma_T; Sc = mu^2 / (5.954e-15 * T_k * rho); 
    c_T = 0.072^(4/3) / Pr; c_S = 0.072^(4/3) / Sc; c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    R_rho = 2.6e-4 * abs(H_ratio) / 7.6e-4;
    if R_rho >= 1, d_r = R_rho + sqrt(R_rho)*sqrt(R_rho - 1); elseif R_rho >= 0.5, d_r = 1.85*R_rho - 0.85; else, d_r = 0.15*R_rho; end
    chi_S = chi_T * d_r / (H_ratio^2); chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    
    nu_kinematic = mu / rho; eta_physical = (nu_kinematic^3 / epsilon)^(0.25); 
    Phi_Hill = @(K, chi_M, c_M) (0.72 / (4 * pi)) * chi_M * epsilon^(-1/3) .* (K.^2 + 1e-25).^(-11/6) .* exp(-176.90 * (K * eta_physical).^2 .* c_M^(0.96)) .* (1 + 21.61 * (K * eta_physical).^(0.61) .* c_M^(0.02) - 18.18 * (K * eta_physical).^(0.55) .* c_M^(0.04));
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + B^2 * Phi_Hill(K, chi_S, c_S) + 2 * A * B * Phi_Hill(K, chi_TS, c_TS));
end

function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2*pi/D; dx = D/N; kx = (-N/2 : N/2-1)*dk; [KX, KY] = meshgrid(kx, kx); K_grid = sqrt(KX.^2 + KY.^2); K_grid(N/2+1, N/2+1) = 1e-10;
    C_nm = (randn(N)+1i*randn(N)).*sqrt(2*pi*k_wave^2*delta_z*Phi_n_func(K_grid))*dk; phase_high = real(ifftshift(ifft2(ifftshift(C_nm)))) * N^2;
    phase_low = zeros(N, N); [xx, yy] = meshgrid((-N/2 : N/2-1)*dx);
    for p = 1:3, dk_p = dk/(3^p);
        for m = -1:1, for n = -1:1, if m==0&&n==0, continue; end
            k_p = sqrt((m*dk_p)^2+(n*dk_p)^2); phase_low = phase_low + real((randn(1)+1i*randn(1))*sqrt(2*pi^2*k_wave^2*delta_z*Phi_n_func(k_p))*dk_p * exp(1i*(m*dk_p*xx+n*dk_p*yy)));
        end; end
    end
    phase_screen = phase_high + phase_low;
end

function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L, eta_physical)
    kappa_eval = 0.01 / eta_physical; Cn2_eq = Phi_n_func(kappa_eval) / (0.033 * kappa_eval^(-11/3));
    try, opts = optimset('Display', 'off'); rho0_exact = fzero(@(rho) 8*pi^2*k_wave^2*L*integral2(@(K, xi) K.*Phi_n_func(K).*(1-besselj(0, K.*rho.*xi)), 1e-2, 1e5, 0, 1, 'Method', 'iterated', 'RelTol', 1e-3) - 2, (0.545*k_wave^2*Cn2_eq*L)^(-3/5), opts); catch, rho0_exact = (0.545*k_wave^2*Cn2_eq*L)^(-3/5); end
end