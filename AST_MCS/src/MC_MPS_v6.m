%% 水下光通信仿真: 纯几何步进 + 直射路径湍流解析惩罚 (方案A)
%  核心算法: 方法一(b-决定步长, a-解析吸收) + HG-PIS(g=0.97) + 强制接收
%  物理模型: 无相位屏物理偏折(纯直线步进)，但保留直射光的湍流展宽解析权重补偿
%  直射路径: 采用方案A (Importance Sampling PDF Likelihood Ratio) 模拟湍流展宽
clc; 
clear; 
close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- 核心模式选择 (相函数开关) ---
param.phase_func = 'Empirical';  
param.n_max = 200;          % 最大散射阶数

% --- 光源与波长参数 ---
lambda_nm = 514;
lambda = lambda_nm * 1e-9;
k_wave = 2 * pi / lambda;

% --- 初始波束参数 ---
w0 = 0.005;                                % 初始束腰 5mm，维持准直性
div_angle = 0.1 * pi / 180;                % 极窄波束：全发散角 0.1 度
theta_half_div = div_angle / 2; 

% --- 3D 空间链路布局 ---
Tx_Pos = [0, 0, 0];         
Rx_Pos = [0, 60, 0];         % 距离设定为 40m
Link_Vec = Rx_Pos - Tx_Pos;
Link_Dist = norm(Link_Vec);
Link_Dir = Link_Vec / Link_Dist;

% --- 光轴对准偏差 ---
theta_tx_error = deg2rad(0); phi_tx_error = 0;
theta_rx_error = deg2rad(0); phi_rx_error = pi;
mu_T = rotate_direction(Link_Dir, theta_tx_error, phi_tx_error);
Rx_Normal = rotate_direction(-Link_Dir, theta_rx_error, phi_rx_error);

if abs(mu_T(3)) < 0.9, up_temp_Tx = [0, 0, 1]; else, up_temp_Tx = [1, 0, 0]; end
u_vec_Tx = cross(up_temp_Tx, mu_T); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(mu_T, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

% --- 接收端物理参数 ---
Rx_Aperture = 0.05;                        % 接收孔径缩减至 5cm
Rx_Radius = Rx_Aperture / 2;
Rx_FOV = 5 * pi / 180;                     % 视场角缩紧至 5 度
Rx_Area = pi * (Rx_Radius)^2;

param.c_water = 2.237e8;

% Petzold清澈海水
param.coef_c = 0.1514;                % 衰减系数 c (1/m)
param.coef_a = 0.114;                 % 吸收系数 a (1/m)
param.coef_b = 0.0374;                % 散射系数 b (0.0374 1/m)

% --- 经验相函数拟合系数计算 (Haltrin) ---
param.c_e = param.coef_c; 
param.a_e = param.coef_a; 
b_e = param.coef_b; 
albedo_e = param.coef_b / param.coef_c;

param.q_e = 2.598 + 17.748 * sqrt(b_e) - 16.722 * b_e + 5.932 * b_e * sqrt(b_e);
param.k1 = 1.188 - 0.688 * albedo_e; 
param.k2 = 0.1 * (3.07 - 1.90 * albedo_e);
param.k3 = 0.01 * (4.58 - 3.02 * albedo_e); 
param.k4 = 0.001 * (3.24 - 2.25 * albedo_e);
param.k5 = 0.0001 * (0.84 - 0.61 * albedo_e);

%% ================= 2. 预计算与提议分布参数设置 =================
fprintf('正在预计算相函数归一化系数...\n');
th_test = linspace(0, pi, 2000); 
val_test_emp = zeros(size(th_test));

for i = 1:length(th_test)
    t_deg = max(th_test(i) * 180 / pi, 1e-6); 
    term = 1 - param.k1 * t_deg^0.5 + param.k2 * t_deg^1.0 - param.k3 * t_deg^1.5 + param.k4 * t_deg^2.0 - param.k5 * t_deg^2.5;
    val_test_emp(i) = exp(param.q_e * term);
end
param.b_emp_norm = 2 * pi * trapz(th_test, val_test_emp .* sin(th_test));

% --- 强制指定 HG-PIS 的不对称因子 ---
param.g_prop = 0.97;

% --- 弱湍流参数 (OTOPS 海洋湍流谱参数) ---
T_avg = 20;             % 平均温度 20°C 
S_avg = 35;             % 平均盐度 35 ppt 
H_ratio = -1;          % 温盐梯度比 -20 °C ppt^(-1) 
epsilon = 1e-10;         % 湍流动能耗散率 10^(-9) m^2 s^(-3) 
chi_T = 1e-5;           % 均方温度耗散率 10^(-7) K^2 s^(-1)

% --- 动态时间轴设置 (CIR) ---
t_LOS = Link_Dist / param.c_water;
delta_t_max = 10e-9;      % 限制探测时间窗口为 10 ns
dt = 0.01e-9;             % 时间分辨率
t_min = t_LOS; 
t_max = t_LOS + delta_t_max;
param.T_bins = t_min : dt : t_max;
param.Delta_T_bins = param.T_bins - t_LOS;
N_bins = length(param.T_bins);

% --- 仿真控制 ---
N_packets = 1e5; % 样本量 10万
scenarios = {'With Turbulence', 'Without Turbulence'};
results = struct();

%% ================= 3. 计算湍流理论解析相干长度 =================
% 仅保留 OTOPS 谱积分计算 rho_0，不再生成物理相位屏
[Phi_func, ~, eta_physical] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, H_ratio);

fprintf('预计算 OTOPS 湍流相干参数...\n');
[rho0_Link, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, Link_Dist, eta_physical);
fprintf('  内尺度 eta = %.4f mm\n', eta_physical * 1000);
fprintf('  等效 Cn2 = %.2e m^(-2/3)\n', Cn2_eq);
fprintf('  相干长度 rho_0 = %.4f m\n\n', rho0_Link);

%% ================= 4. 对比仿真主循环 =================
for s_idx = 1:2
    scenario_name = scenarios{s_idx};
    fprintf('=== Running Scenario: %s ===\n', scenario_name);
    
    % 绝对对齐种子，保证光线随机游走轨迹一致
    rng(123456, 'twister'); 
    h_time = zeros(1, N_bins); 
    tic;
    
    for p = 1:N_packets
        
        r0 = w0 * sqrt(-0.5 * log(rand())); 
        phi0 = 2 * pi * rand();
        pos_local = r0 * cos(phi0) * u_vec_Tx + r0 * sin(phi0) * v_vec_Tx;
        pos_init = Tx_Pos + pos_local;
        
        U_init = theta_half_div * sqrt(-0.5 * log(rand()));
        dir_init = rotate_direction(mu_T, U_init, 2 * pi * rand());
        
        weight_init = 1.0; 
        Huge_Aperture = 1e5; 
        
        % --- 0阶直射路径 (纯几何直线穿透) ---
        [pos_end, dir_end, plane_hit, path_len_ballistic] = ray_march_flat(pos_init, dir_init, 1e9, Rx_Pos, Huge_Aperture, pi, Rx_Normal, true);
        
        if plane_hit
            cos_rx_tilt = abs(dot(dir_end, Rx_Normal));
            if acos(cos_rx_tilt) <= Rx_FOV / 2
                r_wander_perp = norm(cross(pos_end - Rx_Pos, dir_end));
                r_eff = Rx_Radius * sqrt(cos_rx_tilt);
                
                % 方案A: 几何截断判定 + PDF重要性采样权重补偿
                W_geo = w0 + path_len_ballistic * tan(theta_half_div);
                
                if strcmp(scenario_name, 'With Turbulence')
                    % 1. 计算湍流引起的短期展宽 W_ST
                    rho0_ballistic = rho0_Link * (Link_Dist / path_len_ballistic)^(3/5);
                    W_turb_LT = 2 * path_len_ballistic / (k_wave * rho0_ballistic);
                    cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*w0))^(1/3));
                    W_ST = W_turb_LT * cf;
                    
                    % 2. 计算包含湍流的有效光斑大小 W_eff
                    W_eff = sqrt(W_geo^2 + W_ST^2);
                    
                    % 3. 核心：计算重要性采样权重因子 (Likelihood Ratio)
                    exponent = -2 * (r_wander_perp^2) * (1/(W_eff^2) - 1/(W_geo^2));
                    weight_factor = (W_geo^2 / W_eff^2) * exp(exponent);
                else
                    weight_factor = 1.0;
                end
                
                % 仅当光子几何上真正落入接收孔径时，才算被接收
                is_received = double(r_wander_perp <= r_eff); 
                
                % 将弹道能量填入时域仓 (附加重要性采样权重因子)
                base_energy = exp(-param.coef_c * path_len_ballistic) * is_received * weight_factor;
                
                bin_idx = floor((path_len_ballistic / param.c_water - t_min) / dt) + 1;
                if bin_idx >= 1 && bin_idx <= N_bins
                    h_time(bin_idx) = h_time(bin_idx) + weight_init * base_energy;
                end
            end
        end
        
        % --- 多重散射分支 PIS 追踪 (纯几何步进) ---
        pos = pos_init; 
        dir = dir_init; 
        current_dist = 0;
        weight_pis_likelihood = weight_init; 
        
        % 初始物理游走
        d_step = -log(rand()) / param.coef_b;
        [pos, dir, ~, step_len] = ray_march_flat(pos, dir, d_step, Rx_Pos, Huge_Aperture, Rx_FOV, Rx_Normal, false);
        current_dist = current_dist + step_len;
        
        if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist, continue; end
        
        for order = 1 : param.n_max
            % ---- 虚拟强制接收 (Local Estimation) ----
            vec_rx = Rx_Pos - pos; dist_rx = norm(vec_rx); dir_rx = vec_rx / dist_rx;
            
            if acos(dot(-dir_rx, Rx_Normal)) <= Rx_FOV/2
                p_phase = pdf_Empirical(dot(dir, dir_rx), param);
                omega = Rx_Area / (dist_rx^2) * abs(dot(dir_rx, Rx_Normal));
                
                current_physical_weight = weight_pis_likelihood * exp(-param.coef_a * current_dist);
                base_energy = current_physical_weight * min(1, p_phase * omega) * exp(-param.coef_c * dist_rx);
                
                if base_energy > 1e-15
                    % 虚拟光子纯几何直达接收面
                    [pos_v_end, ~, v_hit, v_len] = ray_march_flat(pos, dir_rx, dist_rx + 1e-1, Rx_Pos, Huge_Aperture, pi, Rx_Normal, true);
                    
                    if v_hit
                        cos_theta_inc = abs(dot(dir_rx, Rx_Normal));
                        r_eff_v = Rx_Radius * sqrt(cos_theta_inc);
                        r_wander_perp = norm(cross(pos_v_end - Rx_Pos, dir_rx));
                        
                        % 散射分量保留纯几何截断判定
                        f_spread = double(r_wander_perp <= r_eff_v); 
                        
                        base_energy = base_energy * f_spread;
                        bin_idx = floor(((current_dist + v_len) / param.c_water - t_min) / dt) + 1;
                        if bin_idx >= 1 && bin_idx <= N_bins
                            h_time(bin_idx) = h_time(bin_idx) + base_energy; 
                        end
                    end
                end
            end
            
            current_physical_weight = weight_pis_likelihood * exp(-param.coef_a * current_dist);
            if current_physical_weight < 1e-9
                if rand() > 0.1, break; else, weight_pis_likelihood = weight_pis_likelihood / 0.1; end
            end
            
            % ---- HG-PIS 生成提议散射角 ----
            xi_rand = rand();
            term = (1 - param.g_prop^2) / (1 - param.g_prop + 2 * param.g_prop * xi_rand);
            cos_t = max(min((1 + param.g_prop^2 - term^2) / (2 * param.g_prop), 1), -1);
            theta_i = acos(cos_t); 
            
            q_HG = (1 - param.g_prop^2) / (4 * pi * (1 + param.g_prop^2 - 2 * param.g_prop * cos_t)^1.5);
            p_val = pdf_Empirical(cos_t, param);
            
            weight_factor = min(1e3, p_val / q_HG);
            weight_pis_likelihood = weight_pis_likelihood * weight_factor;
            
            dir = rotate_direction(dir, theta_i, 2 * pi * rand());
            d_step = -log(rand()) / param.coef_b;
            
            % 几何光线偏折
            [pos, dir, ~, step_len] = ray_march_flat(pos, dir, d_step, Rx_Pos, Huge_Aperture, Rx_FOV, Rx_Normal, false);
            
            current_dist = current_dist + step_len;
            if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist, break; end
        end
    end
    t_run = toc; 
    fprintf('   Completed in %.2f s\n', t_run);
    
    h_time_norm = h_time / N_packets;      
    results(s_idx).h_time_cir = h_time_norm; 
    P_rx_total = sum(h_time_norm);
    results(s_idx).P_rx = P_rx_total;
    results(s_idx).Path_Loss_dB = -10 * log10(P_rx_total + eps);
    if P_rx_total > 0
        tau_mean = sum(param.T_bins .* h_time_norm) / P_rx_total;
        results(s_idx).tau_rms = sqrt( sum( ((param.T_bins - tau_mean).^2) .* h_time_norm ) / P_rx_total );
    else
        results(s_idx).tau_rms = 0;
    end
end

%% ================= 5. 对比输出与绘图 =================
fprintf('\n=== 仿真对比结果 (Tx-Rx: %.1fm) ===\n', Link_Dist);
fprintf('%-20s | %-14s | %-14s\n', 'Scenario', 'Path Loss (dB)', 'RMS Delay (ns)');
fprintf('------------------------------------------------------\n');
for i = 1:2
    fprintf('%-20s | %14.4f | %14.4f\n', scenarios{i}, results(i).Path_Loss_dB, results(i).tau_rms * 1e9);
end

% 轻度平滑以提升视觉曲线质量
smooth_window = 10;
cir_turb_smooth = movmean(results(1).h_time_cir, smooth_window);
cir_clean_smooth = movmean(results(2).h_time_cir, smooth_window);

figure('Name', 'Impulse Response Comparison (Optimized)', 'Color', 'w');
plot(param.Delta_T_bins * 1e9, cir_turb_smooth, 'r-', 'LineWidth', 1.5, 'DisplayName', 'With Turb'); 
hold on;
plot(param.Delta_T_bins * 1e9, cir_clean_smooth, 'b--', 'LineWidth', 1.5, 'DisplayName', 'No Turb');
grid on; 
legend('Location', 'best'); 
xlabel('Time Increment \Delta t (ns)', 'FontWeight', 'bold'); 
ylabel('Intensity (Smoothed)', 'FontWeight', 'bold');
title(sprintf('CIR (Narrow Beam + High Sensitivity, L=%.1fm)', Link_Dist));

%% ================= 6. 核心辅助函数区域 =================

function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(max(min(cos_theta, 1), -1)) * 180 / pi, 1e-6); 
    term = 1 - param.k1*t_deg^0.5 + param.k2*t_deg^1.0 - param.k3*t_deg^1.5 + param.k4*t_deg^2.0 - param.k5*t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end

function [pos, dir, hit_flag, total_len] = ray_march_flat(pos, dir, dist_limit, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check)
    % 彻底剔除相位屏梯度的精简射线步进
    hit_flag = false; 
    total_len = 0; 
    rem_dist = dist_limit; 
    
    while rem_dist > 1e-9
        t_rx = inf;
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx) > 1e-10
                t_temp = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
                if t_temp > -1e-7
                    t_rx = max(0, t_temp);
                end
            end
        end
        
        [min_dist, event_idx] = min([rem_dist, t_rx]);
        
        pos = pos + dir * min_dist;
        rem_dist = rem_dist - min_dist;
        total_len = total_len + min_dist;
        
        if event_idx == 2 % 击中接收平面
            if norm(pos - Rx_Pos) <= Rx_Aperture/2 && acos(dot(-dir, Rx_Normal)) <= Rx_FOV/2
                hit_flag = true;
            end
            break; 
        elseif event_idx == 1 % 耗尽当前设定步长
            break;
        end
    end
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

function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L, eta_physical)
    kappa_eval = 0.01 / eta_physical; Cn2_eq = Phi_n_func(kappa_eval) / (0.033 * kappa_eval^(-11/3));
    try, opts = optimset('Display', 'off'); rho0_exact = fzero(@(rho) 8*pi^2*k_wave^2*L*integral2(@(K, xi) K.*Phi_n_func(K).*(1-besselj(0, K.*rho.*xi)), 1e-2, 1e5, 0, 1, 'Method', 'iterated', 'RelTol', 1e-3) - 2, (0.545*k_wave^2*Cn2_eq*L)^(-3/5), opts); catch, rho0_exact = (0.545*k_wave^2*Cn2_eq*L)^(-3/5); end
end

function nd = rotate_direction(d, t, p)
    denom = sqrt(1 - d(3)^2);
    if denom < 1e-10, nd = [sin(t)*cos(p), sin(t)*sin(p), sign(d(3))*cos(t)]; else, nd = [sin(t)/denom*(d(1)*d(3)*cos(p)-d(2)*sin(p))+d(1)*cos(t), sin(t)/denom*(d(2)*d(3)*cos(p)+d(1)*sin(p))+d(2)*cos(t), -sin(t)*cos(p)*denom+d(3)*cos(t)]; end
    nd = nd / norm(nd);
end