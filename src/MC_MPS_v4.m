%% 水下/大气光通信混合仿真: MC-MPS (Turbulence Comparison)
%  核心算法: 方法一(b-决定步长, a-解析吸收) + HG-PIS(g=0.97 提议分布) + 强制接收
%  相函数模: Fournier-Forand (FF) / Haltrin 经验拟合 双模切换
clc; 
clear; 
close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- 核心模式选择 (相函数开关) ---
% 可选: 'FFPF' (Fournier-Forand模型, 论文指定) 或 'Empirical' (Haltrin经验拟合)
param.phase_func = 'Empirical';  
param.n_max = 50;          % 蒙特卡洛追踪的最大多重散射阶数

% --- 光源与波长参数 ---
lambda_nm = 514;
lambda = lambda_nm * 1e-9;
k_wave = 2 * pi / lambda;
w0 = 0.002;                      
div_angle = 6*pi/180; 

% --- 3D 空间链路布局 ---
Tx_Pos = [0, 0, 0];         
Rx_Pos = [0, 50, 0];       
Link_Vec = Rx_Pos - Tx_Pos;
Link_Dist = norm(Link_Vec);
Link_Dir = Link_Vec / Link_Dist;

% --- 光轴对准偏差 ---
theta_tx_error = deg2rad(0); 
phi_tx_error = 0;
theta_rx_error = deg2rad(0); 
phi_rx_error = pi;
mu_T = rotate_direction(Link_Dir, theta_tx_error, phi_tx_error);
Rx_Normal = rotate_direction(-Link_Dir, theta_rx_error, phi_rx_error);

if abs(mu_T(3)) < 0.9
    up_temp_Tx = [0, 0, 1]; 
else
    up_temp_Tx = [1, 0, 0]; 
end
u_vec_Tx = cross(up_temp_Tx, mu_T); 
u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(mu_T, u_vec_Tx);   
v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

% --- 接收端物理参数 ---
Rx_Aperture = 0.05;          
Rx_Radius = Rx_Aperture / 2;
Rx_FOV = 40 * pi / 180;        
Rx_Area = pi * (Rx_Radius)^2;

param.c_water = 2.237e8;

% Petzold清澈海水
param.coef_c = 0.1514;                % 衰减系数 c (1/m)
param.coef_a = 0.114;                 % 吸收系数 a (1/m)
param.coef_b = 0.0374; % 散射系数 b (0.0374 1/m)

% %Petzold近岸海水
% param.coef_c = 0.398;
% param.coef_a = 0.179;
% param.coef_b=0.219;

% Petzold港口海水
% param.coef_c = 2.190;
% param.coef_a = 0.366;
% param.coef_b=1.824;

% --- Fournier-Forand (FF) 相函数物理参数 ---
param.n_water = 1.1549;  
param.mu_water = 3.5688; 

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
fprintf('正在预计算真实相函数 %s 的归一化系数...\n', param.phase_func);
th_test = logspace(log10(1e-5), log10(pi), 5000); 
val_test_emp = zeros(size(th_test));
val_test_ff = zeros(size(th_test));

for i = 1:length(th_test)
    % 经验公式
    t_deg = max(th_test(i) * 180 / pi, 1e-6); 
    term = 1 + (-1)^1 * param.k1 * t_deg^0.5 + ...
           (-1)^2 * param.k2 * t_deg^1.0 + ...
           (-1)^3 * param.k3 * t_deg^1.5 + ...
           (-1)^4 * param.k4 * t_deg^2.0 + ...
           (-1)^5 * param.k5 * t_deg^2.5;
    val_test_emp(i) = exp(param.q_e * term);
    
    % FF 公式
    val_test_ff(i) = calc_FF_raw(cos(th_test(i)), param);
end

param.b_emp_norm = 2 * pi * trapz(th_test, val_test_emp .* sin(th_test));
param.b_ff_norm = 2 * pi * trapz(th_test, val_test_ff .* sin(th_test));

% --- 强制指定 HG-PIS 的不对称因子 ---
param.g_prop = 0.97;
fprintf('已设定 HG 提议分布不对称因子 g = %.4f\n', param.g_prop);

% --- 强湍流参数 (OTOPS 海洋湍流谱参数, 参考论文 Table 1) ---
T_avg = 20;       % 平均温度 (保持不变)
S_avg = 35;       % 平均盐度 (保持不变)
H_ratio = -1;     % 温盐贡献比 (强湍流)
epsilon = 1e-10;  % 湍动能耗散率 (强湍流)
chi_T = 1e-5;     % 均方温度耗散率 (强湍流)
eta = 1e-3;       % 湍流内尺度 (保持不变)

%弱湍流
% T_avg = 20; 
% S_avg = 35; 
% epsilon = 1e-6; 
% chi_T = 1e-8; 
% eta = 1e-3; 
% H_ratio = -20;

N_screens = 20; 
D_screen = 1; 
N_grid = 2^8;
delta_z_screen = Link_Dist / N_screens;

% --- 动态时间轴设置  ---
t_LOS = Link_Dist / param.c_water;
delta_t_max = 2e-9;
dt = 0.01e-9; % 时间分辨率 
t_min = t_LOS; 
t_max = t_LOS + delta_t_max;
param.T_bins = t_min : dt : t_max;
param.Delta_T_bins = param.T_bins - t_LOS;
N_bins = length(param.T_bins);

% --- 仿真控制 ---
N_packets = 1e6; % 若需曲线完全平滑，建议在测试通过后增加至 1e6 或 1e7
scenarios = {'With Turbulence', 'Without Turbulence'};
results = struct();

% 计算水体固有等效发散角 (用于湍流展宽)
th_test_w = linspace(1e-6, pi/4, 2000); 
p_test = zeros(size(th_test_w));
for idx = 1:length(th_test_w)
    p_test(idx) = pdf_phase_func(cos(th_test_w(idx)), param); 
end
param.theta_water = sqrt(2 * pi * trapz(th_test_w, (th_test_w.^2) .* p_test .* sin(th_test_w)));

%% ================= 3. 预生成冻结湍流场 (移出主循环) =================
dx = D_screen / N_grid; 
x_axis = (-N_grid/2 : N_grid/2-1) * dx;
[Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);

fprintf('预计算 OTOPS 湍流相干参数...\n');
[rho0_Link, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, Link_Dist);
fprintf('  基于 OTOPS 提取的等效 Cn2 = %.2e m^(-2/3)\n', Cn2_eq);
fprintf('  链路终点精确相干长度 rho_0 = %.4f m\n\n', rho0_Link);

n_vec = Link_Dir;
if abs(n_vec(3)) < 0.9
    up_temp = [0, 0, 1]; 
else
    up_temp = [1, 0, 0]; 
end
u_vec = cross(up_temp, n_vec); 
u_vec = u_vec / norm(u_vec);
v_vec = cross(n_vec, u_vec);   
v_vec = v_vec / norm(v_vec);

fprintf('正在预生成单一微观冻结湍流相位屏列...\n');
rng('shuffle'); % 使用系统时间随机化该次微观湍流的形貌

Screen_Chain_Turb = repmat(struct('Center', [0,0,0], 'Normal', [0,0,1], ...
    'u_vec', u_vec, 'v_vec', v_vec, 'grad_x', [], 'grad_y', []), 1, N_screens);
Screen_Chain_Clean = Screen_Chain_Turb;

for i = 1:N_screens
    Screen_Chain_Turb(i).Center = Tx_Pos + i * delta_z_screen * Link_Dir; 
    Screen_Chain_Turb(i).Normal = Link_Dir;
    Screen_Chain_Clean(i).Center = Screen_Chain_Turb(i).Center;
    Screen_Chain_Clean(i).Normal = Link_Dir;
    
    % 生成具有湍流的相位屏梯度
    phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
    [gx, gy] = gradient(phi, dx);
    Screen_Chain_Turb(i).grad_x = gx; 
    Screen_Chain_Turb(i).grad_y = gy;
    
    % 生成无湍流时的零梯度矩阵
    Screen_Chain_Clean(i).grad_x = zeros(N_grid, N_grid); 
    Screen_Chain_Clean(i).grad_y = zeros(N_grid, N_grid);
end

%% ================= 4. 对比仿真主循环 =================
for s_idx = 1:2
    scenario_name = scenarios{s_idx};
    fprintf('=== Running Scenario: %s ===\n', scenario_name);
    
    % 根据场景装载相位屏
    if strcmp(scenario_name, 'With Turbulence')
        Screen_Chain_Current = Screen_Chain_Turb;
    else
        Screen_Chain_Current = Screen_Chain_Clean;
    end
    
    % [核心重点]：绝对对齐 PRNG 种子，保证两组光子抽取的步长与散射角完全相同
    rng(2026, 'twister'); 
    
    h_time = zeros(1, N_bins); 
    tic;
    for p = 1:N_packets
        
        % --- 4.0 光源空间与角度高斯分布抽样 (Gaussian Beam Profile) ---
        r0 = w0 * sqrt(-0.5 * log(rand())); 
        phi0 = 2 * pi * rand();
        pos_local = r0 * cos(phi0) * u_vec_Tx + r0 * sin(phi0) * v_vec_Tx;
        pos_init = Tx_Pos + pos_local;
        
        theta_half_div_physics = div_angle/2; 
        U = theta_half_div_physics * sqrt(-0.5 * log(rand()));
        xi_psi = rand();
        psi_ini = 2 * pi * xi_psi;
        dir_init = rotate_direction(mu_T, U, psi_ini);
        
        weight_init = 1.0; 
        Huge_Aperture = 1e5; 
        
        % --- 4.1 [Ballistic Branch] 弹道分量 ---
        [pos_end, ~, plane_hit, path_len_ballistic] = ray_march_generic(pos_init, dir_init, 1e9, Rx_Pos, ...
            Huge_Aperture, pi, Rx_Normal, true, Screen_Chain_Current, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);

        if plane_hit
            % 【修正1】获取入射角的余弦投影值
            cos_theta_inc = abs(dot(dir_init, Rx_Normal));

            % 仅当满足接收视场角(FOV)限制时才进行能量计算
            if acos(cos_theta_inc) <= Rx_FOV / 2

                % 【修正2】计算正交截面上的真实离轴偏心距 (Orthogonal Beam Wander)
                % 几何原理：利用空间向量叉乘，求取 Rx_Pos 到射线 (pos_end, dir_init) 的垂直距离
                r_wander_perp = norm(cross(pos_end - Rx_Pos, dir_init));

                % 【修正3】计算接收孔径在光束正交截面上的等效投影半径
                r_eff = Rx_Radius * sqrt(cos_theta_inc);

                % 束散半径 W_spot 仍基于传播距离在正交截面上计算，维持不变
                W_spot = calc_beam_spot_size(w0, lambda, path_len_ballistic, strcmp(scenario_name, 'With Turbulence')*Cn2_eq, theta_half_div_physics);

                % 利用 Marcum-Q 函数在正交截面上进行二维高斯能量积分
                if W_spot > 1e-6
                    received_fraction = 1 - marcumq(2 * r_wander_perp / W_spot, 2 * r_eff / W_spot);
                else
                    % 若无展宽效应，则蜕化为等效投影面积内的点包容测试
                    received_fraction = double(r_wander_perp <= r_eff);
                end

                attenuation = exp(-param.coef_c * path_len_ballistic);
                bin_idx = floor((path_len_ballistic / param.c_water - t_min) / dt) + 1;

                if bin_idx >= 1 && bin_idx <= N_bins
                    h_time(bin_idx) = h_time(bin_idx) + weight_init * attenuation * received_fraction; 
                end
            end
        end

        % --- 4.2 [Scattering Branch] PIS 追踪 (基于物理 b 游走) ---
        pos = pos_init; 
        dir = dir_init; 
        current_dist = 0;
        weight_pis_likelihood = weight_init; 
        
        % 初始物理游走
        d_step = -log(rand()) / param.coef_b;
        [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, Rx_Pos, Huge_Aperture, Rx_FOV, ...
            Rx_Normal, false, Screen_Chain_Current, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
        current_dist = current_dist + step_len;
        
        if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist
            continue; 
        end
        
        for order = 1 : param.n_max
            % ---- 强制接收评估 (Local Estimation) ----
            vec_rx = Rx_Pos - pos; 
            dist_rx = norm(vec_rx); 
            dir_rx = vec_rx / dist_rx;
            
            if acos(dot(-dir_rx, Rx_Normal)) <= Rx_FOV/2
                p_phase = pdf_phase_func(dot(dir, dir_rx), param);
                omega = Rx_Area / (dist_rx^2) * abs(dot(dir_rx, Rx_Normal));
                
                % 当前物理真实权重 = 累积似然比 * 物理吸收
                current_physical_weight = weight_pis_likelihood * exp(-param.coef_a * current_dist);
                base_energy = current_physical_weight * min(1, p_phase * omega) * exp(-param.coef_c * dist_rx);
                
                if base_energy > 1e-15
                    if strcmp(scenario_name, 'With Turbulence')
                        % 1. 利用相位屏求取真实的质心偏折落点
                        [pos_v_end, ~, v_hit, v_len] = ray_march_generic(pos, dir_rx, dist_rx + 1e-1, Rx_Pos, ...
                            Huge_Aperture, pi, Rx_Normal, true, Screen_Chain_Current, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);

                        if v_hit
                            % 2. 考虑斜入射的余弦投影，计算等效接收半径
                            cos_theta_inc = abs(dot(dir_rx, Rx_Normal));
                            r_eff = Rx_Radius * sqrt(cos_theta_inc);

                            % 3. 【修正项1】计算正交截面上的真实离轴偏心距 (Orthogonal Beam Wander)
                            % 物理映射：将落点误差投影到垂直于局部光子传播方向的二维截面上
                            r_wander_perp = norm(cross(pos_v_end - Rx_Pos, dir_rx));

                            % 4. 【修正项2】计算球面波短曝光湍流展宽 (W_ST)
                            % 剔除无物理意义的 W_geo_PIS，纯粹考量次级点光源的湍流弥散
                            rho_0_spherical = rho0_Link * (Link_Dist / v_len)^(3/5); 
                            
                            % 【修正项3】利用等效孔径直径 (2*r_eff) 替代物理孔径，参与倾斜剥除 (Tilt Removal) 计算
                            W_ST = (2 * v_len / (k_wave * rho_0_spherical)) * max(0, 1 - 0.37*(rho_0_spherical/(2*r_eff))^(1/3));

                            % 5. 离轴偏心高斯积分 (二维空间 PSF 能量截获率计算)
                            if W_ST > 1e-6
                                f_spread = 1 - marcumq(2 * r_wander_perp / W_ST, 2 * r_eff / W_ST);
                            else
                                f_spread = double(r_wander_perp <= r_eff); % 若无湍流展宽，蜕化为等效投影面积内的点包容测试
                            end

                            % 6. 更新基础能量权重与时间仓
                            base_energy = base_energy * f_spread;
                            bin_idx = floor(((current_dist + v_len) / param.c_water - t_min) / dt) + 1;
                        end
                    else
                        bin_idx = floor(((current_dist + dist_rx) / param.c_water - t_min) / dt) + 1;
                    end
                
                    if exist('bin_idx','var') && bin_idx >= 1 && bin_idx <= N_bins
                        h_time(bin_idx) = h_time(bin_idx) + base_energy; 
                    end
                    clear bin_idx; 
                end
            end
            
            % 俄罗斯轮盘赌控制算力 (使用物理期望权重)
            current_physical_weight = weight_pis_likelihood * exp(-param.coef_a * current_dist);
            if current_physical_weight < 1e-9
                if rand() > 0.1
                    break; 
                else
                    weight_pis_likelihood = weight_pis_likelihood / 0.1; 
                end
            end
            
            % ---- HG-PIS 生成提议散射角 ----
            g_prop = param.g_prop;
            xi_rand = rand();
            term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi_rand);
            cos_t = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1);
            theta_i = acos(cos_t); 
            phi_scat = 2 * pi * rand();
            
            % 计算提议概率与真实概率
            q_HG = (1 - g_prop^2) / (4 * pi * (1 + g_prop^2 - 2 * g_prop * cos_t)^1.5);
            p_val = pdf_phase_func(cos_t, param);
            
            % 更新代数似然比补偿 (温和截断)
            weight_factor = min(1e3, p_val / q_HG);
            weight_pis_likelihood = weight_pis_likelihood * weight_factor;
            
            % 坐标旋转与下一步游走
            dir = rotate_direction(dir, theta_i, phi_scat);
            d_step = -log(rand()) / param.coef_b;
            [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, Rx_Pos, Huge_Aperture, Rx_FOV, ...
                Rx_Normal, false, Screen_Chain_Current, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
            
            current_dist = current_dist + step_len;
            if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist
                break; 
            end
        end
    end
    t_run = toc; 
    fprintf('   Completed in %.2f s\n', t_run);
    
    % --- 4.3 计算统计参数 ---
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

figure('Name', 'Impulse Response Comparison (Aligned to Lit.)', 'Color', 'w');
semilogy(param.Delta_T_bins * 1e9, results(1).h_time_cir, 'r-', 'LineWidth', 1.5, 'DisplayName', 'With Turb'); 
hold on;
semilogy(param.Delta_T_bins * 1e9, results(2).h_time_cir, 'b.-', 'LineWidth', 1.0, 'MarkerSize', 8, 'DisplayName', 'No Turb');
grid on; 
legend('Location', 'best'); 
xlabel('Time Increment \Delta t (ns)', 'FontWeight', 'bold'); 
ylabel('Intensity (Discrete)', 'FontWeight', 'bold');
title(sprintf('Impulse Response (%s, HG-PIS g=%.2f)', param.phase_func, param.g_prop));

%% ================= 6. 核心辅助函数区域 =================
% 统一相函数分配器
function p = pdf_phase_func(cos_theta, param)
    if strcmp(param.phase_func, 'Empirical')
        p = pdf_Empirical(cos_theta, param);
    elseif strcmp(param.phase_func, 'FFPF')
        p = pdf_FFPF(cos_theta, param);
    else
        error('未定义的相函数模型');
    end
end

function p = pdf_Empirical(cos_theta, param)
    cos_theta = max(min(cos_theta, 1), -1); 
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6);
    term = 1 + (-1)^1 * param.k1 * t_deg^0.5 + ...
           (-1)^2 * param.k2 * t_deg^1.0 + ...
           (-1)^3 * param.k3 * t_deg^1.5 + ...
           (-1)^4 * param.k4 * t_deg^2.0 + ...
           (-1)^5 * param.k5 * t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end

function p = pdf_FFPF(cos_theta, param)
    p = calc_FF_raw(cos_theta, param) / param.b_ff_norm;
end

function P_raw = calc_FF_raw(cos_theta, param)
    cos_theta = max(min(cos_theta, 1), -1); 
    theta = acos(cos_theta);
    nu = (3 - param.mu_water) / 2;
    delta = max((4 / (3 * (param.n_water - 1)^2)) .* sin(theta/2).^2, 1e-12);
    delta_180 = (4 / (3 * (param.n_water - 1)^2)); 
    
    term1 = nu .* (1 - delta) - (1 - delta.^nu) + ...
            (delta .* (1 - delta.^nu) - nu .* (1 - delta)) .* (sin(theta/2).^(-2));
            
    P_raw = max(0, (1 ./ (4 * pi .* (1 - delta).^2 .* delta.^nu)) .* term1 + ...
            ((1 - delta_180^nu) / (16 * pi * (delta_180 - 1) * delta_180^nu)) .* (3 .* cos(theta).^2 - 1));
end

function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, Rx_Pos, Rx_Aperture, Rx_FOV, ...
    Rx_Normal, enable_hit_check, Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen) 
    
    hit_flag = false; 
    step_actual = dist_limit; 
    hit_rx_plane = false; 
    N_screens = length(Screen_Chain);
    
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
            scr = Screen_Chain(i);
            t_i = (dot(scr.Center - Tx_Pos, Link_Dir) - z_start) / dir_z;
            pos_i = pos + dir_old * t_i;
            
            loc_u = dot(pos_i - scr.Center, scr.u_vec); 
            loc_v = dot(pos_i - scr.Center, scr.v_vec);
            
            idx_x = mod(round((loc_u - x_axis(1)) / dx), N_grid) + 1; 
            idx_y = mod(round((loc_v - x_axis(1)) / dx), N_grid) + 1;
            
            delta_dir = delta_dir + (scr.grad_x(idx_y, idx_x) * scr.u_vec + scr.grad_y(idx_y, idx_x) * scr.v_vec);
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

function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    lambda_nm = lambda * 1e9; 
    k_wave = 2 * pi / lambda; 
    
    A = -1.05e-6 * S + 2 * 1.6e-8 * T * S + 2 * -2.02e-6 * T + -4.23e-3 / lambda_nm; 
    B = 1.779e-4 + -1.05e-6 * T + 1.6e-8 * T^2 + 1.155e-2 / lambda_nm; 
    T_k = T + 273.15; 
    s_frac = S * 1e-3; 
    
    cp = 1000 * ((5.328 - 9.76e-2 * S + 4.04e-4 * S^2) + ...
         (-6.913e-3 + 7.351e-4 * S - 3.15e-6 * S^2) * T + ...
         (9.6e-6 - 1.927e-6 * S + 8.23e-9 * S^2) * T^2 + ...
         (2.5e-9 + 1.666e-9 * S - 7.125e-12 * S^2) * T^3); 
         
    rho = (9.9992293295e2 + 2.0341179217e-2 * T - 6.1624591598e-3 * T^2 + ...
          2.2614664708e-5 * T^3 - 4.6570659168e-8 * T^4) + ...
          s_frac * (8.0200240891e2 - 2.0005183488 * T + 1.6771024982e-2 * T^2 - ...
          3.0600536746e-5 * T^3 - 1.6132224742e-5 * T * S);
          
    mu = ((0.15700386464 * (T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5) * ...
         (1 + (1.5409136040 + 1.9981117208e-2 * T - 9.5203865864e-5 * T^2) * s_frac + ...
         (7.9739318223 - 7.561456881e-2 * T + 4.7237011074e-4 * T^2) * s_frac^2);
         
    sigma_T = 10^(log10(240 + 0.0002 * (S / 1.00472)) - 3 + 0.434 * ...
              (2.3 - (343.5 + 0.037 * (S / 1.00472)) / (1.00024 * T + 273.15)) * ...
              (1 - (1.00024 * T + 273.15) / (647.3 + 0.03 * (S / 1.00472)))^(1/3)); 
              
    Pr = mu * cp / sigma_T; 
    Sc = mu^2 / (5.954e-15 * T_k * rho); 
    
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1); 
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    R_rho = 2.6e-4 * abs(H_ratio) / 7.6e-4;
    
    if R_rho >= 1
        d_r = R_rho + sqrt(R_rho) * sqrt(R_rho - 1); 
    elseif R_rho >= 0.5
        d_r = 1.85 * R_rho - 0.85; 
    else
        d_r = 0.15 * R_rho; 
    end
    
    chi_S = chi_T * d_r / (H_ratio^2); 
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    
    Phi_Hill = @(K, chi_M, c_M) (0.72 / (4 * pi)) * chi_M * epsilon^(-1/3) .* ...
               (K.^2 + 1e-25).^(-11/6) .* exp(-176.90 * (K * eta).^2 .* c_M^(0.96)) .* ...
               (1 + 21.61 * (K * eta).^(0.61) .* c_M^(0.02) - 18.18 * (K * eta).^(0.55) .* c_M^(0.04));
               
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + B^2 * Phi_Hill(K, chi_S, c_S) + 2 * A * B * Phi_Hill(K, chi_TS, c_TS));
end

function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2 * pi / D; 
    dx = D / N; 
    kx = (-N/2 : N/2-1) * dk;
    
    [KX, KY] = meshgrid(kx, kx); 
    K_grid = sqrt(KX.^2 + KY.^2); 
    K_grid(N/2+1, N/2+1) = 1e-10; 
    
    C_nm = (randn(N) + 1i * randn(N)) .* sqrt(2 * pi * k_wave^2 * delta_z * Phi_n_func(K_grid)) * dk;
    phase_high = real(ifftshift(ifft2(ifftshift(C_nm)))) * N^2; 
    phase_low = zeros(N, N); 
    [xx, yy] = meshgrid((-N/2 : N/2-1) * dx); 
    
    for p = 1:3
        dk_p = dk / (3^p); 
        for m = -1:1
            for n = -1:1
                if (m == 0 && n == 0)
                    continue; 
                end
                k_p = sqrt((m * dk_p)^2 + (n * dk_p)^2);
                phase_low = phase_low + real((randn(1) + 1i * randn(1)) * ...
                            sqrt(2 * pi^2 * k_wave^2 * delta_z * Phi_n_func(k_p)) * ...
                            dk_p * exp(1i * (m * dk_p * xx + n * dk_p * yy)));
            end
        end
    end
    phase_screen = phase_high + phase_low;
end

function new_dir = rotate_direction(dir, theta_s, psi_s)
    denom = sqrt(1 - dir(3)^2);
    if denom < 1e-10
        if dir(3) > 0
            new_dir = [sin(theta_s) * cos(psi_s), sin(theta_s) * sin(psi_s), cos(theta_s)]; 
        else
            new_dir = [sin(theta_s) * cos(psi_s), sin(theta_s) * sin(psi_s), -cos(theta_s)]; 
        end
    else
        new_dir = [sin(theta_s)/denom*(dir(1)*dir(3)*cos(psi_s)-dir(2)*sin(psi_s))+dir(1)*cos(theta_s), ...
                   sin(theta_s)/denom*(dir(2)*dir(3)*cos(psi_s)+dir(1)*sin(psi_s))+dir(2)*cos(theta_s), ...
                   -sin(theta_s)*cos(psi_s)*denom+dir(3)*cos(theta_s)];
    end
    new_dir = new_dir / norm(new_dir);
end

function W_L = calc_beam_spot_size(w0, lambda, L, Cn2,halfdiv)
    if w0 <= 1e-6
        W_L = 0; 
        return; 
    end 
    
    k = 2 * pi / lambda; 
    D = 2 * w0; 
    % 用几何投影半径 W_geo 替代原来的 W_diff
    W_geo = w0 + L * tan(halfdiv);
    % z_R = (pi * w0^2) / lambda;
    % W_diff = w0 * sqrt(1 + (L / z_R)^2);
    
    if Cn2 > 1e-15
        rho_0 = (0.545 * k^2 * Cn2 * L)^(-3/5);
        W_turb_ST = (2 * L / (k * rho_0)) * max(1 - 0.37 * (rho_0 / D)^(1/3), 0);
        W_L = sqrt(W_geo^2 + W_turb_ST^2);
        %W_L = sqrt(W_diff^2 + W_turb_ST^2);
    else
        W_L = W_geo; 
    end
end

function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L)
    kappa_eval = 100; 
    Cn2_eq = Phi_n_func(kappa_eval) / (0.033 * kappa_eval^(-11/3));
    
    try
        opts = optimset('Display', 'off'); 
        rho0_exact = fzero(@(rho) 8 * pi^2 * k_wave^2 * L * ...
                     integral2(@(K, xi) K .* Phi_n_func(K) .* (1 - besselj(0, K .* rho .* xi)), ...
                     1e-1, 1e4, 0, 1, 'Method', 'iterated', 'RelTol', 1e-3) - 2, ...
                     (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5), opts); 
    catch
        rho0_exact = (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5); 
    end
end