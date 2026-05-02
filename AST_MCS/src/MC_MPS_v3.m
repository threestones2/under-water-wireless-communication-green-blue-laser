%% 水下/大气光通信混合仿真: MC-MPS (Turbulence Comparison)
%  Fournier-Forand (FF) 相函数 / Haltrin 经验拟合切换
%  物理散射步进解耦 + O(1)散射角查表法 (LUT)
clc; 
clear; 
close all;

%% ================= 1. 参数初始化 =================
param = struct();
% --- 核心模式选择 (相函数开关) ---
% 可选: 'FFPF' (Fournier-Forand模型) 或 'Empirical' (Haltrin经验拟合)
param.phase_func = 'FFPF';  
% 蒙特卡洛追踪的最大多重散射阶数 (对齐论文图5的长尾要求)
param.n_max = 20;

% --- 光源与波长参数 ---
lambda_nm = 514;
lambda = lambda_nm * 1e-9;
k_wave = 2 * pi / lambda;
% 为严格对齐论文冲激响应基准测试，光源退化为理想准直点光源
w0 = 0.0;                       
div_angle = 0.0;           
z_R = realmin; % 避免除零

% --- 3D 空间链路布局 ---
Tx_Pos = [0, 0, 0];         
Rx_Pos = [0, 10.93, 0];       
Link_Vec = Rx_Pos - Tx_Pos;
Link_Dist = norm(Link_Vec);
Link_Dir = Link_Vec / Link_Dist;

% --- 光轴对准偏差 ---
theta_tx_error = 0;      
phi_tx_error = 0;           
theta_rx_error = 0;      
phi_rx_error = pi;          
mu_T = rotate_direction(Link_Dir, theta_tx_error, phi_tx_error);            
Rx_Normal = rotate_direction(-Link_Dir, theta_rx_error, phi_rx_error);

if abs(mu_T(3)) < 0.9
    up_temp_Tx = [0, 0, 1]; 
else
    up_temp_Tx = [1, 0, 0]; 
end
u_vec_Tx = cross(up_temp_Tx, mu_T); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(mu_T, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

% --- 接收端物理参数 ---
Rx_Aperture = 0.4; % 对齐论文 40 cm 孔径
Rx_Radius = Rx_Aperture / 2;
Rx_FOV = 20 * pi / 180;
Rx_Area = pi * (Rx_Radius)^2;

% --- 介质参数 (港口水质 Harbor Water) ---
param.c_water = 2.237e8; % 对齐论文水中光速  
param.coef_a = 0.366;    % 吸收系数
param.coef_b = 1.829;    % 散射系数
param.coef_c = param.coef_a + param.coef_b; % 衰减系数

% --- 弱湍流参数 (OTOPS 海洋湍流谱参数) ---
T_avg = 20; S_avg = 35; H_ratio = -20;                  
epsilon = 1e-9; chi_T = 1e-7; eta = 1e-3;                     
N_screens = 20; D_screen = 10; N_grid = 2^8;                   
delta_z_screen = Link_Dist / N_screens; 

% --- 动态时间轴设置 (增量到达时间 \Delta t) ---
t_LOS = Link_Dist / param.c_water; 
delta_t_max = 2.5e-9;  % 观察窗口设为 2.5 ns 以对齐图5                       
dt = 2.5e-11;          % 时间分辨率 (0.025 ns，提高直方图分辨率)
t_min = t_LOS;                  
t_max = t_LOS + delta_t_max;    
param.T_bins = t_min : dt : t_max;
param.Delta_T_bins = param.T_bins - t_LOS;  
N_bins = length(param.T_bins);

% --- 仿真控制 ---
N_packets = 1e6; % 调大光子数以弥补关闭 PIS 后的统计波动                
scenarios = {'With Turbulence', 'Without Turbulence'};
results = struct(); 

%% ================= 2. 核心：相函数预计算与 LUT 查表法构建 =================
fprintf('正在预计算相函数 %s 并构建 O(1) 查找表 (LUT)...\n', param.phase_func);
% 离散化天顶角 (1e-5 避免奇异点)
theta_array = logspace(log10(1e-5), log10(pi), 5000); 
if strcmp(param.phase_func, 'FFPF')
    n_best = 1.05;   % 典型悬浮颗粒相对水折射率
    mu_best = 3.65;  % 典型粒径分布斜率
    P_raw = calc_FF_Core(theta_array, n_best, mu_best);
else
    P_raw = calc_Empirical_Core(theta_array, param.coef_b, param.coef_b/param.coef_c);
end

% 严格归一化并计算累积分布函数 (CDF)
b_integral = 2 * pi * trapz(theta_array, P_raw .* sin(theta_array));
P_norm = 4 * pi * P_raw ./ b_integral; 
pdf_theta = 0.5 * P_norm .* sin(theta_array); % 几何概率密度
cdf_theta = cumtrapz(theta_array, pdf_theta);
cdf_theta = cdf_theta / cdf_theta(end); % 强制末端归一化为1

% 构建 O(1) 反变换查找表 (LUT)
M_LUT = 1e6; 
P_grid = linspace(0, 1, M_LUT); 
param.theta_LUT = interp1(cdf_theta, theta_array, P_grid, 'linear', 'extrap');

%% ================= 3. 对比仿真主循环 =================
dx = D_screen / N_grid;
x_axis = (-N_grid/2 : N_grid/2-1) * dx;
[Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
fprintf('预计算 OTOPS 湍流相干参数...\n');
[rho0_Link, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, Link_Dist);
fprintf('  基于 OTOPS 提取的等效 Cn2 = %.2e m^(-2/3)\n', Cn2_eq);
fprintf('  链路终点精确相干长度 rho_0 = %.4f m\n\n', rho0_Link);

n_vec = Link_Dir;
if abs(n_vec(3)) < 0.9; up_temp = [0, 0, 1]; else; up_temp = [1, 0, 0]; end
u_vec = cross(up_temp, n_vec); u_vec = u_vec / norm(u_vec);
v_vec = cross(n_vec, u_vec);   v_vec = v_vec / norm(v_vec);

for s_idx = 1:2
    scenario_name = scenarios{s_idx};
    fprintf('=== Running Scenario: %s ===\n', scenario_name);
    Screen_Chain = repmat(struct('Center', [0,0,0], 'Normal', [0,0,1], ...
                             'u_vec', [], 'v_vec', [], 'grad_x', [], 'grad_y', []), 1, N_screens);
    for i = 1:N_screens
        Screen_Chain(i).Center = Tx_Pos + i * delta_z_screen * Link_Dir;
        Screen_Chain(i).Normal = Link_Dir;
        Screen_Chain(i).u_vec = u_vec; Screen_Chain(i).v_vec = v_vec;
        if strcmp(scenario_name, 'With Turbulence')
            phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
            [gx, gy] = gradient(phi, dx);
            Screen_Chain(i).grad_x = gx; Screen_Chain(i).grad_y = gy;
        else
            Screen_Chain(i).grad_x = zeros(N_grid, N_grid); Screen_Chain(i).grad_y = zeros(N_grid, N_grid);
        end
    end

    h_time = zeros(1, N_bins); 
    tic;
    for p = 1:N_packets
        % --- 3.1 微观粒子状态初始化 ---
        pos_init = Tx_Pos;
        dir_init = mu_T;
        weight_init = 1.0;
        Huge_Aperture = 1e5; 

        % --- 3.2 [Ballistic Branch] 弹道光子宏观质心解耦评估 ---
        pos_centroid = Tx_Pos; 
        dir_centroid = mu_T; 
        [pos_end, ~, plane_hit, path_len_ballistic] = ray_march_generic(pos_centroid, dir_centroid, 1e9, ...
            Rx_Pos, Huge_Aperture, pi, Rx_Normal, true, ...
            Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
        if plane_hit
            r_wander = norm(pos_end - Rx_Pos);
            if strcmp(scenario_name, 'With Turbulence')
                 W_spot = calc_beam_spot_size(w0, lambda, path_len_ballistic, Cn2_eq);
            else
                 W_spot = calc_beam_spot_size(w0, lambda, path_len_ballistic, 0);
            end
            cos_rx_tilt = abs(dot(dir_centroid, Rx_Normal));
            if acos(cos_rx_tilt) <= Rx_FOV / 2
                if W_spot > 1e-6
                    r_eff = Rx_Radius * sqrt(cos_rx_tilt);
                    a_param = 2 * r_wander / W_spot;
                    b_param = 2 * r_eff / W_spot;
                    received_fraction = 1 - marcumq(a_param, b_param);
                else
                    received_fraction = double(r_wander <= Rx_Radius);
                end
                attenuation = exp(-param.coef_c * path_len_ballistic);
                final_weight = weight_init * attenuation * received_fraction;
                t_arrival = path_len_ballistic / param.c_water;
                bin_idx = floor((t_arrival - t_min) / dt) + 1;
                if bin_idx >= 1 && bin_idx <= N_bins
                    h_time(bin_idx) = h_time(bin_idx) + final_weight;
                end
            end
        end

        % --- 3.3 [Scattering Branch] 散射光微观演化 (纯物理追踪) ---
        pos = pos_init; 
        dir = dir_init; 
        path_len_total = 0;
        for order = 1 : param.n_max
            d_step = -log(rand()) / param.coef_b; 
            [pos_new, dir_new, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
                Rx_Pos, Huge_Aperture, Rx_FOV, Rx_Normal, false, ...
                Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);

            dist_y_old = dot(pos - Tx_Pos, Link_Dir);
            dist_y_new = dot(pos_new - Tx_Pos, Link_Dir);
            if dist_y_old < Link_Dist && dist_y_new >= Link_Dist
                d_remain = (Link_Dist - dist_y_old) / dot(dir_new, Link_Dir);
                hit_pos = pos + dir_new * d_remain;
                hit_path_len = path_len_total + d_remain;
                final_weight = weight_init * exp(-param.coef_a * hit_path_len);
                r_hit = norm(hit_pos - Rx_Pos);
                cos_theta_inc = abs(dot(dir_new, Rx_Normal));
                if (r_hit <= Rx_Radius) && (acos(cos_theta_inc) <= Rx_FOV/2)
                    t_arrival = hit_path_len / param.c_water;
                    bin_idx = floor((t_arrival - t_min) / dt) + 1;
                    if bin_idx >= 1 && bin_idx <= N_bins
                        h_time(bin_idx) = h_time(bin_idx) + final_weight;
                    end
                end
                break; 
            end

            pos = pos_new; 
            dir = dir_new;
            path_len_total = path_len_total + step_len;

            current_weight = weight_init * exp(-param.coef_a * path_len_total);
            if current_weight < 1e-7
                if rand() <= 0.1
                    weight_init = weight_init / 0.1;
                else
                    break; 
                end
            end

            idx_LUT = floor(rand() * (M_LUT - 1)) + 1;
            theta_s = param.theta_LUT(idx_LUT);
            phi_scat = 2 * pi * rand();
            dir = rotate_direction(dir, theta_s, phi_scat);
        end
    end
    t_run = toc;
    fprintf('   Completed in %.2f s\n', t_run);

    % --- 3.4 存储并计算信道参数与 CIR ---
    h_time_norm = h_time / N_packets;      
    dt_ns = dt * 1e9;

    % 1. 计算总接收概率并推导路径损耗 (加入 eps 防止接收全0时对数报错)
    P_rx_total = sum(h_time_norm);
    results(s_idx).P_rx = P_rx_total;
    results(s_idx).Path_Loss_dB = -10 * log10(P_rx_total + eps);

    % 2. 计算均方根延迟拓展 (RMS Delay Spread)
    if P_rx_total > 0
        tau_mean = sum(param.T_bins .* h_time_norm) / P_rx_total;
        tau_rms = sqrt( sum( ((param.T_bins - tau_mean).^2) .* h_time_norm ) / P_rx_total );
        results(s_idx).tau_rms = tau_rms;
    else
        results(s_idx).tau_rms = 0;
    end

    % 3. 存储归一化后的密度函数 (为绘图乘回 dt_ns 做准备)
    results(s_idx).h_time_cir = h_time_norm;
end

%% ================= 4. 控制台定量输出与绘图 =================
% 打印终端对比表格
fprintf('\n=== 仿真对比结果 (Tx-Rx: %.1fm) ===\n', Link_Dist);
fprintf('%-20s | %-14s | %-14s\n', 'Scenario', 'Path Loss (dB)', 'RMS Delay (ns)');
fprintf('------------------------------------------------------\n');
for i = 1:2
    fprintf('%-20s | %14.4f | %14.4f\n', scenarios{i}, results(i).Path_Loss_dB, results(i).tau_rms * 1e9);
end

% 绘制冲激响应图
figure('Name', 'Impulse Response Comparison', 'Color', 'w');
% 注意：此处将概率质量除以 dt_ns 转换为物理密度函数 [ns^-1] 进行展示
plot(param.Delta_T_bins * 1e9, results(1).h_time_cir , 'r.-', 'LineWidth', 1, 'DisplayName', 'With Turbulence');
hold on;
plot(param.Delta_T_bins * 1e9, results(2).h_time_cir , 'bo', 'MarkerSize', 4, 'DisplayName', 'Without Turbulence');
set(gca, 'YScale', 'log'); 
grid on;
legend('Location', 'best');
xlim([0, 2.5]);
%ylim([1e-9, 1e-6]);
xlabel('\Delta t [ns]', 'FontWeight', 'bold');
ylabel('Intensity [ns^{-1}]', 'FontWeight', 'bold');
title(sprintf('UOWC Impulse Response (Harbor, %s Model)', param.phase_func));

%% ================= 4.5 冲激响应 (IRF) 依据真实文献公式与 Table 1 理论对比 =================
% 真实拟合公式: h(\Delta t) = C1 * (\Delta t^\alpha) / (\Delta t + C2)^\beta * exp(-a * v * (\Delta t + t0))

fprintf('\n=== IRF 真实文献公式 (前 0.5 ns 截断) 拟合结果 ===\n');

% --- 提取系统物理常量 ---
v_water = param.c_water;      % 水中光速 (m/s)
a_coef = param.coef_a;        % 吸收系数 (m^-1)
t0_s = Link_Dist / v_water;   % 直达光飞行时间 (s)

t0_ns = t0_s * 1e9;               % t0 转换为纳秒
a_v_ns = a_coef * v_water * 1e-9; % a*v 乘积项转换为 ns^-1 尺度

% ----------------- 理论基准曲线生成区 -----------------
% 填入文献 Table 1 中对应水质的真实参数
C1_tab    = 1.677e-5; 
C2_tab    = 0.2730;  
alpha_tab = 0.6577;  
beta_tab  = 3.169;  

% 设定局部时间轴 (1e-4 到 0.5 ns，避开 0 点防止奇异值)
t_ns_theory = linspace(1e-4, 1, 500); 
% 依据真实公式生成理论解析数据
h_theory = C1_tab .* (t_ns_theory.^alpha_tab) ./ ((t_ns_theory + C2_tab).^beta_tab) ...
           .* exp(-a_v_ns .* (t_ns_theory + t0_ns));
% ------------------------------------------------------

% 定义非线性拟合类型 (动态注入物理常数)
eq_str = sprintf('c1 * (x.^alpha) ./ ((x + c2).^beta) .* exp(-%.6e * (x + %.6e))', a_v_ns, t0_ns);
true_model = fittype(eq_str, ...
    'independent', 'x', ...
    'dependent', 'y', ...
    'coefficients', {'c1', 'c2', 'alpha', 'beta'});
    
% 设定非线性最小二乘拟合的优化选项
fit_opts = fitoptions('Method', 'NonlinearLeastSquares');
% 【核心修正】：将 Table 1 的物理参数直接作为寻优算法的全局初始起点
fit_opts.StartPoint = [C1_tab, C2_tab, alpha_tab, beta_tab]; 
fit_opts.Lower = [0, 0, 0, 0]; 
fit_opts.Upper = [Inf, 10, 10, 10]; % 适当增加上界约束防止梯度发散
fit_opts.MaxIter = 2000;
fit_opts.MaxFunEvals = 3000;

figure('Name', 'True Formula Fit (<0.5ns) & Table 1 Reference', 'Color', 'w', 'Position', [150, 150, 900, 400]);
colors = {'r', 'b'};

for s_idx = 1:2
    % 提取 >0 且 <= 0.5 ns 的有效数据区间
    valid_idx = param.Delta_T_bins > 0 & param.Delta_T_bins <= 0.5e-9 & results(s_idx).h_time_cir > 0;
    
    x_fit = param.Delta_T_bins(valid_idx)' * 1e9; 
    y_fit = results(s_idx).h_time_cir(valid_idx)'; 
    
    % （from reasoning）由于能量跨越数个对数量级，设置权重向量以削弱峰值对拟合的绝对统治，
    % 强迫优化器兼顾尾部指数衰落的形状拟合。
    weights = 1 ./ (y_fit + eps); 
    fit_opts.Weights = weights;
    
    try
        % 执行截断区间拟合
        [fit_result, gof] = fit(x_fit, y_fit, true_model, fit_opts);
        
        % 终端输出解析参数
        fprintf('[%s] 前 0.5ns 拟合参数:\n', scenarios{s_idx});
        fprintf('   C1 = %e\n', fit_result.c1);
        fprintf('   C2 = %.4f\n', fit_result.c2);
        fprintf('   alpha = %.4f\n', fit_result.alpha);
        fprintf('   beta = %.4f\n', fit_result.beta);
        fprintf('   拟合优度 R^2 = %.4f\n\n', gof.rsquare);
        
        % 绘制对比子图
        subplot(1, 2, s_idx);
        
        % 1. 蒙特卡洛离散统计散点
        semilogy(x_fit, y_fit, 'o', 'Color', colors{s_idx}, 'MarkerSize', 4, 'DisplayName', 'MC Data (\le 0.5ns)');
        hold on;
        
        % 2. 基于当前数据的拟合曲线
        x_smooth = linspace(min(x_fit), 0.5, 500);
        y_smooth = fit_result(x_smooth);
        y_smooth(y_smooth < 1e-20) = 1e-20; 
        semilogy(x_smooth, y_smooth, 'k-', 'LineWidth', 2, 'DisplayName', 'Fit Curve');
        
        % 3. 文献 Table 1 理论基准曲线
        h_theory_plot = h_theory;
        h_theory_plot(h_theory_plot < 1e-20) = 1e-20;
        semilogy(t_ns_theory, h_theory_plot, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Table 1 Theory');
        
        grid on;
        legend('Location', 'best');
        xlabel('Time Increment \Delta t (ns)', 'FontWeight', 'bold');
        ylabel('Intensity', 'FontWeight', 'bold');
        title(sprintf('%s\n(Fit R^2 = %.4f)', scenarios{s_idx}, gof.rsquare));
        
        xlim([0, 0.5]);
        ylim([min(y_fit)*0.1, max(y_fit)*5]);
        
    catch ME
        fprintf('[%s] 拟合失败: %s\n', scenarios{s_idx}, ME.message);
    end
end







%% ================= 5. 辅助函数区域 =================
% Fournier-Forand 核心公式 (从您的 scatter_test5.m 转化)
function P_raw = calc_FF_Core(theta, n, mu)
    nu = (3 - mu) / 2;
    delta = (2 * sin(theta./2)).^2 ./ (3 * (n - 1)^2);
    delta(delta < 1e-9) = 1e-9;
    term1 = nu .* (1 - delta);
    term2 = (1 - delta.^nu);
    term3 = (4 ./ (sin(theta./2).^2 + 1e-14)) .* (delta .* (1 - delta.^nu) - nu .* (1 - delta));
    numerator = term1 - term2 + term3;
    denominator = (1 - delta).^2 .* delta.^nu;
    P_raw = numerator ./ denominator;
end

% Haltrin 经验拟合核心公式
function P_raw = calc_Empirical_Core(theta, b_e, albedo_e)
    q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    k1 = 1.188 - 0.688*albedo_e; k2 = 0.1 * (3.07 - 1.90*albedo_e);
    k3 = 0.01 * (4.58 - 3.02*albedo_e); k4 = 0.001 * (3.24 - 2.25*albedo_e);
    k5 = 0.0001 * (0.84 - 0.61*albedo_e);
    t_deg = theta * 180 / pi;
    term = 1 + (-1)^1*k1*t_deg.^0.5 + (-1)^2*k2*t_deg.^1.0 + (-1)^3*k3*t_deg.^1.5 + ...
           (-1)^4*k4*t_deg.^2.0 + (-1)^5*k5*t_deg.^2.5;
    P_raw = exp(q_e * term);
end

% 弱偏折批量游走求交
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, ...
    Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, ...
    Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen) 
    
    hit_flag = false; 
    N_screens = length(Screen_Chain);
    step_actual = dist_limit;
    hit_rx_plane = false;
    
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
        idx_start = floor(z_start / delta_z_screen) + 1;
        idx_end = floor(z_end / delta_z_screen);
        crossed_screens = idx_start : idx_end;
    elseif dir_z < -1e-6
        idx_start = ceil(z_start / delta_z_screen) - 1;
        idx_end = ceil(z_end / delta_z_screen);
        crossed_screens = idx_start : -1 : idx_end;
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
    lambda_nm = lambda * 1e9; k_wave = 2 * pi / lambda; 
    a1 = 1.779e-4; a2 = -1.05e-6; a3 = 1.6e-8; a4 = -2.02e-6; a5 = 1.155e-2; a6 = -4.23e-3;
    A = a2 * S + 2 * a3 * T * S + 2 * a4 * T + a6 / lambda_nm; 
    B = a1 + a2 * T + a3 * T^2 + a5 / lambda_nm; 
    T_k = T + 273.15; s_frac = S * 1e-3; 
    a11 = 5.328 - 9.76e-2 * S + 4.04e-4 * S^2; a12 = -6.913e-3 + 7.351e-4 * S - 3.15e-6 * S^2;
    a13 = 9.6e-6 - 1.927e-6 * S + 8.23e-9 * S^2; a14 = 2.5e-9 + 1.666e-9 * S - 7.125e-12 * S^2;
    cp = 1000 * (a11 + a12 * T + a13 * T^2 + a14 * T^3); 
    rho_T = 9.9992293295e2 + 2.0341179217e-2 * T - 6.1624591598e-3 * T^2 + 2.2614664708e-5 * T^3 - 4.6570659168e-8 * T^4;
    rho_S = s_frac * (8.0200240891e2 - 2.0005183488 * T + 1.6771024982e-2 * T^2 - 3.0600536746e-5 * T^3 - 1.6132224742e-5 * T * S);
    rho = rho_T + rho_S;
    mu_0 = (0.15700386464 * (T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5;
    a21 = 1.5409136040 + 1.9981117208e-2 * T - 9.5203865864e-5 * T^2;
    a22 = 7.9739318223 - 7.561456881e-2 * T + 4.7237011074e-4 * T^2;
    mu = mu_0 * (1 + a21 * s_frac + a22 * s_frac^2);
    T_b = 1.00024 * T; S_b = S / 1.00472;
    term1 = log10(240 + 0.0002 * S_b); term2 = 0.434 * (2.3 - (343.5 + 0.037 * S_b) / (T_b + 273.15));
    term3 = (1 - (T_b + 273.15) / (647.3 + 0.03 * S_b))^(1/3);
    sigma_T = 10^(term1 - 3 + term2 * term3); 
    Pr = mu * cp / sigma_T; Sc = mu^2 / (5.954e-15 * T_k * rho); 
    c_T = 0.072^(4/3) * Pr^(-1); c_S = 0.072^(4/3) * Sc^(-1); c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    alpha_c = 2.6e-4; beta_c = 7.6e-4; R_rho = alpha_c * abs(H_ratio) / beta_c;
    if R_rho >= 1; d_r = R_rho + sqrt(R_rho) * sqrt(R_rho - 1); elseif R_rho >= 0.5; d_r = 1.85 * R_rho - 0.85; else; d_r = 0.15 * R_rho; end
    chi_S = chi_T * d_r / (H_ratio^2); chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    coeff_Hill = 0.72 / (4 * pi); 
    Phi_Hill = @(K, chi_M, c_M) coeff_Hill * chi_M * epsilon^(-1/3) .* (K.^2 + 1e-25).^(-11/6) .* exp(-176.90 * (K * eta).^2 .* c_M^(0.96)) .* (1 + 21.61 * (K * eta).^(0.61) .* c_M^(0.02) - 18.18 * (K * eta).^(0.55) .* c_M^(0.04));
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + B^2 * Phi_Hill(K, chi_S, c_S) + 2 * A * B * Phi_Hill(K, chi_TS, c_TS));
end

function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2 * pi / D; dx = D / N; kx = (-N/2 : N/2-1) * dk;
    [KX, KY] = meshgrid(kx, kx); K_grid = sqrt(KX.^2 + KY.^2); K_grid(N/2+1, N/2+1) = 1e-10; 
    Phi_n_val = Phi_n_func(K_grid); Phi_n_val(N/2+1, N/2+1) = 0; 
    F_phi = 2 * pi * k_wave^2 * delta_z * Phi_n_val; noise = (randn(N) + 1i * randn(N));
    C_nm = noise .* sqrt(F_phi) * dk; phase_high = real(ifftshift(ifft2(ifftshift(C_nm)))) * N^2;
    phase_low = zeros(N, N); [xx, yy] = meshgrid((-N/2 : N/2-1) * dx); n_sub = 3; 
    for p = 1:n_sub
        dk_p = dk / (3^p); 
        for m = -1:1
            for n = -1:1
                if (m == 0 && n == 0); continue; end
                kx_p = m * dk_p; ky_p = n * dk_p; k_p = sqrt(kx_p^2 + ky_p^2);
                Phi_n_p = Phi_n_func(k_p); F_phi_p = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_p;
                amp = sqrt(F_phi_p) * dk_p; r_c = (randn(1) + 1i * randn(1));
                phase_low = phase_low + real(r_c * amp * exp(1i * (kx_p * xx + ky_p * yy)));
            end
        end
    end
    phase_screen = phase_high + phase_low;
end

function new_dir = rotate_direction(dir, theta_s, psi_s)
    mu_x = dir(1); mu_y = dir(2); mu_z = dir(3); denom = sqrt(1 - mu_z^2);
    if denom < 1e-10
        if mu_z > 0; new_dir = [sin(theta_s) * cos(psi_s), sin(theta_s) * sin(psi_s), cos(theta_s)]; else; new_dir = [sin(theta_s) * cos(psi_s), sin(theta_s) * sin(psi_s), -cos(theta_s)]; end
    else
        sin_theta = sin(theta_s); cos_theta = cos(theta_s); cos_psi = cos(psi_s); sin_psi = sin(psi_s);
        new_dir_x = sin_theta / denom * (mu_x * mu_z * cos_psi - mu_y * sin_psi) + mu_x * cos_theta;
        new_dir_y = sin_theta / denom * (mu_y * mu_z * cos_psi + mu_x * sin_psi) + mu_y * cos_theta;
        new_dir_z = -sin_theta * cos_psi * denom + mu_z * cos_theta;
        new_dir = [new_dir_x, new_dir_y, new_dir_z];
    end
    new_dir = new_dir / norm(new_dir);
end

function W_L = calc_beam_spot_size(w0, lambda, L, Cn2)
    if w0 <= 1e-6; W_L = 0; return; end 
    k = 2 * pi / lambda; D = 2 * w0; z_R = (pi * w0^2) / lambda;
    W_diff = w0 * sqrt(1 + (L / z_R)^2);
    if Cn2 > 1e-15
        rho_0 = (0.545 * k^2 * Cn2 * L)^(-3/5); W_turb_LT = 2 * L / (k * rho_0);
        if rho_0 < D
            correction_factor = max(1 - 0.37 * (rho_0 / D)^(1/3), 0); W_turb_ST = W_turb_LT * correction_factor;
        else
            W_turb_ST = W_turb_LT; 
        end
        W_L = sqrt(W_diff^2 + W_turb_ST^2);
    else
        W_L = W_diff; 
    end
end

function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L)
    kappa_eval = 100; Phi_val = Phi_n_func(kappa_eval);
    Cn2_eq = Phi_val / (0.033 * kappa_eval^(-11/3)); rho_guess = (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5);
    calc_Dw = @(rho) 8 * pi^2 * k_wave^2 * L * integral2(@(K, xi) K .* Phi_n_func(K) .* (1 - besselj(0, K .* rho .* xi)), 1e-1, 1e4, 0, 1, 'Method', 'iterated', 'RelTol', 1e-3);
    try; opts = optimset('Display', 'off'); rho0_exact = fzero(@(rho) calc_Dw(rho) - 2, rho_guess, opts); catch; rho0_exact = rho_guess; end
end