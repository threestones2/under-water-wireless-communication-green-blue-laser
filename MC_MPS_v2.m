%% 水下/大气光通信混合仿真: MC-MPS (Turbulence Comparison)
%  Comparison: With Turbulence vs. Without Turbulence
%  Unified Petzold Clear Ocean parameters + Adaptive HG-PIS + Wavefront Decoupling
clc; 
clear; 
close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- 模式选择 ---
% 相函数模型选择经验 Petzold 模型
param.phase_func = 'Empirical';  
% 蒙特卡洛追踪的最大多重散射阶数
param.n_max = 300;          

% --- 光源与波长参数 ---
% 载波波长 (nm)，选择 514 nm 处于海水透射窗口（蓝绿波段）
lambda_nm = 514; 
lambda = lambda_nm * 1e-9;
% 波数
k_wave = 2 * pi / lambda;
% % 发射高斯光束的束腰半径 (m)
% w0 = 0.01;                       
% % 光束发散角 (rad)
% div_angle = 0.01 * pi / 180;           
% % 瑞利长度 (Rayleigh range)，表征光束衍射展宽的特征距离
% z_R = (pi * w0^2) / lambda;     


w0 = 0.0;                       
div_angle = 0.0;           
z_R = realmin; % 避免除零

% --- 3D 空间链路布局 ---
% 发射端三维坐标 (m)
Tx_Pos = [0, 0, 0];         
% 接收端三维坐标 (m)，此处设定 Y 轴正向距离为 10m
Rx_Pos = [0, 10.93, 0];       
Link_Vec = Rx_Pos - Tx_Pos;
% 物理链路绝对距离 (m)
Link_Dist = norm(Link_Vec);
% 链路主光轴方向向量 (归一化)
Link_Dir = Link_Vec / Link_Dist;

% --- 光轴对准偏差 (Misalignment / Pointing Errors) ---
% 发射端指向偏差极角 (Tx Pointing Error)
theta_tx_error = 0;      
% 发射端偏离的方位角 (rad)
phi_tx_error = 0;           
% 接收端法线偏离对准方向的极角 (rad)
theta_rx_error = 0;      
% 接收端偏离的方位角 (rad)，pi 表示接收面朝向发射端
phi_rx_error = pi;          

% 利用旋转矩阵生成偏转后的真实发射主轴方向与接收端法向量
mu_T = rotate_direction(Link_Dir, theta_tx_error, phi_tx_error);            
Rx_Normal = rotate_direction(-Link_Dir, theta_rx_error, phi_rx_error);

% 生成发射端局部发射平面的正交基，用于高斯光束截面的二维空间采样
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
% 接收孔径直径 (m)
Rx_Aperture = 0.4;          
% 孔径物理半径 (m)
Rx_Radius = Rx_Aperture / 2;
% 接收视场角 (Field of View, rad)
Rx_FOV = 20 * pi / 180;        
% 接收孔径面积 (m^2)
Rx_Area = pi * (Rx_Radius)^2;

% --- 介质参数 (Petzold 水体光学特性) ---
% 光在水中的群速度 (m/s)
param.c_water = 2.237e8;   

% Petzold清澈海水
% param.coef_c = 0.1514;                % 衰减系数 c (1/m)
% param.coef_a = 0.114;                 % 吸收系数 a (1/m)
% param.coef_b = 0.0374; % 散射系数 b (0.0374 1/m)

% %Petzold近岸海水
% param.coef_c = 0.398;
% param.coef_a = 0.179;
% param.coef_b=0.219;

% Petzold港口海水
param.coef_c = 1.829+0.366;
param.coef_a = 0.366;
param.coef_b=1.829;

% 单次散射反照率 (Albedo)，表征散射事件在总衰减中的占比
param.albedo = param.coef_b / param.coef_c;

% --- 相函数参数 (基于 Petzold 数据的经验公式) ---
param.c_e = param.coef_c; 
param.a_e = param.coef_a; 
b_e = param.coef_b; 
albedo_e = param.albedo;

% 经验相函数拟合系数计算
param.q_e = 2.598 + 17.748 * sqrt(b_e) - 16.722 * b_e + 5.932 * b_e * sqrt(b_e);
param.k1 = 1.188 - 0.688 * albedo_e; 
param.k2 = 0.1 * (3.07 - 1.90 * albedo_e);
param.k3 = 0.01 * (4.58 - 3.02 * albedo_e); 
param.k4 = 0.001 * (3.24 - 2.25 * albedo_e);
param.k5 = 0.0001 * (0.84 - 0.61 * albedo_e);

% 预计算经验相函数的归一化系数，确保在 4*pi 立体角内积分为 1
% 该经验系数可以用param.coef_b代替，误差不大，10%以内，大多情况下都是5%以内
th_test = linspace(0, pi, 2000); 
val_test = zeros(size(th_test));
for i = 1:length(th_test)
    t_deg = th_test(i) * 180 / pi; 
    if t_deg < 1e-6
        t_deg = 1e-6; 
    end
    term = 1 + (-1)^1 * param.k1 * t_deg^0.5 + ...
           (-1)^2 * param.k2 * t_deg^1.0 + ...
           (-1)^3 * param.k3 * t_deg^1.5 + ...
           (-1)^4 * param.k4 * t_deg^2.0 + ...
           (-1)^5 * param.k5 * t_deg^2.5;
    val_test(i) = exp(param.q_e * term);
end
param.b_emp_norm = 2 * pi * trapz(th_test, val_test .* sin(th_test));

% --- 提取自适应最优 HG 提议分布参数 (HG-PIS 方差缩减技术) ---
% 提取前向散射峰值
P_max = pdf_Empirical(1.0, param); 
C_val = 4 * pi * P_max;
% 计算等效的 Henyey-Greenstein 不对称因子 g
g_prop = 1 + (1 - sqrt(8 * C_val + 1)) / (2 * C_val);
param.g_prop = g_prop;
fprintf('自适应最优 HG-PIS 不对称因子 g_prop = %.4f\n', param.g_prop);

% --- 自动计算水体固有相函数的等效发散角 ---
theta_peak = 1e-6; 
p_peak = pdf_Empirical(cos(theta_peak), param);
target_p = p_peak * exp(-2); 
obj_func = @(th) pdf_Empirical(cos(th), param) - target_p;
try
    opts = optimset('Display', 'off');
    theta_water_calc = fzero(obj_func, [1e-6, 0.1], opts); 
catch
    % 若求根失败，使用二阶矩估计等效发散角
    th_test_w = linspace(1e-6, pi/4, 2000); 
    p_test = zeros(size(th_test_w));
    for idx = 1:length(th_test_w)
        p_test(idx) = pdf_Empirical(cos(th_test_w(idx)), param); 
    end
    mean_sq_theta = 2 * pi * trapz(th_test_w, (th_test_w.^2) .* p_test .* sin(th_test_w));
    theta_water_calc = sqrt(mean_sq_theta);
end
param.theta_water = theta_water_calc;
fprintf('等效水体发散角 theta_water = %.4f rad\n', param.theta_water);

% --- 弱湍流参数 (OTOPS 海洋湍流谱参数) ---
T_avg = 20;                     % 平均温度 (摄氏度)
S_avg = 35;                     % 平均盐度 (ppt)
H_ratio = -20;                  % 温度和盐度对湍流贡献的比值 (Wadell 参数)
epsilon = 1e-9;                 % 湍流动能耗散率 (m^2/s^3)
chi_T = 1e-7;                   % 温度方差耗散率 (K^2/s)
eta = 1e-3;                     % 柯尔莫哥洛夫微尺度 (m)
N_screens = 20;                 % 多重相位屏数量
D_screen = 10;                  % 相位屏物理尺寸 (m)
N_grid = 2^8;                   % 相位屏网格分辨率
delta_z_screen = Link_Dist / N_screens; % 相位屏间距

% --- 动态时间轴设置 (增量到达时间 \Delta t) ---
% 计算视距传播的极短到达时间
t_LOS = Link_Dist / param.c_water; 
% 设定散射造成的最大时间延迟拓展 (设为 10 ns)
delta_t_max = 10e-9;                       
dt = 1e-10;                     % 时间分辨率 (0.1 ns)
t_min = t_LOS;                  % 时间窗口起始点设定为 LOS 到达时间
t_max = t_LOS + delta_t_max;    % 时间窗口终止点
% 绝对时间 bins，用于物理运算
param.T_bins = t_min : dt : t_max;
% 相对时间延迟 bins，用于图表展示 (\Delta t)
param.Delta_T_bins = param.T_bins - t_LOS;  
N_bins = length(param.T_bins);

% --- 仿真控制 ---
N_packets = 1e5;                % 发射光子包总数
scenarios = {'With Turbulence', 'Without Turbulence'};
results = struct(); 

%% ================= 2. 对比仿真循环 =================
dx = D_screen / N_grid;
x_axis = (-N_grid/2 : N_grid/2-1) * dx;

% 获取 OTOPS 海洋折射率湍流功率谱句柄
[Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
fprintf('预计算 OTOPS 湍流相干参数...\n');
[rho0_Link, Cn2_eq] = calc_turb_coherence_params(Phi_func, k_wave, Link_Dist);
fprintf('  基于 OTOPS 提取的等效 Cn2 = %.2e m^(-2/3)\n', Cn2_eq);
fprintf('  链路终点精确相干长度 rho_0 = %.4f m\n\n', rho0_Link);

% 预计算链路平面的正交基
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

for s_idx = 1:2
    scenario_name = scenarios{s_idx};
    fprintf('=== Running Scenario: %s ===\n', scenario_name);

    % 初始化相位屏链条结构体
    Screen_Chain = repmat(struct('Center', [0,0,0], 'Normal', [0,0,1], ...
                             'u_vec', [], 'v_vec', [], 'grad_x', [], 'grad_y', []), 1, N_screens);
    for i = 1:N_screens
        Screen_Chain(i).Center = Tx_Pos + i * delta_z_screen * Link_Dir;
        Screen_Chain(i).Normal = Link_Dir;
        Screen_Chain(i).u_vec = u_vec; 
        Screen_Chain(i).v_vec = v_vec;
        
        if strcmp(scenario_name, 'With Turbulence')
            % 频谱反演法生成相位屏
            phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
            % 计算相位梯度用于光线偏转
            [gx, gy] = gradient(phi, dx);
            Screen_Chain(i).grad_x = gx; 
            Screen_Chain(i).grad_y = gy;
        else
            Screen_Chain(i).grad_x = zeros(N_grid, N_grid); 
            Screen_Chain(i).grad_y = zeros(N_grid, N_grid);
        end
    end

    % 初始化时域统计直方图 (加入微小底噪避免对数域错误)
    h_time = 1e-11 * ones(1, N_bins); 
    tic;

    for p = 1:N_packets
        % --- 2.1 微观粒子状态初始化 ---
        % 发射平面内的高斯空间抽样
        r0 = w0 * sqrt(-0.5 * log(rand())); 
        phi0 = 2 * pi * rand();
        pos_local = r0 * cos(phi0) * u_vec_Tx + r0 * sin(phi0) * v_vec_Tx;
        pos_init = Tx_Pos + pos_local;

        % 高斯角分布抽样
        U_init = div_angle * sqrt(-0.5 * log(rand()));
        psi_ini = 2 * pi * rand;
        dir_init = rotate_direction(mu_T, U_init, psi_ini);

        weight_init = 1.0;
        Huge_Aperture = 1e5; % 用于全截面接收检测
        
        % --- 2.2 [Ballistic Branch] 弹道光子宏观质心解耦评估 ---
        % 弹道分量不受水体散射影响，直接经历湍流相位屏和水体吸收
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
                % 基于 Marcum Q 函数的严格二维高斯接收能量积分
                r_eff = Rx_Radius * sqrt(cos_rx_tilt);
                a_param = 2 * r_wander / W_spot;
                b_param = 2 * r_eff / W_spot;
                received_fraction = 1 - marcumq(a_param, b_param);

                % 引入水体总衰减 c
                attenuation = exp(-param.coef_c * path_len_ballistic);
                final_weight = weight_init * attenuation * received_fraction;

                t_arrival = path_len_ballistic / param.c_water;
                bin_idx = floor((t_arrival - t_min) / dt) + 1;
                if bin_idx >= 1 && bin_idx <= N_bins
                    h_time(bin_idx) = h_time(bin_idx) + final_weight;
                end
            end
        end

        % --- 2.3 [Scattering Branch] 散射光微观演化 ---
        pos = pos_init; 
        dir = dir_init; 
        weight = weight_init;
        current_dist_traveled = 0;

        % 第 0 阶初始步长抽样 (Beer-Lambert 定律)
        d_step = -log(rand()) / param.coef_c;
        [pos_new, dir_new, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
            Rx_Pos, Huge_Aperture, Rx_FOV, Rx_Normal, false, ...
            Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
            
        current_dist_traveled = current_dist_traveled + step_len;
        % 若未散射前已经越过接收面，则跳过后续演化
        if dot(pos_new - Tx_Pos, Link_Dir) >= Link_Dist
            continue; 
        end
        pos = pos_new; 
        dir = dir_new;

        % 1~n阶多重散射循环
        for order = 1 : param.n_max
            vec_to_rx = Rx_Pos - pos; 
            dist_to_rx = norm(vec_to_rx); 
            dir_to_rx = vec_to_rx / dist_to_rx;

            % ---------------- 接收估计 (Local Estimation) ----------------
            is_in_FOV_theta = acos(dot(-dir_to_rx, Rx_Normal));
            if is_in_FOV_theta <= Rx_FOV/2
                cos_theta_s = dot(dir, dir_to_rx);
                % 取经验相函数对应角度的概率
                p_phase = pdf_Empirical(cos_theta_s, param);
                % 接收立体角权重
                omega = Rx_Area / (dist_to_rx^2) * abs(dot(dir_to_rx, Rx_Normal));
                % 剩余传播距离的生存概率
                prob_survival = exp(-param.coef_c * dist_to_rx);
                % 基础接收期望能量
                base_energy_pis = weight * param.albedo * min(1, p_phase * omega) * prob_survival;

                if base_energy_pis > 1e-13
                     if strcmp(scenario_name, 'With Turbulence')
                        % 评估该强制接收分支的路径是否穿过相位屏引起额外偏移
                        [pos_v_end, ~, v_hit_flag, v_len] = ray_march_generic(pos, dir_to_rx, dist_to_rx + 1e-1, ...
                            Rx_Pos, Huge_Aperture, pi, Rx_Normal, true, ...
                            Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);

                        if v_hit_flag
                            r_wander_pis = norm(pos_v_end - Rx_Pos);
                            % 基于球面波近似衰减修正相干长度
                            rho_0 = rho0_Link * (Link_Dist / v_len)^(3/5);
                            theta_turb_LT = 2 / (k_wave * rho_0); 
                            correction_factor = max(0, 1 - 0.37 * (rho_0 / (2 * Rx_Radius))^(1/3));
                            theta_turb_ST = theta_turb_LT * correction_factor;
                            theta_water = param.theta_water; 

                            % 综合水体扩散与湍流展宽
                            W_total_turb = v_len * sqrt(theta_water^2 + theta_turb_ST^2);
                            cos_rx_tilt_pis = abs(dot(dir_to_rx, Rx_Normal));
                            r_eff_pis = Rx_Radius * sqrt(cos_rx_tilt_pis);

                            a_turb = 2 * r_wander_pis / W_total_turb;
                            b_turb = 2 * r_eff_pis / W_total_turb;
                            fraction_turb = 1 - marcumq(a_turb, b_turb);

                            W_clean = v_len * theta_water;
                            b_clean = 2 * r_eff_pis / W_clean;
                            fraction_clean = 1 - marcumq(0, b_clean);

                            % 以比值形式施加物理惩罚
                            if fraction_clean > 1e-15
                                final_energy_pis = base_energy_pis * (fraction_turb / fraction_clean);
                            else
                                final_energy_pis = 0;
                            end

                            t_pis_arrival = (current_dist_traveled + v_len) / param.c_water;
                            bin_idx = floor((t_pis_arrival - t_min) / dt) + 1;
                            if bin_idx >= 1 && bin_idx <= N_bins
                                h_time(bin_idx) = h_time(bin_idx) + final_energy_pis;
                            end
                        end
                    else
                        t_pis_arrival = (current_dist_traveled + dist_to_rx) / param.c_water;
                        bin_idx = floor((t_pis_arrival - t_min) / dt) + 1;
                        if bin_idx >= 1 && bin_idx <= N_bins
                            h_time(bin_idx) = h_time(bin_idx) + base_energy_pis;
                        end
                    end
                end
            end

            % 若达到最大散射次数限制则终止
            if order == param.n_max
                break; 
            end
            
            % 更新光子代数权重 (吸收效应)
            weight = weight * param.albedo;
            if weight < 1e-9
                break; 
            end

            % ---------------- 采用自适应 HG-PIS 进行真实的物理偏转抽样 ----------------
            % 使用解析求逆的 Henyey-Greenstein 函数生成提议偏转角
            xi_rand = rand();
            term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi_rand);
            cos_t = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1);
            theta_i = acos(cos_t); 
            phi_scat = 2 * pi * rand();

            % 计算 HG 提议概率与真实经验概率的似然比
            q_HG = (1 - g_prop^2) / (4 * pi * (1 + g_prop^2 - 2 * g_prop * cos_t)^1.5);
            p_val = pdf_Empirical(cos_t, param);

            % 代数权重补偿 (Importance Weighting)
            weight_factor = min(2e2, p_val / q_HG);
            %weight_factor =  p_val / q_HG;
            weight = weight * weight_factor;

            dir = rotate_direction(dir, theta_i, phi_scat);

            % 抽取下一段传播步长并推进物理位置
            d_step = -log(rand()) / param.coef_c;
            [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
                Rx_Pos, Huge_Aperture, Rx_FOV, Rx_Normal, false, ...
                Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
                
            current_dist_traveled = current_dist_traveled + step_len;
            
            % 如果高阶散射后光子已越过接收平面对应的 Z 深度，终止演化
            if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist
                break; 
            end
        end
    end
    t_run = toc;
    fprintf('   Completed in %.2f s\n', t_run);

    % --- 2.4: 存储并计算本场景通道特性 ---
    h_time_norm = h_time / N_packets;      
    P_rx_total = sum(h_time_norm);         

    results(s_idx).h_time_energy = h_time_norm;          
    % 计算物理单位为 s^-1 的 CIR
    %results(s_idx).h_time_cir = h_time_norm / dt; 
    results(s_idx).h_time_cir = h_time_norm ; 
    results(s_idx).P_rx = P_rx_total;
    results(s_idx).Path_Loss_dB = -10 * log10(P_rx_total);

    if P_rx_total > 0
        % 均方根延迟拓展度量 (RMS Delay Spread)
        tau_mean = sum(param.T_bins .* h_time_norm) / P_rx_total;
        tau_rms = sqrt( sum( ((param.T_bins - tau_mean).^2) .* h_time_norm ) / P_rx_total );
        results(s_idx).tau_rms = tau_rms;
    else
        results(s_idx).tau_rms = 0;
    end
end

%% ================= 3. 对比绘图与分析 =================
fprintf('\n=== 仿真对比结果 (Tx-Rx: %.1fm) ===\n', Link_Dist);
fprintf('%-20s | Path Loss (dB) | RMS Delay (ns)\n', 'Scenario');
fprintf('---------------------------------------------------\n');
for i = 1:2
    fprintf('%-20s | %14.4f | %14.4f\n', scenarios{i}, results(i).Path_Loss_dB, results(i).tau_rms * 1e9);
end

% --- Plot 1: CIR 对比 (对数域) ---
figure('Name', 'CIR Comparison', 'Color', 'w');
% 横轴使用时间增量 \Delta t (ns)
plot(param.Delta_T_bins * 1e9, 10 * log10(results(1).h_time_cir), 'r-', 'LineWidth', 1, 'DisplayName', 'With Turbulence');
hold on;
plot(param.Delta_T_bins * 1e9, 10 * log10(results(2).h_time_cir), 'b--', 'LineWidth', 1, 'DisplayName', 'Without Turbulence');
grid on;
legend('Location', 'best');
xlabel('Time Increment \Delta t (ns)');
ylabel('Channel Impulse Response (dB s^{-1})');
title(['CIR Comparison: Turbulence vs. No Turbulence (d=' num2str(Link_Dist) 'm)']);

% --- Plot 2: 归一化连续能量分布对比 (线性域) ---
figure('Name', 'Linear Power Comparison', 'Color', 'w');
area(param.Delta_T_bins * 1e9, results(2).h_time_cir, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'No Turb');
hold on;
plot(param.Delta_T_bins * 1e9, results(1).h_time_cir, 'r-', 'LineWidth', 1.5, 'DisplayName', 'With Turb');
grid on;
legend;
xlabel('Time Increment \Delta t (ns)'); 
ylabel('Channel Impulse Response (s^{-1})');
title('Linear Impulse Response Comparison');


%% ================= 4. 辅助函数 =================

% 光束步进追踪函数（包含与相位屏的求交测试）
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, ...
    Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, ...
    Screen_Chain, k_wave, x_axis, dx, N_grid, ...
    Tx_Pos, Link_Dir, delta_z_screen) 

    hit_flag = false;
    total_len = 0;
    remaining_dist = dist_limit;
    N_screens = length(Screen_Chain);
    
    while remaining_dist > 1e-6
        min_dist = remaining_dist;
        event_type = 'none';
        target_idx = -1;
        
        % 计算在传播主轴上的投影深度，定位当前层
        dist_projected = dot(pos - Tx_Pos, Link_Dir);
        current_layer_idx = floor(dist_projected / delta_z_screen);
        target_screen_idx = -1;
        cos_theta = dot(dir, Link_Dir);
        
        % 确定前进方向并匹配下一个可能相交的相位屏
        if cos_theta > 1e-6       
            target_screen_idx = current_layer_idx + 1;
        elseif cos_theta < -1e-6  
            target_screen_idx = current_layer_idx;
        end
        
        if target_screen_idx >= 1 && target_screen_idx <= N_screens
            scr = Screen_Chain(target_screen_idx);
            denom = dot(dir, scr.Normal);
            if abs(denom) > 1e-6
                % 射线与屏幕平面求交
                t = dot(scr.Center - pos, scr.Normal) / denom;
                if t > 1e-6 && t < min_dist
                    min_dist = t;
                    event_type = 'screen';
                    target_idx = target_screen_idx;
                end
            end
        end
        
        % 检查是否抵达接收平面
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx) > 1e-6
                t_rx = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
                if t_rx > 1e-6 && t_rx <= min_dist
                    min_dist = t_rx;
                    event_type = 'rx';
                end
            end
        end
        
        % 演化位置与距离消耗
        pos = pos + dir * min_dist;
        remaining_dist = remaining_dist - min_dist;
        total_len = total_len + min_dist; 
        
        if strcmp(event_type, 'rx')
            % 抵达接收端，判决几何与角度条件
            dist_to_center = norm(pos - Rx_Pos);
            angle_inc = acos(dot(-dir, Rx_Normal));
            if dist_to_center <= Rx_Aperture/2 && angle_inc <= Rx_FOV/2
                hit_flag = true;
                return;
            end
            
        elseif strcmp(event_type, 'screen')
            % 穿透湍流相位屏，引入波前折射偏转
            scr = Screen_Chain(target_idx);
            vec_on_plane = pos - scr.Center;
            loc_u = dot(vec_on_plane, scr.u_vec);
            loc_v = dot(vec_on_plane, scr.v_vec);
            
            % 将物理坐标映射至网格索引
            idx_x = mod(round((loc_u - x_axis(1)) / dx), N_grid) + 1;
            idx_y = mod(round((loc_v - x_axis(1)) / dx), N_grid) + 1;
            g_u = scr.grad_x(idx_y, idx_x); 
            g_v = scr.grad_y(idx_y, idx_x);
            
            % 波矢偏转量
            delta_vec = (g_u * scr.u_vec + g_v * scr.v_vec) / k_wave;
            dir = dir + delta_vec; 
            dir = dir / norm(dir);
        end
    end
end

% 海洋湍流光学谱模型 (OTOPS) 构建
function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    lambda_nm = lambda * 1e9; 
    k_wave = 2 * pi / lambda; 
    
    % 折射率对温度和盐度的偏导数经验系数
    a1 = 1.779e-4;  
    a2 = -1.05e-6;  
    a3 = 1.6e-8; 
    a4 = -2.02e-6;  
    a5 = 1.155e-2;  
    a6 = -4.23e-3;
    
    A = a2 * S + 2 * a3 * T * S + 2 * a4 * T + a6 / lambda_nm; 
    B = a1 + a2 * T + a3 * T^2 + a5 / lambda_nm; 
    
    T_k = T + 273.15; 
    s_frac = S * 1e-3; 
    
    a11 = 5.328 - 9.76e-2 * S + 4.04e-4 * S^2;
    a12 = -6.913e-3 + 7.351e-4 * S - 3.15e-6 * S^2;
    a13 = 9.6e-6 - 1.927e-6 * S + 8.23e-9 * S^2;
    a14 = 2.5e-9 + 1.666e-9 * S - 7.125e-12 * S^2;
    cp = 1000 * (a11 + a12 * T + a13 * T^2 + a14 * T^3); 
    
    rho_T = 9.9992293295e2 + 2.0341179217e-2 * T - ...
            6.1624591598e-3 * T^2 + 2.2614664708e-5 * T^3 - ...
            4.6570659168e-8 * T^4;
    rho_S = s_frac * (8.0200240891e2 - 2.0005183488 * T + ...
            1.6771024982e-2 * T^2 - 3.0600536746e-5 * T^3 - ...
            1.6132224742e-5 * T * S);
    rho = rho_T + rho_S;
    
    mu_0 = (0.15700386464 * (T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5;
    a21 = 1.5409136040 + 1.9981117208e-2 * T - 9.5203865864e-5 * T^2;
    a22 = 7.9739318223 - 7.561456881e-2 * T + 4.7237011074e-4 * T^2;
    mu = mu_0 * (1 + a21 * s_frac + a22 * s_frac^2);
    
    T_b = 1.00024 * T; 
    S_b = S / 1.00472;
    term1 = log10(240 + 0.0002 * S_b);
    term2 = 0.434 * (2.3 - (343.5 + 0.037 * S_b) / (T_b + 273.15));
    term3 = (1 - (T_b + 273.15) / (647.3 + 0.03 * S_b))^(1/3);
    sigma_T = 10^(term1 - 3 + term2 * term3); 
    
    Pr = mu * cp / sigma_T;  
    Sc = mu^2 / (5.954e-15 * T_k * rho); 
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1);
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    alpha_c = 2.6e-4; 
    beta_c = 7.6e-4;  
    R_rho = alpha_c * abs(H_ratio) / beta_c;
    
    if R_rho >= 1
        d_r = R_rho + sqrt(R_rho) * sqrt(R_rho - 1);
    elseif R_rho >= 0.5
        d_r = 1.85 * R_rho - 0.85;
    else
        d_r = 0.15 * R_rho; 
    end
    
    chi_S = chi_T * d_r / (H_ratio^2);
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    
    coeff_Hill = 0.72 / (4 * pi); 
    % 整合湍流功率谱的通用函数型
    Phi_Hill = @(K, chi_M, c_M) coeff_Hill * chi_M * epsilon^(-1/3) .* ...
        (K.^2 + 1e-25).^(-11/6) .* exp(-176.90 * (K * eta).^2 .* c_M^(0.96)) .* ...
        (1 + 21.61 * (K * eta).^(0.61) .* c_M^(0.02) - 18.18 * (K * eta).^(0.55) .* c_M^(0.04));
    
    % 构建多组分组合后的完整 OTOPS 海水谱
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + ...
                       B^2 * Phi_Hill(K, chi_S, c_S) + ...
                       2 * A * B * Phi_Hill(K, chi_TS, c_TS));
end

% 频谱反演法（次谐波补偿）生成单层随机相位屏
function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2 * pi / D; 
    dx = D / N;
    kx = (-N/2 : N/2-1) * dk;
    [KX, KY] = meshgrid(kx, kx);
    K_grid = sqrt(KX.^2 + KY.^2);
    K_grid(N/2+1, N/2+1) = 1e-10; 
    
    Phi_n_val = Phi_n_func(K_grid);
    Phi_n_val(N/2+1, N/2+1) = 0; 
    F_phi = 2 * pi * k_wave^2 * delta_z * Phi_n_val;
    noise = (randn(N) + 1i * randn(N));
    C_nm = noise .* sqrt(F_phi) * dk;
    % 高频相位屏部分 (FFT)
    phase_high = real(ifftshift(ifft2(ifftshift(C_nm)))) * N^2;
    
    % 添加次谐波以补偿低频能量缺失
    phase_low = zeros(N, N);
    [xx, yy] = meshgrid((-N/2 : N/2-1) * dx);
    n_sub = 3; 
    for p = 1:n_sub
        dk_p = dk / (3^p); 
        for m = -1:1
            for n = -1:1
                if (m == 0 && n == 0)
                    continue; 
                end
                kx_p = m * dk_p; 
                ky_p = n * dk_p;
                k_p = sqrt(kx_p^2 + ky_p^2);
                Phi_n_p = Phi_n_func(k_p);
                F_phi_p = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_p;
                amp = sqrt(F_phi_p) * dk_p;
                r_c = (randn(1) + 1i * randn(1));
                phase_low = phase_low + real(r_c * amp * exp(1i * (kx_p * xx + ky_p * yy)));
            end
        end
    end
    phase_screen = phase_high + phase_low;
end

% 依据给定的极角和方位角旋转基准方向向量
function new_dir = rotate_direction(dir, theta_s, psi_s)
    mu_x = dir(1); 
    mu_y = dir(2); 
    mu_z = dir(3);
    denom = sqrt(1 - mu_z^2);
    
    if denom < 1e-10
        if mu_z > 0
            new_dir = [sin(theta_s) * cos(psi_s), sin(theta_s) * sin(psi_s), cos(theta_s)];
        else
            new_dir = [sin(theta_s) * cos(psi_s), sin(theta_s) * sin(psi_s), -cos(theta_s)]; 
        end
    else
        sin_theta = sin(theta_s); 
        cos_theta = cos(theta_s);
        cos_psi = cos(psi_s); 
        sin_psi = sin(psi_s);
        new_dir_x = sin_theta / denom * (mu_x * mu_z * cos_psi - mu_y * sin_psi) + mu_x * cos_theta;
        new_dir_y = sin_theta / denom * (mu_y * mu_z * cos_psi + mu_x * sin_psi) + mu_y * cos_theta;
        new_dir_z = -sin_theta * cos_psi * denom + mu_z * cos_theta;
        new_dir = [new_dir_x, new_dir_y, new_dir_z];
    end
    new_dir = new_dir / norm(new_dir);
end

% 经验相函数求值（保证物理阈值截断防止越界）
function p = pdf_Empirical(cos_theta, param)
    cos_theta = max(min(cos_theta, 1), -1);
    theta_rad = acos(cos_theta);
    t_deg = theta_rad * 180 / pi;
    if t_deg < 1e-6
        t_deg = 1e-6; 
    end
    term = 1 + (-1)^1 * param.k1 * t_deg^0.5 + ...
           (-1)^2 * param.k2 * t_deg^1.0 + ...
           (-1)^3 * param.k3 * t_deg^1.5 + ...
           (-1)^4 * param.k4 * t_deg^2.0 + ...
           (-1)^5 * param.k5 * t_deg^2.5;
    VSF = exp(param.q_e * term);
    p = VSF / param.b_emp_norm;
end

% 弹道分量波束斑展宽计算（考虑衍射与弱湍流影响）
function W_L = calc_beam_spot_size(w0, lambda, L, Cn2)
    k = 2 * pi / lambda; 
    D = 2 * w0; 
    z_R = (pi * w0^2) / lambda;
    
    % 真空中的衍射展宽
    W_diff = w0 * sqrt(1 + (L / z_R)^2);
    
    if Cn2 > 1e-15
        % 平面波近似下的空间相干长度
        rho_0 = (0.545 * k^2 * Cn2 * L)^(-3/5);
        W_turb_LT = 2 * L / (k * rho_0);
        
        if rho_0 < D
            correction_factor = 1 - 0.37 * (rho_0 / D)^(1/3);
            correction_factor = max(correction_factor, 0); 
            W_turb_ST = W_turb_LT * correction_factor;
        else
            W_turb_ST = W_turb_LT; 
        end
        % 利用方差组合求联合展宽
        W_L = sqrt(W_diff^2 + W_turb_ST^2);
    else
        W_L = W_diff; 
    end
end

% 由 OTOPS 频谱求解系统有效的大气等效折射率结构常数 C_n^2 及相干长度
function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L)
    kappa_eval = 100; 
    Phi_val = Phi_n_func(kappa_eval);
    % 利用标准 Kolmogorov 谱律倒推等效 Cn2
    Cn2_eq = Phi_val / (0.033 * kappa_eval^(-11/3));
    rho_guess = (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5);
    
    % 利用波结构函数 (Wave Structure Function) 求解相干长度
    calc_Dw = @(rho) 8 * pi^2 * k_wave^2 * L * integral2(...
        @(K, xi) K .* Phi_n_func(K) .* (1 - besselj(0, K .* rho .* xi)), ...
        1e-1, 1e4, 0, 1, 'Method', 'iterated', 'RelTol', 1e-3);
    
    try
        opts = optimset('Display', 'off');
        rho0_exact = fzero(@(rho) calc_Dw(rho) - 2, rho_guess, opts);
    catch
        rho0_exact = rho_guess;
    end
end