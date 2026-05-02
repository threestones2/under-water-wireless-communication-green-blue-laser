%% 大气紫外光散射通信仿真: MC-MPS 框架重构版
% 逻辑: 基于 MC_MPS_6 框架，同步 MCI_PIS 大气散射参数，停用相位屏与直射分支
% [修改说明]:
% 1. 物理环境: 统一采用大气光速 (2.997e8 m/s) 与紫外波段消光参数。
% 2. 几何构型: 严格匹配 MCI_PIS 的 NLOS 45°/45° 收发指向。
% 3. 核心算法: 保留 MC_MPS_6 的 PIS 散射循环，停用 ray_march_generic 中的相位屏梯度计算。
% 4. 散射模型: 启用 Ding 混合散射相函数 (Rayleigh + Mie)。

clc; clear; close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- 模式选择 ---
param.phase_func = 'Mix';  
param.n_max = 3;           % 最大散射阶数 (同步 MCI_PIS)

% --- 物理常量 (完全对齐大气环境) ---
param.c_air = 2.997046e8;  % 大气光速 (m/s)
lambda_nm = 266;           % 紫外波段 (典型值)
lambda = lambda_nm * 1e-9;
k_wave = 2*pi/lambda;

% --- 介质参数 (同步自 MCI_PIS.m 大气参数) ---
param.coef_a = 0.802e-3;   % 吸收系数 ka (1/m)
param.coef_kr = 0.266e-3;  % Rayleigh 散射系数 kr
param.coef_km = 0.284e-3;  % Mie 散射系数 km
param.coef_b = param.coef_kr + param.coef_km; 
param.coef_c = param.coef_a + param.coef_b; 
param.albedo = param.coef_b / param.coef_c;

% --- Ding 混合模型参数 ---
param.gamma = 0.017;       % Rayleigh 参数
param.g_mie = 0.72;        % Mie 不对称因子
param.f_mie = 0.5;         % Mie 参数

% --- 几何布局 (同步自 MCI_PIS.m) ---
Link_Dist = 100;           % 节点间距 r = 100m
theta_T = deg2rad(45);     % 发射器天顶角
phi_T = 2*pi-deg2rad(90);  % 发射器方位角
theta_R = deg2rad(45);     % 接收器天顶角
phi_R = deg2rad(90);       % 接收器方位角

% 计算主轴方向矢量
mu_T = [cos(phi_T)*sin(theta_T); sin(phi_T)*sin(theta_T); cos(theta_T)];
mu_R = [cos(phi_R)*sin(theta_R); sin(phi_R)*sin(theta_R); cos(theta_R)];

Tx_Pos = [0; Link_Dist; 0]; % 发射器位置
Rx_Pos = [0; 0; 0];         % 接收器位于原点
Rx_Normal = mu_R;           % 接收主轴朝向

Rx_Area = 1.77e-4;          % 接收面积
Rx_Aperture = sqrt(Rx_Area/pi)*2;
Rx_FOV = deg2rad(30);       % 接收器视场角 (FOV)
param.beta_T = deg2rad(17); % 发射束发散角

% --- 时间轴设置 (同步自 MCI_PIS.m) ---
dt = 1e-8;                  % 时间分辨率 10ns
t_min = 0;
t_max = 2e-6;                       
param.T_bins = t_min : dt : t_max;
N_bins = length(param.T_bins);
h_time = zeros(1, N_bins); 

% --- 仿真控制 ---
N_packets = 1e6;            % 仿真光子数

%% ================= 2. 预计算: 相位屏链 (已停用) =================
% 按照要求停用相位屏功能，保持结构但参数设为空
Screen_Chain = []; 
x_axis = 0; dx = 0; N_grid = 1;
delta_z_screen = 1; 

%% ================= 3. 主仿真循环 (MC_MPS_6 框架) =================
fprintf('开始仿真 (MC_MPS_6 框架 + 大气 PIS 参数)...\n');
tic;

for p = 1:N_packets
    % --- 初始化光子 (在发射器坐标系下) ---
    % 采样初始发射方向 (锥形束内均匀分布)
    theta_init = (param.beta_T / 2) * rand;
    phi_init = 2 * pi * rand;
    dir_init = rotate_direction(mu_T, theta_init, phi_init);
    
    pos_init = Tx_Pos;
    weight_init = 1.0;
    
    % [分支 A]: 停用直射分支 (NLOS 环境下直射概率极低，且不符合 PIS 重点)
    % 注释掉原 MC_MPS_6 的 Ballistic 模块
    
    % [分支 B]: 经历散射后打到接收平面 (Scattering)
    pos = pos_init; 
    dir = dir_init;
    weight = weight_init;
    current_dist_traveled = 0;
    
    % --- 高阶散射循环 (核心框架保留) ---
    for order = 1 : param.n_max
        % 1. 采样传输步长
        d_step = -log(1 - rand) / param.coef_c;
        
        % 2. 推进光子位置 (框架调用 ray_march_generic，内部已屏蔽屏干扰)
        [pos_new, dir_new, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
            Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, ...
            Screen_Chain, k_wave, x_axis, dx, N_grid, ...
            Tx_Pos, mu_T, delta_z_screen);
        
        current_dist_traveled = current_dist_traveled + step_len;
        pos = pos_new;
        dir = dir_new;
        
        % 3. PIS 估计 (强制连接接收机)
        vec_to_rx = Rx_Pos - pos;
        dist_to_rx = norm(vec_to_rx);
        dir_to_rx = vec_to_rx / dist_to_rx;
        
        % 检查接收视场 (FOV) - 逻辑对齐 MCI_PIS
        cos_phi_r = dot(mu_R, pos) / dist_to_rx; 
        if cos_phi_r >= cos(Rx_FOV/2)
            % 计算散射相函数 p(theta)
            theta_s = acos(dot(-dir, pos) / dist_to_rx); 
            cos_theta_s = cos(theta_s);
            
            % Ding 混合模型 PDF
            p_ray = (3 * (1 + 3*param.gamma + (1-param.gamma)*cos_theta_s^2)) / (16 * pi * (1 + 2*param.gamma));
            p_mie = ((1 - param.g_mie^2)/2 * (1/(1 + param.g_mie^2 - 2*param.g_mie*cos_theta_s)^1.5 + ...
                    param.f_mie * 0.5 * (3*cos_theta_s^2 - 1) / (1 + param.g_mie^2)^1.5)) / (2*pi);
            p_phase = (param.coef_kr/param.coef_b) * p_ray + (param.coef_km/param.coef_b) * p_mie;
            
            % PIS 权重计算
            omega = (Rx_Area / dist_to_rx^2) * cos_phi_r;
            prob_survival = exp(-param.coef_c * dist_to_rx);
            
            energy_pis = weight * (param.coef_b/param.coef_c) * p_phase * omega * prob_survival;
            
            % 计入时间槽
            t_arrival = (current_dist_traveled + dist_to_rx) / param.c_air;
            bin_idx = floor(t_arrival / dt) + 1;
            if bin_idx >= 1 && bin_idx <= N_bins
                h_time(bin_idx) = h_time(bin_idx) + energy_pis;
            end
        end
        
        % 4. 准备下一阶散射 (重要性采样)
        weight = weight * (param.coef_b/param.coef_c);
        if weight < 1e-10, break; end
        
        % 采样下一个物理散射方向
        dir = rotate_direction(dir, pi*rand, 2*pi*rand);
    end
    
    if mod(p, N_packets/10) == 0, fprintf('   进度: %.0f%%\n', p/N_packets*100); end
end
toc;

%% ================= 4. 结果处理 =================
h_time_norm = h_time / (N_packets * dt); 
P_rx_total = sum(h_time) / N_packets;
Path_Loss_dB = 10 * log10(1 / P_rx_total);

fprintf('\n=== 仿真结果 (MC_MPS 框架) ===\n');
fprintf('总接收概率 P: %.4e\n', P_rx_total);
fprintf('总路径损耗 PL: %.2f dB\n', Path_Loss_dB);

figure;
plot(param.T_bins * 1e6, h_time_norm, 'b-', 'LineWidth', 1.5);
xlabel('Time (\mus)'); ylabel('h(t) (W/s)');
title(['Atmospheric UV Channel - PL: ' num2str(Path_Loss_dB, '%.2f') ' dB']);
grid on;

%% ================= 5. 辅助函数 =================

% --- [逻辑精简] 移除了相位屏梯度计算的 Ray Marching ---
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, ...
    Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, ...
    Screen_Chain, k_wave, x_axis, dx, N_grid, ...
    Tx_Pos, Link_Dir, delta_z_screen)

    hit_flag = false;
    % 由于停用了相位屏，光子在散射点之间沿直线飞行
    total_len = dist_limit;
    pos = pos + dir * dist_limit;
    
    % 如果是直射检测分支 (enable_hit_check=true)，则判断是否穿过接收机孔径
    if enable_hit_check
        % 此处逻辑在散射主循环中默认不开启，保持结构完整
    end
end

function new_dir = rotate_direction(mu, theta, phi)
    % 采用与 MCI_PIS.m 一致的方向旋转逻辑
    mu_x = mu(1); mu_y = mu(2); mu_z = mu(3);
    if abs(mu_z) < 1 - 1e-10
        denom = sqrt(1 - mu_z^2);
        new_dir = [cos(theta)*mu_x + sin(theta)*(mu_x*mu_z*cos(phi) - mu_y*sin(phi))/denom;
                  cos(theta)*mu_y + sin(theta)*(mu_y*mu_z*cos(phi) + mu_x*sin(phi))/denom;
                  cos(theta)*mu_z - sin(theta)*cos(phi)*denom];
    else
        new_dir = [sin(theta)*cos(phi); sin(theta)*sin(phi); sign(mu_z)*cos(theta)];
    end
    new_dir = new_dir / norm(new_dir);
end