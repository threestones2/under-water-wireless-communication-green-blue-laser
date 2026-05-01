%% 3. 紫外大气光通信仿真: MC-MPS 框架适配 MCI_PIS_2 参数 (无湍流 NLOS)
clc; clear; close all;

%% ================= 参数初始化 (基于 MCI_PIS_2 设定) =================
param = struct();
param.n_max = 3;           % 紫外通信最大散射阶数

% --- 介质参数 (Atmospheric UV) ---
param.ka = 0.802e-3;       % 吸收系数
param.kr = 0.266e-3;       % Rayleigh
param.km = 0.284e-3;       % Mie
param.ks = param.kr + param.km; 
param.ke = param.ka + param.ks; 
param.albedo = param.ks / param.ke;
param.c_air = 2.997046e8;  % 真空/空气光速

param.gamma = 0.017; param.g = 0.72; param.f = 0.5;

% --- 空间布局与方向 (NLOS Intersection) ---
r_dist = 100;
Tx_Pos = [0, r_dist, 0];   
Rx_Pos = [0, 0, 0];        

% 发射端角度：theta_T = 45度, phi_T = 270度
t_T = 45*pi/180; p_T = 270*pi/180;
mu_T = [cos(p_T)*sin(t_T), sin(p_T)*sin(t_T), cos(t_T)]; 

% 接收端角度：theta_R = 45度, phi_R = 90度
t_R = 45*pi/180; p_R = 90*pi/180;
mu_R = [cos(p_R)*sin(t_R), sin(p_R)*sin(t_R), cos(t_R)];
Rx_Normal = mu_R; % 接收平面的法向量即观测朝向

div_angle = 17*pi/180;     
Rx_FOV = 30*pi/180;        
Rx_Area = 1.77e-4;         
Rx_Aperture = 2 * sqrt(Rx_Area / pi);

% 发射局部正交基
if abs(mu_T(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
u_vec_Tx = cross(up_temp, mu_T); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(mu_T, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

dt = 5e-9; t_min = 0; t_max = 2e-6; param.T_bins = t_min : dt : t_max;
N_bins = length(param.T_bins);
N_packets = 1e6; % 针对弱 NLOS 收敛增加样本

%% ================= 仿真核心 =================
h_time = 1e-12*ones(1, N_bins); 
tic;
for p = 1:N_packets
    % 均匀圆锥采样 (匹配原 MCI_PIS 的 beta_T 采样)
    th_emit = (div_angle / 2) * rand();
    phi_emit = 2 * pi * rand();
    dir_init = rotate_direction(mu_T, th_emit, phi_emit);
    
    pos = Tx_Pos; dir = dir_init; 
    % 使用锥内均匀分布补正初始能量权重
    weight = (sin(th_emit) / (1 - cos(div_angle/2))) * (div_angle/2) * pi; 
    dist_travel = 0;
    
    for order = 1 : param.n_max
        d_step = -log(rand()) / param.ke;
        pos = pos + dir * d_step;
        dist_travel = dist_travel + d_step;
        
        v_rx = Rx_Pos - pos; d_rx = norm(v_rx); d_dir = v_rx / d_rx;
        
        % PIS 接收判定 (检查是否落入接收视场)
        if acos(dot(-d_dir, Rx_Normal)) <= Rx_FOV/2
            cos_tilt = abs(dot(d_dir, Rx_Normal));
            omega = (Rx_Area / d_rx^2) * cos_tilt;
            p_phase = pdf_UV(dot(dir, d_dir), param);
            
            % PIS 能量提取公式（匹配 UV NLOS 物理量纲）
            base_eng = weight * param.albedo * min(1, p_phase * omega) * exp(-param.ke * d_rx);
            if base_eng > 1e-15
                idx = floor(((dist_travel + d_rx)/param.c_air - t_min) / dt) + 1;
                if idx >= 1 && idx <= N_bins, h_time(idx) = h_time(idx) + base_eng; end
            end
        end
        
        weight = weight * param.albedo; if weight < 1e-9, break; end
        
        % 产生下一次散射角度并依据 UV 相函数做重要性权重补正
        th_i = pi * rand; p_scat = 2 * pi * rand;
        weight = weight * (2 * pi^2 * pdf_UV(cos(th_i), param) * sin(th_i));
        dir = rotate_direction(dir, th_i, p_scat);
    end
end
fprintf('Completed in %.2f s\n', toc);

% 绘图
figure('Color', 'w'); 
plot(param.T_bins, h_time/(N_packets*dt), 'k-'); grid on;
xlabel('Time (s)'); ylabel('Impulse Response (W/m^2)'); 
title('UV NLOS CIR (MC-MPS adaptation)');

%% ================= 辅助函数 =================
function new_dir = rotate_direction(dir, theta, psi)
    mz = dir(3); denom = sqrt(1 - mz^2);
    if denom < 1e-10, new_dir = [sin(theta)*cos(psi), sin(theta)*sin(psi), sign(mz)*cos(theta)];
    else, new_dir = [sin(theta)/denom*(dir(1)*mz*cos(psi) - dir(2)*sin(psi))+dir(1)*cos(theta), ...
                     sin(theta)/denom*(dir(2)*mz*cos(psi) + dir(1)*sin(psi))+dir(2)*cos(theta), -sin(theta)*cos(psi)*denom+mz*cos(theta)];
    end
    new_dir = new_dir / norm(new_dir);
end

function p = pdf_UV(cos_theta, param)
    % 混合 Rayleigh 与 Mie 散射的三维空间归一化相函数
    p_ray = 3 * (1 + 3*param.gamma + (1-param.gamma)*cos_theta^2) / (16 * pi * (1 + 2*param.gamma));
    p_mie = (1 - param.g^2)/2 * (1./(1 + param.g^2 - 2*param.g*cos_theta).^1.5 + ...
            param.f * 0.5 * (3*cos_theta^2 - 1) / (1 + param.g^2)^1.5) / (2*pi);
    p = (param.kr/param.ks) * p_ray + (param.km/param.ks) * p_mie;
end