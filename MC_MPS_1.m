%% 水下光通信混合仿真: MC-PIS (重要性采样) + 动态相位屏 + Eq.10衍射衰减
% 
% 核心逻辑:
% 1. 粒子性: 光子在水中随机游走 (Monte Carlo)。
% 2. 波动性: 每次移动步长内，根据步长插入若干相位屏，利用梯度偏折光子方向。
% 3. 能量修正: 使用 Model 论文 Eq.10 计算湍流导致的衍射能量损耗。
% 4. 加速: 使用部分重要性采样 (PIS) 在每次散射点估计接收功率。

clc; clear; close all;

%% ================= 1. 参数初始化 =================
% --- 光源与几何 ---
lambda_nm = 532; 
lambda = lambda_nm * 1e-9;
k_wave = 2*pi/lambda;
w0 = 0.01;                  % 束腰半径
div_angle = 1e-3;           % 初始发散角 (rad)
Tx_Pos = [0, 0, 0];         % 发射机位置
Rx_Pos = [0, 0, 200];        % 接收机位置 
Rx_Aperture = 0.1;          % 接收孔径直径 (m)
Rx_FOV = 30 * pi/180;       % 视场角
Rx_Area = pi*(Rx_Aperture/2)^2;

% --- 海水光学参数 ---
% water_type = 'coastal';
coef_a = 0.179;             % 吸收系数 (1/m)
coef_b = 0.219;             % 散射系数 (1/m)
coef_c = coef_a + coef_b;
g_HG = 0.924;               % HG散射各向异性因子
albedo = coef_b / coef_c;   % 反照率

% --- 湍流参数 (OTOPS) ---
T_avg = 20; S_avg = 35; 
epsilon = 1e-5; chi_T = 1e-7; eta = 1e-3; H_ratio = -20;
D_screen = 2.0;             % 相位屏物理尺寸 (要足够大以覆盖散射范围)
N_grid = 256;               % 网格数
delta_z_screen = 2.0;       % 定义相位屏的物理间隔 (即每飞2m遇到一张屏)

% --- 仿真控制 ---
N_packets = 1e4;            % 光子数 (根据性能调整)
max_interactions = 5;      % 最大散射次数

%% ================= 2. 预计算: Eq.10 衰减与相位屏池 =================
fprintf('1. 预生成相位屏池 (模拟无限介质)...\n');
% 为了避免在MC循环中实时做FFT (太慢)，我们预生成一堆典型的相位屏
% 每次光子经过时，随机抽取一张作为当前的湍流切片
Num_Screen_Pool = 20; 
Screen_Pool = struct('grad_x', [], 'grad_y', []);
dx = D_screen / N_grid;
x_axis = (-N_grid/2 : N_grid/2-1) * dx;

for i = 1:Num_Screen_Pool
    % 调用之前的次谐波补偿代码 (假设已保存为函数)
    phi = generate_phase_screen_SH_Corrected(D_screen, N_grid, lambda, delta_z_screen, ...
                                             T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
    [gx, gy] = gradient(phi, dx);
    Screen_Pool(i).grad_x = gx; % 存储梯度
    Screen_Pool(i).grad_y = gy;
end

%% ================= 3. MC-PIS 主循环 =================
fprintf('3. 开始 MC-PIS 仿真...\n');
Rx_Power_Total = 0; % 归一化接收功率

tic;
for p = 1:N_packets
    % --- Step A: 初始化光子 ---
    % 初始位置 (高斯光束分布)
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos = [r0*cos(phi0), r0*sin(phi0), 0];
    
    % 初始方向 (高斯发散角)
    theta0 = div_angle * sqrt(-log(rand())); phi_dir0 = 2*pi*rand();
    dir = [sin(theta0)*cos(phi_dir0), sin(theta0)*sin(phi_dir0), cos(theta0)];
    
    weight = 1.0;
    
    % --- Step B: 散射追踪 ---
    for s = 1:max_interactions
        % 1. 采样散射步长 (仅由散射系数 b 决定，吸收和湍流放权重里)
        % 注意: 标准MC用 c 采样，这里用 PIS 方法常用 b 采样防止光子过快灭亡，
        % 或者用 c 采样并在权重中扣除吸收。这里采用经典方式:
        step_dist = -log(rand()) / coef_c; 
        
        % 2. 判断是否会越过接收平面 (简单的 Z 轴判断)
        dist_to_plane = (Rx_Pos(3) - pos(3)) / dir(3);
        if dist_to_plane > 0 && dist_to_plane < step_dist
            % 光子可能会直接打到接收平面 (Ballistic / Low order scattering)
            % 但这部分由 "强制探测" (PIS) 覆盖吗？
            % 传统的 PIS 专注于散射后的贡献。直射光(Ballistic)需要单独算，
            % 或者让光子真的撞上去。这里为了简单，让光子真的撞上去只做判定，不重复计分。
            % (此处逻辑略复杂，为简化，我们假设 PIS 贡献为主，但保留直接撞击)
            step_dist = dist_to_plane;
            is_hit_plane = true;
        else
            is_hit_plane = false;
        end
        
        % 3. --- [核心] 动态相位屏传输过程 ---
        % 在 step_dist 这么长的路径里，光子会穿过多少个相位屏？
        num_screens_in_step = floor(step_dist / delta_z_screen);
        
        % 逐步推进光子，每一步都受到湍流偏折和衰减
        segment_left = step_dist;
        
        for k = 1 : max(1, num_screens_in_step)
            % 本次子步长
            if num_screens_in_step == 0
                d_sub = step_dist; % 不足一个屏间距
            else
                d_sub = min(segment_left, delta_z_screen);
            end
            
            % 3.1 移动光子
            pos = pos + dir * d_sub;
            segment_left = segment_left - d_sub;
            
            % 3.2 衰减 (吸收 + Eq.10 湍流衍射)
            % 散射衰减已经在步长采样中体现了，所以这里只乘吸收和额外湍流损耗
            % 修正: 如果用 c 采样步长，权重乘 albedo。这里再乘 Eq.10。
            loss_factor = (coef_b / coef_c) * exp(-alpha_turb * d_sub);
            % 注意: 如果是第一步(未散射)，不乘 albedo，只乘 exp(-a*z)。
            % 为严谨，这里简化处理：假设步长采样包含了 scattering loss
            % 权重更新：w = w * (b/c) * exp(-alpha * d) 是不对的。
            % 标准做法：步长 -ln(r)/c。到达散射点时，权重 w = w * albedo。
            % 额外衰减：w = w * exp(-alpha * d)
            weight = weight * exp(-alpha_turb * d_sub); 
            
            if weight < 1e-9, break; end
            
            % 3.3 湍流偏折 (几何光学近似)
            % 从池中随机取一个屏
            pool_idx = randi(Num_Screen_Pool);
            
            % 寻找光子在屏上的位置索引
            % (假设相位屏中心对应 (0,0)，且随光子横向位置平铺或周期边界)
            % 这里使用周期边界条件防止越界
            idx_x = mod(round((pos(1) - x_axis(1))/dx), N_grid) + 1;
            idx_y = mod(round((pos(2) - x_axis(1))/dx), N_grid) + 1;
            
            % 获取梯度 (单位: rad/m)
            g_x = Screen_Pool(pool_idx).grad_x(idx_y, idx_x);
            g_y = Screen_Pool(pool_idx).grad_y(idx_y, idx_x);
            
            % 偏折角度变化 d_theta = grad_phi / k
            delta_theta_x = g_x / k_wave;
            delta_theta_y = g_y / k_wave;
            
            % 更新方向向量 (小量近似)
            dir(1) = dir(1) + delta_theta_x;
            dir(2) = dir(2) + delta_theta_y;
            dir = dir / norm(dir); % 归一化
        end
        
        if weight < 1e-9, break; end
        
        % --- Step C: 处理到达事件 ---
        if is_hit_plane
            % 检查是否在孔径和视场内
            r_spot = sqrt((pos(1)-Rx_Pos(1))^2 + (pos(2)-Rx_Pos(2))^2);
            angle_inc = acos(dot(dir, [0,0,1])); % 假设接收机朝向 Z
            
            if r_spot <= Rx_Aperture/2 && angle_inc <= Rx_FOV
                % 这是一个直接打到的光子 (Ballistic or purely geometric wander)
                Rx_Power_Total = Rx_Power_Total + weight;
            end
            break; % 撞到接收平面，停止追踪该光子
        end
        
        % --- Step D: 重要性采样 (PIS) - 强制探测 ---
        % 在发生散射的位置 pos，计算“如果光子散射向接收机”的概率
        
        % 1. 接收机方向向量
        vec_to_rx = Rx_Pos - pos;
        dist_to_rx = norm(vec_to_rx);
        dir_to_rx = vec_to_rx / dist_to_rx;
        
        % 2. HG 相函数概率 p(cos_theta)
        cos_theta_s = dot(dir, dir_to_rx);
        p_phase = (1 - g_HG^2) / (4*pi * (1 + g_HG^2 - 2*g_HG*cos_theta_s)^1.5);
        
        % 3. 接收机立体角 (近似: A / r^2 * cos_incidence)
        cos_inc_rx = dot(-dir_to_rx, [0,0,-1]); % 接收机法线向下? 假设接收机朝向 -Z (面对光源)
        % 修正: Rx_Pos=[0,0,20], Tx=[0,0,0]. 光往+Z传。Rx朝向-Z。
        % dir_to_rx 指向 +Z。 Rx法线应为 [0,0,-1]。
        % dot([0,0,1], [0,0,-1]) = -1。 入射角 0。
        % 简单起见，假设接收机平面法线是 -vec_to_rx (对准光子) 或者固定 Z 轴
        if dist_to_rx > 0
            omega = Rx_Area / (dist_to_rx^2) * abs(dir_to_rx(3)); 
        else
            omega = 0;
        end
        
        % 4. 强制路径上的衰减 (Beer-Lambert + Eq.10 湍流)
        % 这是一条虚拟直线路径
        prob_survival = exp(-(coef_c + alpha_turb) * dist_to_rx);
        
        % 5. 视场角判断
        angle_in_rx = acos(abs(dir_to_rx(3))); 
        if angle_in_rx <= Rx_FOV
            % 累加 PIS 贡献
            % P_scat = W_current * albedo * p_phase * omega * attenuation
            P_pis = weight * albedo * p_phase * omega * prob_survival;
            Rx_Power_Total = Rx_Power_Total + P_pis;
        end
        
        % --- Step E: 真实散射 (更新状态用于下一步) ---
        % 能量被散射衰减
        weight = weight * albedo;
        
        % 采样新的飞行方向 (HG 散射)
        dir = scatter_Henyey_Greenstein(dir, g_HG);
    end
    
    if mod(p, N_packets/10) == 0, fprintf('   进度: %.0f%%\n', p/N_packets*100); end
end
toc;

Rx_Power_Normalized = Rx_Power_Total / N_packets;
fprintf('\n=== 仿真结果 ===\n');
fprintf('接收功率 (归一化): %.6e\n', Rx_Power_Normalized);
fprintf('链路损耗: %.2f dB\n', 10*log10(Rx_Power_Normalized));


%% ================= 辅助函数 =================

% HG 散射方向更新 (标准算法)
function dir_new = scatter_Henyey_Greenstein(dir, g)
    % 1. 采样散射角 theta
    rand_val = rand();
    if g == 0
        cos_theta = 2*rand_val - 1;
    else
        temp = (1 - g^2) / (1 - g + 2*g*rand_val);
        cos_theta = (1 + g^2 - temp^2) / (2*g);
    end
    sin_theta = sqrt(1 - cos_theta^2);
    theta=acos(cos_theta);
    phi = 2*pi*rand();
    
    % 2. 坐标旋转 (Rodrigues 旋转或局部坐标系)
    dir_new=rotate_direction(dir,theta,phi);
end
%%
% 旋转方向向量（根据散射角和方位角）
function new_dir = rotate_direction(dir, theta_s, psi_s)
    % 当前方向余弦
    mu_x = dir(1);
    mu_y = dir(2);
    mu_z = dir(3);
    
    % 处理数值稳定性问题
    denom = sqrt(1 - mu_z^2);
    if denom < 1e-10
        % 当原方向接近z轴时的特殊处理
        if mu_z > 0
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), cos(theta_s)];
        else
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), -cos(theta_s)];
        end
    else
        % 常规计算
        sin_theta = sin(theta_s);
        cos_theta = cos(theta_s);
        cos_psi = cos(psi_s);
        sin_psi = sin(psi_s);
        
        new_dir_x = sin_theta/denom * (mu_x*mu_z*cos_psi - mu_y*sin_psi) + mu_x*cos_theta;
        new_dir_y = sin_theta/denom * (mu_y*mu_z*cos_psi + mu_x*sin_psi) + mu_y*cos_theta;
        new_dir_z = -sin_theta*cos_psi/denom + mu_z*cos_theta;
        
        new_dir = [new_dir_x, new_dir_y, new_dir_z];
    end
    
    % 归一化
    new_dir = new_dir / norm(new_dir);
end
%%
% OTOPS 谱计算
function [Phi_n, Phi_Hill_Handle] = calc_otops_spectrum(k_vals, eps, chi, eta, T, S, lam_nm, Hr)
    % 计算 OTOPS 功率谱密度
    % 系数
    a1=1.779e-4; a2=-1.05e-6; a3=1.6e-8; a4=-2.02e-6; a5=1.155e-2; a6=-4.23e-3;
    A = a2*S + 2*a3*T*S + 2*a4*T + a6/lam_nm; 
    B = a1 + a2*T + a3*T^2 + a5/lam_nm;
    
    rho=1025; mu=1e-3; cp=4182; sigma_T = 0.6; 
    Pr=mu*cp/sigma_T; 
    Sc=mu^2/(5.954e-15*(T+273.15)*rho);
    
    alpha_c=2.6e-4; beta_c=7.6e-4;
    R_rho=alpha_c*abs(Hr)/beta_c;
    if R_rho>=1, dr=R_rho+sqrt(R_rho*(R_rho-1));
    elseif R_rho>=0.5,dr=1.85*R_rho-0.85; 
    else,dr=0.15*R_rho;
    end
    
    chi_S = chi*dr/Hr^2; 
    chi_TS = chi*(1+dr)/(2*Hr);
    
    % Hill Spectrum Function
    f_Hill = @(K, chi_x, c_x) 0.033 * chi_x * eps^(-1/3) .* (K.^2 + 1e-20).^(-11/6) .* ...
        exp(-176.9*(K*eta).^2 .* c_x^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_x^(0.02) - 18.18*(K*eta).^(0.55).*c_x^(0.04));
    
    c_T=0.072^(4/3)/Pr; 
    c_S=0.072^(4/3)/Sc; 
    c_TS=0.072^(4/3)*(Pr+Sc)/(2*Pr*Sc);
    
    Phi_n = A^2*f_Hill(k_vals, chi, c_T) + B^2*f_Hill(k_vals, chi_S, c_S) + ...
            2*A*B*f_Hill(k_vals, chi_TS, c_TS);
    
    Phi_Hill_Handle = f_Hill; % 返回句柄以备用
end
