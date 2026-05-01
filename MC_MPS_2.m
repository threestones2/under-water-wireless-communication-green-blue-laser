%% 水下光通信混合仿真: MC-PIS (修正版) + 动态相位屏
% 
% 核心逻辑修正:
% 1. 粒子性: 光子在水中随机游走，步长由总衰减系数 c 采样。
% 2. 波动性: 相位屏仅提供几何偏折(梯度)，不产生能量衰减(避免与 alpha_turb 重复)。
% 3. 能量修正: 删除 Eq.10 衰减项，能量损耗仅由海水吸收(a)和几何发散决定。
% 4. 加速: 使用 PIS 加速收敛。

clc; clear; close all;

%% ================= 1. 参数初始化 =================
% --- 光源与几何 ---
lambda_nm = 532; %算折射率谱应该用nm，不然量级太大了
lambda = lambda_nm * 1e-9; %m为单位
k_wave = 2*pi/lambda;
w0 = 0.01;                  % 束腰半径
div_angle = 1e-3;           % 初始发散角 (rad)
Tx_Pos = [0, 0, 20];         % 发射机位置(距离改回20m以便测试，200m可能全黑)
Rx_Pos = [0, 0, 0];        % 接收机位置 （注意这里接收机放在原点）
Rx_Aperture = 0.1;          % 接收孔径直径 (m)
Rx_FOV = 30 * pi/180;       % 视场角
Rx_Area = pi*(Rx_Aperture/2)^2;

% --- 海水光学参数 ---
coef_a = 0.179;             % 吸收系数 (1/m)
coef_b = 0.219;             % 散射系数 (1/m)
coef_c = coef_a + coef_b;   % 总衰减系数
g_HG = 0.924;               % HG散射各向异性因子
albedo = coef_b / coef_c;   % 反照率 (每次碰撞后的生存概率)

% --- 湍流参数 (OTOPS) ---
T_avg = 20; S_avg = 35; 
epsilon = 1e-5; chi_T = 1e-7; eta = 1e-3; H_ratio = -20;
D_screen = 2.0;             % 相位屏尺寸
N_grid = 256;               % 网格数
delta_z_screen = 2.0;       % 相位屏物理间隔 (每飞2m遇到一次湍流切片)

% --- 仿真控制 ---
N_packets = 1e4;            % 光子数
max_interactions = 4;      % 最大散射次数

%% ================= 2. 预计算: 相位屏池 =================
% 而 MC 模拟的是能量传输，相位屏的偏折已经包含了这一物理过程。

fprintf('1. 预生成相位屏池 (模拟无限介质)...\n');
Num_Screen_Pool = 20; 
Screen_Pool = struct('grad_x', [], 'grad_y', []);
dx = D_screen / N_grid;
x_axis = (-N_grid/2 : N_grid/2-1) * dx;

for i = 1:Num_Screen_Pool
    % 生成相位屏 (rad)
    phi = generate_phase_screen_SH_Corrected(D_screen, N_grid, lambda, delta_z_screen, ...
                                             T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
    % 计算梯度 (rad/m)
    [gx, gy] = gradient(phi, dx);
    Screen_Pool(i).grad_x = gx;
    Screen_Pool(i).grad_y = gy;
end

%% ================= 3. MC-PIS 主循环 =================
fprintf('2. 开始 MC-PIS 仿真...\n');
Rx_Power_Total = 0; % 归一化接收功率

tic;
for p = 1:N_packets
    % --- Step A: 初始化光子 ---
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos = [r0*cos(phi0), r0*sin(phi0), 0];
    
    theta0 = div_angle * sqrt(-log(rand())); phi_dir0 = 2*pi*rand();
    dir = [sin(theta0)*cos(phi_dir0), sin(theta0)*sin(phi_dir0), cos(theta0)];
    dir = dir / norm(dir);
    
    weight = 1.0;
    
    % --- Step B: 散射追踪 ---
    for s = 1:max_interactions
        % 1. 采样散射步长 (Sampling Step Size)
        % 物理含义: 光子在 coef_c 的介质中，按概率能飞多远才撞到粒子。
        step_dist = -log(rand()) / coef_c; 
        
        % 2. 接收平面判定
        dist_to_plane = (Rx_Pos(3) - pos(3)) / dir(3);
        if dist_to_plane > 0 && dist_to_plane < step_dist
            % 如果步长内穿过了接收平面，截断步长
            step_dist = dist_to_plane;
            is_hit_plane = true;
        else
            is_hit_plane = false;
        end
        
        % 3. --- [核心] 动态相位屏传输过程 (Ray Marching) ---
        % 计算这段路程中会经过几个相位屏
        num_screens_in_step = floor(step_dist / delta_z_screen);
        segment_left = step_dist;
        
        % 在到达下一个散射点(或接收面)之前，分段飞行并与相位屏交互
        for k = 1 : max(1, num_screens_in_step)
            % 确定本次子步长
            if num_screens_in_step == 0
                d_sub = step_dist; 
            else
                d_sub = min(segment_left, delta_z_screen);
            end
            
            % [3.1] 移动光子
            pos = pos + dir * d_sub;
            segment_left = segment_left - d_sub;
            
            % [3.2] 权重更新 (修正点!)
            % 在飞行过程中，权重不变！
            % 因为 step_dist 是根据 c 采样的，意味着光子"存活"并飞过这段距离的概率已经被采样过程包含了。
            % 吸收损耗 (a) 和散射损耗 (b) 将在碰撞发生时统一处理。
            
            if weight < 1e-9, break; end
            
            % [3.3] 湍流偏折 (仅当确实穿过了屏时)
            % 如果这只是最后一段不足 delta_z 的路程(且没穿屏)，就不偏折
            if num_screens_in_step > 0 && segment_left >= 0
                % 随机抽取相位屏
                pool_idx = randi(Num_Screen_Pool);
                
                % 周期性边界条件查表
                idx_x = mod(round((pos(1) - x_axis(1))/dx), N_grid) + 1;
                idx_y = mod(round((pos(2) - x_axis(1))/dx), N_grid) + 1;
                
                % 读取梯度
                g_x = Screen_Pool(pool_idx).grad_x(idx_y, idx_x);
                g_y = Screen_Pool(pool_idx).grad_y(idx_y, idx_x);
                
                % 几何光学偏折: theta = gradient / k
                delta_theta_x = g_x / k_wave;
                delta_theta_y = g_y / k_wave;
                
                % 更新方向 (小角度近似叠加)
                dir(1) = dir(1) + delta_theta_x;
                dir(2) = dir(2) + delta_theta_y;
                dir = dir / norm(dir);
            end
        end
        
        if weight < 1e-9, break; end
        
        % --- Step C: 处理到达接收平面 ---
        if is_hit_plane
            r_spot = sqrt((pos(1)-Rx_Pos(1))^2 + (pos(2)-Rx_Pos(2))^2);
            % 计算入射角
            angle_inc = acos(dot(dir, [0,0,1])); 
            
            if r_spot <= Rx_Aperture/2 && angle_inc <= Rx_FOV
                % 这是一个"真实"到达的光子 (可能是直射，也可能是经过几次散射后撞上来的)
                % 注意: 如果使用 PIS，通常只统计散射分量 PIS，直射分量由这里统计。
                % 为避免重复，PIS 通常计算的是"散射后"的贡献。
                % 这里为了简单，我们认为 Ballistic 光子在这里被捕获，Scattered 光子由 PIS 捕获。
                if s == 1 % 只有未发生散射的光子才在这里计入 (Ballistic)
                     Rx_Power_Total = Rx_Power_Total + weight;
                end
            end
            break; % 光子终止
        end
        
        % --- Step D: 重要性采样 (PIS) - 估计散射贡献 ---
        % 光子在这里(pos)发生散射，我们强制计算一部分能量散射到接收机的概率
        
        % 1. 几何向量
        vec_to_rx = Rx_Pos - pos;
        dist_to_rx = norm(vec_to_rx);
        dir_to_rx = vec_to_rx / dist_to_rx;
        
        % 2. 判断是否被接收 (视场角 + 遮挡)
        % 接收机朝向 -Z (假设)
        cos_inc_rx = dot(-dir_to_rx, [0,0,-1]); 
        angle_in_rx = acos(abs(dir_to_rx(3))); % 简化：与Z轴夹角
        
        if angle_in_rx <= Rx_FOV
            % 3. HG 相函数概率 p(cos_theta)
            cos_theta_s = dot(dir, dir_to_rx);
            p_phase = (1 - g_HG^2) / (4*pi * (1 + g_HG^2 - 2*g_HG*cos_theta_s)^1.5);
            
            % 4. 立体角近似
            omega = Rx_Area / (dist_to_rx^2) * abs(dir_to_rx(3));
            
            % 5. 虚拟路径衰减 (修正点!)
            % 这里只乘 exp(-c * dist)，不再乘 alpha_turb。
            % 含义: 光子在这条虚拟直线上不被任何东西阻挡的概率。
            prob_survival = exp(-coef_c * dist_to_rx);
            
            % 6. 累加贡献
            % 公式: P_rx = Weight * Albedo * Phase_Prob * Solid_Angle * Transmittance
            P_pis = weight * albedo * p_phase * omega * prob_survival;
            Rx_Power_Total = Rx_Power_Total + P_pis;
        end
        
        % --- Step E: 真实散射 (Update Status) ---
        % 1. 能量衰减: 只有在这里，权重才因为吸收而减少
        weight = weight * albedo; 
        
        % 2. 改变方向: 为下一次追踪做准备
        dir = scatter_Henyey_Greenstein(dir, g_HG);
    end
    
    if mod(p, N_packets/10) == 0, fprintf('   进度: %.0f%%\n', p/N_packets*100); end
end
toc;

Rx_Power_Normalized = Rx_Power_Total / N_packets;
fprintf('\n=== 仿真结果 ===\n');
fprintf('接收功率 (归一化): %.6e\n', Rx_Power_Normalized);
fprintf('链路损耗: %.2f dB\n', 10*log10(Rx_Power_Normalized));

%% ================= 辅助函数 (保持不变) =================
function dir_new = scatter_Henyey_Greenstein(dir, g)
    rand_val = rand();
    if g == 0
        cos_theta = 2*rand_val - 1;
    else
        temp = (1 - g^2) / (1 - g + 2*g*rand_val);
        cos_theta = (1 + g^2 - temp^2) / (2*g);
    end
    sin_theta = sqrt(1 - cos_theta^2);
    theta = acos(cos_theta);
    phi = 2*pi*rand();
    dir_new = rotate_direction(dir, theta, phi);
end

function new_dir = rotate_direction(dir, theta_s, psi_s)
    mu_x = dir(1); mu_y = dir(2); mu_z = dir(3);
    denom = sqrt(1 - mu_z^2);
    if denom < 1e-10
        if mu_z > 0
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), cos(theta_s)];
        else
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), -cos(theta_s)];
        end
    else
        sin_theta = sin(theta_s); cos_theta = cos(theta_s);
        cos_psi = cos(psi_s); sin_psi = sin(psi_s);
        new_dir_x = sin_theta/denom * (mu_x*mu_z*cos_psi - mu_y*sin_psi) + mu_x*cos_theta;
        new_dir_y = sin_theta/denom * (mu_y*mu_z*cos_psi + mu_x*sin_psi) + mu_y*cos_theta;
        new_dir_z = -sin_theta*cos_psi/denom + mu_z*cos_theta;
        new_dir = [new_dir_x, new_dir_y, new_dir_z];
    end
    new_dir = new_dir / norm(new_dir);
end

% 请确保 generate_phase_screen_SH_Corrected 函数在路径中