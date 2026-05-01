%% 水下/大气光通信通用混合仿真: MC-PIS (双模式版) + 动态相位屏
% 
% 支持模式 (param.phase_func):
%   1. 'HG'  : Henyey-Greenstein 相函数 (典型水下)
%   2. 'Mix' : Rayleigh + Mie 混合相函数 (典型大气/紫外或特定水体)
%
% 核心逻辑:
% 1. 粒子性: 蒙特卡洛光子追踪，步长由总衰减系数 c 采样。
% 2. 波动性: 动态相位屏提供几何偏折(梯度)，模拟湍流引起的波前畸变。
% 3. PIS加速: 在每个散射点进行部分重要性采样(PIS)估计接收概率。

clc; clear; close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- 模式选择 ---
% 可选: 'HG' (水下默认) 或 'Mix' (大气/混合)
param.phase_func = 'HG'; 

% --- 光源与几何 ---
lambda_nm = 532; 
lambda = lambda_nm * 1e-9;
k_wave = 2*pi/lambda;
w0 = 0.01;                  % 束腰半径 (m)
div_angle = 1e-3;           % 初始发散角 (rad)
Tx_Pos = [0, 0, 0];         % 发射机位置
Rx_Pos = [0, 0, 20];        % 接收机位置 (m)
Rx_Aperture = 0.1;          % 接收孔径直径 (m)
Rx_FOV = 30 * pi/180;       % 视场角 (rad)
Rx_Area = pi*(Rx_Aperture/2)^2;

% --- 介质光学参数 (根据模式自动计算) ---
param.coef_a = 0.179;       % 吸收系数 (1/m)

if strcmp(param.phase_func, 'HG')
    % [HG 模式参数]
    param.coef_b = 0.219;   % 散射系数 (1/m)
    param.g_HG = 0.924;     % HG各向异性因子
    
    fprintf('模式: HG (水下), g=%.3f\n', param.g_HG);
    
elseif strcmp(param.phase_func, 'Mix')
    % [Mix 模式参数 - Rayleigh + Mie]
    param.coef_kr = 0.05;   % Rayleigh 散射系数 (1/m)
    param.coef_km = 0.169;  % Mie 散射系数 (1/m)
    param.coef_b = param.coef_kr + param.coef_km; % 总散射
    
    % Mix 相函数形状参数
    param.gamma = 0.017;    % Rayleigh 参数
    param.g_mie = 0.72;     % Mie g参数
    param.f_mie = 0.5;      % Mie f参数 (Generalized HG)
    
    fprintf('模式: Mix (Rayleigh+Mie)\n');
else
    error('未知的相函数模式');
end

param.coef_c = param.coef_a + param.coef_b; % 总衰减
param.albedo = param.coef_b / param.coef_c; % 反照率

% --- 湍流参数 (OTOPS) ---
T_avg = 20; S_avg = 35; 
epsilon = 1e-5; chi_T = 1e-7; eta = 1e-3; H_ratio = -20;
D_screen = 2.0;             % 相位屏尺寸
N_grid = 256;               % 网格数
delta_z_screen = 2.0;       % 相位屏间距

% --- 仿真控制 ---
N_packets = 1e4;            % 光子数 (建议 1e5 以上以获得平滑结果)
max_interactions = 10;      % 最大散射阶数

%% ================= 2. 预计算: 相位屏池 =================
fprintf('1. 预生成相位屏池 (模拟无限介质)...\n');
Num_Screen_Pool = 10; % 示例用10张
Screen_Pool = struct('grad_x', [], 'grad_y', []);
dx = D_screen / N_grid;
x_axis = (-N_grid/2 : N_grid/2-1) * dx;

% 这里假设您已有 generate_phase_screen_SH_Corrected 函数
% 为防止报错，若无该函数请注释掉下行并使用全零梯度
try
    for i = 1:Num_Screen_Pool
        phi = generate_phase_screen_SH_Corrected(D_screen, N_grid, lambda, delta_z_screen, ...
                                                 T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
        [gx, gy] = gradient(phi, dx);
        Screen_Pool(i).grad_x = gx;
        Screen_Pool(i).grad_y = gy;
    end
catch
    warning('未找到 generate_phase_screen_SH_Corrected，使用零梯度(无湍流)代替。');
    for i = 1:Num_Screen_Pool
        Screen_Pool(i).grad_x = zeros(N_grid);
        Screen_Pool(i).grad_y = zeros(N_grid);
    end
end

%% ================= 3. MC-PIS 主循环 =================
fprintf('2. 开始 MC-PIS 仿真...\n');
Rx_Power_Total = 0; 
tic;

for p = 1:N_packets
    % --- Step A: 初始化光子 ---
    % 空间分布 (高斯)
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos = [r0*cos(phi0), r0*sin(phi0), 0];
    
    % 角度分布 (高斯发散)
    theta0 = div_angle * sqrt(-log(rand())); phi_dir0 = 2*pi*rand();
    dir = [sin(theta0)*cos(phi_dir0), sin(theta0)*sin(phi_dir0), cos(theta0)];
    dir = dir / norm(dir);
    
    weight = 1.0;
    
    % --- Step B: 散射追踪 ---
    for s = 1:max_interactions
        % 1. 采样散射步长 (按总衰减系数 c)
        step_dist = -log(rand()) / param.coef_c; 
        
        % 2. 接收平面判定
        dist_to_plane = (Rx_Pos(3) - pos(3)) / dir(3);
        
        % 防止背向散射导致除以负数(虽然 dir(3) 负数时 dist_to_plane 也是负数，逻辑上不会触发，但加个保护)
        if dir(3) <= 0, dist_to_plane = inf; end
        
        if dist_to_plane > 0 && dist_to_plane < step_dist
            step_dist = dist_to_plane;
            is_hit_plane = true;
        else
            is_hit_plane = false;
        end
        
        % 3. --- 动态相位屏传输 (Ray Marching) ---
        num_screens = floor(step_dist / delta_z_screen);
        seg_left = step_dist;
        
        for k = 1 : max(1, num_screens)
            d_sub = (num_screens == 0) * step_dist + (num_screens > 0) * min(seg_left, delta_z_screen);
            pos = pos + dir * d_sub;
            seg_left = seg_left - d_sub;
            
            if weight < 1e-9, break; end
            
            % 穿过相位屏时发生偏折
            if num_screens > 0 && seg_left >= 0
                pool_idx = randi(Num_Screen_Pool);
                idx_x = mod(round((pos(1) - x_axis(1))/dx), N_grid) + 1;
                idx_y = mod(round((pos(2) - x_axis(1))/dx), N_grid) + 1;
                
                g_x = Screen_Pool(pool_idx).grad_x(idx_y, idx_x);
                g_y = Screen_Pool(pool_idx).grad_y(idx_y, idx_x);
                
                % 更新方向
                dir(1) = dir(1) + g_x / k_wave;
                dir(2) = dir(2) + g_y / k_wave;
                dir = dir / norm(dir);
            end
        end
        
        if weight < 1e-9, break; end
        
        % --- Step C: 处理到达接收平面 (直射/弹道光子) ---
        if is_hit_plane
            r_spot = sqrt((pos(1)-Rx_Pos(1))^2 + (pos(2)-Rx_Pos(2))^2);
            angle_inc = acos(dot(dir, [0,0,1])); 
            
            if r_spot <= Rx_Aperture/2 && angle_inc <= Rx_FOV
                % 仅统计未发生散射的光子(Ballistic)，散射光子由PIS统计
                % (若希望这里也统计散射后的撞击，需去掉 s==1 限制并注意PIS的去重)
                if s == 1 
                     Rx_Power_Total = Rx_Power_Total + weight;
                end
            end
            break; % 撞到平面后终止
        end
        
        % --- Step D: 重要性采样 (PIS) - 估计散射贡献 ---
        vec_to_rx = Rx_Pos - pos;
        dist_to_rx = norm(vec_to_rx);
        dir_to_rx = vec_to_rx / dist_to_rx;
        
        % 视场角判断 (接收机朝向 -Z)
        if acos(dot(-dir_to_rx, [0,0,-1])) <= Rx_FOV
            
            cos_theta_s = dot(dir, dir_to_rx);
            theta_s = acos(cos_theta_s);
            
            % [核心修改] 根据模式计算相函数概率密度 p(theta)
            if strcmp(param.phase_func, 'HG')
                p_phase = pdf_HG(cos_theta_s, param.g_HG);
            else
                % Mix 模式: 加权概率
                p_ray = pdf_Rayleigh(cos_theta_s, param.gamma);
                p_mie = pdf_Mie(cos_theta_s, param.g_mie, param.f_mie);
                ratio_r = param.coef_kr / param.coef_b;
                ratio_m = param.coef_km / param.coef_b;
                p_phase = ratio_r * p_ray + ratio_m * p_mie;
            end
            
            omega = Rx_Area / (dist_to_rx^2) * abs(dir_to_rx(3));
            prob_survival = exp(-param.coef_c * dist_to_rx);
            
            Rx_Power_Total = Rx_Power_Total + weight * param.albedo * p_phase * omega * prob_survival;
        end
        
        % --- Step E: 真实散射 (采样新方向) ---
        weight = weight * param.albedo; 
        
        % [核心修改] 根据模式采样新的散射方向
        if strcmp(param.phase_func, 'HG')
            dir = sample_dir_HG(dir, param.g_HG);
        else
            % Mix 模式: 先随机决定是 Rayleigh 还是 Mie
            if rand() < (param.coef_kr / param.coef_b)
                % 发生 Rayleigh 散射
                dir = sample_dir_Rayleigh(dir, param.gamma);
            else
                % 发生 Mie 散射
                dir = sample_dir_Mie(dir, param.g_mie, param.f_mie);
            end
        end
    end
    
    if mod(p, N_packets/10) == 0, fprintf('   进度: %.0f%%\n', p/N_packets*100); end
end
toc;

Rx_Power_Normalized = Rx_Power_Total / N_packets;
fprintf('\n=== 仿真结果 [%s] ===\n', param.phase_func);
fprintf('接收功率 (归一化): %.6e\n', Rx_Power_Normalized);
fprintf('链路损耗: %.2f dB\n', 10*log10(Rx_Power_Normalized));


%% ================= 4. 概率密度函数 (PDFs) =================
% 注意: 这里的 PDF 是关于立体角的归一化值，即 int p(theta) dOmega = 1

function p = pdf_HG(cos_theta, g)
    p = (1 - g^2) ./ (4*pi * (1 + g^2 - 2*g*cos_theta).^1.5);
end

function p = pdf_Rayleigh(cos_theta, gamma)
    % 参考文献公式 (3) / 2pi (因为积分 dPhi = 2pi)
    % f_Theta_Rayleigh 返回的是 2*pi*p(theta)*sin(theta)
    % 这里我们需要 p(theta) per steradian
    % p(theta) = [3/(16*pi*(1+2gamma))] * (1+3gamma + (1-gamma)cos^2)
    numerator = 3 * (1 + 3*gamma + (1-gamma)*cos_theta^2);
    denominator = 16 * pi * (1 + 2*gamma);
    p = numerator / denominator;
end

function p = pdf_Mie(cos_theta, g, f)
    % Generalized HG (TTHG similar form)
    % 参考文献公式 (4) / 2pi
    term1 = (1 - g^2) ./ (1 + g^2 - 2*g*cos_theta).^1.5;
    term2 = f * 0.5 * (3*cos_theta^2 - 1) ./ (1 + g^2)^1.5;
    p = (1/(4*pi)) * (1 - g^2)/2 * (term1 + term2) * 2 / (1-g^2); % 修正系数
    % 简化直接引用标准 TTHG 或论文公式:
    % 论文公式(4): f_Theta = ... sin(theta). 
    % p(theta) = f_Theta / (2*pi*sin(theta))
    val = (1 - g^2)/2 * ( 1./(1 + g^2 - 2*g*cos_theta).^1.5 + ...
           f * 0.5 * (3*cos_theta^2 - 1) / (1 + g^2)^1.5 );
    p = val / (2*pi);
end

%% ================= 5. 采样函数 (Samplers) =================

% --- HG 采样 (解析解) ---
function dir_new = sample_dir_HG(dir, g)
    rand_val = rand();
    if abs(g) < 1e-6
        cos_theta = 2*rand_val - 1;
    else
        temp = (1 - g^2) / (1 - g + 2*g*rand_val);
        cos_theta = (1 + g^2 - temp^2) / (2*g);
    end
    
    % 限制范围防止数值误差
    cos_theta = max(min(cos_theta, 1), -1);
    theta = acos(cos_theta);
    phi = 2*pi*rand();
    dir_new = rotate_direction(dir, theta, phi);
end

% --- Rayleigh 采样 (拒绝采样法) ---
function dir_new = sample_dir_Rayleigh(dir, gamma)
    % PDF shape ~ (1 + 3gamma + (1-gamma)cos^2(theta)) * sin(theta)
    % 包络函数 M 可以取最大值
    % 这里简化使用拒绝采样对 cos_theta 进行采样
    while true
        cos_t = 2*rand() - 1; % Uniform [-1, 1]
        % 目标 PDF (未归一化部分)
        p_val = 1 + 3*gamma + (1-gamma)*cos_t^2;
        p_max = 1 + 3*gamma + (1-gamma); % max at cos=1 or -1
        
        if rand() * p_max <= p_val
            theta = acos(cos_t);
            break;
        end
    end
    phi = 2*pi*rand();
    dir_new = rotate_direction(dir, theta, phi);
end

% --- Mie 采样 (拒绝采样法) ---
% 使用 Generalized HG 分布
function dir_new = sample_dir_Mie(dir, g, f)
    % 这是一个复杂的分布，使用拒绝采样
    % 包络函数可以用一个简单的 HG 分布或者均匀分布(效率低)
    % 这里为简单起见，对比 pdf_Mie(theta) 和一个常数 M
    
    % 估算最大值 (通常在 theta=0)
    max_val = pdf_Mie(1, g, f) * 2 * pi; % 去掉2pi因子方便计算
    
    while true
        cos_t = 2*rand() - 1; % Uniform [-1, 1] for cos_theta
        % 计算当前角度的非归一化概率
        % 借用 pdf_Mie 函数计算核心值
        current_val = pdf_Mie(cos_t, g, f) * 2 * pi;
        
        if rand() * max_val <= current_val
            theta = acos(cos_t);
            break;
        end
    end
    phi = 2*pi*rand();
    dir_new = rotate_direction(dir, theta, phi);
end

% --- 通用旋转函数 ---
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