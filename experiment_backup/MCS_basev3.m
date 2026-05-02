% =========================================================================
% 高性能水下光无线通信 (UOWC) 信道冲激响应蒙特卡洛仿真
% 集成: 静态时间窗、散射角 LUT 查表法、俄罗斯轮盘赌 (RR)
% 对应文献: Impulse response modeling for underwater optical wireless channels
% =========================================================================
clear; clc; close all;

%% 1. 系统与环境物理参数设定
% 典型港口水质 (Harbor water) 参数
a = 0.366;                  % 吸收系数 (m^-1) [cite: 239]
b = 1.829;                  % 散射系数 (m^-1) [cite: 239]
v = 2.237e8;                % 水中光速 (m/s) [cite: 236]

L = 10.93;                  % 传输距离 (m) [cite: 239]
Aper_D = 0.4;               % 接收器孔径直径 (m) [cite: 251]
r_0 = Aper_D / 2;           % 接收器半径
FOV_deg = 20;               % 接收器视场角 (度) [cite: 235]
FOV = FOV_deg * pi / 180;   % 视场角 (弧度)
max_scat = 300;             % 最大散射次数截断 [cite: 141]
N_photons = 1e6;            % 仿真发射光子总数 

%% 2. 核心数据结构与预计算 (Pre-computation)
% 2.1 Fournier-Forand 相函数 (FFPF) 预计算
n_water = 1.339;            
mu = 3.5835;                
theta_array = linspace(1e-5, pi, 2000); 
[~, cdf_FF] = generate_FFPF_CDF(theta_array, n_water, mu);

% 2.2 构建 O(1) 复杂度散射角查找表 (LUT) (from reasoning)
M_LUT = 1e6;                % 查找表分辨率
P_grid = linspace(0, 1, M_LUT); 
theta_LUT = interp1(cdf_FF, theta_array, P_grid, 'linear', 'extrap');

% 2.3 范式 A: 静态物理时间窗口预分配 (from reasoning)
max_delta_t_ns = 100;       % 相对时延探测上限设定为 100 ns
num_bins = 500;             % 直方图时间仓数量
intensity_hist = zeros(1, num_bins); % O(1) 静态空间分配
bin_width = max_delta_t_ns / num_bins;
bin_centers = linspace(bin_width/2, max_delta_t_ns - bin_width/2, num_bins);

% 2.4 俄罗斯轮盘赌 (RR) 参数 (from reasoning)
W_th = 1e-5;                % 极低权重触发阈值
p_survive = 0.1;            % 存活概率

%% 3. 光子游走主循环 (并行安全结构)
fprintf('启动高性能蒙特卡洛仿真，发射光子数: %d...\n', N_photons);
tic;

% 若需进一步加速，可直接将 for 替换为 parfor (需预设并行池)
for i = 1:N_photons
    % 初始化：光子位于原点，准直发射 [cite: 115]
    pos = [0, 0, 0]; 
    dir = [0, 0, 1]; 
    
    path_len = 0;     % 累计光程 
    weight = 1;       % 初始权重
    scat_count = 0;   
    
    % 传输阶段
    while scat_count < max_scat
        % 生成服从负指数分布的无碰撞传输步长 [cite: 116, 127]
        rand_U = rand();
        delta_s = -log(1 - rand_U) / b; 
        
        % 越界检测: 是否穿过接收探测平面 Z = L [cite: 153]
        if pos(3) + dir(3)*delta_s >= L
            % 计算刚好击中平面时的精确剩余步长
            d_remain = (L - pos(3)) / dir(3);
            hit_pos = pos + dir * d_remain;
            hit_path_len = path_len + d_remain;
            
            % 吸收与散射的分解: 基于最终总光程独立计算吸收衰减 [cite: 145, 146]
            final_weight = weight * exp(-a * d_remain);
            
            % 接收器几何尺寸与视场角验证 [cite: 159]
            r_hit = sqrt(hit_pos(1)^2 + hit_pos(2)^2);
            theta_inc = acos(dir(3)); 
            
            if (r_hit <= r_0) && (theta_inc <= FOV/2)
                % 计算绝对飞行时间与相对时延
                hit_time = hit_path_len / v;
                delta_t_ns = (hit_time - (L / v)) * 1e9;
                
                % 范式 A: 实时直方图时间仓映射
                if delta_t_ns >= 0 && delta_t_ns < max_delta_t_ns
                    bin_idx = floor(delta_t_ns / bin_width) + 1;
                    % 累加权重 (对于 parfor，此处强度数组需采用特殊同步或规约形式)
                    intensity_hist(bin_idx) = intensity_hist(bin_idx) + final_weight;
                end
            end
            break; % 到达接收平面，无论是否被接收均终止该光子追踪
        end
        
        % 状态更新: 执行步长推移
        pos = pos + dir * delta_s;
        path_len = path_len + delta_s;
        
        % 每步更新吸收权重 (为触发轮盘赌做准备)
        weight = weight * exp(-a * delta_s);
        
        % 无偏方差缩减: 俄罗斯轮盘赌机制 (from reasoning)
        if weight < W_th
            if rand() <= p_survive
                weight = weight / p_survive; % 存活，权重放大以确保期望无偏
            else
                break; % 死亡，提前释放算力
            end
        end
        
        scat_count = scat_count + 1;
        
        % O(1) LUT 查表法生成散射天顶角 (from reasoning)
        rand_X = rand();
        idx_LUT = floor(rand_X * (M_LUT - 1)) + 1;
        theta_s = theta_LUT(idx_LUT);
        phi_s = 2 * pi * rand(); 
        
        % 局部球坐标系向全局笛卡尔坐标系的转移矩阵映射
        dir = update_direction(dir, theta_s, phi_s); 
    end
end
toc;

%% 4. 后处理与对数坐标域绘图
% 强度归一化 (此处代表单位时间仓内的平均光子权重)
intensity_hist = intensity_hist / N_photons;

% 过滤掉未接收到光子的零值区间以避免对数绘图警告
valid_idx = intensity_hist > 0;

figure;
semilogy(bin_centers(valid_idx), intensity_hist(valid_idx), 'k.-', 'LineWidth', 1.5, 'MarkerSize', 10);
grid on;
xlim([0, 2.5]); % 依文献图4设定横坐标范围以便于对照
xlabel('\Delta t [ns]', 'FontWeight', 'bold');
ylabel('Intensity', 'FontWeight', 'bold');
title(sprintf('Monte Carlo UOWC Impulse Response (LUT+RR+Static Bin)\nHarbor Water, L=%.2fm, FOV=%d^{\\circ}', L, FOV_deg));
set(gca, 'FontSize', 12);

%% =========================================================================
% 辅助函数：解析相函数累积分布计算与空间矩阵旋转
% =========================================================================

function [pdf, cdf] = generate_FFPF_CDF(theta, n, mu)
    nu = (3 - mu) / 2; % [cite: 81]
    delta = (4 / (3 * (n - 1)^2)) .* sin(theta/2).^2; % [cite: 81]
    delta_180 = (4 / (3 * (n - 1)^2)); % [cite: 81]
    
    term1 = nu .* (1 - delta) - (1 - delta.^nu) + ...
            (delta .* (1 - delta.^nu) - nu .* (1 - delta)) .* (sin(theta/2).^(-2));
    part1 = (1 ./ (4 * pi .* (1 - delta).^2 .* delta.^nu)) .* term1;
    
    part2 = ((1 - delta_180^nu) / (16 * pi * (delta_180 - 1) * delta_180^nu)) .* (3 .* cos(theta).^2 - 1);
    
    beta_FF = part1 + part2;
    pdf = 2 * pi * beta_FF .* sin(theta);
    cdf = cumtrapz(theta, pdf);
    cdf = cdf / cdf(end);
end

function dir_new = update_direction(dir_old, theta, phi)
    ux = dir_old(1); uy = dir_old(2); uz = dir_old(3);
    costh = cos(theta); sinth = sin(theta);
    cosph = cos(phi); sinph = sin(phi);
    
    if abs(uz) > 0.99999
        dir_new = [sinth * cosph, sinth * sinph, sign(uz) * costh];
    else
        temp = sqrt(1 - uz^2);
        ux_new = sinth * (ux * uz * cosph - uy * sinph) / temp + ux * costh;
        uy_new = sinth * (uy * uz * cosph + ux * sinph) / temp + uy * costh;
        uz_new = -sinth * cosph * temp + uz * costh;
        dir_new = [ux_new, uy_new, uz_new];
    end
    dir_new = dir_new / norm(dir_new);
end