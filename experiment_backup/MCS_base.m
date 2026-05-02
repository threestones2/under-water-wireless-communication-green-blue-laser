% =========================================================================
% 纯物理接收蒙特卡洛仿真 (Analog Monte Carlo)
% 对应 MC_MPS_v4 无湍流参数基准测试
% 集成: 物理高斯光束发射、静态时间窗、散射角 LUT、俄罗斯轮盘赌
% =========================================================================
clear; clc; close all;

%% 1. 系统与环境物理参数设定 (对齐 v4 代码的 Clear Ocean 无湍流参数)
% Petzold 清澈海水 (Clear Ocean) 参数
a = 0.114;                  % 吸收系数 (m^-1) 
b = 0.0374;                 % 散射系数 (m^-1) 
c = a + b;                  % 衰减系数
v = 2.237e8;                % 水中光速 (m/s) 

L = 50;                     % 传输距离 (m) 
Aper_D = 0.05;              % 接收器孔径直径 (m) (对应 5cm)
r_0 = Aper_D / 2;           % 接收器半径
FOV_deg = 5;                % 接收器视场角 (度) 
FOV = FOV_deg * pi / 180;   % 视场角 (弧度)
max_scat = 30;              % 最大散射次数截断
N_photons = 1e7;            % 仿真发射光子总数 (物理接收需要极高光子数)

% 光源参数
lambda = 514e-9;
w0 = 0.002;                 % 束腰半径 2mm

%% 2. 核心数据结构与预计算 (Pre-computation)
% 2.1 Fournier-Forand 相函数 (FFPF) 预计算 (对齐 v4 的 FF 参数)
n_water = 1.1549;            
mu = 3.5688;                
theta_array = linspace(1e-5, pi, 2000); 
[~, cdf_FF] = generate_FFPF_CDF(theta_array, n_water, mu);

% 2.2 构建 O(1) 复杂度散射角查找表 (LUT)
M_LUT = 1e6;                % 查找表分辨率
P_grid = linspace(0, 1, M_LUT); 
theta_LUT = interp1(cdf_FF, theta_array, P_grid, 'linear', 'extrap');

% 2.3 静态物理时间窗口预分配
max_delta_t_ns = 5;         % 清澈水质时延较小，探测上限设定为 5 ns
num_bins = 5000;             % 直方图时间仓数量
intensity_hist = zeros(1, num_bins); 
bin_width = max_delta_t_ns / num_bins;
bin_centers = linspace(bin_width/2, max_delta_t_ns - bin_width/2, num_bins);

% 2.4 俄罗斯轮盘赌 (RR) 参数
W_th = 1e-10;                % 极低权重触发阈值
p_survive = 0.1;            % 存活概率

%% 3. 光子游走主循环
fprintf('启动物理接收蒙特卡洛仿真，发射光子数: %d...\n', N_photons);
tic;

for i = 1:N_photons
    
    % --- 物理高斯光束空间与角度抽样 ---
    r_init = w0 * sqrt(-0.5 * log(rand())); 
    phi_init = 2 * pi * rand();
    pos = [r_init * cos(phi_init), r_init * sin(phi_init), 0]; 
    
    theta_div_physics = lambda / (pi * w0); 
    U = theta_div_physics * sqrt(-0.5 * log(rand()));
    psi_ini = 2 * pi * rand();
    dir = [sin(U)*cos(psi_ini), sin(U)*sin(psi_ini), cos(U)]; 
    
    path_len = 0;     % 累计光程 
    weight = 1;       % 初始权重
    scat_count = 0;   
    
    % 传输阶段
    while scat_count < max_scat
        % 生成服从负指数分布的无碰撞传输步长 (基于散射系数 b)
        rand_U = rand();
        delta_s = -log(1 - rand_U) / b; 
        
        % 越界检测: 是否穿过接收探测平面 Z = L
        if pos(3) + dir(3)*delta_s >= L
            % 计算刚好击中平面时的精确剩余步长
            d_remain = (L - pos(3)) / dir(3);
            hit_pos = pos + dir * d_remain;
            hit_path_len = path_len + d_remain;
            
            % 吸收衰减
            final_weight = weight * exp(-a * d_remain);
            
            % 接收器几何尺寸与视场角验证 (纯物理接收)
            r_hit = sqrt(hit_pos(1)^2 + hit_pos(2)^2);
            theta_inc = acos(dir(3)); 
            
            if (r_hit <= r_0) && (theta_inc <= FOV/2)
                hit_time = hit_path_len / v;
                delta_t_ns = (hit_time - (L / v)) * 1e9;
                
                if delta_t_ns >= 0 && delta_t_ns < max_delta_t_ns
                    bin_idx = floor(delta_t_ns / bin_width) + 1;
                    intensity_hist(bin_idx) = intensity_hist(bin_idx) + final_weight;
                end
            end
            break; % 到达接收平面，终止
        end
        
        % 状态更新
        pos = pos + dir * delta_s;
        path_len = path_len + delta_s;
        weight = weight * exp(-a * delta_s);
        
        % 俄罗斯轮盘赌机制 
        if weight < W_th
            if rand() <= p_survive
                weight = weight / p_survive; 
            else
                break; 
            end
        end
        
        scat_count = scat_count + 1;
        
        % 散射角查表
        rand_X = rand();
        idx_LUT = floor(rand_X * (M_LUT - 1)) + 1;
        theta_s = theta_LUT(idx_LUT);
        phi_s = 2 * pi * rand(); 
        
        % 坐标系旋转
        dir = update_direction(dir, theta_s, phi_s); 
    end
end
t_run = toc;
fprintf('仿真完成，耗时: %.2f 秒\n', t_run);

%% 4. 统计结果与绘图
P_rx_total = sum(intensity_hist) / N_photons;
Path_Loss_dB = -10 * log10(P_rx_total + eps);

fprintf('\n=== 物理接收模式结果 ===\n');
fprintf('总接收比例 (P_rx): %e\n', P_rx_total);
fprintf('理论比尔定律纯弹道损耗: %.2f dB\n', 10 * c * L / log(10)); % -10*log10(exp(-c*L))
fprintf('实际物理接收路径损耗 : %.2f dB\n', Path_Loss_dB);
fprintf('========================\n');

intensity_hist_norm = (intensity_hist / N_photons);
valid_idx = intensity_hist_norm > 0;

figure;
semilogy(bin_centers(valid_idx), intensity_hist_norm(valid_idx), 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Time Increment \Delta t (ns)', 'FontWeight', 'bold');
ylabel('Intensity', 'FontWeight', 'bold');
title(sprintf('Physical Reception MC (Clear Ocean, No Turb)\nPath Loss: %.2f dB', Path_Loss_dB));
set(gca, 'FontSize', 12);

%% =========================================================================
% 辅助函数
% =========================================================================
function [pdf, cdf] = generate_FFPF_CDF(theta, n, mu)
    nu = (3 - mu) / 2; 
    delta = (4 / (3 * (n - 1)^2)) .* sin(theta/2).^2; 
    delta_180 = (4 / (3 * (n - 1)^2)); 
    
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