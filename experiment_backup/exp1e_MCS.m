% =========================================================================
% 水下光通信 2D 空间光斑扩展可视化 (Analog Monte Carlo)
% 目的: 验证大发散角 + 水体多次散射导致的光斑宏观扩散现象
% 距离: 50 m
% =========================================================================
clear; clc; close all;


%% 1. 系统与环境物理参数设定 (严格对齐论文 Jerlov I 水质)
% 介质参数 (520 nm, Jerlov I)
a = 0.046;                  % 吸收系数 (m^-1)
b = 0.00205;                % 散射系数 (m^-1)
c = a + b;                  % 衰减系数 (0.04805 m^-1)
v = 2.237e8;                % 水中光速 (m/s) 

L = 100;                    % 传输距离 (m) [核心修改：对齐论文的 100 m 距离]

% 光源宏观参数
lambda = 520e-9;            % 波长 520 nm
w0 = 0.002;                 % 初始束腰半径 2mm
div_angle_deg = 1.0;        % 宏观全发散角 1 度

% 仿真控制参数
max_scat = 100;             % 允许的最大散射次数
N_photons = 2e6;            % 仿真光子总数 (200万即可获得平滑的空间分布)
W_th = 1e-6;                % 俄罗斯轮盘赌阈值
p_survive = 0.1;            % 存活概率

%% 2. Fournier-Forand (FF) 相函数与 LUT 预计算
fprintf('正在预计算相函数与查找表 (LUT)...\n');
n_water = 1.1549;            
mu = 3.5688;                
theta_array = logspace(log10(1e-6), log10(pi), 5000); 

% 计算 FF 相函数 CDF
nu_ff = (3 - mu) / 2; 
delta = (4 / (3 * (n_water - 1)^2)) .* sin(theta_array/2).^2; 
delta_180 = (4 / (3 * (n_water - 1)^2)); 
term1 = nu_ff .* (1 - delta) - (1 - delta.^nu_ff) + ...
        (delta .* (1 - delta.^nu_ff) - nu_ff .* (1 - delta)) .* (sin(theta_array/2).^(-2));
part1 = (1 ./ (4 * pi .* (1 - delta).^2 .* delta.^nu_ff)) .* term1;
part2 = ((1 - delta_180^nu_ff) / (16 * pi * (delta_180 - 1) * delta_180^nu_ff)) .* (3 .* cos(theta_array).^2 - 1);
beta_FF = max(0, part1 + part2);

pdf = 2 * pi * beta_FF .* sin(theta_array);
cdf_FF = cumtrapz(theta_array, pdf);
cdf_FF = cdf_FF / cdf_FF(end);

% 构建 LUT
M_LUT = 1e6;                
P_grid = linspace(0, 1, M_LUT); 
[unique_cdf, unique_idx] = unique(cdf_FF); % 去重以支持插值
theta_LUT = interp1(unique_cdf, theta_array(unique_idx), P_grid, 'linear', 'extrap');

%% 3. 数据记录结构初始化
% 记录所有成功穿过 Z = L 平面的光子 X, Y 坐标及权重
hit_X = zeros(1, N_photons);
hit_Y = zeros(1, N_photons);
hit_W = zeros(1, N_photons);
hit_count = 0;

%% 4. 物理光子追踪主循环
fprintf('启动物理射线追踪 (距离: %d m, 发散角: %.1f 度)...\n', L, div_angle_deg);
tic;
theta_half_macro = (div_angle_deg / 2) * (pi / 180); % 转换为弧度半角

for i = 1:N_photons
    % --- 初始物理发射 (宏观几何发散抽样) ---
    r_init = w0 * sqrt(-0.5 * log(rand())); 
    phi_init = 2 * pi * rand();
    pos = [r_init * cos(phi_init), r_init * sin(phi_init), 0]; 
    
    % 使用宏观发散半角进行高斯抽样
    U = theta_half_macro * sqrt(-0.5 * log(rand()));
    psi_ini = 2 * pi * rand();
    dir = [sin(U)*cos(psi_ini), sin(U)*sin(psi_ini), cos(U)]; 
    
    weight = 1.0;       
    scat_count = 0;   
    
    % --- 传输阶段 ---
    while scat_count < max_scat
        % 生成游走步长 (基于散射系数 b)
        delta_s = -log(rand()) / b; 
        
        % 越界检测: 是否穿过接收探测平面 Z = L
        if pos(3) + dir(3)*delta_s >= L
            % 计算刚好击中平面时的精确剩余步长
            d_remain = (L - pos(3)) / dir(3);
            hit_pos = pos + dir * d_remain;
            
            % 更新最终吸收权重并记录
            final_weight = weight * exp(-a * d_remain);
            
            hit_count = hit_count + 1;
            hit_X(hit_count) = hit_pos(1);
            hit_Y(hit_count) = hit_pos(2);
            hit_W(hit_count) = final_weight;
            break; 
        end
        
        % 状态更新
        pos = pos + dir * delta_s;
        weight = weight * exp(-a * delta_s);
        
        % 俄罗斯轮盘赌 
        if weight < W_th
            if rand() <= p_survive
                weight = weight / p_survive; 
            else
                break; 
            end
        end
        
        scat_count = scat_count + 1;
        
        % 散射角抽样与旋转
        idx_LUT = floor(rand() * (M_LUT - 1)) + 1;
        theta_s = theta_LUT(idx_LUT);
        phi_s = 2 * pi * rand(); 
        
        % 坐标旋转 (优化版)
        ux = dir(1); uy = dir(2); uz = dir(3);
        costh = cos(theta_s); sinth = sin(theta_s);
        cosph = cos(phi_s); sinph = sin(phi_s);
        
        if abs(uz) > 0.99999
            dir = [sinth * cosph, sinth * sinph, sign(uz) * costh];
        else
            temp = sqrt(1 - uz^2);
            ux_new = sinth * (ux * uz * cosph - uy * sinph) / temp + ux * costh;
            uy_new = sinth * (uy * uz * cosph + ux * sinph) / temp + uy * costh;
            uz_new = -sinth * cosph * temp + uz * costh;
            dir = [ux_new, uy_new, uz_new];
        end
        dir = dir / norm(dir);
    end
end
t_run = toc;
fprintf('追踪完成！有效击中光子数: %d, 耗时: %.2f 秒\n', hit_count, t_run);

% 截断有效数据
hit_X = hit_X(1:hit_count);
hit_Y = hit_Y(1:hit_count);
hit_W = hit_W(1:hit_count);

%% 5. 空间数据统计与光斑半径计算
r_hit = sqrt(hit_X.^2 + hit_Y.^2);
R_geo_ideal = w0 + L * tan(theta_half_macro);

% =========================================================================
% [方法 A]: 论文的 1m x 1m 粗网格平滑法 (模拟论文的 Peak 1/e^2 机制)
% =========================================================================
% 强制使用大网格(例如 dr = 0.5m) 来平滑尖峰
r_edges_coarse = linspace(0, 5, 11); 
r_centers_coarse = (r_edges_coarse(1:end-1) + r_edges_coarse(2:end)) / 2;
irradiance_coarse = zeros(1, length(r_centers_coarse));
for k = 1:length(r_centers_coarse)
    idx = (r_hit >= r_edges_coarse(k)) & (r_hit < r_edges_coarse(k+1));
    ring_area = pi * (r_edges_coarse(k+1)^2 - r_edges_coarse(k)^2);
    irradiance_coarse(k) = sum(hit_W(idx)) / ring_area;
end
I_norm_coarse = irradiance_coarse / max(irradiance_coarse);
boundary_idx_coarse = find(I_norm_coarse < exp(-2), 1);
if ~isempty(boundary_idx_coarse)
    R_scatter_coarse = r_centers_coarse(boundary_idx_coarse);
else
    R_scatter_coarse = r_centers_coarse(end);
end

% =========================================================================
% [方法 B]: 学术界通用的 86.5% 能量包围度法 (Encircled Energy, 极高稳定性)
% =========================================================================
num_r_bins = 500; % 恢复高精度网格
r_edges = linspace(0, 10, num_r_bins+1);
r_centers = (r_edges(1:end-1) + r_edges(2:end)) / 2;
energy_bins = zeros(1, num_r_bins);

for k = 1:num_r_bins
    idx = (r_hit >= r_edges(k)) & (r_hit < r_edges(k+1));
    energy_bins(k) = sum(hit_W(idx));
end

% 计算累积能量占比 (CDF)
cumulative_energy = cumsum(energy_bins) / sum(energy_bins);

% 寻找包含 1 - 1/e^2 (约 86.47%) 能量的物理半径
ee_threshold = 1 - exp(-2);
ee_idx = find(cumulative_energy >= ee_threshold, 1);
if ~isempty(ee_idx)
    R_scatter_EE = r_centers(ee_idx);
else
    R_scatter_EE = r_centers(end);
end

fprintf('\n=== 光斑展宽分析 (Z = %d m) ===\n', L);
fprintf('理想纯几何投影半径: %.3f m\n', R_geo_ideal);
fprintf('实际水体散射半径 (模拟论文粗网格 1/e^2): %.3f m (放大 %.1f 倍)\n', R_scatter_coarse, R_scatter_coarse / R_geo_ideal);
fprintf('实际水体散射半径 (高精度能量包围度 86.5%%): %.3f m (放大 %.1f 倍)\n', R_scatter_EE, R_scatter_EE / R_geo_ideal);
% ===== 接在第5部分的 fprintf 之后补充高精度光强计算用于绘图 =====
% 计算高精度网格下的一维径向光强 (Irradiance) 用于子图2显示
irradiance_fine = zeros(1, num_r_bins);
for k = 1:num_r_bins
    ring_area = pi * (r_edges(k+1)^2 - r_edges(k)^2);
    irradiance_fine(k) = energy_bins(k) / ring_area;
end
irradiance_norm = irradiance_fine / max(irradiance_fine);
threshold = exp(-2); % 传统的 1/e^2 阈值

%% 6. 绘图可视化
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

% 子图 1: 2D 散点密度图 (随机抽取几万个点展示即可，避免卡顿)
subplot(1, 2, 1);
plot_pts = min(hit_count, 30000); 
scatter(hit_X(1:plot_pts), hit_Y(1:plot_pts), 2, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
hold on;

% 绘制几何边界和 能量包围度 (EE) 实际边界
viscircles([0, 0], R_geo_ideal, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
viscircles([0, 0], R_scatter_EE, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
axis equal; grid on;

% 自适应坐标轴范围
plot_limit = max(2, R_scatter_EE * 1.2);
xlim([-plot_limit, plot_limit]); ylim([-plot_limit, plot_limit]);

xlabel('X Position (m)', 'FontWeight', 'bold');
ylabel('Y Position (m)', 'FontWeight', 'bold');
title('2D Spatial Photon Distribution at Rx Plane');
legend('Photon Hits', 'Ideal Geometric Beam', '86.5% Encircled Energy Boundary', 'Location', 'northeast');

% 子图 2: 1D 径向光强剖面图
subplot(1, 2, 2);
plot(r_centers, irradiance_norm, 'b-', 'LineWidth', 2);
hold on;
xline(R_geo_ideal, 'k--', 'LineWidth', 1.5);
xline(R_scatter_EE, 'r-', 'LineWidth', 2);
yline(threshold, 'r:', 'LineWidth', 1.5);
grid on;
xlim([0, max(2, R_scatter_EE * 1.5)]);
xlabel('Radial Distance r (m)', 'FontWeight', 'bold');
ylabel('Normalized Spatial Irradiance I(r)', 'FontWeight', 'bold');
title('Radial Irradiance Profile I(r)');
legend('Simulated Irradiance (Fine Grid)', 'Geometric Radius', ...
    sprintf('Encircled Energy Radius (%.2fm)', R_scatter_EE), '1/e^2 Peak Threshold', ...
    'Location', 'northeast');