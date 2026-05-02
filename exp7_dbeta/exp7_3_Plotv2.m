%% Exp7: Combined Plot (4 Algorithms x 3 Water Types at Fixed Off-axis Angle)
clc; clear; close all;

% =========================================================
% 1. 数据加载与容错处理
% =========================================================
has_data = false;
try
    d_WCI = load('data_exp7_WCIMC.mat');
    d_MCS = load('data_exp7_MCS.mat');
    water_types = d_WCI.water_types;
    dist_cell = d_WCI.dist_cell;
    has_data = true;
catch
    warning('未能找到 data_exp7_*.mat 文件，将采用物理近似数据进行视觉效果演示。');
    water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
    dist_cell = {10:5:40, 10:5:40, 10:5:40};
end

num_W = length(water_types);
fixed_a = 1; % 固定展示第一组偏轴角度 (通常为 0 度)

% =========================================================
% 2. 视觉映射规范
% =========================================================
% 颜色区分水质 (宏观聚类)
c_water = {'#0072BD', '#D95319', '#EDB120'}; 

% 线型与标记区分算法 (微观解析)
lines   = {'-', '--', ':', '-.'};
markers = {'o', 'd', 's', '^'};
alg_names = {'WCI-MC (Strong)', 'WCI-MC (Weak)', 'WCI-MC (No Turb)', 'SA-MCS'};

figure('Name', 'Exp7: Algorithms vs Water Types', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);

h_leg_alg = gobjects(1, 4); % 预分配算法图例句柄

% =========================================================
% 3. 曲线绘制与空间排版
% =========================================================
for w = 1:num_W
    dist_axis = dist_cell{w};
    c_val = c_water{w};
    m_size = 6;
    
    % --- 数据提取或生成 ---
    if has_data
        y_strong = -d_WCI.PL_Cell_WCI_Strong{w}(fixed_a, :);
        y_weak   = -d_WCI.PL_Cell_WCI_Weak{w}(fixed_a, :);
        y_none   = -d_WCI.PL_Cell_WCI_None{w}(fixed_a, :);
        y_mcs    = -d_MCS.PL_Cell_MCS{w}(fixed_a, :);
    else
        % 用于预览的近似衰减斜率模拟
        base_slope = [0.15, 0.8, 1.5]; base_att = base_slope(w) * dist_axis + (w-1)*20;
        y_strong = base_att + 2.5; y_weak = base_att + 1.0; 
        y_none   = base_att;       y_mcs  = base_att - 0.5;
    end
    
    % --- 绘制 4 种算法 ---
    p1 = plot(dist_axis, y_strong, lines{1}, 'Color', c_val, 'Marker', markers{1}, 'LineWidth', 1.5, 'MarkerSize', m_size, 'MarkerFaceColor', c_val);
    p2 = plot(dist_axis, y_weak,   lines{2}, 'Color', c_val, 'Marker', markers{2}, 'LineWidth', 1.2, 'MarkerSize', m_size, 'MarkerFaceColor', 'w');
    p3 = plot(dist_axis, y_none,   lines{3}, 'Color', c_val, 'Marker', markers{3}, 'LineWidth', 1.5, 'MarkerSize', m_size, 'MarkerFaceColor', 'w');
    p4 = plot(dist_axis, y_mcs,    lines{4}, 'Color', c_val, 'Marker', markers{4}, 'LineWidth', 1.0, 'MarkerSize', m_size, 'MarkerFaceColor', 'w');
    
    % 提取第一种水质的句柄供图例使用
    if w == 1
        h_leg_alg = [p1, p2, p3, p4];
    end
    
    % --- 空间文本锚定 (标注水质) ---
    % 将文本锚定在 'WCI-MC (Strong)' 曲线的倒数第一个数据点附近
    idx_txt = length(dist_axis);
    x_txt = dist_axis(idx_txt);
    y_txt = y_strong(idx_txt);
    
    % Y轴方向微调，避免与端点重叠
    text(x_txt, y_txt + 3, water_types{w}, 'FontSize', 12, ...
        'FontName', 'Times New Roman', 'FontWeight', 'bold', 'Color', c_val, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

% =========================================================
% 4. 坐标系与图例精调
% =========================================================
xlabel('Transmission Distance \it{L} \rm{(m)}', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Path Loss (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Performance Comparison across Water Types (\theta_{err} = 0^\circ)', 'FontSize', 15, 'FontName', 'Times New Roman');

% 由于颜色图例被文本替代，此处仅需保留算法图例
lgd = legend(h_leg_alg, alg_names, 'Location', 'northwest', 'FontSize', 11, 'FontName', 'Times New Roman');
title(lgd, 'Algorithm Models');

set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');

% 自动适配全局 X 轴范围
min_dist = min(cellfun(@min, dist_cell));
max_dist = max(cellfun(@max, dist_cell));
xlim([min_dist, max_dist + 5]); % 右侧预留一定空间供文本显示

% exportgraphics(gcf, 'Figure_Exp7_Algorithms_vs_Water.png', 'Resolution', 600);