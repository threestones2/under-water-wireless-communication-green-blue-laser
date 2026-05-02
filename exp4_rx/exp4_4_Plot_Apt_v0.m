%% Exp4: Data Visualization (Path Loss & Increments vs Aperture)
clc; clear; close all;

try
    d_WCI = load('data_exp4_apt_WCIMC.mat');
    d_HG  = load('data_exp4_apt_HGPIS.mat');
    d_MCS = load('data_exp4_apt_MCS.mat');
catch
    error('数据加载失败，请务必先运行前三个仿真脚本生成 .mat 数据文件。');
end

dist_axis = d_WCI.dist_axis;
apertures_m = d_WCI.apertures_m;
num_apt = length(apertures_m);
c_val = d_HG.coef_c; 
w0 = 0.002; 
theta_half_div = 0.1 * pi / 180 / 2;

% 视觉色系规范
c_Strong = '#D95319'; c_Weak   = '#EDB120'; 
c_MCS    = '#77AC30'; c_Beer   = '#555555'; 

% =========================================================================
% 线型设计：
% 0.05m -> '--' (虚线)
% 0.10m -> '-'  (实线)
% 0.20m -> '-.' (点划线)
% 线宽与标记大小已统一一致
% =========================================================================
line_styles  = {'--', '-', '-.'};

fig_width = 640; fig_height = 480;

% =========================================================================
% 图 1：绝对路径衰减绘图 (含局部放大图)
% =========================================================================
figure('Name', 'Path Loss vs Aperture', 'Color', 'w', 'Position', [100, 100, fig_width, fig_height]);
hold on;

for i = 1:num_apt
    r_rx = apertures_m(i) / 2;
    W_geo = w0 + dist_axis .* tan(theta_half_div);
    geom_loss = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo.^2));
    PL_Beer = -10 * log10(max(exp(-c_val .* dist_axis) .* geom_loss, 1e-300));
    
    ls = line_styles{i}; 
    
    plot(dist_axis, -d_WCI.PL_Cell_WCI_Strong{i}, ls, 'Color', c_Strong, ...
        'Marker', 'o', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Strong);
    plot(dist_axis, -d_WCI.PL_Cell_WCI_Weak{i}, ls, 'Color', c_Weak, ...
        'Marker', 'd', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
        
    % 注意：已按照要求彻底移除了 HG-PIS 的绘制
    
    plot(dist_axis, -d_MCS.PL_Cell_MCS{i}, ls, 'Color', c_MCS, ...
        'Marker', '^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_MCS);
    plot(dist_axis, PL_Beer, ls, 'Color', c_Beer, ...
        'Marker', '*', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Beer);
end

grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Path Loss (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Path Loss under Different Receiver Apertures', 'FontSize', 15, 'FontName', 'Times New Roman');
xlim([min(dist_axis), max(dist_axis)]);

% 主图图例 (采用分离式构建，剔除 HG-PIS)
h_alg1 = plot(nan, nan, '-', 'Color', c_Strong, 'Marker', 'o', 'LineWidth', 1.5, 'MarkerFaceColor', c_Strong);
h_alg2 = plot(nan, nan, '-', 'Color', c_Weak, 'Marker', 'd', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
h_alg4 = plot(nan, nan, '-', 'Color', c_MCS, 'Marker', '^', 'LineWidth', 1.5, 'MarkerFaceColor', c_MCS);
h_alg5 = plot(nan, nan, '-', 'Color', c_Beer, 'Marker', '*', 'LineWidth', 1.5, 'MarkerFaceColor', c_Beer);

h_apt1 = plot(nan, nan, line_styles{1}, 'Color', 'k', 'LineWidth', 1.5);
h_apt2 = plot(nan, nan, line_styles{2}, 'Color', 'k', 'LineWidth', 1.5);
h_apt3 = plot(nan, nan, line_styles{3}, 'Color', 'k', 'LineWidth', 1.5);

legend([h_alg1, h_alg2, h_alg4, h_alg5, h_apt1, h_apt2, h_apt3], ...
    {'WCI-MC (Strong)', 'WCI-MC (Weak)', 'SA-MCS (w/o Turb)', 'Beer Law', ...
     sprintf('D_{rx} = %.2f m', apertures_m(1)), ...
     sprintf('D_{rx} = %.2f m', apertures_m(2)), ...
     sprintf('D_{rx} = %.2f m', apertures_m(3))}, ...
    'Location', 'southwest', 'FontSize', 11, 'NumColumns', 2, 'FontName', 'Times New Roman');

% -------------------------------------------------------------------------
% 2. 添加局部放大图 (Inset Plot)
% -------------------------------------------------------------------------
x_zoom = [45, 60];
y_zoom = [25, 50];

ax_inset = axes('Position', [0.55, 0.22, 0.32, 0.28]); 
hold(ax_inset, 'on');
box(ax_inset, 'on');

for i = 1:num_apt
    r_rx = apertures_m(i) / 2;
    W_geo = w0 + dist_axis .* tan(theta_half_div);
    geom_loss = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo.^2));
    PL_Beer = -10 * log10(max(exp(-c_val .* dist_axis) .* geom_loss, 1e-300));
    
    ls = line_styles{i}; 
    
    plot(ax_inset, dist_axis, -d_WCI.PL_Cell_WCI_Strong{i}, ls, 'Color', c_Strong, ...
        'Marker', 'o', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', c_Strong);
    plot(ax_inset, dist_axis, -d_WCI.PL_Cell_WCI_Weak{i}, ls, 'Color', c_Weak, ...
        'Marker', 'd', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
        
    % 移除 HG-PIS
    
    plot(ax_inset, dist_axis, -d_MCS.PL_Cell_MCS{i}, ls, 'Color', c_MCS, ...
        'Marker', '^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', c_MCS);
    plot(ax_inset, dist_axis, PL_Beer, ls, 'Color', c_Beer, ...
        'Marker', '*', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', c_Beer);
end

xlim(ax_inset, x_zoom);
ylim(ax_inset, y_zoom);
grid(ax_inset, 'on');
set(ax_inset, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 10, 'LineWidth', 1, 'FontName', 'Times New Roman');


% =========================================================================
% 图 2：核心机制对比：以 SA-MCS 为基准的绝对衰减增量演化图
% =========================================================================
figure('Name', 'Increment relative to SA-MCS', 'Color', 'w', 'Position', [400, 200, fig_width, fig_height]);
hold on;

color_apertures = {'#7E2F8E', '#D95319', '#EDB120'};
h_strong = zeros(1, num_apt);
h_weak   = zeros(1, num_apt);
leg_labels_inc = {};
h_leg_all = [];

for i = 1:num_apt
    PL_MCS    = -d_MCS.PL_Cell_MCS{i};         
    PL_Strong = -d_WCI.PL_Cell_WCI_Strong{i};  
    PL_Weak   = -d_WCI.PL_Cell_WCI_Weak{i};    
    
    inc_strong = PL_Strong - PL_MCS;
    inc_weak   = PL_Weak - PL_MCS;
    
    % 强湍流增量 (实线，实心圆)
    h_strong(i) = plot(dist_axis, inc_strong, '-', 'Color', color_apertures{i}, ...
        'Marker', 'o', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', color_apertures{i});
    
    % 弱湍流增量 (长虚线，空心菱形)
    h_weak(i)   = plot(dist_axis, inc_weak, '--', 'Color', color_apertures{i}, ...
        'Marker', 'd', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
        
    leg_labels_inc{end+1} = sprintf('D_{rx}=%.2fm (Strong)', apertures_m(i));
    leg_labels_inc{end+1} = sprintf('D_{rx}=%.2fm (Weak)', apertures_m(i));
    h_leg_all = [h_leg_all, h_strong(i), h_weak(i)];
end

plot([min(dist_axis), max(dist_axis)], [0, 0], 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Path Loss Increment relative to SA-MCS (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Aperture Averaging Effect on Turbulence Penalty', 'FontSize', 15, 'FontName', 'Times New Roman');
xlim([min(dist_axis), max(dist_axis)]);

legend(h_leg_all, leg_labels_inc, 'Location', 'northwest', 'FontSize', 11, 'NumColumns', 2, 'FontName', 'Times New Roman');