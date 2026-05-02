%% Exp2: Plotting Macroscopic Channel Characteristics (Final Integrated Version)
% 包含: WCI-MC(强湍流/弱湍流/无湍流), SA-MCS(无湍流基准), Beer-Lambert解析解
clc; clear; close all;
        
% 1. 加载所有分离的数据文件
try
    d_WCI_S = load('data_exp2_WCIMC_Strong.mat');
    d_WCI_W = load('data_exp2_WCIMC_Weak.mat');
    d_None  = load('data_exp2_WCIMC_None.mat'); % 无湍流修正基准
    d_MCS   = load('data_exp2_MCS.mat');        % 传统 O(1) MCS 基准
catch
    error('未能找到数据文件，请确保已运行 Strong, Weak, None 及 MCS 对应的仿真脚本。');
end

% 2. 提取公共参数
dist_cell = d_WCI_S.dist_cell;
num_W = length(dist_cell);

coef_c_arr = [0.1514, 0.398, 2.190]; 
w0 = 0.002;                          
theta_half_div = 10 * pi / 180 / 2; 
Rx_Aperture = 0.05;                  
r_rx = Rx_Aperture / 2;              

% 统一配色方案
c_Strong = '#D95319'; 
c_Weak   = '#EDB120'; 
c_NoTurb = '#0072BD'; 
c_MCS    = '#77AC30'; 
c_Beer   = '#555555'; 

short_wt_names = {'Clear', 'Coastal', 'Harbor'}; 

% =========================================================
% 图 1: 路径损耗随距离变化全局对比图
% =========================================================
figure('Name', 'Path Loss vs Distance', 'Color', 'w', 'Position', [100, 100, 640, 480]);
hold on;

h_leg = []; 

for w = 1:num_W
    dist_axis = dist_cell{w};
    c_val = coef_c_arr(w);
    wt_name = short_wt_names{w};
    
    % Beer-Lambert 解析解 (加入几何束散衰减)
    W_geo = sqrt(w0^2 + (dist_axis .* tan(theta_half_div)).^2);
    geom_loss = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo.^2));
    PL_Beer = -10 * log10(max(exp(-c_val .* dist_axis) .* geom_loss, 1e-300));
    
    % 绘制各算法曲线
    h1 = plot(dist_axis, -d_WCI_S.PL_Cell_WCI_Strong{w}, '-', 'Color', c_Strong, ...
                   'Marker', 'o', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Strong);
    h2 = plot(dist_axis, -d_WCI_W.PL_Cell_WCI_Weak{w}, '-', 'Color', c_Weak, ...
                   'Marker', 'd', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    h3 = plot(dist_axis, -d_None.PL_Cell_WCI_None{w}, '--', 'Color', c_NoTurb, ...
                   'Marker', 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_NoTurb);
    h4 = plot(dist_axis, -d_MCS.PL_Cell_MCS{w}, '-.', 'Color', c_MCS, ...
                   'Marker', '^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_MCS);
    h5 = plot(dist_axis, PL_Beer, ':', 'Color', c_Beer, ...
                   'Marker', '*', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Beer);
    
    % 仅提取第一种水质的句柄用于图例
    if w == 1
        h_leg = [h1, h2, h3, h4, h5];
    end
    
    % 添加水质标签
    x_end = dist_axis(end);
    y_end = -d_WCI_S.PL_Cell_WCI_Strong{w}(end);
    text(x_end - 10, y_end, wt_name, 'FontSize', 12, 'FontWeight', 'normal', ...
         'Color', '#222222', 'HorizontalAlignment', 'left', 'FontName', 'Times New Roman');
end

grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Path Loss (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Path Loss Evaluation under Different Turbulence Strengths', 'FontSize', 14, 'FontName', 'Times New Roman');
xlim([5, 63]); 

leg_labels = {'WCI-MC (w/ Strong Turb)', 'WCI-MC (w/ Weak Turb)', 'WCI-MC (w/o Turb)', 'SA-MCS (w/o Turb)', 'Beer-Lambert'};
lgd = legend(h_leg, leg_labels, 'Location', 'northeast', 'FontSize', 12, 'FontName', 'Times New Roman');
title(lgd, 'Algorithms', 'FontWeight', 'normal', 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');

% =========================================================
% 附加：为纯净海水 (Clear) 绘制局部放大图 (Inset Plot)
% =========================================================
x_zoom = [50, 80];   
y_zoom = [30, 60]; 

ax_inset = axes('Position', [0.6, 0.30, 0.28, 0.25]); 
hold(ax_inset, 'on');
box(ax_inset, 'on');

w_zoom = 1;
dist_zoom = dist_cell{w_zoom};
c_val_zoom = coef_c_arr(w_zoom);

W_geo_zoom = sqrt(w0^2 + (dist_zoom .* tan(theta_half_div)).^2);
geom_loss_zoom = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo_zoom.^2));
PL_Beer_zoom = -10 * log10(max(exp(-c_val_zoom .* dist_zoom) .* geom_loss_zoom, 1e-300));

plot(ax_inset, dist_zoom, -d_WCI_S.PL_Cell_WCI_Strong{w_zoom}, '-', 'Color', c_Strong, 'Marker', 'o', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Strong);
plot(ax_inset, dist_zoom, -d_WCI_W.PL_Cell_WCI_Weak{w_zoom}, '-', 'Color', c_Weak, 'Marker', 'd', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
plot(ax_inset, dist_zoom, -d_None.PL_Cell_WCI_None{w_zoom}, '--', 'Color', c_NoTurb, 'Marker', 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_NoTurb);
plot(ax_inset, dist_zoom, -d_MCS.PL_Cell_MCS{w_zoom}, '-.', 'Color', c_MCS, 'Marker', '^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_MCS);
plot(ax_inset, dist_zoom, PL_Beer_zoom, ':', 'Color', c_Beer, 'Marker', '*', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Beer);

xlim(ax_inset, x_zoom);
ylim(ax_inset, y_zoom);
grid(ax_inset, 'on');
set(ax_inset, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 10, 'LineWidth', 1, 'FontName', 'Times New Roman');

% =========================================================
% 衰减量及衰减率计算与绘图 (仅针对 Clear Ocean)
% =========================================================
w_clear = 1; 
dist_clear = dist_cell{w_clear};

PL_Strong = -d_WCI_S.PL_Cell_WCI_Strong{w_clear};
PL_Weak   = -d_WCI_W.PL_Cell_WCI_Weak{w_clear};
PL_NoTurb = -d_None.PL_Cell_WCI_None{w_clear};
PL_MCS    = -d_MCS.PL_Cell_MCS{w_clear};

% 计算湍流引发的附加惩罚
diff_Strong_NoTurb = PL_Strong - PL_NoTurb;
diff_Weak_NoTurb   = PL_Weak - PL_NoTurb;  
diff_Strong_Weak   = PL_Strong - PL_Weak;  

W_geo_clear = sqrt(w0^2 + (dist_clear .* tan(theta_half_div)).^2);
geom_loss_clear = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo_clear.^2));
PL_Beer_clear = -10 * log10(max(exp(-coef_c_arr(w_clear) .* dist_clear) .* geom_loss_clear, 1e-300));

% 计算每 10m 的衰减增量
delta_L = diff(dist_clear);
inc10_Strong = (diff(PL_Strong) ./ delta_L) * 10;
inc10_Weak   = (diff(PL_Weak) ./ delta_L) * 10;
inc10_NoTurb = (diff(PL_NoTurb) ./ delta_L) * 10;
inc10_MCS    = (diff(PL_MCS) ./ delta_L) * 10;
inc10_Beer   = (diff(PL_Beer_clear) ./ delta_L) * 10;

% 终端数据输出
fprintf('\n========================================================================\n');
fprintf('纯净海水 (Clear) 场景下各数据点路径衰减及差值统计 (单位: dB)\n');
fprintf('========================================================================\n');
fprintf('%-8s | %-14s | %-14s | %-14s\n', '距离 (m)', '强湍流-无湍流', '弱湍流-无湍流', '强湍流-弱湍流');
fprintf('------------------------------------------------------------------------\n');
for i = 1:length(dist_clear)
    fprintf('%-8.1f | %-14.4f | %-14.4f | %-14.4f\n', ...
        dist_clear(i), diff_Strong_NoTurb(i), diff_Weak_NoTurb(i), diff_Strong_Weak(i));
end
fprintf('========================================================================\n');

fprintf('\n========================================================================================\n');
fprintf('纯净海水 (Clear) 路径衰减每10m增量统计 (单位: dB/10m)\n');
fprintf('========================================================================================\n');
fprintf('%-12s | %-10s | %-10s | %-10s | %-10s | %-10s\n', '距离段 (m)', '强湍流', '弱湍流', '无湍流', 'SA-MCS', 'Beer解析');
fprintf('----------------------------------------------------------------------------------------\n');
for i = 1:length(delta_L)
    segment_str = sprintf('%.0f-%.0f', dist_clear(i), dist_clear(i+1));
    fprintf('%-12s | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f\n', ...
        segment_str, inc10_Strong(i), inc10_Weak(i), inc10_NoTurb(i), inc10_MCS(i), inc10_Beer(i));
end
fprintf('========================================================================================\n');

% =========================================================
% 图 2: 湍流惩罚随距离变化图
% =========================================================
figure('Name', 'Turbulence Penalty', 'Color', 'w', 'Position', [150, 150, 640, 480]);
hold on;

plot(dist_clear, diff_Strong_NoTurb, '-o', 'Color', c_Strong, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_Strong);
plot(dist_clear, diff_Weak_NoTurb, '-d', 'Color', c_Weak, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
plot(dist_clear, diff_Strong_Weak, '-s', 'Color', '#7E2F8E', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', '#7E2F8E');

grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Path Loss Penalty (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Path Loss Penalty induced by Turbulence', 'FontSize', 14, 'FontName', 'Times New Roman');
legend({'Strong Turb vs No Turb (WCI-MC)', 'Weak Turb vs No Turb (WCI-MC)', 'Strong Turb vs Weak Turb'}, ...
       'Location', 'northwest', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');

% =========================================================
% 图 3: 衰减增率 (dB/10m) 柱状折线图
% =========================================================
figure('Name', 'Path Loss Increment per 10m', 'Color', 'w', 'Position', [200, 200, 640, 480]);
hold on;

L_mid = (dist_clear(1:end-1) + dist_clear(2:end)) / 2;

plot(L_mid, inc10_Strong, '-o', 'Color', c_Strong, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_Strong);
plot(L_mid, inc10_Weak, '-d', 'Color', c_Weak, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
plot(L_mid, inc10_NoTurb, '--s', 'Color', c_NoTurb, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_NoTurb);
plot(L_mid, inc10_MCS, '-.^', 'Color', c_MCS, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_MCS);
plot(L_mid, inc10_Beer, ':*', 'Color', c_Beer, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_Beer);

grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);

xtick_labels = cell(1, length(delta_L));
for i = 1:length(delta_L)
    xtick_labels{i} = sprintf('%.0f-%.0f', dist_clear(i), dist_clear(i+1));
end

set(gca, 'XTick', L_mid, 'XTickLabel', xtick_labels);
xlim([L_mid(1)-4, L_mid(end)+4]);

xlabel('Distance Segments (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Path Loss Increment (dB / 10m)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Path Loss Attenuation Rate across Segments', 'FontSize', 14, 'FontName', 'Times New Roman');
legend({'WCI-MC (Strong Turb)', 'WCI-MC (Weak Turb)', 'WCI-MC (w/o Turb)', 'SA-MCS', 'Beer-Lambert'}, ...
       'Location', 'northwest', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');