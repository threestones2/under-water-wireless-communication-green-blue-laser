%% Exp2: Plotting Macroscopic Channel Characteristics (Raw Value + Log Scale)
clc; clear; close all;
        
try
    d_WCI  = load('data_exp2_WCIMC.mat');
    d_None = load('data_exp2_WCIMC_None.mat'); 
    d_MCS  = load('data_exp2_MCS.mat');
catch
    error('未能找到数据文件，请确保先运行 WCIMC, WCIMC_None 及 MCS 仿真脚本。');
end

dist_cell = d_WCI.dist_cell;
num_W = length(dist_cell);

coef_c_arr = [0.1514, 0.398, 2.190]; 
w0 = 0.002;                          
theta_half_div = 0.1 * pi / 180 / 2; 
Rx_Aperture = 0.05;                  
r_rx = Rx_Aperture / 2;              

c_Strong = '#D95319'; 
c_Weak   = '#EDB120'; 
c_NoTurb = '#0072BD'; 
c_MCS    = '#77AC30'; 
c_Beer   = '#555555'; 

short_wt_names = {'Clear', 'Coastal', 'Harbor'}; 

% =========================================================================
% 图 1：接收功率衰减绘图 (真实大小，对数坐标轴)
% =========================================================================
figure('Name', 'Received Power vs Distance', 'Color', 'w', 'Position', [100, 100, 640, 480]);
hold on;

h_leg = []; 

for w = 1:num_W
    dist_axis = dist_cell{w};
    c_val = coef_c_arr(w);
    wt_name = short_wt_names{w};
    
    W_geo = sqrt(w0^2 + (dist_axis .* tan(theta_half_div)).^2);
    geom_loss = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo.^2));
    
    % 从 dB 还原为原始大小(真实接收能量比例)
    P_Beer   = max(exp(-c_val .* dist_axis) .* geom_loss, 1e-300);
    P_Strong = 10.^(d_WCI.PL_Cell_WCI_Strong{w} / 10);
    P_Weak   = 10.^(d_WCI.PL_Cell_WCI_Weak{w} / 10);
    P_NoTurb = 10.^(d_None.PL_Cell_WCI_None{w} / 10);
    P_MCS    = 10.^(d_MCS.PL_Cell_MCS{w} / 10);
    
    h1 = plot(dist_axis, P_Strong, '-', 'Color', c_Strong, ...
                   'Marker', 'o', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Strong);
    h2 = plot(dist_axis, P_Weak, '-', 'Color', c_Weak, ...
                   'Marker', 'd', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    h3 = plot(dist_axis, P_NoTurb, '--', 'Color', c_NoTurb, ...
                   'Marker', 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_NoTurb);
    h4 = plot(dist_axis, P_MCS, '-.', 'Color', c_MCS, ...
                   'Marker', '^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_MCS);
    h5 = plot(dist_axis, P_Beer, ':', 'Color', c_Beer, ...
                   'Marker', '*', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Beer);
    
    if w == 1
        h_leg = [h1, h2, h3, h4, h5];
    end
    
    x_end = dist_axis(end);
    y_end = P_Strong(end);
    text(x_end - 10, y_end, wt_name, 'FontSize', 12, 'FontWeight', 'normal', ...
         'Color', '#222222', 'HorizontalAlignment', 'left', 'FontName', 'Times New Roman');
end

grid on; 
% 核心设置：将 Y 轴设置为对数刻度 (Log Scale)
set(gca, 'YScale', 'log', 'GridLineStyle', ':', 'GridAlpha', 0.6);
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
ylabel('Normalized Received Power', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
title('Received Power Evaluation under Different Turbulence Strengths', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
xlim([5, 63]); 
leg_labels = {'WCI-MC (w/ Strong)', 'WCI-MC (w/ Weak)', 'WCI-MC (w/o)', 'SA-MCS', 'Beer Law'};
lgd = legend(h_leg, leg_labels, 'Location', 'northeast', 'FontSize', 12, 'FontName', 'Times New Roman');
title(lgd, 'Algorithms', 'FontWeight', 'normal', 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');

%%
% =========================================================
% 附加：为纯净海水 (Clear) 绘制局部放大图 (Inset Plot)
% =========================================================
x_zoom = [50, 80];   
% 针对线性数值的 Y 轴放缩范围 (依据之前 30~60dB 映射约为 1e-3 到 1e-6)
y_zoom = [1e-6, 1e-3]; 

ax_inset = axes('Position', [0.6, 0.30, 0.28, 0.25]); 
hold(ax_inset, 'on');
box(ax_inset, 'on');

w_zoom = 1;
dist_zoom = dist_cell{w_zoom};
c_val_zoom = coef_c_arr(w_zoom);

W_geo_zoom = sqrt(w0^2 + (dist_zoom .* tan(theta_half_div)).^2);
geom_loss_zoom = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo_zoom.^2));
P_Beer_zoom = max(exp(-c_val_zoom .* dist_zoom) .* geom_loss_zoom, 1e-300);

P_Strong_zoom = 10.^(d_WCI.PL_Cell_WCI_Strong{w_zoom} / 10);
P_Weak_zoom   = 10.^(d_WCI.PL_Cell_WCI_Weak{w_zoom} / 10);
P_NoTurb_zoom = 10.^(d_None.PL_Cell_WCI_None{w_zoom} / 10);
P_MCS_zoom    = 10.^(d_MCS.PL_Cell_MCS{w_zoom} / 10);

plot(ax_inset, dist_zoom, P_Strong_zoom, '-', 'Color', c_Strong, 'Marker', 'o', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Strong);
plot(ax_inset, dist_zoom, P_Weak_zoom, '-', 'Color', c_Weak, 'Marker', 'd', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
plot(ax_inset, dist_zoom, P_NoTurb_zoom, '--', 'Color', c_NoTurb, 'Marker', 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_NoTurb);
plot(ax_inset, dist_zoom, P_MCS_zoom, '-.', 'Color', c_MCS, 'Marker', '^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_MCS);
plot(ax_inset, dist_zoom, P_Beer_zoom, ':', 'Color', c_Beer, 'Marker', '*', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Beer);

xlim(ax_inset, x_zoom);
ylim(ax_inset, y_zoom);
grid(ax_inset, 'on');
% 同样为局部放大图设置 Log Scale
set(ax_inset, 'YScale', 'log', 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 10, 'LineWidth', 1, 'FontName', 'Times New Roman');

%%
% =========================================================
% 衰减量及衰减率计算与绘图 (真实比例)
% =========================================================
w_clear = 1; 
dist_clear = dist_cell{w_clear};

P_Strong = 10.^(d_WCI.PL_Cell_WCI_Strong{w_clear} / 10);
P_Weak   = 10.^(d_WCI.PL_Cell_WCI_Weak{w_clear} / 10);
P_NoTurb = 10.^(d_None.PL_Cell_WCI_None{w_clear} / 10);
P_MCS    = 10.^(d_MCS.PL_Cell_MCS{w_clear} / 10);

W_geo_clear = sqrt(w0^2 + (dist_clear .* tan(theta_half_div)).^2);
geom_loss_clear = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo_clear.^2));
P_Beer_clear = max(exp(-coef_c_arr(w_clear) .* dist_clear) .* geom_loss_clear, 1e-300);

% 计算比例 (Ratio)，替代原先的减法 dB 差值
Ratio_Strong_NoTurb = P_Strong ./ P_NoTurb;
Ratio_Weak_NoTurb   = P_Weak ./ P_NoTurb;  
Ratio_Strong_Weak   = P_Strong ./ P_Weak;  

delta_L = diff(dist_clear);
% 计算每 10m 的能量保留比例 (Retention Ratio per 10m)
retention10_Strong = (P_Strong(2:end) ./ P_Strong(1:end-1)) .^ (10 ./ delta_L);
retention10_Weak   = (P_Weak(2:end) ./ P_Weak(1:end-1)) .^ (10 ./ delta_L);
retention10_NoTurb = (P_NoTurb(2:end) ./ P_NoTurb(1:end-1)) .^ (10 ./ delta_L);
retention10_MCS    = (P_MCS(2:end) ./ P_MCS(1:end-1)) .^ (10 ./ delta_L);
retention10_Beer   = (P_Beer_clear(2:end) ./ P_Beer_clear(1:end-1)) .^ (10 ./ delta_L);

% 更新控制台输出
fprintf('\n========================================================================\n');
fprintf('纯净海水 (Clear) 场景下各数据点接收功率比例 (Turb / Ref)\n');
fprintf('========================================================================\n');
fprintf('%-8s | %-14s | %-14s | %-14s\n', '距离 (m)', '强湍流/无湍流', '弱湍流/无湍流', '强湍流/弱湍流');
fprintf('------------------------------------------------------------------------\n');
for i = 1:length(dist_clear)
    fprintf('%-8.1f | %-14.4f | %-14.4f | %-14.4f\n', ...
        dist_clear(i), Ratio_Strong_NoTurb(i), Ratio_Weak_NoTurb(i), Ratio_Strong_Weak(i));
end
fprintf('========================================================================\n');

fprintf('\n========================================================================================\n');
fprintf('纯净海水 (Clear) 接收功率每 10m 衰减保留率\n');
fprintf('========================================================================================\n');
fprintf('%-12s | %-10s | %-10s | %-10s | %-10s | %-10s\n', '距离段 (m)', '强湍流', '弱湍流', '无湍流', 'SA-MCS', 'Beer解析');
fprintf('----------------------------------------------------------------------------------------\n');
for i = 1:length(delta_L)
    segment_str = sprintf('%.0f-%.0f', dist_clear(i), dist_clear(i+1));
    fprintf('%-12s | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f\n', ...
        segment_str, retention10_Strong(i), retention10_Weak(i), retention10_NoTurb(i), retention10_MCS(i), retention10_Beer(i));
end
fprintf('========================================================================================\n');

% 图 2：接收能量比例图 (替换 Penalty)
figure('Name', 'Turbulence Ratio', 'Color', 'w', 'Position', [150, 150, 640, 480]);
hold on;

plot(dist_clear, Ratio_Strong_NoTurb, '-o', 'Color', c_Strong, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_Strong);
plot(dist_clear, Ratio_Weak_NoTurb, '-d', 'Color', c_Weak, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
plot(dist_clear, Ratio_Strong_Weak, '-s', 'Color', '#7E2F8E', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', '#7E2F8E');

grid on; set(gca, 'YScale', 'log', 'GridLineStyle', ':', 'GridAlpha', 0.6);
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
ylabel('Received Power Ratio', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
title('Received Power Ratio induced by Turbulence', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
legend({'Strong Turb / No Turb', 'Weak Turb / No Turb', 'Strong Turb / Weak Turb'}, ...
       'Location', 'southwest', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');

% 图 3：每 10m 能量保留率
figure('Name', 'Power Retention Ratio per 10m', 'Color', 'w', 'Position', [200, 200, 640, 480]);
hold on;

L_mid = (dist_clear(1:end-1) + dist_clear(2:end)) / 2;

plot(L_mid, retention10_Strong, '-o', 'Color', c_Strong, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_Strong);
plot(L_mid, retention10_Weak, '-d', 'Color', c_Weak, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
plot(L_mid, retention10_NoTurb, '--s', 'Color', c_NoTurb, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_NoTurb);
plot(L_mid, retention10_MCS, '-.^', 'Color', c_MCS, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_MCS);
plot(L_mid, retention10_Beer, ':*', 'Color', c_Beer, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', c_Beer);

grid on; set(gca, 'YScale', 'log', 'GridLineStyle', ':', 'GridAlpha', 0.6);

xtick_labels = cell(1, length(delta_L));
for i = 1:length(delta_L)
    xtick_labels{i} = sprintf('%.0f-%.0f', dist_clear(i), dist_clear(i+1));
end

set(gca, 'XTick', L_mid, 'XTickLabel', xtick_labels);
xlim([L_mid(1)-4, L_mid(end)+4]);

xlabel('Distance Segments (m)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
ylabel('Power Retention Ratio (per 10m)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
title('Power Retention Rate across Segments', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
legend({'WCI-MC (Strong)', 'WCI-MC (Weak)', 'WCI-MC (w/o Turb)', 'SA-MCS', 'Beer-Lambert'}, ...
       'Location', 'southwest', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');