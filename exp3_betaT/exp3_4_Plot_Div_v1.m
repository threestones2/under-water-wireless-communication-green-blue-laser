%% Exp3.3.1: Data Visualization (Path Loss & Turbulence Penalty vs Divergence)
clc; clear; close all;

% --- 1. 加载所有产生的数据 ---
try
    d_WCI_Strong = load('data_exp3_div_WCIMC_Strong.mat');
    d_WCI_Weak   = load('data_exp3_div_WCIMC_Weak.mat');
    d_WCI_None   = load('data_exp3_div_WCIMC_None.mat'); 
    d_MCS        = load('data_exp3_div_MCS.mat');
catch
    error('数据加载失败，请务必先运行对应的仿真脚本生成 .mat 数据文件。');
end

dist_axis = d_WCI_Strong.dist_axis;
div_angles_deg = d_WCI_Strong.div_angles_deg;
num_div = length(div_angles_deg);
c_val = d_WCI_Strong.coef_c; 
w0 = 0.002; r_rx = 0.005;              

% 视觉色系规范 (固定算法映射颜色)
c_Strong = '#D95319'; c_Weak   = '#EDB120'; c_NoTurb = '#0072BD'; 
c_MCS    = '#77AC30'; c_Beer   = '#555555'; 

% 严格设置 4:3 物理尺寸参数
fig_width = 640; fig_height = 480;

% =========================================================================
% 图 1：绝对路径衰减绘图 (所有束散角合并至同一张图)
% =========================================================================
figure('Color', 'w', 'Position', [100, 100, fig_width, fig_height]);
hold on;
h_leg = []; % 用于存储图例句柄

for i = 1:num_div
    theta_half_rad = (div_angles_deg(i) / 2) * (pi / 180);
    
    % Beer-Lambert 几何衰减基准解析解
    W_geo = sqrt(w0^2 + (dist_axis .* tan(theta_half_rad)).^2);
    geom_loss = 1 - exp(-2 .* (r_rx.^2) ./ (W_geo.^2));
    PL_Beer = -10 * log10(max(exp(-c_val .* dist_axis) .* geom_loss, 1e-300));
    
    % 绘制性能曲线
    h1 = plot(dist_axis, -d_WCI_Strong.PL_Cell_WCI_Strong{i}, '-', 'Color', c_Strong, ...
                   'Marker', 'o', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Strong);
    h2 = plot(dist_axis, -d_WCI_Weak.PL_Cell_WCI_Weak{i}, '-', 'Color', c_Weak, ...
                   'Marker', 'd', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
                   
    h3 = plot(dist_axis, -d_WCI_None.PL_Cell_WCI_None{i}, '--', 'Color', c_NoTurb, ...
                   'Marker', 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_NoTurb);
                   
    h4 = plot(dist_axis, -d_MCS.PL_Cell_MCS{i}, '-.', 'Color', c_MCS, ...
                   'Marker', '^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_MCS);
                   
    h5 = plot(dist_axis, PL_Beer, ':', 'Color', c_Beer, ...
                   'Marker', '*', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', c_Beer);
    
    % 收集第一组束散角的曲线句柄，用于生成全局图例
    if i == 1
        h_leg = [h1, h2, h3, h4, h5];
    end
    
    % 使用 text 动态标注当前束散角，定位于最远距离数据的右侧
    x_text = dist_axis(end);
    y_text = -d_WCI_Strong.PL_Cell_WCI_Strong{i}(end);
    
    % 为防止文本重叠，将文本放置在点位的稍微偏右位置
    text(x_text + 1.5, y_text, sprintf('\\theta_{div} = %.2f^\\circ', div_angles_deg(i)), ...
        'FontSize', 12, 'FontName', 'Times New Roman', 'Color', 'k', 'HorizontalAlignment', 'left');
end

grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1);
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Path Loss (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Path Loss Evaluation under Different Divergence Angles', 'FontSize', 15, 'FontName', 'Times New Roman');

% 适度扩展 X 轴范围，以容纳右侧附加的文本标注
xlim([min(dist_axis), max(dist_axis) + 12]); 

% 仅使用第一组句柄生成代表算法的图例，更新标签为 WCI-MC (w/o Turb)
lgd = legend(h_leg, {'WCI-MC (w/ Strong Turb)', 'WCI-MC (w/ Weak Turb)', 'WCI-MC (w/o Turb)', 'MCS (w/o Turb)', 'Beer Law'}, ...
    'Location', 'southwest', 'FontSize', 11, 'FontName', 'Times New Roman');

% =========================================================================
% 图 2：核心机制对比：强/弱湍流附加损耗演化图
% =========================================================================
figure('Color', 'w', 'Position', [400, 200, fig_width, fig_height]);
hold on;

% 为不同的束散角分配递进色系
color_penalties = {'#7E2F8E', '#D95319', '#EDB120'};
markers = {'o', 's', '^'};

h_pen_strong = zeros(1, num_div);
h_pen_weak   = zeros(1, num_div);

for i = 1:num_div
    % 计算强/弱湍流引起的相对惩罚 (相对无湍流 WCI-MC 基准)
    penalty_strong = (-d_WCI_Strong.PL_Cell_WCI_Strong{i}) - (-d_WCI_None.PL_Cell_WCI_None{i});
    penalty_weak   = (-d_WCI_Weak.PL_Cell_WCI_Weak{i}) - (-d_WCI_None.PL_Cell_WCI_None{i});
    
    % 绘制强湍流惩罚（实线，实心标记）
    h_pen_strong(i) = plot(dist_axis, penalty_strong, '-', 'Color', color_penalties{i}, ...
        'Marker', markers{i}, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', color_penalties{i});
    
    % 绘制弱湍流惩罚（虚线，空心标记）
    h_pen_weak(i)   = plot(dist_axis, penalty_weak, '--', 'Color', color_penalties{i}, ...
        'Marker', markers{i}, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
end

grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1);
xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Turbulence Penalty (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Turbulence Penalty vs. Transmission Distance', 'FontSize', 15, 'FontName', 'Times New Roman');
xlim([min(dist_axis), max(dist_axis)]);

% 构建图例标签 (结合束散角与湍流强度)
leg_labels_pen = {};
h_leg_all = [];
for i = 1:num_div
    leg_labels_pen{end+1} = sprintf('\\theta_{div} = %.2f^\\circ (Strong)', div_angles_deg(i));
    leg_labels_pen{end+1} = sprintf('\\theta_{div} = %.2f^\\circ (Weak)', div_angles_deg(i));
    h_leg_all = [h_leg_all, h_pen_strong(i), h_pen_weak(i)];
end

% 将图例分两列显示以提升空间利用率和排版美观度
legend(h_leg_all, leg_labels_pen, 'Location', 'northwest', 'FontSize', 11, 'FontName', 'Times New Roman', 'NumColumns', 1);

disp('数据可视化完毕。所有图像均按照要求 4:3 物理比例出图，已包含强/弱双重湍流独立解耦加载。');