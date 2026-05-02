%% Exp1: Plotting Statistical Convergence (Bar Chart / Updated Order & Names)
clc; clear; close all;

try
    d_HG  = load('data_exp1_HGPIS.mat');
    d_MCI = load('data_exp1_MCIPIS.mat');
    d_MCS = load('data_exp1_MCS.mat');
catch
    error('未能找到数据文件，请确保先运行相应的仿真脚本。');
end

N_axis = d_HG.N_arr;

figure('Name', 'Path Loss Convergence (Bar)', 'Color', 'w', 'Position', [100, 100, 640, 480]);

% =========================================================
% 视觉映射定义
% =========================================================
c_MCI    = '#EDB120'; % 金黄 (MCI-PIS)
c_NoTurb = '#0072BD'; % 蓝   (WCI-MC w/o Turb)
c_MCS    = '#77AC30'; % 绿   (SA-MCS)

% 组合数据矩阵 (行：N_packets, 列：不同算法)
% 数据顺序已对调：MCI, HG(WCI-MC), SA-MCS
PL_matrix = [-d_MCI.PL_arr(:), -d_HG.PL_arr(:), -d_MCS.PL_arr(:)];

% 绘制分组柱状图
b = bar(1:length(N_axis), PL_matrix, 'grouped', 'BarWidth', 0.8);

% 单独设置每组柱子的颜色（需匹配上面矩阵的排列顺序）
b(1).FaceColor = c_MCI;
b(2).FaceColor = c_NoTurb;
b(3).FaceColor = c_MCS;

hold on;
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7);

% 为了更清晰地展示微小的收敛差异，限制 y 轴范围 (建议保留)
ylim([31, 33.5]); 

% 文本设置
xlabel('Number of Simulated Photons (N_{packets})', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
ylabel('Path Loss (dB)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
title('Statistical Convergence of Channel Path Loss', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');

% 格式化 X 轴刻度标签为 10^n 形式
xtick_labels = arrayfun(@(n) sprintf('10^{%d}', log10(n)), N_axis, 'UniformOutput', false);
set(gca, 'XTick', 1:length(N_axis), 'XTickLabel', xtick_labels);

% 图例名称与字体指定 (注意顺序已与数据矩阵对应，并且更新为 SA-MCS)
legend('MCI-PIS', ...
       'WCI-MC (w/o Turb)', ...
       'SA-MCS', ...
       'Location', 'northeast', ...
       'FontSize', 12, 'FontName', 'Times New Roman');

% 坐标轴整体格式
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');

% --- 学术图表导出指令 ---
% exportgraphics(gcf, 'Figure1_Convergence_Bar.pdf', 'ContentType', 'vector');
% exportgraphics(gcf, 'Figure1_Convergence_Bar.png', 'Resolution', 600);