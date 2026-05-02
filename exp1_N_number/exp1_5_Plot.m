%% Exp1: Plotting Statistical Convergence (Unified with Exp2 Style - Hollow Diamond)
clc; clear; close all;

try
    d_HG  = load('data_exp1_HGPIS.mat');
    d_MCI = load('data_exp1_MCIPIS.mat');
    d_MCS = load('data_exp1_MCS.mat');
catch
    error('未能找到数据文件，请确保先运行相应的仿真脚本。');
end

N_axis = d_HG.N_arr;

figure('Name', 'Path Loss Convergence', 'Color', 'w', 'Position', [100, 100, 640, 480]);

% =========================================================
% 视觉映射定义 (严格对齐 Exp2)
% =========================================================
c_NoTurb = '#0072BD'; % 蓝   (对齐 Exp2: WCI-MC w/o Turb)
c_MCI    = '#EDB120'; % 金黄 (借用 Exp2 色系，用于区分 MCI-PIS)
c_MCS    = '#77AC30'; % 绿   (对齐 Exp2: MCS)

% 1. WCI-MC (w/o Turb) - 对应 d_HG: 蓝色，虚线，实心正方形
semilogx(N_axis, -d_HG.PL_arr, '--', 'Color', c_NoTurb, ...
    'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', 7, 'MarkerFaceColor', c_NoTurb); 
hold on;

% 2. MCI-PIS - 对应 d_MCI: 金黄色，实线，空心菱形 (修改了 MarkerFaceColor 为 'w')
semilogx(N_axis, -d_MCI.PL_arr, '-', 'Color', c_MCI, ...
    'LineWidth', 1.5, 'Marker', 'd', 'MarkerSize', 7, 'MarkerFaceColor', 'w');

% 3. MCS - 对应 d_MCS: 绿色，点划线，实心三角形
semilogx(N_axis, -d_MCS.PL_arr, '-.', 'Color', c_MCS, ...
    'LineWidth', 1.5, 'Marker', '^', 'MarkerSize', 7, 'MarkerFaceColor', c_MCS);

grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7);

% 文本设置对齐 14pt normal，并指定 Times New Roman 字体
xlabel('Number of Simulated Photons (N_{packets})', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
ylabel('Path Loss (dB)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
title('Statistical Convergence of Channel Path Loss', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Times New Roman');

xlim([N_axis(1), N_axis(end)]);

% 图例名称与 Exp2 保持术语统一，并指定字体
legend('WCI-MC (w/o Turb)', ...
       'MCI-PIS', ...
       'MCS', ...
       'Location', 'northeast', ...
       'FontSize', 12, 'FontName', 'Times New Roman');

% 坐标轴刻度对齐 12pt，并指定字体
set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');

% --- 学术图表导出指令 ---
% exportgraphics(gcf, 'Figure1_Convergence_Unified.pdf', 'ContentType', 'vector');
% exportgraphics(gcf, 'Figure1_Convergence_Unified.png', 'Resolution', 600);