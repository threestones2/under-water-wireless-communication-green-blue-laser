%% Plot_Exp8_FOV_Comparison.m
% 实验目的: 可视化并对比 WCI-MC 与 MCS 算法在不同 FOV 和水质下的路径损耗
clc; clear; close all;

% ================= 1. 加载仿真数据 =================
file_wci = 'data_exp8_WCIMC_Strong_FOV.mat';
file_mcs = 'data_exp8_MCS_FOV.mat';

if ~exist(file_wci, 'file') || ~exist(file_mcs, 'file')
    error('未找到数据文件，请确保 Exp8 的两份仿真代码已成功运行并保存数据。');
end

data_wci = load(file_wci);
data_mcs = load(file_mcs);

water_types = data_wci.water_types;
dist_cell = data_wci.dist_cell;
fov_angles = data_wci.fov_angles_deg;

num_W = length(water_types);
num_F = length(fov_angles);

% ================= 2. 绘图参数设置 =================
% 为不同的 FOV 分配不同的颜色
colors = lines(num_F); 
% 定义两种算法的线型和标记
% WCI-MC: 实线 + 实心标记
% MCS: 虚线 + 空心标记
markers_wci = {'o', 's', '^', 'd', 'v'};
markers_mcs = {'o', 's', '^', 'd', 'v'};

% ================= 3. 循环生成独立图表 =================
for w_idx = 1:num_W
    figure('Name', sprintf('FOV Comparison - %s', water_types{w_idx}), ...
           'Position', [100+w_idx*50, 100+w_idx*50, 700, 500]);
    hold on; box on; grid on;
    
    dist_arr = dist_cell{w_idx};
    
    % 提取当前水质下的路径损耗矩阵 (转换为正值的 Path Loss)
    % 注意：此前代码保存的 PL_matrix 为 10*log10(P_rx)，通常 Path Loss 取其相反数
    PL_WCI = -data_wci.PL_Cell_WCI_Strong_FOV{w_idx};
    PL_MCS = -data_mcs.PL_Cell_MCS_FOV{w_idx};
    
    h_plots = zeros(1, num_F * 2); % 用于存储句柄以生成图例
    legend_str = cell(1, num_F * 2);
    
    plot_idx = 1;
    for f_idx = 1:num_F
        % 绘制 WCI-MC 结果
        h_plots(plot_idx) = plot(dist_arr, PL_WCI(f_idx, :), ...
            'LineStyle', '-', 'LineWidth', 1.5, ...
            'Color', colors(f_idx, :), ...
            'Marker', markers_wci{f_idx}, 'MarkerFaceColor', colors(f_idx, :), 'MarkerSize', 6);
        legend_str{plot_idx} = sprintf('WCI-MC (FOV=%.0f\\circ)', fov_angles(f_idx));
        plot_idx = plot_idx + 1;
        
        % 绘制 MCS 结果
        h_plots(plot_idx) = plot(dist_arr, PL_MCS(f_idx, :), ...
            'LineStyle', '--', 'LineWidth', 1.5, ...
            'Color', colors(f_idx, :), ...
            'Marker', markers_mcs{f_idx}, 'MarkerFaceColor', 'none', 'MarkerSize', 6);
        legend_str{plot_idx} = sprintf('MCS (FOV=%.0f\\circ)', fov_angles(f_idx));
        plot_idx = plot_idx + 1;
    end
    
    % ================= 4. 图表格式化 =================
    title(sprintf('Path Loss vs Distance (%s)', water_types{w_idx}), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Transmission Distance (m)', 'FontSize', 11);
    ylabel('Path Loss (dB)', 'FontSize', 11);
    
    % 根据水质调整坐标轴范围 (可选)
    xlim([min(dist_arr) max(dist_arr)]);
    
    % 生成图例
    legend(h_plots, legend_str, 'Location', 'best', 'FontSize', 10, 'NumColumns', 2);
    
    % 提升整体视觉质量
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end

fprintf('绘图完成，共生成 3 张对比图表。\n');