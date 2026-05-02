%% Exp7: Plotting MCS Pointing Error Simulation Only
% 功能: 独立加载并绘制 Exp7 中 MCS 算法的偏轴误差路径损耗曲线，方便参数调试
clc; clear; close all;

try
    % 尝试加载 MCS 仿真产生的数据文件
    load('data_exp7_MCS.mat');
catch
    error('未能找到 data_exp7_MCS.mat，请确保已成功运行 exp7_2_MCS.m。');
end

% --- 补充仿真中硬编码的物理元数据 ---
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
off_axis_angles = [0, 1, 3]; % 偏轴角度 (度)
num_W = length(water_types);
num_A = length(off_axis_angles);

% --- 绘图视觉规范 ---
% 采用标准学术配色与不同的标记区分多组偏轴误差
colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'};
markers = {'o', 's', '^', 'd'};

fprintf('\n=========================================================================\n');
fprintf('>>> MCS 偏轴误差物理基准独立绘图 <<<\n');
fprintf('=========================================================================\n');

for w = 1:num_W
    dist_axis = dist_cell{w};
    PL_matrix = PL_Cell_MCS{w}; % 提取当前水质下的路径损耗矩阵 (尺寸: num_A x num_D)
    
    figure('Name', sprintf('MCS Pointing Error - %s', water_types{w}), ...
           'Color', 'w', 'Position', [30 + w*50, 30 + w*50, 700, 500]);
    hold on;
    
    % 预分配图例句柄数组
    h_leg = gobjects(1, num_A);
    leg_str = cell(1, num_A);
    
    for a = 1:num_A
        c_val = colors{a};
        m_val = markers{a};
        
        % 提取损耗数据。仿真中计算的 PL 为负值，绘图时取相反数转换为正向衰减 dB
        h_leg(a) = plot(dist_axis, -PL_matrix(a, :), '-.', 'Color', c_val, ...
            'Marker', m_val, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
            
        leg_str{a} = sprintf('\\theta_{err} = %.1f^\\circ', off_axis_angles(a));
    end
    
    % --- 图表坐标轴与网格排版 ---
    grid on; 
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
    xlabel('Transmission Distance L (m)', 'FontSize', 14, 'FontName', 'Times New Roman');
    ylabel('Path Loss (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
    title(sprintf('MCS Path Loss under Pointing Errors (%s)', water_types{w}), 'FontSize', 15, 'FontName', 'Times New Roman');
    
    % --- 图例配置 ---
    lgd = legend(h_leg, leg_str, 'Location', 'bestoutside', 'FontSize', 12, 'FontName', 'Times New Roman');
    title(lgd, 'Off-axis Angle');
    
    set(gca, 'FontSize', 12, 'LineWidth', 1, 'FontName', 'Times New Roman');
    xlim([min(dist_axis)-2, max(dist_axis)+2]);
end