%% Exp7: 逐水质独立绘图 (仅对比 强湍流 与 SA-MCS)
% 语言风格：学术、严谨、客观
clc; clear; close all;

try
    % 仅加载强湍流与 SA-MCS 数据
    d_Strong = load('data_exp7_WCIMC_Strong.mat');
    %d_Strong = load('data_exp7_WCIMC_Strong_Ensemble.mat');
    d_MCS    = load('data_exp7_MCS.mat');
catch
    error('未能找到数据文件，请确保路径正确且已运行 Strong 和 MCS 仿真脚本。');
end

water_types = d_Strong.water_types;
dist_cell = d_Strong.dist_cell;
off_axis_angles = d_Strong.off_axis_angles;

num_W = length(water_types);
num_A = length(off_axis_angles);

% --- 视觉规范定义 ---
% 颜色区分偏轴角度 (θ_err)
colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'}; 
% 线型区分算法模型
% 实线: Strong, 点划线: SA-MCS
lines  = {'-', '-.'};
markers = {'o', '^'};

for w = 1:num_W
    dist_axis = dist_cell{w};
    
    % 为每种水质创建独立画布
    figure('Name', sprintf('Exp7: Off-axis Analysis - %s', water_types{w}), ...
           'Color', 'w', 'Position', [100 + w*50, 100 + w*50, 640, 480]);
    hold on;
    
    % --- 绘制数据 ---
    for a = 1:num_A
        c_color = colors{a};
        
        % 1. 强湍流 (WCI-MC Strong)
        plot(dist_axis, -d_Strong.PL_Cell_WCI_Strong{w}(a, :), lines{1}, 'Color', c_color, ...
            'LineWidth', 1.5, 'Marker', markers{1}, 'MarkerSize', 6, 'MarkerFaceColor', c_color);
            
        % 2. SA-MCS (w/o Turb)
        plot(dist_axis, -d_MCS.PL_Cell_MCS{w}(a, :), lines{2}, 'Color', c_color, ...
            'LineWidth', 1.5, 'Marker', markers{2}, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    end
    
    % --- 生成矩阵图例 ---
    % 记录算法图例
    h_alg = zeros(1, 2);
    alg_str = {'WCI-MC (Strong)', 'SA-MCS (w/o Turb)'};
    for i = 1:2
        % 根据算法设置占位线条的标记填充状态
        if i == 1, face_col = 'k'; else, face_col = 'w'; end
        h_alg(i) = plot(NaN, NaN, lines{i}, 'Color', 'k', 'LineWidth', 1.5, ...
            'Marker', markers{i}, 'MarkerFaceColor', face_col);
    end
    
    % 记录偏轴角度图例
    h_color = zeros(1, num_A);
    angle_str = cell(1, num_A);
    for a = 1:num_A
        h_color(a) = plot(NaN, NaN, '-', 'Color', colors{a}, 'LineWidth', 2, 'Marker', 'none'); 
        angle_str{a} = sprintf('\\theta_{err} = %.1f^\\circ', off_axis_angles(a));
    end
    
    % --- 坐标轴与网格优化 ---
    grid on; 
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
    xlabel('Transmission Distance \it{L} \rm{(m)}', 'FontSize', 14, 'FontName', 'Times New Roman');
    ylabel('Path Loss (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
    title(sprintf('Pointing Error Analysis: %s', water_types{w}), ...
        'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
    
    % 组合图例
    h_combined = [h_alg, h_color];
    str_combined = [alg_str, angle_str];
    
    legend(h_combined, str_combined, ...
        'Location', 'northwest', ...    
        'FontSize', 10, ...            
        'FontName', 'Times New Roman', ...
        'NumColumns', 2);
end