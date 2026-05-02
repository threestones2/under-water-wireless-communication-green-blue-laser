%% Exp7: 逐水质独立绘图 (指向误差与不同湍流算法对比)
clc; clear; close all;

try
    % 独立加载分离后的仿真数据
    d_Strong = load('data_exp7_WCIMC_Strong.mat');
    d_Weak   = load('data_exp7_WCIMC_Weak.mat');
    d_None   = load('data_exp7_WCIMC_None.mat');
    d_MCS    = load('data_exp7_MCS.mat');
catch
    error('未能找到数据文件，请确保路径正确且已运行 Strong/Weak/None/MCS 仿真脚本。');
end

water_types = d_Strong.water_types;
dist_cell = d_Strong.dist_cell;
off_axis_angles = d_Strong.off_axis_angles;

num_W = length(water_types);
num_A = length(off_axis_angles);

% 颜色区分偏轴角度 (θ_err)
colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'}; 
% 线型区分算法模型
lines  = {'-', '--', ':', '-.'};
markers = {'o', 'd', 's', '^'};

for w = 1:num_W
    dist_axis = dist_cell{w};
    
    figure('Name', sprintf('Exp7: Off-axis Analysis - %s', water_types{w}), ...
           'Color', 'w', 'Position', [100 + w*50, 100 + w*50, 640, 480]);
    hold on;
    
    % --- 绘制数据 ---
    for a = 1:num_A
        c_val = colors{a};
        
        plot(dist_axis, -d_Strong.PL_Cell_WCI_Strong{w}(a, :), lines{1}, 'Color', c_val, ...
            'LineWidth', 1.5, 'Marker', markers{1}, 'MarkerSize', 6, 'MarkerFaceColor', c_val);
            
        plot(dist_axis, -d_Weak.PL_Cell_WCI_Weak{w}(a, :), lines{2}, 'Color', c_val, ...
            'LineWidth', 1.5, 'Marker', markers{2}, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
            
        plot(dist_axis, -d_None.PL_Cell_WCI_None{w}(a, :), lines{3}, 'Color', c_val, ...
            'LineWidth', 1.5, 'Marker', markers{3}, 'MarkerSize', 6, 'MarkerFaceColor', c_val);
            
        plot(dist_axis, -d_MCS.PL_Cell_MCS{w}(a, :), lines{4}, 'Color', c_val, ...
            'LineWidth', 1.5, 'Marker', markers{4}, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    end
    
    % --- 生成矩阵图例 ---
    h_alg = zeros(1, 4);
    alg_str = {'WCI-MC (Strong)', 'WCI-MC (Weak)', 'WCI-MC (w/o Turb)', 'SA-MCS'};
    for i = 1:4
        h_alg(i) = plot(NaN, NaN, lines{i}, 'Color', 'k', 'LineWidth', 1.5, ...
            'Marker', markers{i}, 'MarkerFaceColor', 'k');
    end
    
    h_color = zeros(1, num_A);
    angle_str = cell(1, num_A);
    for a = 1:num_A
        h_color(a) = plot(NaN, NaN, '-', 'Color', colors{a}, 'LineWidth', 2, 'Marker', 'none'); 
        angle_str{a} = sprintf('\\theta_{err} = %.1f^\\circ', off_axis_angles(a));
    end
    
    grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
    xlabel('Transmission Distance \it{L} \rm{(m)}', 'FontSize', 14, 'FontName', 'Times New Roman');
    ylabel('Path Loss (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
    title(sprintf('Pointing Error Analysis: %s', water_types{w}), ...
        'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
    
    h_combined = [h_alg, h_color];
    str_combined = [alg_str, angle_str];
    
    lgd = legend(h_combined, str_combined, 'Location', 'northwest', ...
                 'FontSize', 10, 'FontName', 'Times New Roman', 'NumColumns', 2);
end