%% 绘图代码：最大通信速率联合分析与归一化 CIR 波形抽取 (640x480 版)
clc; clear; close all;

try
    % 加载多距离 WCI-MC CIR 仿真数据
    load('data_exp8_CIR_Dist.mat'); 
catch
    error('未找到 data_exp8_CIR_Dist.mat，请先运行仿真脚本。');
end

% ==================== 参数定制配置区 ====================
% 1. CIR 波形抽取距离设置 (单位: m)
sel_distances = [50, 20]; 

% 2. CIR 波形显示的时间截断上限 (单位: ns)
time_display_limits = [0.1, 0.1]; 

% 3. 图像尺寸固定为 640x480
fig_size = [640, 480];

% 4. 视觉配置
colors = [0, 0.27, 0.53;   % 深蓝色 (Oxford Blue)
          0.85, 0.33, 0.1]; % 砖红色 (Brick Red)
line_width_main = 1.2; 
% ========================================================

num_water = length(water_params);

fprintf('\n=========================================================================\n');
fprintf('>>> 实验数据绘图：通信速率联合分析与归一化 CIR 波形抽取 <<<\\n');
fprintf('=========================================================================\n');

%% 任务 1: 计算最大通信速率并在同一张图中联合绘制 (右上角图例)
figure('Name', 'Combined Max Bit Rate', 'Color', 'w', 'Position', [100, 200, fig_size(1), fig_size(2)]);
hold on; box on;

h_plots = zeros(1, num_water); 
labels = cell(1, num_water);

for w = 1:num_water
    wp = water_params(w);
    L_arr = wp.L_arr;
    num_d = length(L_arr);
    rate_arr = zeros(1, num_d);
    
    for d_idx = 1:num_d
        E_bal = E_Ballistic_Total{w, d_idx};
        E_scat = CIR_WCIMC_Scattered{w, d_idx};
        t_full_ns = [0, (Time_Axis{w, d_idx} + dt) * 1e9];
        E_full = [E_bal, E_scat];
        
        % 严格基于原始能量计算 RMS Delay
        total_e = sum(E_full);
        if total_e > 0
            tau_mean = sum(t_full_ns .* E_full) / total_e;
            rms_val = sqrt(sum(((t_full_ns - tau_mean).^2) .* E_full) / total_e);
        else
            rms_val = 0;
        end
        
        % 计算最大码速率 R <= 1 / (10 * D_rms)
        if rms_val > 0
            rate_arr(d_idx) = 1 / (10 * rms_val * 1e-9) / 1e6;
        else
            rate_arr(d_idx) = NaN; 
        end
    end
    
    % 绘制速率-距离曲线
    h_plots(w) = plot(L_arr, rate_arr, '-o', 'Color', colors(w,:), ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    labels{w} = wp.name;
end

% 图表 1 格式化
grid on; 
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('Transmission Distance (m)', 'FontSize', 12);
ylabel('Maximum Bit Rate (Mbps)', 'FontSize', 12);
title('Maximum Bit Rate vs Distance', 'FontSize', 13, 'FontWeight', 'bold');
legend(h_plots, labels, 'Location', 'northeast', 'FontSize', 10);
xlim([0, max(cellfun(@max, {water_params.L_arr})) + 5]); 

%% 任务 2: 抽取选定距离的原始 CIR 波形并进行归一化绘图
for w = 1:num_water
    wp = water_params(w);
    t_limit = time_display_limits(w);
    
    [~, d_idx] = min(abs(wp.L_arr - sel_distances(w)));
    actual_L = wp.L_arr(d_idx);
    
    E_bal = E_Ballistic_Total{w, d_idx};
    E_scat = CIR_WCIMC_Scattered{w, d_idx};
    E_tot = E_bal + sum(E_scat);
    direct_ratio = E_bal / E_tot;
    
    t_full_ns = [0, (Time_Axis{w, d_idx} + dt) * 1e9];
    h_full_dens = [E_bal, E_scat] / (E_tot * dt * 1e9);
    
    % 按照时间截断上限筛选数据
    mask = (t_full_ns <= t_limit);
    t_plot = t_full_ns(mask);
    h_raw = h_full_dens(mask);
    
    % --- 核心：进行峰值归一化 ---
    if ~isempty(h_raw) && max(h_raw) > 0
        h_plot = h_raw / max(h_raw);
    else
        h_plot = h_raw;
    end
    
    % 计算物理指标 (基于原始能量，不随绘图归一化改变)
    total_e_tmp = sum([E_bal, E_scat]);
    tau_mean_tmp = sum(t_full_ns .* [E_bal, E_scat]) / total_e_tmp;
    rms_val = sqrt(sum(((t_full_ns - tau_mean_tmp).^2) .* [E_bal, E_scat]) / total_e_tmp);
    
    % 视觉连接用的保形插值曲线
    t_curve = linspace(min(t_plot), max(t_plot), 1000);
    h_curve = pchip(t_plot, h_plot, t_curve);
    
    % 绘制独立 CIR 波形图
    figure('Name', sprintf('Normalized CIR - %s (%dm)', wp.name, actual_L), ...
        'Color', 'w', 'Position', [150+w*50, 150, fig_size(1), fig_size(2)]);
    hold on; box on;
    
    plot(t_curve, h_curve, '-', 'Color', colors(w,:), 'LineWidth', line_width_main, 'HandleVisibility', 'off');
    plot(t_plot, h_plot, 'o', 'Color', colors(w,:), ...
        'LineWidth', 1, 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', ...
        'DisplayName', 'WCI-MC Simulation');
    
    % 坐标轴设置
    xlim([0, t_limit]); 
    ylim([0, 1.1]); % 归一化后 y 轴锁定在 1.1 左右即可
    drawnow;
    
    ax = gca;
    current_xticks = xticks; 
    if length(current_xticks) > 1
        exp_val = floor(log10(current_xticks(2) - current_xticks(1))); 
    else
        exp_val = floor(log10(t_limit));
    end
    ax.XAxis.Exponent = min(exp_val, 0);
    ax.XAxis.TickLabelFormat = '%g'; 
    
    grid on; set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontSize', 11, 'FontName', 'Times New Roman');
    xlabel('Time (ns)', 'FontSize', 12); 
    ylabel('Normalized CIR h(t)/max[h(t)]', 'FontSize', 12);
    title(sprintf('Impulse Response: %s (%dm)', wp.name, actual_L), 'FontSize', 13, 'FontWeight', 'bold');
    
    legend({'WCI-MC Simulation'}, 'Location', 'northeast', 'FontSize', 10);
    
    metrics_str = {
        sprintf('Distance: %d m', actual_L), ...
        sprintf('Direct Ratio: %.4f', direct_ratio), ...
        sprintf('RMS Delay: %.4f ns', rms_val)
    };
    text(0.52, 0.55, metrics_str, 'Units', 'normalized', 'FontSize', 10, 'FontName', 'Times New Roman', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.8 0.8 0.8], 'Margin', 6);
end