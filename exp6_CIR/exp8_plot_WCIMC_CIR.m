%% 绘图代码：WCI-MC CIR 原始数据美化 (科学计数法坐标轴 + 学术配色)
clc; clear; close all;

try
    load('data_exp8_CIR_WCIMC.mat'); 
catch
    error('未找到数据文件，请先运行仿真脚本生成 data_exp8_CIR_WCIMC.mat。');
end

% ================= 绘图参数配置 =================
% 针对不同水质独立设定显示的显示时间上限 (单位: ns)
t_display_limits = [0.1, 0.1]; 

% 学术配色方案
color_main = [0, 0.27, 0.53];      % 深蓝色 (Oxford Blue)
line_width_main = 1.2; 
% ================================================

num_water = length(water_params);
if length(t_display_limits) < num_water
    error('t_display_limits 数组长度必须大于等于水质种类数量。');
end

fig_width = 500; fig_height = 400;

for w = 1:num_water
    wp = water_params(w);
    dt_ns = dt * 1e9;
    t_limit = t_display_limits(w);
    
    % 1. 提取原始数据并进行全信道归一化密度计算
    E_bal = E_Ballistic_Total{w};
    E_scat_array = CIR_WCIMC_Scattered{w};
    E_total = E_bal + sum(E_scat_array);
    
    direct_ratio = E_bal / E_total;
    
    % 构造全信道时间轴与密度轴
    t_full_ns = [0, (Time_Axis{w} + dt) * 1e9];
    h_full_dens = [E_bal, E_scat_array] / (E_total * dt_ns);
    
    % 2. 量化指标计算 (-40dB 寻优内部包含平滑)
    [rms_val, dt40_val, ~, ~] = analyze_mixed_cir_40dB(t_full_ns, [E_bal, E_scat_array], dt);
    
    % 3. 筛选显示数据与曲线插值
    display_mask = (t_full_ns <= t_limit);
    t_plot = t_full_ns(display_mask);
    h_plot = h_full_dens(display_mask);
    
    t_curve = linspace(min(t_plot), max(t_plot), 1000);
    h_curve = pchip(t_plot, h_plot, t_curve); % 保形插值
    
    % 4. 绘图
    figure('Name', sprintf('CIR Beautified - %s', wp.name), 'Color', 'w', 'Position', [200+w*50, 200, fig_width, fig_height]);
    hold on; box on;
    
    % 绘制平滑连接线
    plot(t_curve, h_curve, '-', 'Color', color_main, 'LineWidth', line_width_main, 'HandleVisibility', 'off');
    
    % 绘制空心圆数据点
    plot(t_plot, h_plot, 'o', 'Color', color_main, ...
        'LineWidth', 1, 'MarkerSize', 4.5, ...
        'MarkerFaceColor', 'w', ...
        'DisplayName', 'WCI-MC Simulation');
    
    % 5. 坐标轴美化与科学计数法设置
    ax = gca;
    ax.XAxis.Exponent = 0;              % 若单位已为 ns，保持 0；若需显示 10^-9，可在此调整
    xtickformat('%.1f');                % 刻度保留一位小数或整数
    
    % 如果您希望横坐标显示为 10^-9 数量级下的整数：
    % ax.XAxis.Exponent = -9;           % 需先将 t_plot 换算回秒
    
    grid on; 
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 0.8);
    
    xlabel('Time (ns)', 'FontSize', 12, 'FontWeight', 'normal'); 
    ylabel('Normalized CIR h(t)', 'FontSize', 12, 'FontWeight', 'normal');
    title(sprintf('Impulse Response: %s', wp.name), 'FontSize', 13, 'FontWeight', 'bold');
    
    % --- 先强制刷新并设置坐标轴范围 ---
    xlim([0, t_limit]); 
    ylim([0, max(h_plot)*1.1]);
    drawnow; % 强制渲染引擎更新刻度

    % --- 后计算并应用科学计数法 ---
    ax = gca;
    current_xticks = xticks; % 获取当前实际生成的刻度 (如 0, 0.02, 0.04...)
    
    if length(current_xticks) > 1
        % 根据刻度间距计算精确的数量级 (0.02 的间距将算出 -2)
        tick_step = current_xticks(2) - current_xticks(1);
        exp_val = floor(log10(tick_step)); 
    else
        exp_val = floor(log10(t_limit));
    end
    
    % 将指数提取至坐标轴末端
    if exp_val < 0
        ax.XAxis.Exponent = exp_val;
    else
        ax.XAxis.Exponent = 0;
    end
    
    % 使用 '%g' 格式，自动丢弃多余的零，保证显示为最简整数
    ax.XAxis.TickLabelFormat = '%g';
    
    % 6. 图例与标注美化
    legend({'WCI-MC Simulation'}, 'Location', 'northeast', 'FontSize', 10);
    
    metrics_str = {
        sprintf('Direct Ratio: %.4f', direct_ratio), ...
        sprintf('RMS Delay: %.4f ns', rms_val), ...
        sprintf('-40dB Spread: %.4f ns', dt40_val)
    };
    text(0.52, 0.65, metrics_str, 'Units', 'normalized', 'FontSize', 10, 'FontName', 'Times New Roman', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.8 0.8 0.8], 'Margin', 6);
end

% --- 辅助量化分析函数 ---
function [rms_delay, dt_40dB, idx_peak, t_peak_ns] = analyze_mixed_cir_40dB(t_array_ns, E_array, dt_sec)
    total_e = sum(E_array);
    if total_e == 0, [rms_delay, dt_40dB, idx_peak, t_peak_ns] = deal(0,0,1,0); return; end
    
    tau_mean = sum(t_array_ns .* E_array) / total_e;
    rms_delay = sqrt(sum(((t_array_ns - tau_mean).^2) .* E_array) / total_e);
    
    h_dens = E_array / (total_e * dt_sec); 
    h_smooth = movmean(h_dens, 5); 
    [peak_val, idx_peak] = max(h_smooth);
    threshold_40dB = peak_val * 1e-4;
    
    idx_40dB = length(h_smooth);
    for i = idx_peak : length(h_smooth)
        if h_smooth(i) < threshold_40dB
            if (i+1 > length(h_smooth)) || (h_smooth(i+1) < threshold_40dB)
                idx_40dB = i; break;
            end
        end
    end
    t_peak_ns = t_array_ns(idx_peak);
    dt_40dB = t_array_ns(idx_40dB) - t_peak_ns;
end