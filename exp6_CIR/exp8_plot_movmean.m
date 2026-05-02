%% 绘图代码 1：包含平滑处理与完整指标 (RMS Delay + 20dB Cut-off)
% 特性：无阴影上色，纯散射分量从 (0,0) 起步
clc; clear; close all;

try
    % 请根据您的实际仿真加载数据
    load('data_exp8_CIR_MCS.mat'); 
catch
    error('未找到数据文件，请先运行相应的 CIR 仿真脚本。');
end

num_water = length(water_params);
fig_width = 850; fig_height = 420;

fprintf('\n=========================================================================\n');
fprintf('>>> 脉冲响应量化分析 (平滑处理与完整指标版) <<<\n');
fprintf('=========================================================================\n');

for w = 1:num_water
    wp = water_params(w);
    
    % 散射光子的时间轴对齐到 Bin 的右边界 (dt, 2dt, 3dt...)
    t_scat_ns = (Time_Axis{w} + dt) * 1e9; 
    E_scat_array = CIR_MCS_Scattered{w};
    
    E_scat_total = sum(E_scat_array);
    E_bal = E_Ballistic_Total{w};
    E_total = E_scat_total + E_bal;
    
    Path_Loss_Total = -10 * log10(E_total + eps);
    impulsive_weight = E_bal / E_total;
    
    % --- 构造隔离版的数据对 ---
    % 全信道 (包含 t=0 的弹道能量)
    t_full_ns = [0, t_scat_ns];
    E_full    = [E_bal, E_scat_array];
    
    % 纯散射 (t=0 强制为 0，仅保留散射拖尾)
    t_scat_only_ns = [0, t_scat_ns];
    E_scat_only    = [0, E_scat_array];
    
    % --- 量化分析 ---
    [rms_full, dt20_full, p_idx_full, t_peak_full] = analyze_mixed_cir_smoothed(t_full_ns, E_full, dt);
    [rms_scat, dt20_scat, p_idx_scat, t_peak_scat] = analyze_mixed_cir_smoothed(t_scat_only_ns, E_scat_only, dt);
    
    fprintf('【 %s 】\n', wp.name);
    fprintf('  全信道 RMS Delay: %.6f ns | 20dB 展宽: %.4f ns\n', rms_full, dt20_full);
    fprintf('  纯散射 RMS Delay: %.6f ns | 20dB 展宽: %.4f ns\n\n', rms_scat, dt20_scat);
    
    % =====================================================================
    % 绘图 A：全信道 (包含 LOS 路径)
    % =====================================================================
    figure('Name', sprintf('Version A: Full Channel - %s', wp.name), 'Color', 'w', 'Position', [50, 100 , fig_width, fig_height]);
    
    h_dens_full_norm = (E_full / dt) / E_total;
    valid_mask = h_dens_full_norm > 0;
    t_plot_A = t_full_ns(valid_mask);
    h_plot_A = h_dens_full_norm(valid_mask);
    
    % [左图: 线性坐标] - 移除 area 阴影
    subplot(1, 2, 1); hold on;
    plot(t_plot_A, h_plot_A, '-', 'Color', '#D95319', 'LineWidth', 1.5);
    
    grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1, 'FontName', 'Times New Roman');
    xlabel('Time (ns)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Normalized CIR h(t)', 'FontSize', 12, 'FontName', 'Times New Roman');
    title('Version A: Full Channel (Linear)', 'FontSize', 13, 'FontName', 'Times New Roman');
    %if w == 1, xlim([0, 4.0]); else, xlim([0, 10.0]); end 
    xlim([0, 1.0]); % 【修改此处】统一将 X 轴截断在 3 ns
    
    % [右图: 对数坐标]
    subplot(1, 2, 2); hold on;
    h_dB_A = 10*log10(h_plot_A);
    plot(t_plot_A, h_dB_A, '-', 'Color', '#D95319', 'LineWidth', 1.5);
    
    peak_dB_A = 10*log10(h_dens_full_norm(p_idx_full));
    plot(t_peak_full, peak_dB_A, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    plot(t_full_ns(p_idx_full + round(dt20_full/(dt*1e9))), peak_dB_A-20, 'rs', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
    plot([min(t_plot_A), max(t_plot_A)], [peak_dB_A-20, peak_dB_A-20], 'r:', 'LineWidth', 1.5);
    
    metrics_str = {
        sprintf('Path Loss: %.2f dB', Path_Loss_Total), ...
        sprintf('Full Channel RMS: %.4f ns', rms_full), ...
        sprintf('20dB Cut-off: %.4f ns', dt20_full)
    };
    text(0.40, 0.85, metrics_str, 'Units', 'normalized', 'FontSize', 11, 'FontName', 'Times New Roman', ...
        'BackgroundColor', [0.95 0.95 0.95], 'EdgeColor', 'k', 'Margin', 5, 'LineWidth', 0.8);
    
    grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1, 'FontName', 'Times New Roman');
    xlabel('Time (ns)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('10\cdotlog_{10}[h(t)] (dB)', 'FontSize', 12, 'FontName', 'Times New Roman');
    title('Version A: Full Channel (Log)', 'FontSize', 13, 'FontName', 'Times New Roman');
    ylim([max(h_dB_A) - 50, max(h_dB_A) + 5]);
    % if w == 1, xlim([0, 4.0]); else, xlim([0, 10.0]); end
    xlim([0, 1.0]); % 【修改此处】统一将 X 轴截断在 3 ns
    legend({'Full CIR', 'Peak', '-20dB Point'}, 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');

    % =====================================================================
    % 绘图 B：纯散射 (剥离 LOS 路径)
    % =====================================================================
    figure('Name', sprintf('Version B: Dispersive Only - %s', wp.name), 'Color', 'w', 'Position', [860, 100, fig_width, fig_height]);
    
    h_dens_scat_norm = (E_scat_only / dt) / E_scat_total;
    
    % [左图: 线性坐标] - 强制从 (0,0) 绘制，移除 area 阴影
    subplot(1, 2, 1); hold on;
    h_smooth_scat = movmean(h_dens_scat_norm, 1);
    % 确保起步点被强制锁定在 (0,0)
    h_smooth_scat(1) = 0; 
    
    plot(t_scat_only_ns, h_smooth_scat, '-', 'Color', '#0072BD', 'LineWidth', 2);
    
    grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1, 'FontName', 'Times New Roman');
    xlabel('Time (ns)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Normalized Dispersive CIR h_{BP}(t)', 'FontSize', 12, 'FontName', 'Times New Roman');
    title('Version B: Dispersive Component (Linear)', 'FontSize', 13, 'FontName', 'Times New Roman');
    % if w == 1, xlim([0, 4.0]); else, xlim([0, 10.0]); end 
    xlim([0, 1.0]); % 【修改此处】统一将 X 轴截断在 3 ns
    
    % [右图: 对数坐标]
    subplot(1, 2, 2); hold on;
    valid_idx_B = h_dens_scat_norm > 0;
    h_dB_B = 10*log10(h_dens_scat_norm(valid_idx_B));
    t_plot_B = t_scat_only_ns(valid_idx_B);
    plot(t_plot_B, h_dB_B, '-', 'Color', '#0072BD', 'LineWidth', 1.5);
    
    peak_dB_B = 10*log10(h_dens_scat_norm(p_idx_scat));
    plot(t_peak_scat, peak_dB_B, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    plot(t_scat_only_ns(p_idx_scat + round(dt20_scat/(dt*1e9))), peak_dB_B-20, 'rs', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
    plot([min(t_plot_B), max(t_plot_B)], [peak_dB_B-20, peak_dB_B-20], 'r:', 'LineWidth', 1.5);
    
    metrics_str_scat = {
        sprintf('Impulsive Weight (1-k): %.3f', impulsive_weight), ...
        sprintf('Pure Scat. RMS: %.4f ns', rms_scat), ...
        sprintf('20dB Cut-off: %.4f ns', dt20_scat)
    };
    text(0.40, 0.85, metrics_str_scat, 'Units', 'normalized', 'FontSize', 11, 'FontName', 'Times New Roman', ...
        'BackgroundColor', [0.95 0.95 0.95], 'EdgeColor', 'k', 'Margin', 5, 'LineWidth', 0.8);
    
    grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1, 'FontName', 'Times New Roman');
    xlabel('Time (ns)', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('10\cdotlog_{10}[h_{BP}(t)] (dB)', 'FontSize', 12, 'FontName', 'Times New Roman');
    title('Version B: Dispersive Component (Log)', 'FontSize', 13, 'FontName', 'Times New Roman');
    ylim([max(h_dB_B) - 50, max(h_dB_B) + 5]);
    % if w == 1, xlim([0, 4.0]); else, xlim([0, 10.0]); end
    xlim([0, 1.0]); % 【修改此处】统一将 X 轴截断在 3 ns
    legend({'Norm. Dispersive CIR', 'Peak', '-20dB Point'}, 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');
end
fprintf('=========================================================================\n');

function [rms_delay, dt_20dB, idx_peak, t_peak_ns] = analyze_mixed_cir_smoothed(t_array_ns, E_array, dt_sec)
    total_e = sum(E_array);
    if total_e == 0
        rms_delay = 0; dt_20dB = 0; idx_peak = 1; t_peak_ns = 0; return;
    end
    
    tau_mean = sum(t_array_ns .* E_array) / total_e;
    rms_delay = sqrt(sum(((t_array_ns - tau_mean).^2) .* E_array) / total_e);
    
    h_dens = E_array / dt_sec; 
    % 寻找峰值与截止点时采用平滑处理
    h_smooth = movmean(h_dens, 5);
    [peak_val, idx_peak] = max(h_smooth);
    threshold = peak_val / 100; 
    
    idx_20dB = length(h_smooth);
    for i = idx_peak : length(h_smooth)
        if h_smooth(i) < threshold
            if (i+2 > length(h_smooth)) || (h_smooth(i+1) < threshold && h_smooth(i+2) < threshold)
                idx_20dB = i;
                break;
            end
        end
    end
    t_peak_ns = t_array_ns(idx_peak);
    dt_20dB = t_array_ns(idx_20dB) - t_peak_ns;
end