%% 绘图代码：MCS 与 WCI-MC CIR 联合对比 (空心数据点 + 平滑曲线)
clc; clear; close all;

% ================= 数据加载 =================
file_mcs = 'data_exp8_CIR_MCS.mat';
file_wci = 'data_exp8_CIR_WCIMC.mat';

if ~exist(file_mcs, 'file') || ~exist(file_wci, 'file')
    error('缺少数据文件。请确保 data_exp8_CIR_MCS.mat 和 data_exp8_CIR_WCIMC.mat 均存在于当前目录。');
end

data_mcs = load(file_mcs);
data_wci = load(file_wci);

% 验证水质参数一致性
num_water = length(data_wci.water_params);

% ================= 绘图参数配置 =================
% 针对不同水质独立设定显示的时间上限 (单位: ns)
t_display_limits = [0.9, 0.2]; 

line_width_thin = 0.8; 
fig_width = 600; fig_height = 450;
% ================================================

fprintf('\n=========================================================================\n');
fprintf('>>> CIR 联合绘图分析 (MCS vs WCI-MC) <<<\n');
fprintf('=========================================================================\n');

for w = 1:num_water
    wp = data_wci.water_params(w);
    
    % --- 提取 MCS 数据 ---
    dt_mcs_ns = data_mcs.dt * 1e9;
    E_bal_mcs = data_mcs.E_Ballistic_Total{w};
    E_scat_mcs = data_mcs.CIR_MCS_Scattered{w};
    E_tot_mcs = E_bal_mcs + sum(E_scat_mcs);
    
    t_full_mcs = [0, (data_mcs.Time_Axis{w} + data_mcs.dt) * 1e9];
    h_full_mcs = [E_bal_mcs, E_scat_mcs] / (E_tot_mcs * dt_mcs_ns);
    
    % --- 提取 WCI-MC 数据 ---
    dt_wci_ns = data_wci.dt * 1e9;
    E_bal_wci = data_wci.E_Ballistic_Total{w};
    E_scat_wci = data_wci.CIR_WCIMC_Scattered{w};
    E_tot_wci = E_bal_wci + sum(E_scat_wci);
    
    t_full_wci = [0, (data_wci.Time_Axis{w} + data_wci.dt) * 1e9];
    h_full_wci = [E_bal_wci, E_scat_wci] / (E_tot_wci * dt_wci_ns);
    
    % --- 量化指标计算 ---
    [rms_mcs, dt40_mcs, ~, ~] = analyze_mixed_cir_40dB(t_full_mcs, [E_bal_mcs, E_scat_mcs], data_mcs.dt);
    [rms_wci, dt40_wci, ~, ~] = analyze_mixed_cir_40dB(t_full_wci, [E_bal_wci, E_scat_wci], data_wci.dt);
    
    % --- 时间窗口截断 ---
    t_limit = t_display_limits(w);
    
    mask_mcs = (t_full_mcs <= t_limit);
    t_plot_mcs = t_full_mcs(mask_mcs);
    h_plot_mcs = h_full_mcs(mask_mcs);
    
    mask_wci = (t_full_wci <= t_limit);
    t_plot_wci = t_full_wci(mask_wci);
    h_plot_wci = h_full_wci(mask_wci);
    
    % --- 构造平滑插值曲线 (pchip) ---
    t_curve_mcs = linspace(min(t_plot_mcs), max(t_plot_mcs), 1000);
    h_curve_mcs = pchip(t_plot_mcs, h_plot_mcs, t_curve_mcs);
    
    t_curve_wci = linspace(min(t_plot_wci), max(t_plot_wci), 1000);
    h_curve_wci = pchip(t_plot_wci, h_plot_wci, t_curve_wci);
    
    % ================= 绘图环节 =================
    figure('Name', sprintf('Combined CIR Data - %s', wp.name), 'Color', 'w', 'Position', [100+w*50, 100, fig_width, fig_height]);
    hold on; box on;
    
    % 1. 绘制平滑连接线 (不参与图例)
    plot(t_curve_mcs, h_curve_mcs, '-', 'Color', '#D95319', 'LineWidth', line_width_thin, 'HandleVisibility', 'off');
    plot(t_curve_wci, h_curve_wci, '-', 'Color', '#0072BD', 'LineWidth', line_width_thin, 'HandleVisibility', 'off');
    
    % 2. 绘制空心数据点
    % MCS 采用橙色空心圆
    h_m1 = plot(t_plot_mcs, h_plot_mcs, 'o', 'Color', '#D95319', ...
        'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'w', ...
        'DisplayName', 'MCS Simulation');
        
    % WCI-MC 采用蓝色空心方块
    h_m2 = plot(t_plot_wci, h_plot_wci, 's', 'Color', '#0072BD', ...
        'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'w', ...
        'DisplayName', 'WCI-MC Simulation');
    
    % 3. 指标标注框
    metrics_str = {
        sprintf('[MCS] RMS: %.4f ns, -40dB: %.4f ns', rms_mcs, dt40_mcs), ...
        sprintf('[WCI] RMS: %.4f ns, -40dB: %.4f ns', rms_wci, dt40_wci)
    };
    text(0.35, 0.70, metrics_str, 'Units', 'normalized', 'FontSize', 10, 'FontName', 'Times New Roman', ...
        'BackgroundColor', [0.98 0.98 0.98], 'EdgeColor', '#A0A0A0', 'Margin', 5);
    
    % 4. 图表格式化
    grid on; 
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'FontName', 'Times New Roman');
    xlabel('Time (ns)', 'FontSize', 12); 
    ylabel('Normalized CIR h(t)', 'FontSize', 12);
    title(sprintf('CIR Comparison: %s', wp.name), 'FontSize', 13);
    
    xlim([0, t_limit]); 
    ylim_max = max([max(h_plot_mcs), max(h_plot_wci)]);
    ylim([0, ylim_max * 1.1]);
    
    legend([h_m1, h_m2], 'Location', 'northeast', 'FontSize', 10);
    
    % 控制台输出
    fprintf('水质: %s | 显示范围: 0-%.2f ns\n', wp.name, t_limit);
    fprintf('  [MCS] RMS: %.4f ns | -40dB: %.4f ns\n', rms_mcs, dt40_mcs);
    fprintf('  [WCI] RMS: %.4f ns | -40dB: %.4f ns\n\n', rms_wci, dt40_wci);
end
fprintf('=========================================================================\n');

% --- 辅助量化分析函数 ---
function [rms_delay, dt_40dB, idx_peak, t_peak_ns] = analyze_mixed_cir_40dB(t_array_ns, E_array, dt_sec)
    total_e = sum(E_array);
    if total_e == 0
        rms_delay = 0; dt_40dB = 0; idx_peak = 1; t_peak_ns = 0; return;
    end
    
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