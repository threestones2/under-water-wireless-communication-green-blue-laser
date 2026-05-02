%% Exp8: WCI-MC (Strong Turbulence) 归一化 CIR 绘图程序 (原点对齐版)
% 参考文献: Boluda-Ruiz, "Impulse Response Modeling of Underwater Optical Scattering Channels..."
clc; clear; close all;

try
    % 加载强湍流仿真生成的原始能量分布数据
    d_WCI = load('data_exp8_CIR_WCIMC_Strong.mat');
catch
    error('未找到数据文件 data_exp8_CIR_WCIMC_Strong.mat，请先运行相应的仿真脚本。');
end

water_params = d_WCI.water_params;
num_water = length(water_params);
dt = d_WCI.dt;

window_size = 30; % 滑动平均窗口，用于抑制蒙特卡洛统计噪声
fig_width = 1000; 
fig_height = 450;

fprintf('--- 绘制 WCI-MC (强湍流) 归一化 CIR (散射曲线从原点起始) ---\n');

for w = 1:num_water
    wp = water_params(w);
    
    % 原始时间轴 (从 dt 开始)
    t_raw_ns = (d_WCI.Time_Axis{w} + dt) * 1e9; 
    E_scat = d_WCI.CIR_WCIMC_Scattered{w};
    E_bal  = d_WCI.E_Ballistic_Total{w};
    
    % --- 1. 执行全局归一化 (参照 Boluda-Ruiz 理论) ---
    E_total = sum(E_scat) + E_bal; 
    if E_total > 0
        % 计算归一化功率密度 h_diff(t) = E_bin / (E_total * dt)
        h_diff = E_scat / (E_total * dt);
        bal_ratio_pct = (E_bal / E_total) * 100;
    else
        h_diff = E_scat * 0;
        bal_ratio_pct = 0;
    end
    
    % --- 2. 预处理：滑动平均平滑 ---
    h_smooth = movmean(h_diff, window_size);
    
    % --- 3. 核心修正：构造从 (0,0) 起始的绘图向量 ---
    t_plot_ns = [0, t_raw_ns];     % 在时间轴最前方插入 0
    h_plot_linear = [0, h_smooth]; % 在功率密度最前方插入 0，确保曲线从原点连出
    
    % 转换为 dB 域 (仅对非零部分进行对数计算)
    h_plot_dB = 10 * log10(max(h_plot_linear, 1e-300));
    
    % 打印控制台信息
    fprintf('[%s] 直射能量占比 (Energy Fraction): %.2f%%\n', wp.name, bal_ratio_pct);
    
    % ================= 创 建 图 窗 =================
    figure('Position', [100+w*20, 100+w*20, fig_width, fig_height], 'Color', 'w');
    sgtitle(sprintf('Normalized Scattered CIR (Strong Turb.) | %s', wp.name), ...
        'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % ---------- 子图 1: 线性坐标 (Linear Scale) ----------
    subplot(1, 2, 1);
    hold on; box on; grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'FontName', 'Times New Roman');
    
    % 绘制从 (0,0) 起始的连续曲线
    plot(t_plot_ns, h_plot_linear, '-', 'Color', '#D95319', 'LineWidth', 2.5);
    
    xlabel('Time relative to LOS (ns)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Normalized CIR (s^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
    title('Linear Scale', 'FontSize', 12, 'FontName', 'Times New Roman');
    
    xlim([0, 3.0]); 
    legend(sprintf('Scattered Component\n(Ballistic: %.2f%%)', bal_ratio_pct), ...
        'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');

    % ---------- 子图 2: 对数坐标 (dB Scale) ----------
    subplot(1, 2, 2);
    hold on; box on; grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'FontName', 'Times New Roman');
    
    plot(t_plot_ns, h_plot_dB, '-', 'Color', '#0072BD', 'LineWidth', 2.5);
    
    xlabel('Time relative to LOS (ns)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Normalized CIR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Logarithmic Scale (dB)', 'FontSize', 12, 'FontName', 'Times New Roman');
    
    xlim([0, 3.0]);
    % 自动锚定散射分量峰值，并提供 50dB 的动态观测范围
    y_max = max(h_plot_dB) + 5;
    ylim([y_max - 50, y_max]); 
    
    legend(sprintf('Scattered Component\n(Ballistic: %.2f%%)', bal_ratio_pct), ...
        'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');
end