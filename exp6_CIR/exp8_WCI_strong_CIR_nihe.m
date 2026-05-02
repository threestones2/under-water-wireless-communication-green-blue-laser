%% 绘图代码：WCI-MC 全信道 CIR 拟合 (双 Gamma 模型，含直射分量)
clc; clear; close all;

try
    load('data_exp8_CIR_WCIMC.mat'); 
catch
    error('未找到数据文件，请先运行仿真脚本。');
end

% ================= 拟合参数配置 =================
% 针对不同水质独立设定参与拟合的最大时间阈值 (单位: ns)
t_fit_limits = [0.5, 0.5]; 
% ================================================

num_water = length(water_params);
if length(t_fit_limits) < num_water
    error('t_fit_limits 数组长度必须大于等于水质种类数量。');
end

fig_width = 500; fig_height = 400;

fprintf('\n=========================================================================\n');
fprintf('>>> 全信道 CIR 拟合分析 (模型: Double Gamma | 含直射分量) <<<\n');
fprintf('=========================================================================\n');

for w = 1:num_water
    wp = water_params(w);
    dt_ns = dt * 1e9;
    t_fit_limit = t_fit_limits(w);
    
    % 1. 提取原始数据并进行全信道归一化密度计算
    E_bal = E_Ballistic_Total{w};
    E_scat_array = CIR_WCIMC_Scattered{w};
    E_total = E_bal + sum(E_scat_array);
    
    direct_ratio = E_bal / E_total;
    
    % 构造全信道时间轴与密度轴 (t=0 处包含直射分量)
    t_full_ns = [0, (Time_Axis{w} + dt) * 1e9];
    h_full_dens = [E_bal, E_scat_array] / (E_total * dt_ns);
    
    % 2. 量化指标计算 (RMS/-40dB)
    [rms_val, dt40_val, ~, ~] = analyze_mixed_cir_40dB(t_full_ns, [E_bal, E_scat_array], dt);
    
    % 3. 应用独立时间窗口进行拟合截断
    valid_mask = (t_full_ns <= t_fit_limit);
    t_fit_window = t_full_ns(valid_mask);
    h_fit_window = h_full_dens(valid_mask);
    
    % 4. Double Gamma 参数寻优
    sse_func = @(b) double_gamma_sse_with_penalty(b, t_fit_window, h_fit_window);
    % 提升迭代次数以应对 6 参数的搜索空间
    options = optimset('Display', 'off', 'MaxIter', 10000, 'MaxFunEvals', 10000);
    
    % 初始猜测值 [C1, a1, b1, C2, a2, b2]
    % 组1负责直射脉冲拟合，组2负责长拖尾拟合
    beta_init = [h_full_dens(1), 0.01, 10.0, max(h_full_dens)*0.1, 2.0, 5.0]; 
    beta_opt = fminsearch(sse_func, beta_init, options);
    
    C1_opt = beta_opt(1); a1_opt = beta_opt(2); b1_opt = beta_opt(3);
    C2_opt = beta_opt(4); a2_opt = beta_opt(5); b2_opt = beta_opt(6);
    
    % 5. 生成拟合平滑曲线 (外推至 1.0 ns)
    t_fine = linspace(0, 1.0, 1000);
    t_safe = max(t_fine, 1e-12); % 避免 0^0 引发数值问题
    h_theory = C1_opt .* (t_safe.^a1_opt) .* exp(-b1_opt .* t_fine) + ...
               C2_opt .* (t_safe.^a2_opt) .* exp(-b2_opt .* t_fine);
    
    % 6. 绘图
    figure('Name', sprintf('Double Gamma Fit - %s', wp.name), 'Color', 'w', 'Position', [200+w*50, 200, fig_width, fig_height]);
    hold on; box on;
    
    % 绘制完整仿真数据点 (包含 t=0 的直射峰)
    h_scat = plot(t_full_ns, h_full_dens, 'o', 'Color', '#0072BD', 'MarkerSize', 4, ...
        'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'none');
    
    % 绘制 Double Gamma 拟合曲线
    h_line = plot(t_fine, h_theory, '-', 'Color', '#D95319', 'LineWidth', 2);
    
    % 标注量化指标
    metrics_str = {
        sprintf('Direct Ratio: %.4f', direct_ratio), ...
        sprintf('RMS Delay: %.4f ns', rms_val), ...
        sprintf('-40dB Spread: %.4f ns', dt40_val)
    };
    text(0.50, 0.50, metrics_str, 'Units', 'normalized', 'FontSize', 11, 'FontName', 'Times New Roman', ...
        'BackgroundColor', [0.95 0.95 0.95], 'EdgeColor', 'k', 'Margin', 5);
    
    grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'FontName', 'Times New Roman');
    xlabel('Time (ns)'); ylabel('Normalized CIR h(t)');
    title(sprintf('Full CIR DG Modeling: %s', wp.name));
    
    xlim([0, 1.0]); 
    ylim([0, max(h_full_dens)*1.1]);
    
    legend([h_scat, h_line], {'Simulation (Ballistic+Scat)', 'Double Gamma Fit'}, 'Location', 'northeast');
    
    fprintf('水质: %s | 截断: %.2f ns | RMS: %.4f ns | -40dB: %.4f ns\n', wp.name, t_fit_limit, rms_val, dt40_val);
    fprintf('  拟合参数 -> C1: %.4e | a1: %.4f | b1: %.4f\n', C1_opt, a1_opt, b1_opt);
    fprintf('           -> C2: %.4e | a2: %.4f | b2: %.4f\n\n', C2_opt, a2_opt, b2_opt);
end
fprintf('=========================================================================\n');

% --- 辅助局部函数 ---
function sse = double_gamma_sse_with_penalty(b, t, h)
    C1 = b(1); a1 = b(2); b1 = b(3);
    C2 = b(4); a2 = b(5); b2 = b(6);
    
    % 参数边界约束：权重及衰减系数需大于0，幂次需非负
    if C1 <= 0 || a1 < 0 || b1 <= 0 || C2 <= 0 || a2 < 0 || b2 <= 0
        sse = 1e15; return;
    end
    
    t_safe = max(t, 1e-12); % 数值保护
    h_model = C1 .* (t_safe.^a1) .* exp(-b1 .* t) + C2 .* (t_safe.^a2) .* exp(-b2 .* t);
    sse = sum((h_model - h).^2);
end

function [rms_delay, dt_40dB, idx_peak, t_peak_ns] = analyze_mixed_cir_40dB(t_array_ns, E_array, dt_sec)
    total_e = sum(E_array);
    tau_mean = sum(t_array_ns .* E_array) / total_e;
    rms_delay = sqrt(sum(((t_array_ns - tau_mean).^2) .* E_array) / total_e);
    h_dens = E_array / (total_e * dt_sec); 
    h_smooth = movmean(h_dens, 5);
    [peak_val, idx_peak] = max(h_smooth);
    threshold_40dB = peak_val * 1e-4;
    idx_40dB = length(h_smooth);
    for i = idx_peak : length(h_smooth)
        if h_smooth(i) < threshold_40dB
            if (i+2 > length(h_smooth)) || (h_smooth(i+1) < threshold_40dB && h_smooth(i+2) < threshold_40dB)
                idx_40dB = i; break;
            end
        end
    end
    t_peak_ns = t_array_ns(idx_peak);
    dt_40dB = t_array_ns(idx_40dB) - t_peak_ns;
end