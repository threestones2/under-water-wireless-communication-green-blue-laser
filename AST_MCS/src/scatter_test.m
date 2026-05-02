%% Check_Petzold_Clear.m
% 脚本功能：检查 Petzold Clear Water 的相函数数据及生成逻辑
% 目的：验证角度制/弧度制混用是否导致 PDF 异常

clc; clear; close all;

% ==========================================
% 1. 直接加载原始数据 (.mat) 进行检查
% ==========================================
fprintf('1. Loading raw data from petzold_ocean.mat...\n');
try
    load('petzold_ocean.mat'); % 确保文件在路径中
    % petzold_ocean 列定义推测: [Angle, Harbor, Coastal, Clear]
    raw_angle = petzold_ocean(:,1); 
    raw_vsf_clear = petzold_ocean(:,4);
    
    figure('Name', 'Petzold Clear Water Analysis');
    
    % --- 子图 1: 原始 VSF (相函数) ---
    subplot(2,2,1);
    semilogy(raw_angle, raw_vsf_clear, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Angle (Raw Value)'); 
    ylabel('Volume Scattering Function (VSF)');
    title('Raw Data: VSF vs Angle');
    legend('干净海水');
    
    % --- 子图 2: 尝试计算 PDF (假设 Angle 是角度制) ---
    % PDF \propto VSF * sin(theta)
    % 如果 raw_angle 是角度，必须用 sind()；如果是弧度，用 sin()
    pdf_deg = raw_vsf_clear .* sin(raw_angle); % （如果是角度）
    %pdf_rad_err = raw_vsf_clear .* sin(raw_angle); % 错误做法（如果是弧度函数）
    
    subplot(2,2,2);
    plot(raw_angle, pdf_deg, 'b-', 'LineWidth', 1); hold on;
    %plot(raw_angle, pdf_rad_err, 'r--', 'LineWidth', 1);
    grid on;
    xlabel('Angle (Raw Value)'); ylabel('Unnormalized PDF');
    title('PDF Check (Integrand)');
    legend( 'Buggy (using sin)', 'Location', 'Best');
    
    fprintf('   -> 观察子图2：红色虚线是否在震荡或为负值？如果是，说明 sin(角度) 导致了错误。\n');
    
catch ME
    fprintf('Error loading .mat file: %s\n', ME.message);
end

% ==========================================
% 2. 调用 generate_scatter 函数查看输出
% ==========================================
fprintf('\n2. Calling generate_scatter(''measured'', ''petzold_clear'')...\n');

% 注意：如果代码未修复，这里的输出可能已经包含错误
[cdf_out, angle_out] = generate_scatter('measured', 'petzold_clear');

% --- 子图 3: 输出的 CDF ---
subplot(2,2,3);
plot(angle_out, cdf_out, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Output Angle (Check Unit!)'); ylabel('CDF');
title('Output CDF from function');

% --- 子图 4: 从 CDF 反推 PDF ---
% PDF = d(CDF) / d(Angle)
d_cdf = diff(cdf_out);
d_ang = diff(angle_out);
pdf_derived = d_cdf ./ d_ang;
mid_angle = angle_out(1:end-1) + d_ang/2;

subplot(2,2,4);
plot(mid_angle, pdf_derived, 'm-', 'LineWidth', 1.5);
grid on;
xlabel('Output Angle'); ylabel('Derived PDF');
title('PDF derived from Output CDF');

fprintf('检查完成。请查看生成的图表。\n');