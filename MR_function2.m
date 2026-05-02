% =========================================================================
% 散射相函数PDF可视化与修正脚本
%
% 重点：处理并展示米氏散射模型中的非物理负值问题
% =========================================================================
clear; clc; close all;

%% 1. 定义模型参数 (与之前相同)
k_s_ray = 0.266e-3;
k_s_mie = 0.284e-3;
k_s = k_s_ray + k_s_mie;
gamma = 0.017;
g_mie = 0.72;
f_param = 0.5;
g_hg = 0.924;

%% 2. 定义角度范围
theta = linspace(0, pi, 10001);
mu = cos(theta);

%% 3. 计算原始的PDF（可能包含负值）

% --- 米氏散射PDF (原始版本) ---
pdf_mie_original = pdf_mie_func(theta, g_mie, f_param);

% 检查是否存在负值
min_pdf_mie = min(pdf_mie_original);
if min_pdf_mie < 0
    fprintf('!!! 警告: 原始米氏散射PDF中检测到负值。\n');
    fprintf('    最小值为: %.4e\n', min_pdf_mie);
end

%% 4. [新增] 修正米氏散射PDF (钳位 + 重新归一化)
pdf_mie_corrected = pdf_mie_original;
pdf_mie_corrected(pdf_mie_corrected < 0) = 0; % 步骤1: 将所有负值设为0

% 步骤2: 重新归一化
integral_after_clamp = trapz(theta, pdf_mie_corrected);
pdf_mie_corrected = pdf_mie_corrected / integral_after_clamp;

%% 5. 重新计算组合函数 (使用修正后的米氏PDF)
w_ray = k_s_ray / k_s;
w_mie = k_s_mie / k_s;

% 瑞利PDF (不变)
pdf_rayleigh = pdf_rayleigh_func(theta, gamma);

% [新增] 使用修正后的米氏PDF来组合
pdf_combined_corrected = w_ray * pdf_rayleigh + w_mie * pdf_mie_corrected;


%% 6. 验证积分
integral_mie_original = trapz(theta, pdf_mie_original);
integral_mie_corrected = trapz(theta, pdf_mie_corrected);
integral_combined_corrected = trapz(theta, pdf_combined_corrected);

fprintf('\n--- 积分验证 ---\n');
fprintf('原始米氏PDF积分值:      %.6f (由于负值存在，可能不为1)\n', integral_mie_original);
fprintf('修正后米氏PDF积分值:    %.6f (应为1)\n', integral_mie_corrected);
fprintf('修正后组合PDF积分值:    %.6f (应为1)\n', integral_combined_corrected);


%% 7. 可视化对比

% --- 图1: 对比原始和修正后的米氏PDF ---
figure('Name', '米氏散射PDF的修正', 'Position', [50, 100, 1000, 600]);
hold on;
plot(rad2deg(theta), pdf_mie_original, 'r--', 'LineWidth', 1.5);
plot(rad2deg(theta), pdf_mie_corrected, 'r-', 'LineWidth', 0.5);
plot([0 180], [0 0], 'k-'); % 画一条y=0的黑线
hold off;
title('米氏散射PDF的负值问题与修正');
xlabel('散射角 \theta (度)');
ylabel('概率密度 f(\theta) (1/rad)');
legend( ...
    sprintf('原始模型 (积分=%.4f)', integral_mie_original), ...
    sprintf('修正后模型 (钳位+归一化, 积分=%.4f)', integral_mie_corrected)...
);
grid on;
set(gca, 'FontSize', 12);
xlim([0, 180]);

% --- 图2: 最终的组合函数对比图 (使用修正后的版本) ---
figure('Name', '最终组合散射PDF','Position', [100, 100, 1000, 600]);
hold on;
plot(rad2deg(theta), pdf_rayleigh, 'b--', 'LineWidth', 1.5);
plot(rad2deg(theta), pdf_mie_corrected, 'r--', 'LineWidth', 1.5); % 使用修正后的米氏
plot(rad2deg(theta), pdf_combined_corrected, 'k-', 'LineWidth', 3); % 使用修正后的组合
hold off;
title('修正后的组合散射角度概率密度函数 f(\theta)');
xlabel('散射角 \theta (度)');
ylabel('概率密度 f(\theta) (1/rad)');
legend( ...
    '纯瑞利', ...
    '纯米氏 (已修正)', ...
    sprintf('最终组合函数 (积分=%.2f)', integral_combined_corrected),'Location','best'...
);
grid on;
set(gca, 'FontSize', 12);
xlim([0, 180]);
ylim([0, max(pdf_mie_corrected)*1.1]);


%% 辅助函数区域 (保持不变)
function f = pdf_rayleigh_func(theta, gamma)
    f = 3 * (1 + 3*gamma + (1-gamma)*cos(theta).^2) .* sin(theta) / (8 * (1 + 2*gamma));
end

function f = pdf_mie_func(theta, g, f_param)
    cos_theta = cos(theta);
    term1 = (1 - g^2)/2 * (1 ./ (1 + g^2 - 2*g*cos_theta).^(3/2));
    term2 = f_param * (1 - g^2)/2*0.5 * (3*cos_theta.^2 - 1) ./ (1 + g^2)^(3/2);
    f = (term1 + term2) .* sin(theta);
end
