% =========================================================================
% 散射相函数概率密度函数(PDF)的可视化与积分验证脚本
% =========================================================================
clear; clc; close all;

%% 1. 定义模型参数 (与论文一致)
gamma = 0.017;  % 瑞利散射参数
g_mie = 0.72;   % 米氏散射不对称因子
f_param = 0.5;  % 米氏散射修正参数
g_hg = 0.924;   % 一个典型的强前向散射HG不对称因子

%% 2. 定义角度范围
theta = linspace(0, pi, 10001); % 使用足够多的点以保证积分精度

%% 3. 计算三种散射模型的天顶角PDF f(theta)

% --- 瑞利散射PDF ---
pdf_val_rayleigh = pdf_rayleigh(theta, gamma);

% --- 米氏散射PDF (根据论文公式) ---
pdf_val_mie = pdf_mie(theta, g_mie, f_param);

% --- 纯Henyey-Greenstein (HG) PDF ---
% f(theta) = 2*pi * P(theta) * sin(theta)
P_hg = (1 - g_hg^2) / (4*pi) ./ (1 + g_hg^2 - 2*g_hg*cos(theta)).^(3/2);
pdf_val_hg = 2 * pi * P_hg .* sin(theta);


%% 4. 验证积分是否为1
% 使用梯形法则(trapz)对f(theta)在[0, pi]上积分
integral_rayleigh = trapz(theta, pdf_val_rayleigh);
integral_mie = trapz(theta, pdf_val_mie);
integral_hg = trapz(theta, pdf_val_hg);

fprintf('瑞利散射PDF积分值: %.6f\n', integral_rayleigh);
fprintf('米氏散射PDF积分值: %.6f\n', integral_mie);
fprintf('纯HG散射PDF积分值: %.6f\n', integral_hg);


%% 5. 可视化PDF图像
figure('Name', '散射相函数概率密度函数 (PDF)', 'Position', [100, 100, 1000, 600]);
hold on;

plot(rad2deg(theta), pdf_val_rayleigh, 'b-', 'LineWidth', 2);
plot(rad2deg(theta), pdf_val_mie, 'r-', 'LineWidth', 2);
plot(rad2deg(theta), pdf_val_hg, 'g--', 'LineWidth', 2);

hold off;

% --- 图像美化 ---
title('不同散射模型的角度概率密度函数 f(\theta)');
xlabel('散射角 \theta (度)');
ylabel('概率密度 f(\theta) (1/rad)');
legend( ...
    sprintf('瑞利散射 (积分=%.4f)', integral_rayleigh), ...
    sprintf('米氏散射 (g=%.2f, f=%.1f, 积分=%.4f)', g_mie, f_param, integral_mie), ...
    sprintf('纯HG散射 (g=%.3f, 积分=%.4f)', g_hg, integral_hg) ...
);
grid on;
set(gca, 'FontSize', 12);
xlim([0, 180]);
ylim([0, inf]); % 自动调整y轴


%% 辅助函数区域

% 瑞利散射PDF
function f = pdf_rayleigh(theta, gamma)
    f = 3 * (1 + 3*gamma + (1-gamma)*cos(theta).^2) .* sin(theta) / (8 * (1 + 2*gamma));
end

% 米氏散射PDF (论文版本)
function f = pdf_mie(theta, g, f_param)
    % 避免在theta=0时出现NaN
    cos_theta = cos(theta);
    cos_theta(1) = 1 - 1e-12; % 微小扰动

    term1 = (1 - g^2)/2 * (1 ./ (1 + g^2 - 2*g*cos_theta).^(3/2));
    term2 = f_param * 0.5 * (3*cos_theta.^2 - 1) ./ (1 + g^2)^(3/2);
    
    f = (term1 + term2) .* sin(theta);
end
