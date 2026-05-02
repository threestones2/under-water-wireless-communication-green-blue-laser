% =========================================================================
% 散射相函数概率密度函数(PDF)的可视化与积分验证脚本
%
% 对比模型:
% 1. 纯瑞利散射 (Rayleigh)
% 2. 纯米氏散射 (Mie)
% 3. 纯Henyey-Greenstein (HG)
% 4. 瑞利+米氏组合散射 (Combined)
% =========================================================================
clear; clc; close all;

%% 1. 定义模型参数
% 散射系数 (来自论文Table I, λ=260nm, 单位从 km^-1 转为 m^-1)
k_s_ray = 0.266e-3; % Rayleigh散射系数 (m^-1)
k_s_mie = 0.284e-3; % Mie散射系数 (m^-1)
k_s = k_s_ray + k_s_mie; % 总散射系数

% 相函数模型参数
gamma = 0.017;      % 瑞利散射参数
g_mie = 0.72;       % 米氏散射不对称因子
f_param = 0.4;      % 米氏散射修正参数
g_hg = 0.924;       % 一个典型的强前向散射HG不对称因子

%% 2. 定义角度范围
theta = linspace(0, pi, 10001); % 散射天顶角 (rad)
mu = cos(theta); % 对应余弦值

%% 3. 计算各个散射模型的天顶角PDF f(theta)

% --- a. 纯瑞利散射PDF ---
pdf_val_rayleigh = pdf_rayleigh_func(theta, gamma);

% --- b. 纯米氏散射PDF (根据论文公式) ---
pdf_val_mie = pdf_mie_func(theta, g_mie, f_param);

% --- c. 纯Henyey-Greenstein (HG) PDF ---
P_hg = (1 - g_hg^2) / (4*pi) ./ (1 + g_hg^2 - 2*g_hg*cos(theta)).^(3/2);
pdf_val_hg = 2 * pi * P_hg .* sin(theta);

% --- d. 瑞利+米氏组合散射PDF ---
% 首先，定义单位立体角的相函数 p(μ)
p_ray_unit_solid_angle = @(mu_val) (3 * (1 + 3*gamma + (1-gamma)*mu_val.^2)) / (16*pi*(1 + 2*gamma));
p_mie_unit_solid_angle = @(mu_val) (1-g_mie^2)/(4*pi) * ( ...
    (1 ./ (1 + g_mie^2 - 2*g_mie*mu_val).^(3/2)) + ...
    f_param * (0.5*(3*mu_val.^2 - 1)) ./ ((1+g_mie^2).^(3/2)) ...
);
% 然后，计算权重
w_ray = k_s_ray / k_s;
w_mie = k_s_mie / k_s;
% 计算加权组合相函数 P(μ)
P_combined = w_ray * p_ray_unit_solid_angle(mu) + w_mie * p_mie_unit_solid_angle(mu);
% 最后，转换为天顶角PDF f(θ)
pdf_combined = 2 * pi * P_combined .* sin(theta);


%% 4. 验证所有PDF的积分是否为1
integral_rayleigh = trapz(theta, pdf_val_rayleigh);
integral_mie = trapz(theta, pdf_val_mie);
integral_hg = trapz(theta, pdf_val_hg);
integral_combined = trapz(theta, pdf_combined);

fprintf('--- 积分验证 ---\n');
fprintf('瑞利散射PDF积分值:      %.6f\n', integral_rayleigh);
fprintf('米氏散射PDF积分值:      %.6f\n', integral_mie);
fprintf('纯HG散射PDF积分值:        %.6f\n', integral_hg);
fprintf('组合散射PDF积分值:      %.6f\n', integral_combined);
fprintf('\n组合散射PDF的权重:\n');
fprintf('  - 瑞利权重 (w_ray): %.4f\n', w_ray);
fprintf('  - 米氏权重 (w_mie): %.4f\n', w_mie);


%% 5. 可视化PDF图像

% --- [新增] 让图窗在屏幕中间显示 ---
% 获取屏幕尺寸
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);
% 定义图窗尺寸
fig_width = 1000;
fig_height = 500;
% 计算图窗左下角的位置
fig_pos_x = (screen_width - fig_width) / 2;
fig_pos_y = (screen_height - fig_height) / 2;

% 创建图窗并设置位置
figure('Name', '散射相函数概率密度函数 (PDF) 对比', 'Position', [fig_pos_x, fig_pos_y, fig_width, fig_height]);
hold on;

% 绘制各个分量
plot(rad2deg(theta), pdf_val_rayleigh, 'b--', 'LineWidth', 1.5);
plot(rad2deg(theta), pdf_val_mie, 'r--', 'LineWidth', 1.5);
plot(rad2deg(theta), pdf_val_hg, 'g-.', 'LineWidth', 1.5); % [新增] 绘制HG曲线

% 绘制最终的组合函数（加粗显示）
plot(rad2deg(theta), pdf_combined, 'k-', 'LineWidth', 3);

hold off;

% --- 图像美化 ---
title('不同散射模型的角度概率密度函数 f(\theta) 对比');
xlabel('散射角 \theta (度)');
ylabel('概率密度 f(\theta) (1/rad)');
legend( ...
    sprintf('纯瑞利 (积分=%.2f)', integral_rayleigh), ...
    sprintf('纯米氏 (g=%.2f, 积分=%.2f)', g_mie, integral_mie), ...
    sprintf('纯HG (g=%.3f, 积分=%.2f)', g_hg, integral_hg), ... % [新增] HG图例
    sprintf('组合函数 (w_{ray}=%.2f, w_{mie}=%.2f, 积分=%.2f)', w_ray, w_mie, integral_combined) ...
    , 'Location', 'northeast' ...
);
grid on;
set(gca, 'FontSize', 12);
xlim([0, 180]);
% 动态调整Y轴范围，以包含所有曲线的峰值
max_y_val = max([max(pdf_val_rayleigh), max(pdf_val_mie), max(pdf_val_hg), max(pdf_combined)]);
ylim([0, max_y_val * 1.1]);


%% 辅助函数区域

% 瑞利散射PDF f(θ)
function f = pdf_rayleigh_func(theta, gamma)
    f = 3 * (1 + 3*gamma + (1-gamma)*cos(theta).^2) .* sin(theta) / (8 * (1 + 2*gamma));
end

% 米氏散射PDF f(θ) (论文版本)
function f = pdf_mie_func(theta, g, f_param)
    cos_theta = cos(theta);
    
    term1 = (1 - g^2)/2 * (1 ./ (1 + g^2 - 2*g*cos_theta).^(3/2));
    term2 = f_param * 0.5 * (3*cos_theta.^2 - 1) ./ (1 + g^2)^(3/2);
    
    f = (term1 + term2) .* sin(theta);
end
