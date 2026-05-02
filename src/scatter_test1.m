%% Compare_Models_Full_Range.m
% 功能：对比 Petzold 实测数据 vs 三种散射模型 (全角度 0.1°-180° 最佳拟合)
% 数据源：Petzold Clear Ocean (Empirical)
% 优化方法：最小二乘法 (Log-MSE)

clc; clear; close all;

%% 1. 加载数据
if exist('petzold_ocean.mat', 'file')
    load('petzold_ocean.mat');
    theta_rad_meas = petzold_ocean(:, 1); 
    vsf_meas = petzold_ocean(:, 4); 
    theta_deg_meas = theta_rad_meas * 180 / pi;
else
    % 模拟数据 (Petzold Empirical Model)
    theta_deg_meas = logspace(log10(0.1), log10(180), 500)';
    theta_rad_meas = theta_deg_meas * pi / 180;
    
    % 使用经验公式生成 VSF
    c=0.1514; a=0.114; b=c-a; albedo=b/c;
    q = 2.598 + 17.748*sqrt(b) - 16.722*b + 5.932*b*sqrt(b);
    k1=1.188-0.688*albedo; k2=0.1*(3.07-1.90*albedo);
    k3=0.01*(4.58-3.02*albedo); k4=0.001*(3.24-2.25*albedo);
    k5=0.0001*(0.84-0.61*albedo);
    t = theta_deg_meas;
    term = 1 + (-1)^1*k1*t.^0.5 + (-1)^2*k2*t.^1.0 + (-1)^3*k3*t.^1.5 + ...
           (-1)^4*k4*t.^2.0 + (-1)^5*k5*t.^2.5;
    vsf_meas = exp(q * term);
end

%% 2. 归一化计算 P_meas
b_meas = 2 * pi * trapz(theta_rad_meas, vsf_meas .* sin(theta_rad_meas));
P_meas = 4 * pi * vsf_meas / b_meas;
fprintf('实测数据散射系数 b = %.4f m^-1\n', b_meas);

%% 3. 定义模型函数
% (A) Generalized Rayleigh
P_Ray_Gen = @(mu, gamma) (3 ./ (4*(1+2*gamma))) .* ( (1+3*gamma) + (1-gamma).*mu.^2 );

% (B) Standard HG
P_HG = @(mu, g) (1 - g^2) ./ (1 + g^2 - 2*g*mu).^(1.5);

% (C) Generalized HG (Ding et al. Eq.6)
P_GHG = @(mu, g, f) (1 - g^2) .* ( ...
    1 ./ (1 + g^2 - 2*g*mu).^(1.5) + ...
    f .* 0.5 .* (3*mu.^2 - 1) ./ (1 + g^2).^(1.5) ...
);

% --- 模型组装 ---
Model_Ding = @(theta, w_mie, gamma, g_mie, f_mie) ...
    (1 - w_mie) * P_Ray_Gen(cos(theta), gamma) + ...
    w_mie * P_GHG(cos(theta), g_mie, f_mie);

Model_TTHG = @(theta, alpha, g1, g2) ...
    alpha * P_HG(cos(theta), g1) + (1 - alpha) * P_HG(cos(theta), g2);

Model_SingleHG = @(theta, g) P_HG(cos(theta), g);

%% 4. 参数设置 (全角度 0.1°-180° 最佳拟合结果)

% --- Ding et al. 混合模型 ---
w_mie_ding = 0.8820;  % Mie 权重
gamma_ding = 0.0000;  % 瑞利参数 (接近0)
g_ding     = 0.9801;  % 极强前向峰
f_ding     = 0.7846;  % 侧向/背向修正

% --- TTHG 双项模型 ---
alpha_tthg = 0.4435;  
g1_tthg    = 0.9900;  
g2_tthg    = 0.8238;

% --- Single HG 单项模型 ---
g_single   = 0.9707; 

%% 5. 作图
theta_deg_plot = logspace(log10(0.1), log10(180), 1000);
theta_rad_plot = theta_deg_plot * pi / 180;

P_Ding   = Model_Ding(theta_rad_plot, w_mie_ding, gamma_ding, g_ding, f_ding);
P_TTHG   = Model_TTHG(theta_rad_plot, alpha_tthg, g1_tthg, g2_tthg);
P_Single = Model_SingleHG(theta_rad_plot, g_single);

figure('Position', [200, 100, 900, 600], 'Color', 'w');

% 绘制实测点
loglog(theta_deg_meas, P_meas, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'DisplayName', 'Petzold 实测'); 
hold on;

% 绘制模型曲线
loglog(theta_deg_plot, P_Ding, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Ding混合模型');
loglog(theta_deg_plot, P_TTHG, 'b--', 'LineWidth', 2, 'DisplayName', 'TTHG 双项模型');
loglog(theta_deg_plot, P_Single, 'g-.', 'LineWidth', 1.5, 'DisplayName', sprintf('Single HG (g=%.4f)', g_single));

% 装饰
grid on;
xlabel('散射角 \theta (度)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('相位函数 P(\theta) [sr^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
title('水下散射相函数全角度模型对比 (0.1^{\circ}-180^{\circ})', 'FontSize', 14);

% 标注参数
dim = [0.15 0.15 0.35 0.25];
str = {['\color{red}Ding: w_{mie}=' num2str(w_mie_ding,'%.3f') ', g=' num2str(g_ding,'%.3f') ', f=' num2str(f_ding,'%.3f')], ...
       ['\color{blue}TTHG: \alpha=' num2str(alpha_tthg,'%.3f') ', g_1=' num2str(g1_tthg,'%.3f') ', g_2=' num2str(g2_tthg,'%.3f')], ...
       ['\color[rgb]{0,0.5,0}Single: g=' num2str(g_single,'%.4f')]};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 10, 'EdgeColor', 'k');

legend('Location', 'SouthWest', 'FontSize', 11);
xlim([0.1, 180]); 
ylim([1e-4, 1e6]);

fprintf('=== 优化结果 (全角度) ===\n');
fprintf('Ding Model: w=%.4f, gamma=%.4f, g=%.4f, f=%.4f\n', w_mie_ding, gamma_ding, g_ding, f_ding);
fprintf('TTHG Model: alpha=%.4f, g1=%.4f, g2=%.4f\n', alpha_tthg, g1_tthg, g2_tthg);