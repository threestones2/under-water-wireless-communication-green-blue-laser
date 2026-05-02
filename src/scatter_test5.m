%% Compare_All_Models_Full_Analysis.m
% 功能：
% Figure 1: 展示 Petzold 三类典型海水的 VSF 实测数据
% Figure 2: 展示相函数 P(theta) 的形状对比 (归一化为 4pi)
% Figure 3: 展示散射角 theta 的概率密度函数 PDF (考虑几何权重 sin(theta))

% Volume scattering functions for selected ocean waters petzold原始数据出处论文
% 作者: Vladimir I. Haltrin (Naval Research Laboratory)
% 机构: Scripps Institution of Oceanography (SIO), San Diego, CA.
% 编号: SIO Ref. 72-78.
% 年份: 1972

clc; clear; close all;
%% ================= [新增] Petzold 经验公式生成模块 =================
% 1. 生成角度向量 (0.1度到180度)
theta_deg_emp = logspace(log10(0.1), log10(180), 1000); 
theta_rad_emp = theta_deg_emp * pi / 180;

% 2.1 定义 Petzold 经验参数(请选择一种水质，注释掉其他)
% --- 选项 A: Clear Ocean (清洁大洋, Station 2040) ---
c_e = 0.1514; 
a_e = 0.114; 

% ==========================================================
% 3.1 petzold自动计算多项式系数
b_e = c_e - a_e;      
albedo_e = b_e / c_e; 

% Haltrin (1997) 经验参数计算
q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
k1 = 1.188 - 0.688*albedo_e;
k2 = 0.1 * (3.07 - 1.90*albedo_e);
k3 = 0.01 * (4.58 - 3.02*albedo_e);
k4 = 0.001 * (3.24 - 2.25*albedo_e);
k5 = 0.0001 * (0.84 - 0.61*albedo_e);

% 4. 计算 VSF (体散射函数)
t = theta_deg_emp; % 角度制输入
term = 1 + (-1)^1*k1*t.^0.5 + (-1)^2*k2*t.^1.0 + (-1)^3*k3*t.^1.5 + ...
       (-1)^4*k4*t.^2.0 + (-1)^5*k5*t.^2.5;
VSF_emp = exp(q_e * term);

% 5. 计算 P (相函数) 和 PDF
b_calc_emp = 2 * pi * trapz(theta_rad_emp, VSF_emp .* sin(theta_rad_emp));
P_emp = 4 * pi * VSF_emp / b_calc_emp; 
PDF_emp = 0.5 * P_emp .* sin(theta_rad_emp);

theta_deg_emp_ext = [0, theta_deg_emp];
PDF_emp_ext = [0, PDF_emp];

fprintf('[系统信息] 经验数据已生成 | 模式: c=%.3f, b=%.3f\n', c_e, b_e);

%% ================== 第一部分：加载与预处理数据 ==================
if exist('petzold_ocean.mat', 'file')
    load('petzold_ocean.mat');
    % Col 1: Angle(rad), Col 2: Harbor, Col 3: Coastal, Col 4: Clear
    theta_rad_meas = petzold_ocean(:, 1); 
    theta_deg_meas = theta_rad_meas * 180 / pi;
    
    vsf_harbor  = petzold_ocean(:, 2);
    vsf_coastal = petzold_ocean(:, 3);
    vsf_clear   = petzold_ocean(:, 4);
    
    data_loaded = true;
else
    warning('未找到 petzold_ocean.mat，使用模拟数据代替。');
    theta_deg_meas = logspace(log10(0.1), log10(180), 500)';
    theta_rad_meas = theta_deg_meas * pi / 180;
    
    % [修改]: 补充所有三种水质的模拟基准以便在缺失数据时能够顺利执行
    vsf_clear   = 0.1 ./ (theta_rad_meas + 0.01).^1.5; 
    vsf_coastal = 0.5 ./ (theta_rad_meas + 0.01).^1.6; 
    vsf_harbor  = 2.0 ./ (theta_rad_meas + 0.01).^1.7; 
    data_loaded = true;
end

%% ================== 第二部分：绘制图1 - VSF 对比 ==================
if data_loaded
    figure('Name', 'Fig1: Petzold VSF', 'Position', [50, 100, 600, 500], 'Color', 'w');
    loglog(theta_deg_meas, vsf_harbor, 's-', 'Color', [0.8 0.4 0], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Turbid Harbor'); hold on;
    loglog(theta_deg_meas, vsf_coastal, '^-', 'Color', [0 0.6 0.2], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Coastal Ocean');
    loglog(theta_deg_meas, vsf_clear,   'o-', 'Color', [0 0.4 0.8], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Clear Ocean');
    loglog(theta_deg_emp, VSF_emp, 'm--', 'LineWidth', 2, 'DisplayName', 'Empirical Formula (Ref)');
    grid on;
    xlabel('Scattering Angle \theta (deg)', 'FontWeight', 'bold');
    ylabel('VSF \beta(\theta) [m^{-1} sr^{-1}]', 'FontWeight', 'bold');
    title('Fig 1: Volume Scattering Function (VSF)');
    legend('Location', 'SouthWest');
    xlim([0.1, 180]); axis tight;
end

%% ================== 第三部分：模型定义与计算 ==================

% 1. 实测数据归一化
% 计算 Clear Ocean 的散射系数 b 并转换为相函数 P
b_meas = 2 * pi * trapz(theta_rad_meas, vsf_clear .* sin(theta_rad_meas));
P_meas = 4 * pi * vsf_clear / b_meas; 

% [新增]: 计算 Harbor 和 Coastal Ocean 的散射系数 b 并转换为相函数 P
b_harbor = 2 * pi * trapz(theta_rad_meas, vsf_harbor .* sin(theta_rad_meas));
P_harbor = 4 * pi * vsf_harbor / b_harbor;

b_coastal = 2 * pi * trapz(theta_rad_meas, vsf_coastal .* sin(theta_rad_meas));
P_coastal = 4 * pi * vsf_coastal / b_coastal;

% 2. 模型函数定义
P_HG = @(mu, g) (1 - g^2) ./ (1 + g^2 - 2*g*mu).^(1.5);
P_Ray_Gen = @(mu, gamma) (3 ./ (4*(1+2*gamma))) .* ( (1+3*gamma) + (1-gamma).*mu.^2 );
P_GHG = @(mu, g, f) (1 - g^2) .* ( ...
    1 ./ (1 + g^2 - 2*g*mu).^(1.5) + ...
    f .* 0.5 .* (3*mu.^2 - 1) ./ (1 + g^2).^(1.5) ...
);

Model_Ding = @(theta, w_mie, gamma, g_mie, f_mie) ...
    (1 - w_mie) * P_Ray_Gen(cos(theta), gamma) + ...
    w_mie * P_GHG(cos(theta), g_mie, f_mie);

Model_TTHG = @(theta, alpha, g1, g2) ...
    alpha * P_HG(cos(theta), g1) + (1 - alpha) * P_HG(cos(theta), g2);

Model_SingleHG = @(theta, g) P_HG(cos(theta), g);

% 3. 参数设置
w_mie_ding = 0.8440;
gamma_ding = 0.017;
g_ding     = 0.9814;
f_ding     = 0.49;  

g_single   = 0.9707;

alpha_tthg = 0.4437;
g1_tthg    = 0.9900;
g2_tthg    = 0.8232;

% 4. 计算绘图数据
theta_deg_plot = logspace(log10(0.01), log10(180), 1000); 
theta_rad_plot = theta_deg_plot * pi / 180;

P_fit_Ding   = Model_Ding(theta_rad_plot, w_mie_ding, gamma_ding, g_ding, f_ding);
P_fit_TTHG   = Model_TTHG(theta_rad_plot, alpha_tthg, g1_tthg, g2_tthg);
P_fit_Single = Model_SingleHG(theta_rad_plot, g_single);

%% ================= FF 模型最小二乘法自动拟合 =================
fprintf('------------------------------------------------------\n');
fprintf('[自动拟合] 正在以 Haltrin 经验公式为参考目标优化 FF 模型参数...\n');

valid_mask = (P_emp > 0) & isfinite(P_emp);
theta_fit_rad = theta_rad_emp(valid_mask);
P_ref_log = log10(P_emp(valid_mask)); 

calc_FF_optimized = @(t, n, mu) FF_Core_Function(t, n, mu);
cost_func = @(params) calculate_cost(params, theta_fit_rad, P_ref_log, calc_FF_optimized);

initial_guesses = [
    1.05, 3.65;
    1.10, 3.80;
    1.15, 4.20;
    1.02, 3.20
];

options = optimset('Display', 'off', 'TolX', 1e-5, 'TolFun', 1e-5, 'MaxFunEvals', 1000);
best_overall_params = initial_guesses(1, :);
min_overall_error = inf;

for i = 1:size(initial_guesses, 1)
    [temp_params, temp_error] = fminsearch(cost_func, initial_guesses(i, :), options);
    if temp_error < min_overall_error
        min_overall_error = temp_error;
        best_overall_params = temp_params;
    end
end

n_best = best_overall_params(1);
mu_best = best_overall_params(2);
fprintf('[优化结果] 全局最佳折射率 n = %.4f\n', n_best);
fprintf('[优化结果] 全局最佳斜率 mu  = %.4f\n', mu_best);
fprintf('[优化结果] 最小对数均方误差 = %.4f\n', min_overall_error);
fprintf('------------------------------------------------------\n');

P_fit_FF = calc_FF_optimized(theta_rad_plot, n_best, mu_best);

%% ================== 第四部分：绘制图2 - 相函数 P(theta) ==================
figure('Name', 'Fig2: Phase Function', 'Position', [660, 100, 600, 500], 'Color', 'w');

% [修改]: 增加 Harbor 与 Coastal 数据散点的绘制
loglog(theta_deg_meas, P_harbor, 's', 'Color', [0.8 0.4 0], 'MarkerSize', 4, 'DisplayName', 'Meas (Harbor)'); hold on;
loglog(theta_deg_meas, P_coastal, '^', 'Color', [0 0.6 0.2], 'MarkerSize', 4, 'DisplayName', 'Meas (Coastal)');
loglog(theta_deg_meas, P_meas, 'ko', 'MarkerSize', 4, 'DisplayName', 'Meas (Clear)');

%loglog(theta_deg_plot, P_fit_Ding, 'r-', 'LineWidth', 2, 'DisplayName', 'Ding Mix');
loglog(theta_deg_plot, P_fit_TTHG, 'b--', 'LineWidth', 2, 'DisplayName', 'TTHG');
loglog(theta_deg_plot, P_fit_Single, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Single HG');
loglog(theta_deg_emp, P_emp, 'm--', 'LineWidth', 2, 'DisplayName', 'Empirical Formula');
loglog(theta_deg_plot, P_fit_FF, 'c-', 'LineWidth', 2, 'DisplayName', 'Fournier-Forand'); 

grid on;
xlabel('Scattering Angle \theta (deg)', 'FontWeight', 'bold');
ylabel('Phase Function P(\theta) [sr^{-1}]', 'FontWeight', 'bold');
title('Fig 2: Phase Function P(\theta) (Normalized to 4\pi)');
legend('Location', 'SouthWest');
xlim([0.1, 180]); ylim([1e-4, 1e6]);

%% ================== 第五部分：绘制图3 - 散射角 PDF f(theta) ==================
PDF_meas   = 0.5 * P_meas .* sin(theta_rad_meas);
PDF_Ding   = 0.5 * P_fit_Ding .* sin(theta_rad_plot);
PDF_TTHG   = 0.5 * P_fit_TTHG .* sin(theta_rad_plot);
PDF_Single = 0.5 * P_fit_Single .* sin(theta_rad_plot);
PDF_FF     = 0.5 * P_fit_FF .* sin(theta_rad_plot);

figure('Name', 'Fig3: Scattering Angle PDF', 'Position', [360, 50, 600, 500], 'Color', 'w');
semilogx(theta_deg_meas, PDF_meas, 'ko', 'MarkerSize', 4, 'DisplayName', 'Meas (Clear)'); hold on;
semilogx(theta_deg_plot, PDF_Ding, 'r-', 'LineWidth', 2, 'DisplayName', 'Ding Mix');
semilogx(theta_deg_plot, PDF_TTHG, 'b--', 'LineWidth', 2, 'DisplayName', 'TTHG');
semilogx(theta_deg_plot, PDF_Single, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Single HG(g=0.9707)');
semilogx(theta_deg_emp_ext, PDF_emp_ext, 'm--', 'LineWidth', 2, 'DisplayName', 'Empirical Formula (Ref)');
semilogx(theta_deg_plot, PDF_FF, 'c-', 'LineWidth', 2, 'DisplayName', 'Fournier-Forand');
grid on;
xlabel('Scattering Angle \theta (deg)', 'FontWeight', 'bold');
ylabel('Probability Density f(\theta) [rad^{-1}]', 'FontWeight', 'bold');
title({'Fig 3: PDF of Plane Scattering Angle \theta', 'f(\theta) = 0.5 \cdot P(\theta) \cdot sin(\theta)'});
legend('Location', 'NorthEast');
xlim([0.1, 180]); 

integral_check = trapz(theta_rad_plot, PDF_Ding);
fprintf('Ding 模型 PDF 积分检查 (应接近1): %.4f\n', integral_check);

%% === 辅助局部函数 ===
function P_norm = FF_Core_Function(theta, n, mu)
    theta = max(theta, 1e-9); 
    nu = (3 - mu) / 2;
    delta = (4 ./ (3 * (n - 1)^2)) .* sin(theta/2).^2;
    delta(abs(delta - 1) < 1e-7) = 1.000001; 
    delta_180 = (4 / (3 * (n - 1)^2)); 
    
    term1 = nu .* (1 - delta) - (1 - delta.^nu) + ...
            (delta .* (1 - delta.^nu) - nu .* (1 - delta)) .* (sin(theta/2).^(-2));
    part1 = (1 ./ (4 * pi .* (1 - delta).^2 .* delta.^nu)) .* term1;
    part2 = ((1 - delta_180^nu) / (16 * pi * (delta_180 - 1) * delta_180^nu)) .* (3 .* cos(theta).^2 - 1);
    
    beta_tilde = max(1e-12, part1 + part2);
    P_norm = 4 * pi * beta_tilde; 
end

function cost = calculate_cost(params, theta_meas, P_meas_log, model_func)
    n = params(1);
    mu = params(2);
    
    penalty = 0;
    if n <= 1.01; penalty = penalty + (1.01 - n) * 1000; n = 1.01; end
    if n >= 1.35; penalty = penalty + (n - 1.35) * 1000; n = 1.35; end
    if mu <= 2.5; penalty = penalty + (2.5 - mu) * 1000; mu = 2.5; end
    if mu >= 5.5; penalty = penalty + (mu - 5.5) * 1000; mu = 5.5; end
    
    P_model = model_func(theta_meas, n, mu);
    P_model_log = log10(max(P_model, 1e-12));
    
    mse = mean((P_meas_log(:) - P_model_log(:)).^2);
    cost = mse + penalty;
    
    if ~isfinite(cost); cost = 1e6; end
end