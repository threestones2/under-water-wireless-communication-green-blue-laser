%% Compare_All_Models_Full_Analysis.m
% 功能：
% Figure 1: 展示 Petzold 三类典型海水的 VSF 实测数据
% Figure 2: 展示相函数 P(theta) 的形状对比 (归一化为 4pi)
% Figure 3: 展示散射角 theta 的概率密度函数 PDF (考虑几何权重 sin(theta))

% Volume scattering functions for selected ocean waters petzold原始数据出处论文
% 作者: Vladimir I. Haltrin (美国海军研究实验室, Naval Research Laboratory)
% 机构: Scripps Institution of Oceanography (SIO), San Diego, CA.
% 编号: SIO Ref. 72-78.
% 年份: 1972

% 论文标题: "Theoretical and empirical phase functions for Monte Carlo calculations of light scattering in seawater"
% 发表出处: Proceedings of SPIE, Vol. 2963, Ocean Optics XIII (1997).
% 拟合公式出处

clc; clear; close all;
%% ================= [新增] Petzold 经验公式生成模块 =================
% 说明：基于 Mobley/Haltrin 经验公式生成 Petzold 标准水质数据
% 目的：作为标准参考曲线 (Ref)，不依赖外部 .mat 文件

% 1. 生成角度向量 (0.1度到180度)
theta_deg_emp = logspace(log10(0.1), log10(180), 1000); 
theta_rad_emp = theta_deg_emp * pi / 180;

% 2.1 定义 Petzold 经验参数(请选择一种水质，注释掉其他)
% petzold在514nm下的实测值的经验拟合

% --- 选项 A: Clear Ocean (清洁大洋, Station 2040) ---
% 特征：衰减小，前向峰值极尖锐，适合验证深海通信
c_e = 0.1514; 
a_e = 0.114; 

% --- 选项 B: Coastal Ocean (近岸海水) ---
% 特征：散射增强，浑浊度中等
% c_e = 0.398;
% a_e = 0.179;

% --- 选项 C: Turbid Harbor (浑浊港口) ---
% 特征：极度浑浊，散射系数 b 很大 (>1.8)，通信距离极短
% c_e = 2.190;
% a_e = 0.366;

% 2.2 利用 Solonenko (2015) 的 Jerlov 参数 + Haltrin (1999) 的经验公式
% 目的：生成特定 Jerlov 水质在 532nm 下的参考相函数

% 2.2.1. 定义 Jerlov 水质参数 (532 nm, 插值结果)
% 请根据需要取消注释某一行

% --- Oceanic Types (大洋水) ---
% Jerlov_Type = 'I';   c_j = 0.0524; b_j = 0.0020; % 极纯净
% Jerlov_Type = 'IA';  c_j = 0.0550; b_j = 0.0039;
% Jerlov_Type = 'IB';  c_j = 0.1069; b_j = 0.0557; % 普通大洋
% Jerlov_Type = 'II';  c_j = 0.4301; b_j = 0.3789;
% Jerlov_Type = 'III'; c_j = 1.0922; b_j = 1.0376;

% --- Coastal Types (近岸水) ---
% Jerlov_Type = '1C';  c_j = 0.4568; b_j = 0.3866;
% Jerlov_Type = '5C';  c_j = 1.5334; b_j = 1.4092; % 典型浑浊近岸
% Jerlov_Type = '9C';  c_j = 3.7142; b_j = 3.3072; % 极浑浊

% ==========================================================
% % 3.1 petzold自动计算多项式系数 (无需修改)
% b_e = c_e - a_e;      % 计算散射系数 b
% albedo_e = b_e / c_e; % 计算单次散射反照率 w0

% 3.2 jerlov经验拟合
b_e = b_j;      
c_e = c_j;
a_e=c_e-b_e;
albedo_e = b_e / c_e; % 单次散射反照率 omega_0

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
% 归一化积分
b_calc_emp = 2 * pi * trapz(theta_rad_emp, VSF_emp .* sin(theta_rad_emp));
P_emp = 4 * pi * VSF_emp / b_calc_emp; % 相函数 P(theta)

% 计算 PDF (用于图3)
PDF_emp = 0.5 * P_emp .* sin(theta_rad_emp);

% 补全 0 点 (优化绘图效果)
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
    % 简单模拟 Clear Ocean 数据
    vsf_clear = 0.1 ./ (theta_rad_meas + 0.01).^1.5; 
    vsf_harbor = []; vsf_coastal = [];
    data_loaded = false;
end

%% ================== 第二部分：绘制图1 - VSF 对比 ==================
if data_loaded
    figure('Name', 'Fig1: Petzold VSF', 'Position', [50, 100, 600, 500], 'Color', 'w');
    loglog(theta_deg_meas, vsf_harbor, 's-', 'Color', [0.8 0.4 0], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Turbid Harbor');
    hold on;
    loglog(theta_deg_meas, vsf_coastal, '^-', 'Color', [0 0.6 0.2], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Coastal Ocean');
    loglog(theta_deg_meas, vsf_clear,   'o-', 'Color', [0 0.4 0.8], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Clear Ocean');
    hold on;
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
% 计算散射系数 b (用于将 VSF 转为 P)
b_meas = 2 * pi * trapz(theta_rad_meas, vsf_clear .* sin(theta_rad_meas));
P_meas = 4 * pi * vsf_clear / b_meas; % 相函数 P(theta)

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

% 3. 参数设置 (Clear Ocean 0-180度拟合参数)
% Ding et al. (0.1-180度 修正)
w_mie_ding = 0.8758;
gamma_ding = 0.0000;
g_ding     = 0.9803;
f_ding     = 4.6721;  % 更大的 f 值以修正背向

% Single HG
g_single   = 0.9807;

% TTHG (0.1-180度)
alpha_tthg = 0.4435;
g1_tthg    = 0.9900;
g2_tthg    = 0.8238;

% 4. 计算绘图数据
theta_deg_plot = logspace(log10(0.01), log10(180), 1000); % 0.1 到 180度
theta_rad_plot = theta_deg_plot * pi / 180;

P_fit_Ding   = Model_Ding(theta_rad_plot, w_mie_ding, gamma_ding, g_ding, f_ding);
P_fit_TTHG   = Model_TTHG(theta_rad_plot, alpha_tthg, g1_tthg, g2_tthg);
P_fit_Single = Model_SingleHG(theta_rad_plot, g_single);

%% ================== 第四部分：绘制图2 - 相函数 P(theta) ==================
figure('Name', 'Fig2: Phase Function', 'Position', [660, 100, 600, 500], 'Color', 'w');

% loglog(theta_deg_meas, P_meas, 'ko', 'MarkerSize', 4, 'DisplayName', 'Meas (Clear)');

loglog(theta_deg_plot, P_fit_Ding, 'r-', 'LineWidth', 2, 'DisplayName', 'Ding Mix');
hold on;
loglog(theta_deg_plot, P_fit_TTHG, 'b--', 'LineWidth', 2, 'DisplayName', 'TTHG');
loglog(theta_deg_plot, P_fit_Single, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Single HG');
loglog(theta_deg_emp, P_emp, 'm--', 'LineWidth', 2, 'DisplayName', 'Empirical Formula (Ref)');
% xline(90, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Fit Range');
grid on;
xlabel('Scattering Angle \theta (deg)', 'FontWeight', 'bold');
ylabel('Phase Function P(\theta) [sr^{-1}]', 'FontWeight', 'bold');
title('Fig 2: Phase Function P(\theta) (Normalized to 4\pi)');
legend('Location', 'SouthWest');
xlim([0.1, 180]); ylim([1e-4, 1e6]);

%% ================== 第五部分：绘制图3 - 散射角 PDF f(theta) ==================
% 公式: f(theta) = 0.5 * P(theta) * sin(theta)

% 计算 PDF
PDF_meas   = 0.5 * P_meas .* sin(theta_rad_meas);
PDF_Ding   = 0.5 * P_fit_Ding .* sin(theta_rad_plot);
PDF_TTHG   = 0.5 * P_fit_TTHG .* sin(theta_rad_plot);
PDF_Single = 0.5 * P_fit_Single .* sin(theta_rad_plot);

figure('Name', 'Fig3: Scattering Angle PDF', 'Position', [360, 50, 600, 500], 'Color', 'w');

% 使用半对数坐标 (semilogx) 或 双对数 (loglog) 
% 这里推荐 semilogx，因为 PDF 值本身跨度不像 P 那么大 (sin压制了极值)
%semilogx(theta_deg_meas, PDF_meas, 'ko', 'MarkerSize', 4, 'DisplayName', 'Meas (Clear)'); hold on;
semilogx(theta_deg_plot, PDF_Ding, 'r-', 'LineWidth', 2, 'DisplayName', 'Ding Mix');
hold on;
semilogx(theta_deg_plot, PDF_TTHG, 'b--', 'LineWidth', 2, 'DisplayName', 'TTHG');
semilogx(theta_deg_plot, PDF_Single, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Single HG');
semilogx(theta_deg_emp_ext, PDF_emp_ext, 'm--', 'LineWidth', 2, 'DisplayName', 'Empirical Formula (Ref)');
grid on;
xlabel('Scattering Angle \theta (deg)', 'FontWeight', 'bold');
ylabel('Probability Density f(\theta) [rad^{-1}]', 'FontWeight', 'bold');
title({'Fig 3: PDF of Plane Scattering Angle \theta', 'f(\theta) = 0.5 \cdot P(\theta) \cdot sin(\theta)'});
legend('Location', 'NorthEast');
xlim([0.1, 180]); 

% 检查积分是否归一化 (Debug用)
integral_check = trapz(theta_rad_plot, PDF_Ding);
fprintf('Ding 模型 PDF 积分检查 (应接近1): %.4f\n', integral_check);
