%% 该函数是水下湍流功率谱生成函数

function [Phi_n, kappa] = water_turbulence_spectrum(T_avg, S_avg, lambda, epsilon, chi_T, H, kappa_range)
% 计算水湍流的空间功率谱
% 输入参数:
%   T_avg: 平均温度 [°C]
%   S_avg: 平均盐浓度 [ppt]
%   lambda: 光波长 [nm]
%   epsilon: 能量耗散率 [m^2/s^3]
%   chi_T: 温度方差耗散率 [K^2/s]
%   H: 温度-盐度梯度比 [°C/ppt]
%   kappa_range: 空间波数范围 [m^-1] (可选，默认为1e-3:1e3:1e4)
%
% 输出参数:
%   Phi_n: 水湍流的空间功率谱
%   kappa: 空间波数

% 如果没有提供kappa_range，则使用默认值
if nargin < 7
    kappa = logspace(-3, 3, 1000); % 从10^-3到10^3，对数间隔的1000个点
else
    kappa = kappa_range;
end

% 常数定义
a0 = 1.31405;
a1 = 1.779e-4;
a2 = -1.05e-6;
a3 = 1.6e-8;
a4 = -2.02e-6;
a5 = 15.868;
a6 = 0.01155;
a7 = -0.00423;
a8 = -4382;
a9 = 1.1455e6;
beta0 = 0.72;

% 1. 计算Prandtl数和Schmidt数
Pr = calculate_Prandtl(T_avg, S_avg);
Sc = calculate_Schmidt(T_avg, S_avg);

% 2. 计算涡流扩散率比
d_r = calculate_eddy_diffusivity_ratio(T_avg, S_avg, H);

% 3. 计算线性系数A和B
A = calculate_coefficient_A(T_avg, S_avg, lambda, a2, a3, a4, a7);
B = calculate_coefficient_B(T_avg, S_avg, lambda, a1, a2, a3, a6);

% 4. 计算动态粘度和密度
mu = calculate_dynamic_viscosity(T_avg, S_avg);
rho = calculate_density(T_avg, S_avg);

% 5. 计算Kolmogorov微尺度
nu = mu ./ rho;  % 运动粘度
eta = (nu.^(3/4)) .* (epsilon.^(-1/4));  % Kolmogorov微尺度

% 6. 计算无量纲参数c_i
c_T = (0.072.^(4/3)) .* beta0 ./ Pr;
c_S = (0.072.^(4/3)) .* beta0 ./ Sc;
c_TS = (0.072.^(4/3)) .* beta0 .* (Pr + Sc) ./ (2 .* Pr .* Sc);

% 7. 计算方差耗散率
chi_S = d_r ./ (H.^2) .* chi_T;
chi_TS = (1 + d_r) ./ (2 .* H) .* chi_T;

% 8. 计算温度谱、盐度谱和互谱
Phi_T = calculate_spectrum_component(kappa, eta, epsilon, chi_T, c_T);
Phi_S = calculate_spectrum_component(kappa, eta, epsilon, chi_S, c_S);
Phi_TS = calculate_spectrum_component(kappa, eta, epsilon, chi_TS, c_TS);

% 9. 计算总功率谱
Phi_n = (A.^2) .* Phi_T + (B.^2) .* Phi_S + 2 .* A .* B .* Phi_TS;

end

function Pr = calculate_Prandtl(T_avg, S_avg)
% 计算Prandtl数
% 公式: Pr = c_p * mu / sigma_T

% 计算比热容 c_p [J/kg/K]
a11 = 5.328 - 9.76e-2 * S_avg + 4.04e-4 * S_avg^2;
a12 = -6.913e-3 + 7.351e-4 * S_avg - 3.15e-6 * S_avg^2;
a13 = 9.6e-6 - 1.927e-6 * S_avg + 8.23e-9 * S_avg^2;
a14 = 2.5e-9 + 1.666e-9 * S_avg - 7.125e-12 * S_avg^2;
c_p = 1000 * (a11 + a12 * T_avg + a13 * T_avg^2 + a14 * T_avg^2);

% 计算动态粘度 mu [N·s/m^2]
mu = calculate_dynamic_viscosity(T_avg, S_avg);

% 计算热导率 sigma_T [W/m/K]
T_h = 1.00024 * T_avg;
S_h = S_avg / 1.00472;
log_sigma_T = log(240 + 0.0002 * S_h) - 3 + 0.434 * ...
    (2.3 - (343.5 + 0.037 * S_h) / (T_h + 273.15)) * ...
    (1 - (T_h + 273.15) / (647.3 + 0.03 * S_h)).^(1/3);
sigma_T = 10.^log_sigma_T;

% 计算Prandtl数
Pr = c_p .* mu ./ sigma_T;
end

function Sc = calculate_Schmidt(T_avg, S_avg)
% 计算Schmidt数
% 公式: Sc ≈ μ^2 / (5.954e-15 * T_avg * ρ)

% 计算动态粘度 mu [N·s/m^2]
mu = calculate_dynamic_viscosity(T_avg, S_avg);

% 计算密度 rho [kg/m^3]
rho = calculate_density(T_avg, S_avg);

% 计算Schmidt数
Sc = mu.^2 ./ (5.954e-15 * T_avg * rho);
end

function d_r = calculate_eddy_diffusivity_ratio(T_avg, S_avg, H)
% 计算涡流扩散率比
% 简化模型，基于论文中的公式(9)-(11)

% 计算热膨胀系数alpha和盐收缩系数beta
% 这里使用简化值，实际应用中可以使用TEOS-10工具箱
% 推荐值: alpha ≈ 2e-4 [1/°C], beta ≈ 7.6e-4 [1/ppt]
alpha = 2e-4;  % 热膨胀系数 [1/°C]
beta = 7.6e-4; % 盐收缩系数 [1/ppt]

% 计算密度比R_rho
R_rho = alpha ./ abs(H) ./ beta;

% 计算涡流扩散率比d_r
d_r = zeros(size(R_rho));
idx1 = (R_rho >= 1);
idx2 = (R_rho >= 0.5) & (R_rho < 1);
idx3 = (R_rho < 0.5);

d_r(idx1) = (R_rho(idx1) + sqrt(R_rho(idx1) .* (R_rho(idx1) - 1))) ./ (R_rho(idx1) - 1);
d_r(idx2) = 1.85 .* R_rho(idx2) - 0.85;
d_r(idx3) = 0.15 .* R_rho(idx3);
end

function A = calculate_coefficient_A(T_avg, S_avg, lambda, a2, a3, a4, a7)
% 计算线性系数A
% 公式: A = a2*S_avg + 2*a3*T_avg*S_avg + 2*a4*T_avg + a7/lambda

A = a2 .* S_avg + 2 .* a3 .* T_avg .* S_avg + 2 .* a4 .* T_avg + a7 ./ lambda;
end

function B = calculate_coefficient_B(T_avg, S_avg, lambda, a1, a2, a3, a6)
% 计算线性系数B
% 公式: B = a1 + a2*T_avg + a3*T_avg^2 + a6/lambda

B = a1 + a2 .* T_avg + a3 .* T_avg.^2 + a6 ./ lambda;
end

function mu = calculate_dynamic_viscosity(T_avg, S_avg)
% 计算动态粘度 [N·s/m^2]

s = S_avg * 1e-3;
a21 = 1.5409136040 + 1.9981117208e-2 * T_avg - 9.5203865864e-5 * T_avg.^2;
a22 = 7.9739318223 - 7.5614568881e-2 * T_avg + 4.7237011074e-4 * T_avg.^2;
mu0 = 1 ./ (0.15700386464 .* (T_avg + 64.992620050).^2 - 91.296496657) + 4.2844324477e-5;
mu = mu0 .* (a21 .* s + a22 .* s.^2);
end

function rho = calculate_density(T_avg, S_avg)
% 计算水的密度 [kg/m^3]

s = S_avg * 1e-3;
% 温度对密度的贡献
rho_T = 9.9992293295e2 + 2.0341179217e-2 * T_avg - 6.1624591598e-3 * T_avg.^2 + ...
    2.2614664708e-5 * T_avg.^3 - 4.6570659168e-8 * T_avg.^4;
% 盐度对密度的调整
rho_TS = s .* (8.0200240891e2 - 2.0005183488 * T_avg + 1.6771024982e-2 * T_avg.^2 - ...
    3.0600536746e-5 * T_avg.^3 - 1.6132224742e-5 * T_avg.^2 .* s);
% 总密度
rho = rho_T + rho_TS;
end

function Phi = calculate_spectrum_component(kappa, eta, epsilon, chi, c)
% 计算谱分量(温度谱、盐度谱或互谱)
% 公式(20)

beta0 = 0.72;
kappa_eta = kappa .* eta;

% 计算谱分量
Phi = (1 + 21.61 .* kappa_eta.^0.61 .* c.^0.02 - 18.18 .* kappa_eta.^0.55 .* c.^0.04) .* ...
    (1 ./ (4 * pi)) .* beta0 .* epsilon.^(-1/3) .* kappa.^(-11/3) .* chi .* ...
    exp(-176.90 .* kappa_eta.^2 .* c.^0.96);
end



% 示例：计算水湍流的空间功率谱
% 参数设置
T_avg = 15;      % 平均温度 [°C]
S_avg = 35;      % 平均盐浓度 [ppt]
lambda = 532;    % 光波长 [nm]
epsilon = 1e-2;  % 能量耗散率 [m^2/s^3]
chi_T = 1e-5;    % 温度方差耗散率 [K^2/s]
H = -20;         % 温度-盐度梯度比 [°C/ppt]

% 计算功率谱
[Phi_n, kappa] = water_turbulence_spectrum(T_avg, S_avg, lambda, epsilon, chi_T, H);

% 绘制结果
figure;
loglog(kappa, Phi_n);
xlabel('空间波数 \kappa [m^{-1}]');
ylabel('功率谱 \Phi_n(\kappa) [m^3]');
title(sprintf('水湍流空间功率谱 (T=%.1f^{\\circ}C, S=%.1fppt, \\lambda=%.0fnm)', T_avg, S_avg, lambda));
grid on;
