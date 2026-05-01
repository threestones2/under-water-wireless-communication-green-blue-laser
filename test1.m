% =========================================================================
% 湍流相位屏“低频漂移(Beam Wander)”统计特性验证脚本
% 验证目的: 检查粗网格相位屏的几何折射是否精准满足解析的光束漂移宏观方差演化
% =========================================================================
clc; clear; close all;

%% 1. 参数设置 (完全沿用您的 MC_MPS_v1.m 设定)
lambda = 514e-9;
k_wave = 2*pi/lambda;
L = 100;                    % 传输距离 100m
N_screens = 20;             % 20张相位屏
delta_z = L / N_screens;    % 屏间距 5m

D_screen = 10;             % 屏幕尺寸 0.2m
N_grid = 2^10;               % 网格数 256
dx = D_screen / N_grid;     % 空间分辨率 ~0.78mm

% OTOPS 弱湍流参数
T_avg = 20; S_avg = 35; H_ratio = -20;
epsilon = 1e-9; chi_T = 1e-7; eta = 1e-3;

% 射线数量 (Monte Carlo 样本数)
N_rays = 50000;

%% 2. 获取谱函数
[Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);

%% 3. 预先生成相位屏并提取理论宏观角方差
fprintf('正在生成 %d 张湍流相位屏并提取空间梯度统计特征...\n', N_screens);
gx_all = zeros(N_grid, N_grid, N_screens);
gy_all = zeros(N_grid, N_grid, N_screens);
var_theta_k_list = zeros(N_screens, 1);

for i = 1:N_screens
    phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z);
    [gx, gy] = gradient(phi, dx);
    gx_all(:,:,i) = gx;
    gy_all(:,:,i) = gy;
    
    % [理论核心] 单张屏幕的二维偏折角方差 = 梯度方差 / k^2
    var_theta_x = var(gx(:)) / k_wave^2;
    var_theta_y = var(gy(:)) / k_wave^2;
    var_theta_k_list(i) = var_theta_x + var_theta_y;
end

% 计算单屏平均角方差
mean_var_theta = mean(var_theta_k_list);

% [解析推导]: 假设各屏相互独立 (Markov近似)，计算接收面上总的宏观漂移方差
% 理论方差 = Σ (单屏角方差 * 剩余距离^2)
z_screens = (1:N_screens) * delta_z;
theoretical_wander_var = sum(mean_var_theta .* (L - z_screens).^2);

%% 4. 向量化射线追踪 (Vectorized Ray Tracing)
fprintf('正在进行 %d 条射线的 Monte Carlo 宏观几何追踪...\n', N_rays);
% 为避免所有射线都在同一个网格点上，在中心微小区域内均匀撒点
pos_x = (rand(N_rays, 1) - 0.5) * 0.05; 
pos_y = (rand(N_rays, 1) - 0.5) * 0.05;
pos_x_init = pos_x;
pos_y_init = pos_y;

dir_x = zeros(N_rays, 1); % 初始方向严格沿 Z 轴
dir_y = zeros(N_rays, 1);

for i = 1:N_screens
    % 1. 飞到当前屏
    pos_x = pos_x + dir_x * delta_z;
    pos_y = pos_y + dir_y * delta_z;
    
    % 2. 计算在当前屏上的网格索引 (周期性边界映射)
    idx_x = mod(round((pos_x - (-D_screen/2)) / dx), N_grid) + 1;
    idx_y = mod(round((pos_y - (-D_screen/2)) / dx), N_grid) + 1;
    
    idx_x(idx_x < 1) = 1; idx_x(idx_x > N_grid) = N_grid;
    idx_y(idx_y < 1) = 1; idx_y(idx_y > N_grid) = N_grid;
    
    % 3. 提取相位屏梯度并更新飞行方向
    linear_indices = sub2ind([N_grid, N_grid], idx_y, idx_x);
    gx_current = gx_all(:,:,i);
    gy_current = gy_all(:,:,i);
    
    delta_dir_x = gx_current(linear_indices) / k_wave;
    delta_dir_y = gy_current(linear_indices) / k_wave;
    
    dir_x = dir_x + delta_dir_x;
    dir_y = dir_y + delta_dir_y;
end

%% 5. 统计结果与精度比对
% 计算相对偏移量
wander_x = pos_x - pos_x_init; 
wander_y = pos_y - pos_y_init;
simulated_wander_var = var(wander_x) + var(wander_y);

fprintf('\n=== 漂移方差 (Beam Wander Variance) 严谨验证 ===\n');
fprintf('相空间理论解析方差: %.6e m^2\n', theoretical_wander_var);
fprintf('MC射线追踪统计方差: %.6e m^2\n', simulated_wander_var);
fprintf('宏观漂移相对误差  : %.2f %%\n', abs(simulated_wander_var - theoretical_wander_var)/theoretical_wander_var * 100);

% 绘制光斑漂移散点图
figure('Color', 'w');
scatter(wander_x*1000, wander_y*1000, 5, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
axis equal; grid on;
xlabel('X Wander / 漂移量 (mm)'); ylabel('Y Wander / 漂移量 (mm)');
title(sprintf('Monte Carlo Beam Wander at Receiver Plane (L = 100m)\nSimulated Var: %.2e | Theoretical Var: %.2e', ...
    simulated_wander_var, theoretical_wander_var));

%% 辅助函数 (保持与您代码完全一致)
function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    lambda_nm = lambda * 1e9; k_wave = 2 * pi / lambda; 
    a1 = 1.779e-4;  a2 = -1.05e-6;  a3 = 1.6e-8; 
    a4 = -2.02e-6;  a5 = 1.155e-2;  a6 = -4.23e-3;
    A = a2*S + 2*a3*T*S + 2*a4*T + a6/lambda_nm; 
    B = a1 + a2*T + a3*T^2 + a5/lambda_nm; 
    
    T_k = T + 273.15; s_frac = S * 1e-3; 
    a11 = 5.328 - 9.76e-2*S + 4.04e-4*S^2; a12 = -6.913e-3 + 7.351e-4*S - 3.15e-6*S^2;
    a13 = 9.6e-6 - 1.927e-6*S + 8.23e-9*S^2; a14 = 2.5e-9 + 1.666e-9*S - 7.125e-12*S^2;
    cp = 1000 * (a11 + a12*T + a13*T^2 + a14*T^3); 
    
    rho_T = 9.9992293295e2 + 2.0341179217e-2*T - 6.1624591598e-3*T^2 + 2.2614664708e-5*T^3 - 4.6570659168e-8*T^4;
    rho_S = s_frac * (8.0200240891e2 - 2.0005183488*T + 1.6771024982e-2*T^2 - 3.0600536746e-5*T^3 - 1.6132224742e-5*T*S);
    rho = rho_T + rho_S;
    
    mu_0 = (0.15700386464*(T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5;
    a21 = 1.5409136040 + 1.9981117208e-2*T - 9.5203865864e-5*T^2;
    a22 = 7.9739318223 - 7.561456881e-2*T + 4.7237011074e-4*T^2;
    mu = mu_0 * (1 + a21*s_frac + a22*s_frac^2);
    
    T_b = 1.00024 * T; S_b = S / 1.00472;
    log_sigma = log10(240 + 0.0002*S_b) - 3 + 0.434 * (2.3 - (343.5 + 0.037*S_b)/(T_b + 273.15)) * (1 - (T_b + 273.15)/(647.3 + 0.03*S_b))^(1/3);
    sigma_T = 10^log_sigma; 
    
    Pr = mu * cp / sigma_T;  Sc = mu^2 / (5.954e-15 * T_k * rho); 
    c_T = 0.072^(4/3) * Pr^(-1); c_S = 0.072^(4/3) * Sc^(-1); c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    alpha_c = 2.6e-4; beta_c = 7.6e-4; R_rho = alpha_c * abs(H_ratio) / beta_c;
    if R_rho >= 1, d_r = R_rho + sqrt(R_rho)*sqrt(R_rho-1);
    elseif R_rho >= 0.5, d_r = 1.85*R_rho - 0.85;
    else, d_r = 0.15*R_rho; end
    
    chi_S = chi_T * d_r / (H_ratio^2); chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    coeff_Hill = 0.72 / (4 * pi); 
    Phi_Hill = @(K, chi_M, c_M) coeff_Hill * chi_M * epsilon^(-1/3) .* (K.^2 + 1e-25).^(-11/6) .* exp(-176.90 * (K*eta).^2 .* c_M^(0.96)) .* (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + B^2 * Phi_Hill(K, chi_S, c_S) + 2*A*B * Phi_Hill(K, chi_TS, c_TS));
end

function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2 * pi / D; dx = D / N;
    kx = (-N/2 : N/2-1) * dk; [KX, KY] = meshgrid(kx, kx);
    K_grid = sqrt(KX.^2 + KY.^2); K_grid(N/2+1, N/2+1) = 1e-10; 
    Phi_n_val = Phi_n_func(K_grid); Phi_n_val(N/2+1, N/2+1) = 0; 
    F_phi = 2 * pi * k_wave^2 * delta_z * Phi_n_val;
    noise = (randn(N) + 1i * randn(N));
    C_nm = noise .* sqrt(F_phi) * dk;
    phase_high = real(ifftshift(ifft2(ifftshift(C_nm))) )*N^2;
    phase_low = zeros(N, N); [xx, yy] = meshgrid((-N/2 : N/2-1) * dx);
    for p = 1:3
        dk_p = dk / (3^p); 
        for m = -1:1
            for n = -1:1
                if (m==0 && n==0), continue; end
                kx_p = m * dk_p; ky_p = n * dk_p; k_p = sqrt(kx_p^2 + ky_p^2);
                F_phi_p = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_func(k_p);
                amp = sqrt(F_phi_p) * dk_p; r_c = (randn(1) + 1i * randn(1));
                phase_low = phase_low + real( r_c * amp * exp(1i * (kx_p * xx + ky_p * yy)) );
            end
        end
    end
    phase_screen = phase_high + phase_low;
end