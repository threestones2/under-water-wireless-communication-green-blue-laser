%% 相位函数归一化验证脚本
% 功能：验证 Rayleigh 和 Mie 混合相位函数在全立体角上的积分是否为 1
% 理论基础：Integral(P(mu) * 2*pi * dmu) from -1 to 1 should be 1.

clc; clear;

% 1. 参数设置 (保持与主程序一致)
params = struct();
params.k_s_ray = 0.266e-3;  % Rayleigh散射系数
params.k_s_mie = 0.284e-3;  % Mie散射系数
params.gamma = 0.017;       % Rayleigh参数
params.g = 0.72;            % Mie不对称参数
params.f = 0.5;             % Mie修正参数

% 2. 计算权重
w_ray = params.k_s_ray / (params.k_s_ray + params.k_s_mie);
w_mie = 1 - w_ray;

fprintf('--- 参数信息 ---\n');
fprintf('Rayleigh 权重: %.4f\n', w_ray);
fprintf('Mie 权重:      %.4f\n', w_mie);
fprintf('----------------\n');

% 3. 定义相位函数句柄 (P_total)
% 输入: mu (cos_theta), 标量或向量
% 输出: 相位函数值

% Rayleigh 部分 P_ray(mu)
P_ray = @(mu) (3 * (1 + 3*params.gamma + (1-params.gamma).*mu.^2)) ./ ...
              (16 * pi * (1 + 2*params.gamma));

% Mie 部分 P_mie(mu)
% Term 1: 标准 HG 相位函数
term1 = @(mu) (1 - params.g^2) ./ (4 * pi) .* ...
              (1 ./ (1 + params.g^2 - 2*params.g.*mu).^(1.5));
          
% Term 2: 修正项 (Correction Term)
% 注意：代码中分母是 ((1+g^2)^(3/2))，这是一个常数
term2 = @(mu) (1 - params.g^2) / (4 * pi) * params.f .* ...
              (0.5 * (3*mu.^2 - 1)) ./ ((1 + params.g^2)^(1.5));

P_mie = @(mu) term1(mu) + term2(mu);

% 总相位函数
P_total = @(mu) w_ray .* P_ray(mu) + w_mie .* P_mie(mu);

% 4. 执行数值积分验证
% 积分区间：mu 从 -1 到 1 (对应 theta 从 pi 到 0)
% 积分公式：2*pi * Integral(P(mu) dmu)

% 定义被积函数 (包含 2*pi 因子)
integrand = @(mu) 2 * pi * P_total(mu);

% 使用 integral 函数进行高精度积分
[integral_val] = integral(integrand, -1, 1);

% 5. 输出结果
fprintf('积分验证结果:\n');
fprintf('Rayleigh 部分积分: %.6f\n', integral(@(mu) 2*pi*P_ray(mu), -1, 1));
fprintf('Mie 部分积分:      %.6f\n', integral(@(mu) 2*pi*P_mie(mu), -1, 1));
fprintf('总相位函数积分:    %.15f\n', integral_val);

if abs(integral_val - 1) < 1e-6
    fprintf('\n[结论] 验证通过！相位函数已正确归一化。\n');
else
    fprintf('\n[警告] 积分不等于1，请检查公式系数！\n');
end

% 6. 可视化相位函数形状 (附加功能)
mu_vals = linspace(-1, 1, 1000);
theta_vals = acos(mu_vals) * 180/pi; % 转换为角度
P_vals = P_total(mu_vals);

figure;
semilogy(theta_vals, P_vals, 'LineWidth', 2);
xlabel('散射角 \theta (degrees)');
ylabel('相位函数 P(\theta)');
title('混合相位函数分布 (对数坐标)');
grid on;
xlim([0 180]);