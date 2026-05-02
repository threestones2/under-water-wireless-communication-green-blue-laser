% 理想高斯光束受控扩束与水体散射展宽对比测试
clc; clear; close all;

% ================= 1. 物理参数定义 =================
w0 = 0.002;               % 初始束腰半径: 2 mm
L = 100;                  % 传输距离: 100 m
div_angle_deg = 1;        % 硬件控制的宏观全发散角: 1 度

% ================= 2. 纯几何展宽计算 (真空中) =================
% 计算半角(弧度)
theta_half_rad = (div_angle_deg / 2) * (pi / 180);

% 100m 处的纯几何光斑半径 (到达 1/e^2 光强处的边界)
w_geometric = w0 + L * tan(theta_half_rad);

fprintf('--- 纯几何展宽评估 (无水体散射) ---\n');
fprintf('初始束腰: %.3f m\n', w0);
fprintf('全发散角: %.1f 度\n', div_angle_deg);
fprintf('传输距离: %d m\n', L);
fprintf('纯几何到达半径 (1/e^2): %.3f m\n\n', w_geometric);

% ================= 3. 论文实际仿真数据 (深海信道中) =================
w_paper_450nm = 5.4;      % 450nm 激光在 Jerlov I 水质下的光斑半径
w_paper_520nm = 2.1;      % 520nm 激光在 Jerlov I 水质下的光斑半径

fprintf('--- 论文仿真展宽评估 (含水体多次散射) ---\n');
fprintf('450nm 接收光斑半径: %.3f m (散射放大倍数: %.1f 倍)\n', w_paper_450nm, w_paper_450nm / w_geometric);
fprintf('520nm 接收光斑半径: %.3f m (散射放大倍数: %.1f 倍)\n', w_paper_520nm, w_paper_520nm / w_geometric);

% ================= 4. 光强径向分布可视化 =================
r = linspace(0, 8, 1000); % 考察径向距离 0 到 8m

% 归一化光强公式: I(r) = I_center * exp(-2 * (r / w)^2)
I_geometric = exp(-2 * (r / w_geometric).^2);
I_scattered_520 = exp(-2 * (r / w_paper_520nm).^2);
I_scattered_450 = exp(-2 * (r / w_paper_450nm).^2);

figure('Name', 'Spot Size Comparison', 'Color', 'w', 'Position', [100, 100, 800, 500]);
plot(r, I_geometric, 'k-', 'LineWidth', 2.5, 'DisplayName', sprintf('Ideal Geometric (w = %.2fm)', w_geometric));
hold on;
plot(r, I_scattered_520, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('520nm with Scattering (w = %.1fm)', w_paper_520nm));
plot(r, I_scattered_450, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('450nm with Scattering (w = %.1fm)', w_paper_450nm));

% 标注 1/e^2 边界 (光斑半径的物理定义线)
yline(exp(-2), 'r:', 'LineWidth', 1.5, 'DisplayName', '1/e^2 Energy Boundary (~13.5%)');

grid on;
xlabel('Radial Distance r from Optical Axis (m)', 'FontWeight', 'bold');
ylabel('Normalized Spatial Irradiance I(r)/I_0', 'FontWeight', 'bold');
title('Spatial Irradiance Profile at 100 m: Geometric vs. Scattered Spread');
legend('Location', 'northeast', 'FontSize', 11);
xlim([0, 8]);
ylim([0, 1.05]);