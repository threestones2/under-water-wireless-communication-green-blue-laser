% =========================================================================
% 示例：调用单步传播函数 ang_spec_prop_circular 并显示相位
% =========================================================================

clear; clc; close all;

% --- 1. 基础参数设置 ---
N = 1024;               % 网格点数 (建议 512 或 1024)
L = 0.1;                % 源平面总宽度 [m] (10 cm)
delta1 = L / N;         % 源平面网格间距
wvl = 1.55e-6;          % 波长 [m] (比如通信波段 1550nm)
k = 2*pi/wvl;           % 波数

% --- 2. 传播几何参数 ---
Dz = 2000;              % 传播距离 [m] (2 km)
% 设置观察平面的网格间距 delta2
% 选项 A: 保持与源平面一致 (适用于准直光或近场)
% delta2 = delta1; 
% 选项 B: 让网格随衍射扩大 (适用于远场，防止光束跑出网格)
delta2 = delta1 * 2;    % 这里演示放大 2 倍网格

% --- 3. 光束参数 (关键) ---
D_ap = 0.05;            % 光阑直径 [m] (5 cm)
w0 = 0.02;              % 高斯束腰半径 [m] (2 cm)
% 也就是光束半径 2cm，通过一个 5cm 的孔径发射

% --- 4. 准备输入场 ---
% 因为光束整形都在函数内部做，这里只需要给一个全 1 的矩阵
Uin_base = ones(N, N); 

% --- 5. 调用函数 ---
% function [x2, y2, Uout] = ang_spec_prop_circular(Uin, wvl, d1, d2, Dz, w0, D_ap)
[x2, y2, Uout] = ang_spec_prop_circular_gaussian(Uin_base, wvl, delta1, delta2, Dz, w0, D_ap);

% --- 6. 数据处理 ---

% (1) 计算光强 (Irradiance)
Iout = abs(Uout).^2;

% (2) 计算相位 (Phase)
% angle 的结果在 [-pi, pi] 之间卷绕 (Wrapped Phase)
raw_phase = angle(Uout);

% (3) 相位掩膜技巧 (Phase Masking)
% 在光强极弱的地方 (比如 < 峰值的 1/1000)，相位计算受数值噪声影响很大，
% 画出来会像雪花一样乱。为了看清主光束相位，我们将暗处的相位设为 NaN。
threshold = max(Iout(:)) * 1e-4; 
masked_phase = raw_phase;
masked_phase(Iout < threshold) = NaN; 

% --- 7. 绘图显示 ---
figure('Name', '单步传播结果：光强与相位', 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% --- 左图：输出光强 ---
subplot(1, 2, 1);
imagesc(x2(1,:), y2(:,1), Iout);
axis image;         % 保持 x,y 比例
colormap('hot');    
colorbar;
title(['输出光强 I (z = ' num2str(Dz) ' m)']);
xlabel('x [m]'); ylabel('y [m]');
% 限制显示范围，聚焦中心
limit_range = max(x2(:)) * 0.8; % 显示中心 80% 的区域
xlim([-limit_range, limit_range]); 
ylim([-limit_range, limit_range]);

% --- 右图：输出相位 ---
subplot(1, 2, 2);
% 使用 AlphaData 让 NaN 区域透明
h_phase = imagesc(x2(1,:), y2(:,1), masked_phase);
set(h_phase, 'AlphaData', ~isnan(masked_phase));
axis image;
colormap(gca, 'jet'); % 相位通常用彩虹色
cb = colorbar;
cb.Label.String = 'Phase [rad]';
title('输出相位 \phi (Wrapped Phase)');
xlabel('x [m]'); ylabel('y [m]');
xlim([-limit_range, limit_range]); 
ylim([-limit_range, limit_range]);

% --- 补充：验证是否为球面波 ---
% 如果是理想传播，相位应该是二次曲面 (抛物面/球面)
% 画一条中心切线来看看相位展开后的样子
figure('Name', '相位中心切线', 'Color', 'w');
center_idx = floor(N/2) + 1;
phase_slice = unwrap(raw_phase(center_idx, :)); % 一维展开
plot(x2(center_idx, :), phase_slice, 'LineWidth', 1.5);
grid on;
title('中心切线处的展开相位 (Unwrapped Phase)');
xlabel('x [m]'); ylabel('Phase [rad]');
xlim([-limit_range, limit_range]);