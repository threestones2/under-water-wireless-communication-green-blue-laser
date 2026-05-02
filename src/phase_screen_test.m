%% 独立运行：基于 OTOPS 谱生成海水相位屏
clc; clear; close all;

%% 1. 参数设置 (提取自 MC_MPS_6.m)
lambda_nm = 514; 
lambda = lambda_nm * 1e-9;  % 波长 (m)

% 环境参数 (Petzold Clear Ocean + OTOPS)
T_avg = 20;        % 温度 (°C)
S_avg = 35;        % 盐度 (ppt)
epsilon = 1e-6;    % 动能耗散率 (m^2/s^3)
chi_T = 1e-8;      % 温度方差耗散率 (K^2/s)
eta = 1e-3;        % Kolmogorov 内尺度 (m)
H_ratio = -20;     % 温盐度梯度比

% 相位屏几何参数
D_screen = 1;             % 相位屏边长 (m)
N_grid = 2^8;               % 网格点数 (256x256)
Link_Dist = 10;            % 总链路距离 (m)
N_screens = 10;             % 屏的数量
delta_z_screen = Link_Dist / N_screens; % 单张屏代表的厚度 (m)

%% 2. 计算 OTOPS 谱句柄
fprintf('正在计算 OTOPS 功率谱系数...\n');
[Phi_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);

%% 3. 生成相位屏 (FFT + 次谐波补偿)
fprintf('正在生成相位屏 (D=%.1fm, N=%d)...\n', D_screen, N_grid);
phase_screen = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);

%% 4. 绘图显示
figure('Name', 'OTOPS Phase Screen', 'Color', 'w');
imagesc(linspace(-D_screen/2, D_screen/2, N_grid), ...
        linspace(-D_screen/2, D_screen/2, N_grid), ...
        phase_screen);
axis xy; axis square;
colormap jet; 
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title({['OTOPS Seawater Phase Screen (\Delta z = ' num2str(delta_z_screen) ' m)']; ...
       ['T=' num2str(T_avg) '^\circC, S=' num2str(S_avg) 'ppt, \epsilon=' num2str(epsilon) ', \chi_T=' num2str(chi_T)]});

fprintf('完成。相位屏最大值: %.4f rad, 最小值: %.4f rad\n', max(phase_screen(:)), min(phase_screen(:)));

%% --- [新增] 相位屏相邻网格梯度核验模块 ---
% 假设您的相位屏变量名为 phase_screen

% 1. 计算相邻相位差
% diff(X, 1, 2) 计算每一行中相邻列的差 (水平梯度)
dx_phase = diff(phase_screen, 1, 2); 

% diff(X, 1, 1) 计算每一列中相邻行的差 (垂直梯度)
dy_phase = diff(phase_screen, 1, 1);

% 将所有梯度值合并到一个数组中以便统计
all_diffs = [dx_phase(:); dy_phase(:)];
abs_diffs = abs(all_diffs);

% 2. 统计计算
max_diff = max(all_diffs);
min_diff = min(all_diffs);
max_abs_diff = max(abs_diffs);
mean_abs_diff = mean(abs_diffs);
std_diff = std(all_diffs);

% 3. 打印结果
fprintf('\n========================================\n');
fprintf('      相位屏相邻网格梯度统计核验      \n');
fprintf('========================================\n');
fprintf('网格尺寸: %d x %d\n', size(phase_screen, 1), size(phase_screen, 2));
fprintf('最大相邻相位差 (+):  %.4f rad\n', max_diff);
fprintf('最小相邻相位差 (-):  %.4f rad\n', min_diff);
fprintf('最大绝对相位差 (|Δ|): %.4f rad\n', max_abs_diff);
fprintf('平均绝对相位差 (|Δ|): %.4f rad\n', mean_abs_diff);
fprintf('相位差标准差 (RMS):   %.4f rad\n', std_diff);
fprintf('----------------------------------------\n');

% 4. 自动判定是否满足采样约束
if max_abs_diff < pi
    fprintf('✅ [通过] 满足 Nyquist 采样约束 (|Δφ| < π)。\n');
    fprintf('   最大梯度 %.2f π，未发生相位卷绕。\n', max_abs_diff/pi);
else
    fprintf('❌ [警告] 不满足采样约束！最大梯度达到 %.2f π。\n', max_abs_diff/pi);
    fprintf('   可能发生相位卷绕 (Aliasing)。建议增加网格数 N 或减小屏间距 delta_z。\n');
end
fprintf('========================================\n');

% 5. (可选) 绘制梯度直方图
figure('Name', 'Phase Gradient Distribution', 'Color', 'w');
histogram(all_diffs, 100, 'Normalization', 'pdf');
xlabel('相邻相位差 (rad)');
ylabel('概率密度');
title(['相位梯度分布 (Max: ' num2str(max_abs_diff, '%.2f') ' rad)']);
xline(pi, 'r--', 'Label', '+\pi');
xline(-pi, 'r--', 'Label', '-\pi');
grid on;