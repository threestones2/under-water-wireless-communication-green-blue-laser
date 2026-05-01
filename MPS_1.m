function phase_screen = generate_phase_screen(D, N, lambda, d, T_avg, S_avg, epsilon, chi_T, eta, H_ratio)
    % MPS: 生成海洋湍流相位屏 (基于 Yao et al. 2020 OTOPS 模型)
    % 修正版: 修复了噪声归一化、FFT缩放系数及物理参数符号
    %
    % 输入参数:
    %   D: 相位屏物理尺寸 (m)
    %   N: 采样点数 (如 256)
    %   lambda: 光波长 (m)
    %   d: 相邻相位屏间距 (m)
    %   T_avg: 平均温度 (°C)
    %   S_avg: 平均盐度 (ppt)
    %   epsilon: 湍流动能耗散率 (m²/s³)
    %   chi_T: 温度方差耗散率 (K²/s)
    %   eta: 湍流内尺度 (m)
    %   H_ratio: 温度盐度变化率 (°C/ppt), 通常为负值 (如 -20)
    %
    % 输出:
    %   phase_screen: 生成的相位屏 (N×N矩阵, 单位: rad)

    % 1. 计算频率域网格
    % 频率分辨率 (波数域)
    dk = 2 * pi / D;          
    
    % 构建频域坐标 (将零频移到中心)
    kx = (-N/2 : N/2-1) * dk;   
    [KX, KY] = meshgrid(kx, kx);
    K = sqrt(KX.^2 + KY.^2);  % 径向波数
    % 这里计算空间分辨率
    dx=1/N/dk;
    
    % 2. 计算OTOPS模型参数 (Yao et al. 2020)
    % 经验系数
    a1 = 1.779e-4;  a2 = -1.05e-6;  a3 = 1.6e-8; 
    a4 = -2.02e-6;  a5 = 1.155e-2;  a6 = -4.23e-3;
    
    % 计算折射率系数 A 和 B
    A = a2*S_avg + 2*a3*T_avg*S_avg + 2*a4*T_avg + a6/(lambda*10^9); 
    B = a1 + a2*T_avg + a4*T_avg^2 + a5/(lambda*10^9);
    
    % 计算 Pr 和 Sc (此处使用典型值，高精度需调用 TEOS-10)
    rho = 1025;         
    mu = 1e-3;          
    cp = 4182;          
    sigma_T = 0.6;      
    Pr = mu * cp / sigma_T;  
    Sc = mu^2 / (5.954e-15 * (T_avg + 273.15) * rho); 
    
    % 计算 Hill 模型参数 c_T, c_S, c_TS
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1);
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    % 计算涡流扩散比 d_r
    alpha_c = 2.6e-4;   
    beta_c = 7.6e-4;  
    R_rho = alpha_c * abs(H_ratio) / beta_c;  % 密度比必须为正
    
    if R_rho >= 1
        d_r = R_rho + sqrt(R_rho) * sqrt(R_rho - 1);
    elseif R_rho >= 0.5
        d_r = 1.85 * R_rho - 0.85;
    else
        d_r = 0.15 * R_rho;
    end
    
    % 计算耗散率
    chi_S = chi_T * d_r / (H_ratio^2);
    % [修正1] chi_TS 符号需与 H_ratio 一致 (通常为负)
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    
    % 3. 计算折射率功率谱 Phi_n (Hill Model)
    % 系数 0.033 * ... (0.033 = 0.072 / (4*pi) * 2pi^2 ? No, standard constant)
    % 标准 Hill 谱公式系数为 0.033
    Phi_M = @(chi_M, c_M) 0.033 * chi_M * epsilon^(-1/3) .* ...
        K.^(-11/3) .* exp(-176.90 * (K*eta).^2 .* c_M^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    
    Phi_T = Phi_M(chi_T, c_T);
    Phi_S = Phi_M(chi_S, c_S);
    Phi_TS = Phi_M(chi_TS, c_TS); 
    
    % 总功率谱 (Eq. 19)
    Phi_n = A^2 * Phi_T + B^2 * Phi_S + 2*A*B*Phi_TS;
    Phi_n(K == 0) = 0; % 移除直流分量奇点
    
    % 4. 计算相位功率谱密度 F_phi (Eq. 12)
    k_wave = 2 * pi / lambda;
    F_phi = 2 * pi^2 * k_wave^2 * d * Phi_n;
    
    % 5. 生成复高斯随机数矩阵 (频域)
    % [修正2] 单位化: 方差应为 1 (实部虚部各 0.5)
    noise = (randn(N) + 1i * randn(N)) / sqrt(2);
    
    % 6. 构建频域相位屏 H
    % [修正3] 缩放系数: 
    % sqrt(F_phi) * dk 是物理幅值密度
    % 乘以 N^2 是为了抵消 ifft2 的 1/N^2
    H = noise .* sqrt(F_phi) * dk;
    
    % 7. 转换到空间域
    %H = ifftshift(H); % 将零频移回左上角以适配 FFT 算法
    phase_screen = real(ift2(H,dk));
    
end


% 设置参数 (根据论文表1)
D = 1;          % 相位屏尺寸 (m)
N = 4096;        % 采样点数
lambda = 532e-9;% 波长 (m)
d = 2;          % 相邻相位屏间距 (m)
T_avg = 20;     % 平均温度 (°C)
S_avg = 35;     % 平均盐度 (ppt)
epsilon = 1e-3; % 湍流动能耗散率 (弱湍流)
chi_T = 1e-7;   % 温度方差耗散率 (弱湍流)
eta = 1e-3;     % 湍流内尺度 (m)
H_ratio = -20;  % 温度盐度变化率 (°C/ppt)

% 生成相位屏
phase_screen = generate_phase_screen(D, N, lambda, d, T_avg, S_avg, ...
                                    epsilon, chi_T, eta, H_ratio);

% 可视化相位屏
% 定义物理坐标向量 (将原点放在中心，范围 -0.5m 到 0.5m)
dx = D / N;  % 网格间距 (约 3.9 mm)
x = (-N/2 : N/2-1) * dx; 
y = x;

% 可视化相位屏 (指定 x 和 y 轴)
figure;
imagesc(x, y, phase_screen); % [关键修改] 加入 x, y 参数
colormap('jet');
colorbar;

% 修正坐标轴标签
xlabel('x (m)');
ylabel('y (m)');
title(['Phase Screen (D = ' num2str(D) ' m)']);
axis square; % 保持横纵比例一致
axis xy;     % 确保 y 轴方向正确（从下到上）

% 计算相邻网格相位差
diff_x = diff(phase_screen, 1, 1); % x方向差分
diff_y = diff(phase_screen, 1, 2); % y方向差分

max_diff = max(max(abs([diff_x(:); diff_y(:)])));
mean_diff = mean(abs([diff_x(:); diff_y(:)]));

fprintf('相位值范围: [%.2f, %.2f] rad\n', min(phase_screen(:)), max(phase_screen(:)));
fprintf('最大相邻相位差: %.4f rad (阈值 pi)\n', max_diff);
fprintf('平均相邻相位差: %.4f rad\n', mean_diff);

if max_diff > pi
    warning('警告: 相位屏存在欠采样！建议增加 N 或减小 D，或者减弱湍流参数。');
else
    disp('检查通过: 相位屏采样良好。');
end
