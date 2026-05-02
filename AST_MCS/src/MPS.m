function phase_screen = generate_phase_screen(D, N, lambda, d, T_avg, S_avg, epsilon, chi_T, eta, H_ratio)
    % 生成海洋湍流相位屏
    % 输入参数:
    %   D: 相位屏物理尺寸 (m)
    %   N: 采样点数
    %   lambda: 光波长 (m)
    %   d: 相邻相位屏间距 (m)
    %   T_avg: 平均温度 (°C)
    %   S_avg: 平均盐度 (ppt)
    %   epsilon: 湍流动能耗散率 (m²/s³)
    %   chi_T: 温度方差耗散率 (K²/s)
    %   eta: 湍流内尺度 (m)
    %   H_ratio: 温度盐度变化率 (°C/ppt)
    % 输出:
    %   phase_screen: 生成的相位屏 (N×N矩阵)

    % 1. 计算空间域和频率域网格
    dx = D / N;  % 空间域采样间隔
    x = (-N/2:N/2-1) * dx;  % 空间域坐标
    [X, Y] = meshgrid(x, x);
    
    dk = 2 * pi / D;  % 频率域采样间隔
    kx = (-N/2:N/2-1) * dk;  % 频率域坐标
    ky = kx;
    [KX, KY] = meshgrid(kx, ky);
    K = sqrt(KX.^2 + KY.^2);  % 径向频率
    
    % 2. 计算OTOPS模型参数 (公式4-9)
    % 系数A和B (公式5)
    a1 = 1.779e-4; a2 = -1.05e-6; a3 = 1.6e-8; 
    a4 = -2.02e-6; a5 = 1.155e-2; a6 = -4.23e-3;
    
    A = a2*S_avg + 2*a3*T_avg*S_avg + 2*a4*T_avg + a6/lambda;
    B = a1 + a2*T_avg + a4*T_avg^2 + a5/lambda;
    
    % 计算Prandtl数(Pr)和Schmidt数(Sc) (公式9)
    % 注意: 实际应用中应使用TEOs-10工具箱计算这些参数
    % 这里使用近似值 (实际应用中需替换为TEOs-10计算结果)
    mu = 1e-3;  % 动态粘度 (Pa·s)
    cp = 4000;  % 特征温度 (J/kg·K)
    sigma_T = 0.6;  % 热导率 (W/m·K)
    rho = 1025;  % 密度 (kg/m³)
    
    Pr = mu * cp / sigma_T;  % Prandtl数
    Sc = mu^2 / (5.954e-15 * (T_avg + 273.15) * rho);  % Schmidt数
    
    % 计算c_T, c_S, c_TS (公式8)
    c_T = 0.072 * (4/3) * 0.72 * Pr^(-1);
    c_S = 0.072 * (4/3) * 0.72 * Sc^(-1);
    c_TS = 0.072 * (4/3) * 0.72 * (Pr + Sc) / (2 * Pr * Sc);
    
    % 计算d_r (公式7)
    alpha_c = 2e-4;  % 热膨胀系数 (1/°C)
    beta_c = 7.6e-4;  % 盐收缩系数 (1/ppt)
    omega = alpha_c * H_ratio / beta_c;  % 温度-盐度扰动相关强度
    
    if abs(omega) >= 1
        d_r = abs(omega) + sqrt(abs(omega)) * (abs(omega) - 1)^0.5;
    elseif abs(omega) >= 0.5
        d_r = 1.85 * abs(omega) - 0.85;
    else
        d_r = 0.15 * abs(omega);
    end
    
    % 计算χ_S和χ_TS
    chi_S = chi_T * d_r / (H_ratio^2);
    %chi_TS = chi_T * (1 + d_r) / (2 * abs(H_ratio));%好像确实应该加绝对值
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio);%好像确实不应该加绝对值，A是负的
    
    % 3. 计算OTOPS功率谱 (公式4-6)
    % 定义计算Φ_M(κ)的函数
    Phi_M = @(chi_M, c_M) 0.72 * chi_M * epsilon^(-1/3) / (4*pi) .* ...
        K.^(-11/3) .* exp(-176.90 * (K*eta).^(2).*c_M^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    
    Phi_T = Phi_M(chi_T, c_T);
    Phi_S = Phi_M(chi_S, c_S);
    Phi_TS = Phi_M(chi_TS, c_TS);
    
    % 计算OTOPS模型Φₙ(κ) (公式4)
    Phi_n = A^2 * Phi_T + B^2 * Phi_S + 2*A*B*Phi_TS;
    
    % 处理κ=0处的奇点
    Phi_n(K == 0) = 0;
    
    % 4. 计算相位功率谱密度Fφ(κ) (公式12)
    k = 2 * pi / lambda;  % 波数
    F_phi = 2 * pi^2 * k^2 * d * Phi_n;
    
    % 5. 生成复高斯随机数矩阵
    real_part = randn(N) / sqrt(2);
    imag_part = randn(N) / sqrt(2);
    H_random = real_part + 1i * imag_part;
    
    % 6. 构建频率域相位屏 (公式11)
    H = H_random .* sqrt(F_phi) * dk * dk;
    
    % 7. 执行逆傅里叶变换得到空间域相位屏
    H_natural = ifftshift(H);  % 转换为自然顺序
    phase_screen = real(ifft2(H_natural));  % 执行IFFT并取实部
    phase_screen = fftshift(phase_screen);  % 重新排列使中心在原点
end


% 设置参数 (根据论文表1)
D = 1;          % 相位屏尺寸 (m)
N = 256;        % 采样点数
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
figure;
imagesc(phase_screen);
colormap('jet');
colorbar;
title('Generated Phase Screen');
xlabel('x (m)');
ylabel('y (m)');
axis square;
