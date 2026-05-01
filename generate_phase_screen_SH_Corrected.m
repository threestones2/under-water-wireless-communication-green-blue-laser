
function phase_screen = generate_phase_screen_SH_Corrected(D, N, lambda, d, T_avg, S_avg, epsilon, chi_T, eta, H_ratio)
    % GENERATE_PHASE_SCREEN_SH_CORRECTED 生成基于 OTOPS 模型和次谐波补偿的海洋湍流相位屏
    %
    %   本函数利用功率谱反演法 (FFT) 结合次谐波补偿技术，模拟光在海洋湍流中传输时累积的相位畸变。
    %   物理模型基于 Yao et al. (2020) 提出的 OTOPS (Optical Turbulence of Ocean Power Spectrum) 模型，
    %   并针对 FFT 方法固有的低频采样不足问题（会导致低估光束漂移），引入了 Schmidt (2010) 
    %   第9章所述的次谐波补偿 (Subharmonics Compensation) 算法。
    %
    %   主要修正与特性：
    %   1. 量纲修正：修复了原论文/代码中的单位错误。在计算经验系数 (A, B) 时正确将波长转换为
    %      纳米 (nm)，而在计算波动光学参数（波数 k）时保持为米 (m)。
    %   2. 谱模型：内核采用 Hill 谱的解析近似形式，准确描述高频耗散区 (Viscous-Convective Range)
    %      的“凸起”特征，优于简单的 Kolmogorov 谱。
    %   3. 次谐波补偿：通过叠加多级低频正弦分量，修复了低频能量缺失，确保长距离传输仿真
    %      中的光束漂移 (Beam Wander) 统计特性正确。
    %   4. 幅度缩放：显式处理 IFFT 的缩放因子，确保输出相位屏具有正确的物理量纲 (rad)。
    %
    %   输入参数:
    %       D       - 相位屏物理边长 (单位: m)
    %       N       - 网格采样点数 (建议为 2 的幂，如 256, 512)
    %       lambda  - 光波长 (单位: m，例如 532e-9)
    %       d       - 相位屏代表的传输路径长度 / 切片厚度 (单位: m)
    %       T_avg   - 海水平均温度 (单位: °C)
    %       S_avg   - 海水平均盐度 (单位: ppt)
    %       epsilon - 湍流动能耗散率 (单位: m²/s³，典型值 1e-1 ~ 1e-10)
    %       chi_T   - 温度方差耗散率 (单位: K²/s，典型值 1e-4 ~ 1e-10)
    %       eta     - Kolmogorov 内尺度 (单位: m，典型值 1e-3)
    %       H_ratio - 温盐度梯度比 (无量纲，通常为负值，例如 -20)
    %
    %   输出参数:
    %       phase_screen - 生成的二维随机相位屏矩阵 (N x N)，单位: 弧度 (rad)
    %
    %   参考文献:
    %       [1] Yao, J., et al. "Spatial power spectrum of natural water turbulence...", JOSA A (2020).
    %       [2] Schmidt, J. D. "Numerical Simulation of Optical Wave Propagation", SPIE (2010), Chap 9.
    %       [3] Wen, H., et al. "Modeling and performance analysis...", Applied Optics (2023).
    %
    %   Usage Example:
    %       scr = generate_phase_screen_SH_Corrected(1, 256, 532e-9, 5, 20, 35, 1e-3, 1e-7, 1e-3, -20);
    %       imagesc(scr); axis image; colormap jet; colorbar;
    %
    % [版本信息]
    % 1. 修正 lambda 单位 (nm for coeff, m for wave)
    % 2. 添加 Schmidt 第9章 次谐波补偿 (Subharmonics)
    % 3. 显式使用 ifft2 避免缩放混淆

    %% [1] 基础参数
    lambda_nm = lambda * 1e9; % 纳米
    k_wave = 2 * pi / lambda; % 米制波数
    dk = 2 * pi / D;          % 频率分辨率 (rad/m)
    dx = D / N;               % 空间分辨率 (m)
    
    %% [2] 计算 OTOPS 参数 (Yao et al. 2020)
    % 系数
    a1 = 1.779e-4;  a2 = -1.05e-6;  a3 = 1.6e-8; 
    a4 = -2.02e-6;  a5 = 1.155e-2;  a6 = -4.23e-3;
    
    % A, B (使用 nm)
    A = a2*S_avg + 2*a3*T_avg*S_avg + 2*a4*T_avg + a6/lambda_nm; 
    B = a1 + a2*T_avg + a3*T_avg^2 + a5/lambda_nm; % 修正 a4->a3
    
    % 流体参数
    rho = 1025; mu = 1e-3; cp = 4182; sigma_T = 0.6;      
    Pr = mu * cp / sigma_T;  
    Sc = mu^2 / (5.954e-15 * (T_avg + 273.15) * rho); 
    
    % Hill 参数
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1);
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    % d_r 参数
    alpha_c = 2.6e-4; beta_c = 7.6e-4;  
    R_rho = alpha_c * abs(H_ratio) / beta_c;
    if R_rho >= 1, d_r = R_rho + sqrt(R_rho)*sqrt(R_rho-1);
    elseif R_rho >= 0.5, d_r = 1.85*R_rho - 0.85;
    else, d_r = 0.15*R_rho; end
    
    chi_S = chi_T * d_r / (H_ratio^2);
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 

    %% [3] 定义谱函数句柄
    % Hill 谱 (修正了 Batchelor 谱的指数)
    % 注意: 这里 K.^2 对应 -11/6 是为了方便后面计算 Phi_n ~ K^-11/3
    % 实际上 Phi_n = 0.033 * chi * eps^-1/3 * K^-11/3 ...
    Phi_Hill = @(K, chi_M, c_M) 0.033 * chi_M * epsilon^(-1/3) .* ...
        (K.^2 + 1e-25).^(-11/6) .* ... 
        exp(-176.90 * (K*eta).^2 .* c_M^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    
    calc_Phi_n = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + ...
                       B^2 * Phi_Hill(K, chi_S, c_S) + ...
                       2*A*B * Phi_Hill(K, chi_TS, c_TS));

    %% [4] 高频相位屏 (High-frequency Part)
    % 频率网格
    kx = (-N/2 : N/2-1) * dk;
    [KX, KY] = meshgrid(kx, kx);
    K_grid = sqrt(KX.^2 + KY.^2);
    K_grid(N/2+1, N/2+1) = 1e-10; % 避开 0 点
    
    % 计算功率谱 Phi_n (折射率涨落谱)
    Phi_n_val = calc_Phi_n(K_grid);
    Phi_n_val(N/2+1, N/2+1) = 0; % 强制 DC 为 0
    
    % 计算相位功率谱 F_phi (积分厚度 d)
    % F_phi = 2 * pi^2 * k^2 * d * Phi_n
    F_phi = 2 * pi^2 * k_wave^2 * d * Phi_n_val;
    
    % 生成频域随机复数 (均值为0，方差为1的复高斯)
    noise = (randn(N) + 1i * randn(N)) / sqrt(2);
    
    % 构造频域场 C_nm
    % 幅度 = sqrt(PSD * 面积) = sqrt(F_phi * dk^2)
    % 这一步的物理量纲是 [rad]
    C_nm = noise .* sqrt(F_phi) * dk;
    
    % 逆变换到空间域
    % Schmidt Eq 9.53: phi(x) = sum [ C_nm * exp(i k x) ]
    % Matlab ifft2 计算的是: sum [ X ] / N^2
    % 所以我们需要乘以 N^2 来抵消 Matlab 的系数
    phase_high = real(ifft2(ifftshift(C_nm))) * N^2;

    %% [5] 次谐波补偿 (Subharmonics Part)
    phase_low = zeros(N, N);
    [xx, yy] = meshgrid((-N/2 : N/2-1) * dx); % 空间坐标
    
    n_sub = 3; % 级数
    
    for p = 1:n_sub
        dk_p = dk / (3^p); % 频率步长指数衰减
        
        for m = -1:1
            for n = -1:1
                if (m==0 && n==0), continue; end % 跳过直流
                
                kx_p = m * dk_p;
                ky_p = n * dk_p;
                k_p = sqrt(kx_p^2 + ky_p^2);
                
                % 计算该低频点的 PSD
                Phi_n_p = calc_Phi_n(k_p);
                F_phi_p = 2 * pi^2 * k_wave^2 * d * Phi_n_p;
                
                % 生成随机系数
                % 幅度 = sqrt(PSD * 面积)
                amp = sqrt(F_phi_p) * dk_p;
                r_real = randn(1); r_imag = randn(1);
                c_sub = (r_real + 1i * r_imag) / sqrt(2);
                
                % 累加低频相位 (显式求和)
                phase_low = phase_low + real( c_sub * amp * exp(1i * (kx_p * xx + ky_p * yy)) );
            end
        end
    end
    
    %% [6] 总相位屏
    phase_screen = phase_high + phase_low;
end