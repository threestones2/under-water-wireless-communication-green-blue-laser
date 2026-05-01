function phase_screen = generate_phase_screen_SH(D, N, lambda, d, T_avg, S_avg, epsilon, chi_T, eta, H_ratio)
    % MPS_SH: 带次谐波补偿的海洋湍流相位屏生成
    % 基于 Schmidt 书中第 9 章的次谐波方法补偿低频
    %
    % 修改点:
    % 1. 修正了 lambda 单位问题 (nm 用于系数, m 用于波动)
    % 2. 将 OTOPS 谱封装为函数句柄
    % 3. 增加了次谐波 (Subharmonics) 叠加模块

    %% [1] 物理常数与基本参数
    lambda_nm = lambda * 1e9; % 纳米单位，用于经验公式
    k_wave = 2 * pi / lambda; % 波动光学波数 (m^-1)
    
    % 频率分辨率 (FFT 基础分辨率)
    dk = 2 * pi / D;          
    
    %% [2] 计算 OTOPS 模型系数 (Yao et al. 2020)
    % 经验系数
    a1 = 1.779e-4;  a2 = -1.05e-6;  a3 = 1.6e-8; 
    a4 = -2.02e-6;  a5 = 1.155e-2;  a6 = -4.23e-3;
    
    % 计算折射率系数 A 和 B (使用 nm 单位)
    % A = dn/dT
    A = a2*S_avg + 2*a3*T_avg*S_avg + 2*a4*T_avg + a6/lambda_nm; 
    % B = dn/dS (修正了原代码中误用 a4 的问题)
    B = a1 + a2*T_avg + a3*T_avg^2 + a5/lambda_nm;
    
    % 计算流体物理参数
    rho = 1025; mu = 1e-3; cp = 4182; sigma_T = 0.6;      
    Pr = mu * cp / sigma_T;  
    Sc = mu^2 / (5.954e-15 * (T_avg + 273.15) * rho); 
    
    % Hill 模型参数
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1);
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    % 涡流扩散比 d_r
    alpha_c = 2.6e-4; beta_c = 7.6e-4;  
    R_rho = alpha_c * abs(H_ratio) / beta_c;
    if R_rho >= 1, d_r = R_rho + sqrt(R_rho)*sqrt(R_rho-1);
    elseif R_rho >= 0.5, d_r = 1.85*R_rho - 0.85;
    else, d_r = 0.15*R_rho; end
    
    % 耗散率参数
    chi_S = chi_T * d_r / (H_ratio^2);
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 

    %% [3] 封装功率谱计算函数 (Phi_n)
    % 输入 K 是波数(标量或矩阵)，输出是折射率功率谱值
    % 这样既可以用在 FFT 大矩阵上，也可以用在次谐波的几个点上
    Phi_Hill = @(K, chi_M, c_M) 0.033 * chi_M * epsilon^(-1/3) .* ...
        (K.^2 + 1e-20).^(-11/6) .* ... % 加上极小值防止除0 (Batchelor谱通常是 -11/3? 这里需确认公式)
        exp(-176.90 * (K*eta).^2 .* c_M^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    % 注: 原代码写的是 K.^(-11/3)，这是 Phi_n(k)。
    % OTOPS 论文公式 19 是 Phi_n(k)。
    
    calc_Phi_n = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + ...
                       B^2 * Phi_Hill(K, chi_S, c_S) + ...
                       2*A*B * Phi_Hill(K, chi_TS, c_TS));

    %% [4] 生成高频相位屏 (FFT 部分)
    % 网格坐标
    kx = (-N/2 : N/2-1) * dk;   
    [KX, KY] = meshgrid(kx, kx);
    K_grid = sqrt(KX.^2 + KY.^2);
    K_grid(N/2+1, N/2+1) = 1e-10; % 处理直流分量
    
    % 计算谱密度
    Phi_n_val = calc_Phi_n(K_grid);
    Phi_n_val(N/2+1, N/2+1) = 0; % 强制去除 DC
    
    % 相位功率谱密度 F_phi (Eq. 12)
    % 注意: lambda 用米, F_phi = 2*pi^2 * k^2 * L * Phi_n
    F_phi = 2 * pi^2 * k_wave^2 * d * Phi_n_val;
    
    % 生成随机复数 (CN + i*CN)
    noise = (randn(N) + 1i * randn(N)) / sqrt(2);
    
    % 频域相位 (sqrt(F_phi) * dk 是为了获得正确的幅度量纲)
    H = noise .* sqrt(F_phi) * dk;
    
    % 逆变换得到高频部分
    phase_high = real(ift2(H, 1)); % 假设 ift2 已经包含了 ifftshift 和 N缩放

    %% [5] 次谐波补偿 (Subharmonics Compensation)
    % Schmidt 第 9 章方法
    % 在中心低频区域叠加更细密的频率分量
    
    phase_low = zeros(N, N);
    
    % 空间网格
    dx = D / N;
    x = (-N/2 : N/2-1) * dx;
    [xx, yy] = meshgrid(x, x);
    
    n_sub = 3; % 次谐波级数 (通常 3-10 级)
    
    for p = 1:n_sub
        % 当前级的频率步长 (呈 3 的指数衰减)
        dk_p = dk / (3^p); 
        
        % 在 [-1, 0, 1] * dk_p 范围内采样 3x3 网格
        % 除去中心 (0,0) 点
        for m = -1:1
            for n = -1:1
                if (m == 0 && n == 0)
                    continue; % 跳过直流
                end
                
                % 当前次谐波的波数
                kx_p = m * dk_p;
                ky_p = n * dk_p;
                k_p = sqrt(kx_p^2 + ky_p^2);
                
                % 计算该频率点的 PSD
                Phi_n_p = calc_Phi_n(k_p);
                F_phi_p = 2 * pi^2 * k_wave^2 * d * Phi_n_p;
                
                % 生成随机权重
                % 幅度 = sqrt(PSD * 面积)
                % 面积是 (dk_p)^2
                amp = sqrt(F_phi_p) * dk_p;
                
                % 生成随机复高斯噪声 (实部虚部独立)
                r_real = randn(1);
                r_imag = randn(1);
                c_mn = (r_real + 1i * r_imag) / sqrt(2);
                
                % 累加到空间域相位
                % exp(i * (kx*x + ky*y))
                phase_low = phase_low + real( c_mn * amp * exp(1i * (kx_p * xx + ky_p * yy)) );
            end
        end
    end
    
    %% [6] 合并结果
    phase_screen = phase_high + phase_low;

end

% 设置参数 (根据论文表1)
D = 1;          % 相位屏尺寸 (m)
N = 2^8;        % 采样点数
lambda = 532e-9;% 波长 (m)
d = 2;          % 相邻相位屏间距 (m)
T_avg = 20;     % 平均温度 (°C)
S_avg = 35;     % 平均盐度 (ppt)
epsilon = 1e-3; % 湍流动能耗散率 (弱湍流)
chi_T = 1e-7;   % 温度方差耗散率 (弱湍流)
eta = 1e-3;     % 湍流内尺度 (m)
H_ratio = -20;  % 温度盐度变化率 (°C/ppt)

% 生成相位屏
phase_screen = generate_phase_screen_SH(D, N, lambda, d, T_avg, S_avg, ...
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
