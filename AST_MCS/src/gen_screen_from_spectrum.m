function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    % GEN_SCREEN_FROM_SPECTRUM 根据给定的谱函数生成相位屏
    %
    % 输入:
    %   Phi_n_func - 折射率谱函数句柄 @(K)
    %   D          - 相位屏边长 (m)
    %   N          - 网格点数
    %   k_wave     - 光波数 (rad/m)
    %   delta_z    - 该相位屏代表的传输距离/厚度 (m)
    
    dk = 2 * pi / D;
    dx = D / N;
    
    %% 1. 高频部分 (FFT)
    kx = (-N/2 : N/2-1) * dk;
    [KX, KY] = meshgrid(kx, kx);
    K_grid = sqrt(KX.^2 + KY.^2);
    K_grid(N/2+1, N/2+1) = 1e-10; % 避开奇异点
    
    % 调用传入的谱函数计算 Phi_n
    Phi_n_val = Phi_n_func(K_grid);
    Phi_n_val(N/2+1, N/2+1) = 0; 
    
    % 转换为相位谱 F_phi = 2 * pi^2 * k^2 * dz * Phi_n
    F_phi = 2 * pi * k_wave^2 * delta_z * Phi_n_val;
    
    % 生成复高斯噪声
    noise = (randn(N) + 1i * randn(N));
    
    % 频域滤波
    C_nm = noise .* sqrt(F_phi) * dk;
    
    % 逆变换 (注意 ifft2 的缩放)
    phase_high = real(ifftshift(ifft2(ifftshift(C_nm))) )*(N^2);

    
    %% 2. 次谐波补偿 (Subharmonics)
    phase_low = zeros(N, N);
    [xx, yy] = meshgrid((-N/2 : N/2-1) * dx);
    n_sub = 3; % 补偿级数
    
    for p = 1:n_sub
        dk_p = dk / (3^p); % 频率步长指数衰减
        for m = -1:1
            for n = -1:1
                if (m==0 && n==0), continue; end
                
                kx_p = m * dk_p;
                ky_p = n * dk_p;
                k_p = sqrt(kx_p^2 + ky_p^2);
                
                % 计算低频点的谱值
                Phi_n_p = Phi_n_func(k_p);
                F_phi_p = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_p;
                
                amp = sqrt(F_phi_p) * dk_p;
                r_c = (randn(1) + 1i * randn(1)) / sqrt(2);
                
                phase_low = phase_low + real( r_c * amp * exp(1i * (kx_p * xx + ky_p * yy)) );
            end
        end
    end
    
    %% 3. 合成
    phase_screen = phase_high + phase_low;
end