function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    % GET_OTOPS_SPECTRUM_HANDLE 计算 OTOPS 折射率谱的函数句柄
    % 基于 Yao et al., JOSA A, 2020.
    
    %% 1. 基础参数
    lambda_nm = lambda * 1e9; 
    k_wave = 2 * pi / lambda; 
    
    %% 2. 计算 OTOPS 系数 A 和 B
    % 注意：此处变量命名沿用 Wen et al. (2023) 的下标习惯
    % a5 对应 Yao 论文中的 a6 (盐度项), a6 对应 a7 (温度项)
    a1 = 1.779e-4;  a2 = -1.05e-6;  a3 = 1.6e-8; 
    a4 = -2.02e-6;  a5 = 1.155e-2;  a6 = -4.23e-3;
    
    A = a2*S + 2*a3*T*S + 2*a4*T + a6/lambda_nm; 
    B = a1 + a2*T + a3*T^2 + a5/lambda_nm; 
    
    %% 3. [修正] 流体热力学参数 (基于 Yao 2020 Appendix A)
    % 原始代码使用了固定值，这里改为动态计算以支持任意 T, S
    
    % 辅助变量
    T_k = T + 273.15; % 开尔文温度
    s_frac = S * 1e-3; % 盐度比例 (ppt -> ratio)
    
    % (1) 比热容 cp (Eq. A1-A2) [J/(kg K)]
    a11 = 5.328 - 9.76e-2*S + 4.04e-4*S^2;
    a12 = -6.913e-3 + 7.351e-4*S - 3.15e-6*S^2;
    a13 = 9.6e-6 - 1.927e-6*S + 8.23e-9*S^2;
    a14 = 2.5e-9 + 1.666e-9*S - 7.125e-12*S^2;
    cp = 1000 * (a11 + a12*T + a13*T^2 + a14*T^3); 
    
    % (2) 密度 rho (Eq. A8-A10) [kg/m^3]
    rho_T = 9.9992293295e2 + 2.0341179217e-2*T - 6.1624591598e-3*T^2 + ...
            2.2614664708e-5*T^3 - 4.6570659168e-8*T^4;
    rho_S = s_frac * (8.0200240891e2 - 2.0005183488*T + 1.6771024982e-2*T^2 - ...
            3.0600536746e-5*T^3 - 1.6132224742e-5*T*S);
    rho = rho_T + rho_S;
    
    % (3) 动力粘度 mu (Eq. A5-A7) [N s / m^2]
    mu_0 = (0.15700386464*(T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5;
    a21 = 1.5409136040 + 1.9981117208e-2*T - 9.5203865864e-5*T^2;
    a22 = 7.9739318223 - 7.561456881e-2*T + 4.7237011074e-4*T^2;
    mu = mu_0 * (1 + a21*s_frac + a22*s_frac^2);
    
    % (4) 热导率 sigma_T (Eq. A3-A4) [W/(m K)]
    T_b = 1.00024 * T;
    S_b = S / 1.00472;
    term1 = log10(240 + 0.0002*S_b);
    term2 = 0.434 * (2.3 - (343.5 + 0.037*S_b)/(T_b + 273.15));
    term3 = (1 - (T_b + 273.15)/(647.3 + 0.03*S_b))^(1/3);
    log_sigma = term1 - 3 + term2 * term3; % 修正公式结构
    % 注意：Yao论文Eq A3 写法比较复杂，这里简化处理或直接使用近似值
    % 考虑到公式实现的复杂性和排版歧义，如果对精度要求不极高，
    % 也可以保留 sigma_T = 0.6 或使用海水标准值，但 mu, rho, cp 建议用上面的公式
    sigma_T = 10^log_sigma; 
    % 备用方案：若上述公式报错，可回退到 sigma_T = 0.6;
    
    %% 4. 计算 Pr 和 Sc
    Pr = mu * cp / sigma_T;  
    Sc = mu^2 / (5.954e-15 * T_k * rho); 
    
    % Hill 参数 (Eq. 23)
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1);
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    %% 5. Eddy Diffusivity Ratio (d_r)
    % 注：严格来说 alpha_c 和 beta_c 也应通过 TEOS-10 计算
    % 但此处保留您的近似常数以维持代码独立性
    alpha_c = 2.6e-4; beta_c = 7.6e-4;  
    R_rho = alpha_c * abs(H_ratio) / beta_c;
    
    if R_rho >= 1
        d_r = R_rho + sqrt(R_rho)*sqrt(R_rho-1);
    elseif R_rho >= 0.5
        d_r = 1.85*R_rho - 0.85;
    else
        d_r = 0.15*R_rho; 
    end
    
    chi_S = chi_T * d_r / (H_ratio^2);
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    
    %% 6. [修正] 定义谱函数 (Hill Spectrum)
    % Eq. 20: 系数应为 (0.72 / 4pi)，而不是 0.033
    coeff_Hill = 0.72 / (4 * pi); % ≈ 0.0573
    
    Phi_Hill = @(K, chi_M, c_M) coeff_Hill * chi_M * epsilon^(-1/3) .* ...
        (K.^2 + 1e-25).^(-11/6) .* ... 
        exp(-176.90 * (K*eta).^2 .* c_M^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    
    % 返回总谱句柄
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + ...
                       B^2 * Phi_Hill(K, chi_S, c_S) + ...
                       2*A*B * Phi_Hill(K, chi_TS, c_TS));
end