%% 水下/大气光通信混合仿真: MC-MPS
% 
% 加权直射 (Weighted Ballistic) + PIS 散射 + 固定相位屏
% 统一使用 Petzold Clear Ocean 光学参数 (c=0.1514, a=0.114),清洁海水


clc; clear; close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- 模式选择(HG/TTHG/Mix/Empirical) ---
% 可选: 'HG', 'TTHG', 'Mix', 'Empirical'
param.phase_func = 'Empirical';  
param.n_max = 10;          

% --- 波长 ---
lambda_nm = 514; 
lambda = lambda_nm * 1e-9;
k_wave = 2*pi/lambda;

w0 = 0.1;                  % 束腰半径 (m)
div_angle = 3 * pi/180;           % 发散角 (rad)

% --- 3D 空间布局 ---
Tx_Pos = [0, 0, 0];         
Rx_Pos = [0, 10, 0];       % 传输距离
Link_Vec = Rx_Pos - Tx_Pos;
Link_Dist = norm(Link_Vec);
Link_Dir = Link_Vec / Link_Dist;

mu_T = Link_Dir;            
Rx_Normal = -Link_Dir;
Rx_Aperture = 0.2;          % 孔径 20cm
Rx_FOV = 10 * pi/180;       
Rx_Area = pi*(Rx_Aperture/2)^2;

% --- 介质参数 (统一使用 Petzold Clear Ocean) ---
param.c_water = 2.25e8;     % 水中光速 (m/s)

% Petzold清澈海水
% param.coef_c = 0.1514;                % 衰减系数 c (1/m)
% param.coef_a = 0.114;                 % 吸收系数 a (1/m)
% param.coef_b = 0.0374; % 散射系数 b (0.0374 1/m)

% %Petzold近岸海水
% param.coef_c = 0.398;
% param.coef_a = 0.179;
% param.coef_b=0.219;

% Petzold港口海水
param.coef_c = 2.190;
param.coef_a = 0.366;
param.coef_b=1.824;


param.albedo = param.coef_b / param.coef_c; % 单次散射反照率

% --- 相函数参数配置 ---
if strcmp(param.phase_func, 'HG')
    param.g_HG = 0.9707; 
elseif strcmp(param.phase_func, 'Mix')
    w_mie = 0.8440; 
    param.coef_kr = param.coef_b * (1 - w_mie); 
    param.coef_km = param.coef_b * w_mie;       
    param.gamma = 0.017; 
    param.g_mie = 0.9814; 
    param.f_mie = 0.49;    
    w_mie_ding = 0.8440;
elseif strcmp(param.phase_func, 'TTHG')
    param.alpha_TTHG = 0.4437;
    param.g1_TTHG = 0.9900;
    param.g2_TTHG = 0.8232;
elseif strcmp(param.phase_func, 'Empirical')
    param.c_e = param.coef_c; 
    param.a_e = param.coef_a;
    b_e = param.coef_b; 
    albedo_e = param.albedo;
    param.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    param.k1 = 1.188 - 0.688*albedo_e;
    param.k2 = 0.1 * (3.07 - 1.90*albedo_e);
    param.k3 = 0.01 * (4.58 - 3.02*albedo_e);
    param.k4 = 0.001 * (3.24 - 2.25*albedo_e);
    param.k5 = 0.0001 * (0.84 - 0.61*albedo_e);
    
    th_test = linspace(0, pi, 2000);
    val_test = zeros(size(th_test));
    for i=1:length(th_test)
        t_deg = th_test(i) * 180 / pi;
        if t_deg < 1e-6, t_deg = 1e-6; end
        term = 1 + (-1)^1*param.k1*t_deg^0.5 + (-1)^2*param.k2*t_deg^1.0 + ...
               (-1)^3*param.k3*t_deg^1.5 + (-1)^4*param.k4*t_deg^2.0 + ...
               (-1)^5*param.k5*t_deg^2.5;
        val_test(i) = exp(param.q_e * term);
    end
    param.b_emp_norm = 2 * pi * trapz(th_test, val_test .* sin(th_test));
end

% --- 湍流与相位屏 ---
T_avg = 20; 
S_avg = 35; 
epsilon = 1e-6; 
chi_T = 1e-8; 
eta = 1e-3; 
H_ratio = -20;
N_screens = 20;             
D_screen = 2.0;             
N_grid = 256;               
delta_z_screen = Link_Dist / N_screens; 

% --- [核心] 时间轴设置 ---
dt = 1e-10;   
t_min = dt;
t_max= 1e-6;                       
param.T_bins = t_min : dt : t_max;

N_bins = length(param.T_bins);
h_time = 1e-14*ones(1, N_bins); 

% --- 仿真控制 ---
N_packets = 1e5;         % 光子数

%% ================= 2. 预计算: 相位屏链 =================
fprintf('1. 生成相位屏链...\n');
n_vec = Link_Dir;
if abs(n_vec(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
u_vec = cross(up_temp, n_vec); u_vec = u_vec / norm(u_vec);
v_vec = cross(n_vec, u_vec);   v_vec = v_vec / norm(v_vec);

Screen_Chain = repmat(struct('Center', [0,0,0], 'Normal', [0,0,1], ...
                             'u_vec', [], 'v_vec', [], 'grad_x', [], 'grad_y', []), 1, N_screens);
dx = D_screen / N_grid;
x_axis = (-N_grid/2 : N_grid/2-1) * dx;
[Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);

for i = 1:N_screens
    Center_Pos = Tx_Pos + i * delta_z_screen * Link_Dir;
    phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
    [gx, gy] = gradient(phi, dx);
    Screen_Chain(i).Center = Center_Pos;
    Screen_Chain(i).Normal = Link_Dir;
    Screen_Chain(i).u_vec = u_vec; 
    Screen_Chain(i).v_vec = v_vec;
    Screen_Chain(i).grad_x = gx; 
    Screen_Chain(i).grad_y = gy;
end

%% ================= 3. 主仿真循环 =================
fprintf('3. 开始仿真 (已启用优化)...\n');
tic;

for p = 1:N_packets
    % 初始化光子
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos_local = r0*cos(phi0)*u_vec + r0*sin(phi0)*v_vec;
    pos_init = Tx_Pos + pos_local;
    
    theta_div_physics = lambda / (pi * w0); 
    U = theta_div_physics * sqrt(-0.5 * log(rand()));
    xi_psi = rand;
    psi_ini = 2 * pi * xi_psi;
    dir_init = rotate_direction(mu_T, U, psi_ini);
    
    weight_init = 1.0;
    
    % =========================================================
    % [分支 A]: 出射后直接打到接收平面 (Ballistic)
    % =========================================================
    pos = pos_init;
    dir = dir_init;

    [~, ~, hit_flag, path_len_ballistic] = ray_march_generic(pos, dir, 1e9, ...
        Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, true, ...
        Screen_Chain, k_wave, x_axis, dx, N_grid, ...
        Tx_Pos, Link_Dir, delta_z_screen); % <--- 新增参数

    if hit_flag
        % [修正] 能量计算: 不需要平方
        attenuation = exp(-param.coef_c * path_len_ballistic);
        energy_val = weight_init * attenuation; 

        t_arrival = path_len_ballistic / param.c_water;
        bin_idx = floor((t_arrival - t_min) / dt) + 1;
        if bin_idx >= 1 && bin_idx <= N_bins
            h_time(bin_idx) = h_time(bin_idx) + energy_val;
        end
    end
    
    % =========================================================
    % [分支 B]: 经历散射后才打到接收平面 (Scattering)
    % =========================================================
    pos = pos_init; 
    dir = dir_init;
    weight = weight_init;
    current_dist_traveled = 0;
    
    d_step = -log(rand()) / param.coef_c;
    
    [pos_new, dir_new, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
        Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, ...
        Screen_Chain, k_wave, x_axis, dx, N_grid, ...
        Tx_Pos, Link_Dir, delta_z_screen); % <--- 新增参数
    
    current_dist_traveled = current_dist_traveled + step_len;
    
    if dot(pos_new - Tx_Pos, Link_Dir) >= Link_Dist
        continue; 
    end
    
    pos = pos_new;
    dir = dir_new;
    
    % --- 高阶散射循环 ---
    for order = 1 : param.n_max
        vec_to_rx = Rx_Pos - pos;
        dist_to_rx = norm(vec_to_rx);
        dir_to_rx = vec_to_rx / dist_to_rx;
        
        % PIS 估计
        is_in_FOV_theta=acos(dot(-dir_to_rx, Rx_Normal));
        if is_in_FOV_theta <= Rx_FOV/2
            cos_theta_s = dot(dir, dir_to_rx);
            if strcmp(param.phase_func, 'HG')
                p_phase = pdf_HG(cos_theta_s, param.g_HG);
            elseif strcmp(param.phase_func, 'Mix')
                p_ray = pdf_Rayleigh(cos_theta_s, param.gamma);
                p_mie = pdf_Mie(cos_theta_s, param.g_mie, param.f_mie);
                p_phase = (param.coef_kr * p_ray + param.coef_km * p_mie) / param.coef_b;
            elseif strcmp(param.phase_func, 'TTHG')
                p_phase = pdf_TTHG(cos_theta_s, param.alpha_TTHG, param.g1_TTHG, param.g2_TTHG);
            elseif strcmp(param.phase_func, 'Empirical')
                p_phase = pdf_Empirical(cos_theta_s, param);
            end
            
            omega = Rx_Area / (dist_to_rx^2) * abs(dot(dir_to_rx, Rx_Normal));
            prob_survival = exp(-param.coef_c * dist_to_rx);
            
            energy_pis = weight * param.albedo * min(1,p_phase * omega) * prob_survival;
            t_pis_arrival = (current_dist_traveled + dist_to_rx) / param.c_water;
            
            bin_idx = floor((t_pis_arrival - t_min) / dt) + 1;
            if bin_idx >= 1 && bin_idx <= N_bins
                h_time(bin_idx) = h_time(bin_idx) + energy_pis;
            end
        end
        
        if order == param.n_max, break; end
        
        weight = weight * param.albedo;
        if weight < 1e-9, break; end
        
        % 重要性采样
        theta_i = pi * rand;
        phi = 2 * pi * rand;
        cos_t = cos(theta_i);
        
        if strcmp(param.phase_func, 'HG')
            p_val = pdf_HG(cos_t, param.g_HG);
        elseif strcmp(param.phase_func, 'Mix')
            p_ray = pdf_Rayleigh(cos_t, param.gamma);
            p_mie = pdf_Mie(cos_t, param.g_mie, param.f_mie);
            p_val = (param.coef_kr * p_ray + param.coef_km * p_mie) / param.coef_b;
        elseif strcmp(param.phase_func, 'TTHG')
            p_val = pdf_TTHG(cos_t, param.alpha_TTHG, param.g1_TTHG, param.g2_TTHG);
        elseif strcmp(param.phase_func, 'Empirical')
            p_val = pdf_Empirical(cos_t, param);
        end
        
        weight_factor = 2 * pi^2 * p_val * sin(theta_i);
        weight = weight * weight_factor;
        dir = rotate_direction(dir, theta_i, phi);
        
        d_step = -log(rand()) / param.coef_c;
        
        [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
            Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, ...
            Screen_Chain, k_wave, x_axis, dx, N_grid, ...
            Tx_Pos, Link_Dir, delta_z_screen); % <--- 新增参数
        
        current_dist_traveled = current_dist_traveled + step_len;
    end
    
    if mod(p, N_packets/10) == 0, fprintf('   进度: %.0f%%\n', p/N_packets*100); end
end
toc;

%% ================= 4. 结果处理 =================
h_time_norm = h_time / N_packets;
P_rx_total = sum(h_time_norm);
Path_Loss_dB = 10 * log10(P_rx_total);
T = param.T_bins; 
if P_rx_total > 0
    tau_mean = sum(T .* h_time_norm) / P_rx_total;
    tau_rms = sqrt( sum( ((T - tau_mean).^2) .* h_time_norm ) / P_rx_total );
else
    tau_mean = 0; tau_rms = 0;
end

fprintf('\n=== 仿真结果 (Petzold Clear Ocean) ===\n');
fprintf('总接收功率: %.4e (线性)\n', P_rx_total);
fprintf('总路径衰减: %.2f dB\n', Path_Loss_dB);
fprintf('平均时延:   %.4f ns\n', tau_mean * 1e9);
fprintf('RMS时延拓展: %.4f ns\n', tau_rms * 1e9);

figure;
semilogy(T, h_time_norm, 'b-', 'LineWidth', 0.5);
xlabel('Time (s)'); ylabel('Received Power (Linear)');
title(['Channel Impulse Response - RMS: ' num2str(tau_rms*1e9, '%.2f') 'ns']);
grid on; xlim([t_min, t_max]);

%第二种处理方法
figure;
dB_h_time_norm=10*log10(h_time_norm);
plot(T, dB_h_time_norm, 'b-', 'LineWidth', 0.5); % 加粗线条方便观察
xlabel('Time (s)'); 
ylabel('Received Power (Linear)');
title(['Channel Impulse Response (Zoomed: 400-500ns) - RMS: ' num2str(tau_rms*1e9, '%.2f') 'ns']);
grid on; 

% --- [关键修改] 将显示范围锁定在 400ns 到 500ns ---
xlim([400e-9, 700e-9]); 

%% ================= 5. 辅助函数 =================

% --- [优化后] O(1) 复杂度 Ray Marching ---
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, ...
    Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, ...
    Screen_Chain, k_wave, x_axis, dx, N_grid, ...
    Tx_Pos, Link_Dir, delta_z_screen) % <--- 新增几何参数

    hit_flag = false;
    total_len = 0;
    remaining_dist = dist_limit;
    
    N_screens = length(Screen_Chain);
    
    while remaining_dist > 1e-6
        min_dist = remaining_dist;
        event_type = 'none';
        target_idx = -1;
        
        % --- [优化核心]: 直接几何计算目标相位屏，替代循环 ---
        % 计算光子在链路方向上的投影距离
        dist_projected = dot(pos - Tx_Pos, Link_Dir);
        % 计算当前所在的区间索引 (floor向下取整)
        current_layer_idx = floor(dist_projected / delta_z_screen);
        
        target_screen_idx = -1;
        cos_theta = dot(dir, Link_Dir);
        
        % 根据飞行方向决定检测哪一个屏
        if cos_theta > 1e-6       % 向前飞: 检测下一屏
            target_screen_idx = current_layer_idx + 1;
        elseif cos_theta < -1e-6  % 向后飞: 检测当前屏(即身后的屏)
            target_screen_idx = current_layer_idx;
        end
        
        % 仅当目标索引有效时进行相交检测
        if target_screen_idx >= 1 && target_screen_idx <= N_screens
            scr = Screen_Chain(target_screen_idx);
            denom = dot(dir, scr.Normal);
            if abs(denom) > 1e-6
                t = dot(scr.Center - pos, scr.Normal) / denom;
                if t > 1e-6 && t < min_dist
                    min_dist = t;
                    event_type = 'screen';
                    target_idx = target_screen_idx;
                end
            end
        end
        % --------------------------------------------------
        
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx) > 1e-6
                t_rx = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
                if t_rx > 1e-6 && t_rx <= min_dist
                    min_dist = t_rx;
                    event_type = 'rx';
                end
            end
        end
        
        pos = pos + dir * min_dist;
        remaining_dist = remaining_dist - min_dist;
        total_len = total_len + min_dist; 
        
        if strcmp(event_type, 'rx')
            dist_to_center = norm(pos - Rx_Pos);
            angle_inc = acos(dot(-dir, Rx_Normal));
            if dist_to_center <= Rx_Aperture/2 && angle_inc <= Rx_FOV/2
                hit_flag = true;
                return;
            end
            
        elseif strcmp(event_type, 'screen')
            scr = Screen_Chain(target_idx);
            vec_on_plane = pos - scr.Center;
            loc_u = dot(vec_on_plane, scr.u_vec);
            loc_v = dot(vec_on_plane, scr.v_vec);
            idx_x = mod(round((loc_u - x_axis(1))/dx), N_grid) + 1;
            idx_y = mod(round((loc_v - x_axis(1))/dx), N_grid) + 1;
            g_u = scr.grad_x(idx_y, idx_x); 
            g_v = scr.grad_y(idx_y, idx_x);
            delta_vec = (g_u * scr.u_vec + g_v * scr.v_vec) / k_wave;
            dir = dir + delta_vec; 
            dir = dir / norm(dir);
        end
    end
end

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
    F_phi = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_val;
    
    % 生成复高斯噪声
    noise = (randn(N) + 1i * randn(N)) / sqrt(2);
    
    % 频域滤波
    C_nm = noise .* sqrt(F_phi) * dk;
    
    % 逆变换 (注意 ifft2 的缩放)
    %phase_high = real(ifft2(ifftshift(C_nm))) * N^2;
    phase_high = real(ifftshift(ifft2(ifftshift(C_nm))) )*N^2;
    %phase_high = real((ifft2(C_nm)))*N^2;
    
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

function new_dir = rotate_direction(dir, theta_s, psi_s)
    mu_x = dir(1); mu_y = dir(2); mu_z = dir(3);
    denom = sqrt(1 - mu_z^2);
    if denom < 1e-10
        if mu_z > 0
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), cos(theta_s)];
        else
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), -cos(theta_s)];
        end
    else
        sin_theta = sin(theta_s); cos_theta = cos(theta_s);
        cos_psi = cos(psi_s); sin_psi = sin(psi_s);
        new_dir_x = sin_theta/denom * (mu_x*mu_z*cos_psi - mu_y*sin_psi) + mu_x*cos_theta;
        new_dir_y = sin_theta/denom * (mu_y*mu_z*cos_psi + mu_x*sin_psi) + mu_y*cos_theta;
        new_dir_z = -sin_theta*cos_psi*denom + mu_z*cos_theta;
        new_dir = [new_dir_x, new_dir_y, new_dir_z];
    end
    new_dir = new_dir / norm(new_dir);
end

% --- PDF 函数 (1/sr) ---
function p = pdf_HG(cos_theta, g)
    p = (1 - g^2) ./ (4*pi * (1 + g^2 - 2*g*cos_theta).^1.5);
end

function p = pdf_Rayleigh(cos_theta, gamma)
    p = 3 * (1 + 3*gamma + (1-gamma)*cos_theta^2) / (16 * pi * (1 + 2*gamma));
end

function p = pdf_Mie(cos_theta, g, f)
    val = (1 - g^2)/2 * ( 1./(1 + g^2 - 2*g*cos_theta).^1.5 + f * 0.5 * (3*cos_theta^2 - 1) / (1 + g^2)^1.5 );
    p = val / (2*pi);
end

function p = pdf_TTHG(cos_theta, alpha, g1, g2)
    p1 = pdf_HG(cos_theta, g1);
    p2 = pdf_HG(cos_theta, g2);
    p = alpha * p1 + (1 - alpha) * p2;
end

function p = pdf_Empirical(cos_theta, param)
    cos_theta = max(min(cos_theta, 1), -1);
    theta_rad = acos(cos_theta);
    t_deg = theta_rad * 180 / pi;
    if t_deg < 1e-6, t_deg = 1e-6; end
    
    term = 1 + (-1)^1*param.k1*t_deg^0.5 + ...
           (-1)^2*param.k2*t_deg^1.0 + ...
           (-1)^3*param.k3*t_deg^1.5 + ...
           (-1)^4*param.k4*t_deg^2.0 + ...
           (-1)^5*param.k5*t_deg^2.5;
       
    VSF = exp(param.q_e * term);
    p = VSF / param.b_emp_norm;
end


