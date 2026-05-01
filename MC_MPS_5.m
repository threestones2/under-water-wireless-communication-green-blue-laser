%% 水下/大气光通信混合仿真: MC-MPS
% 
% 逻辑: 加权直射 (Weighted Ballistic) + PIS 散射 + 固定相位屏
%
% 该代码采用的是HG相函数，用反变换法求散射角（HG有解析反解），本来支持切换到Ding的混合相函数模型，但相函数部分还没写

clc; clear; close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- 模式选择 ---
param.phase_func = 'HG';  
param.n_max = 3;          

% --- 光源与几何 ---
lambda_nm = 532; 
lambda = lambda_nm * 1e-9;
k_wave = 2*pi/lambda;

w0 = 0.01;                  % 束腰半径 (m)
div_angle = 1e-3;          % 发散角 (rad)

% --- 3D 空间布局 ---
Tx_Pos = [0, 0, 0];         
Rx_Pos = [0, 50, 0];        % 传输距离 50m
Link_Vec = Rx_Pos - Tx_Pos;
Link_Dist = norm(Link_Vec);
Link_Dir = Link_Vec / Link_Dist;

mu_T = Link_Dir;            
Rx_Normal = -Link_Dir;      
Rx_Aperture = 0.2;          % 孔径 20cm
Rx_FOV = 30 * pi/180;       
Rx_Area = pi*(Rx_Aperture/2)^2;

% --- 介质参数 ---
param.c_water = 2.25e8;     % [关键] 水中光速 (m/s)
param.coef_a = 0.1;         
if strcmp(param.phase_func, 'HG')
    param.coef_b = 0.2; 
    param.g_HG = 0.924; 
elseif strcmp(param.phase_func, 'Mix')
    param.coef_kr = 0.05; 
    param.coef_km = 0.169;
    param.coef_b = param.coef_kr + param.coef_km;
    param.gamma = 0.017; 
    param.g_mie = 0.72; 
    param.f_mie = 0.5;      
end
param.coef_c = param.coef_a + param.coef_b; 
param.albedo = param.coef_b / param.coef_c; 

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
dt = 1e-9;   
t_min = dt;
t_max= 1e-6;                       
param.T_bins = t_min : dt : t_max;

N_bins = length(param.T_bins);
h_time = zeros(1, N_bins); % [线性能量] 存放 Impulse Response

% --- 仿真控制 ---
N_packets = 1e4;         % 光子数

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
fprintf('3. 开始仿真...\n');
tic;

for p = 1:N_packets
    % 初始化
    % 高斯光束，光子几何位置按照能量分布确定，高斯光束的公式是振幅，振幅平方才是光强，所以有0.5
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos_local = r0*cos(phi0)*u_vec + r0*sin(phi0)*v_vec;
    pos_init = Tx_Pos + pos_local;
    % 光子位置确定后
    % 确定物理发散角 (基于束腰和波长)
    theta_div_physics = lambda / (pi * w0); % 半角 (1/e^2)
    % 采样发射角 (高斯分布)
    % 注意：这里使用的是瑞利分布生成 radial angle (theta)，
    % 对应于正态分布的 x, y 分量合成
    U = theta_div_physics * sqrt(-0.5 * log(rand()));
    xi_psi = rand;
    psi_ini = 2 * pi * xi_psi;
    dir_init = rotate_direction(mu_T, U, psi_ini);
    
    weight_init = 1.0;
    
    % =========================================================
    % [分支 A]: 出射后直接打到接收平面
    % =========================================================
    pos = pos_init;
    dir = dir_init;
    
    % 强制穿过所有屏，返回实际飞行路径长度
    [~, ~, hit_flag, path_len_ballistic] = ray_march_generic(pos, dir, 1e9, ...
        Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, true, ...
        Screen_Chain, k_wave, x_axis, dx, N_grid);
    
    if hit_flag
        % 1. 计算能量 (Beer-Lambert 衰减),因为情况A没有对距离采样,所以
        attenuation = exp(-param.coef_c * path_len_ballistic);
        energy_val = weight_init * attenuation;
        
        % 2. 计算到达时间
        t_arrival = path_len_ballistic / param.c_water;
        
        % 3. 放入时间桶 (线性累加)
        bin_idx = floor((t_arrival - t_min) / dt) + 1;
        if bin_idx >= 1 && bin_idx <= N_bins
            h_time(bin_idx) = h_time(bin_idx) + energy_val;
        end
    end
    
    % =========================================================
    % [分支 B]: 经历散射后才打到接收平面
    % =========================================================
    pos = pos_init; 
    dir = dir_init;
    weight = weight_init;
    
    % [新增] 记录光子已飞行的总距离
    current_dist_traveled = 0;
    
    % 0阶飞行
    d_step = -log(rand()) / param.coef_c;
    
    [pos_new, dir_new, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
        Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, ...
        Screen_Chain, k_wave, x_axis, dx, N_grid);
    
    % 更新已飞距离
    current_dist_traveled = current_dist_traveled + step_len;
    
    % 检查是否直接穿越接收屏的平面 (如果是，跳过分支B，因为分支A已算)
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
        
        % 如果发射或者几次散射后，还在接收机视场角内，直接往接收机平面打
        if acos(dot(-dir_to_rx, Rx_Normal)) <= Rx_FOV
            cos_theta_s = dot(dir, dir_to_rx);
            if strcmp(param.phase_func, 'HG')
                p_phase = pdf_HG(cos_theta_s, param.g_HG);
            else
                p_ray = pdf_Rayleigh(cos_theta_s, param.gamma);
                p_mie = pdf_Mie(cos_theta_s, param.g_mie, param.f_mie);
                p_phase = (param.coef_kr * p_ray + param.coef_km * p_mie) / param.coef_b;
            end
            
            omega = Rx_Area / (dist_to_rx^2) * abs(dot(dir_to_rx, Rx_Normal));
            prob_survival = exp(-param.coef_c * dist_to_rx);
            
            % 1. 计算 PIS 能量 (线性)
            energy_pis = weight * param.albedo * p_phase * omega * prob_survival;
            
            % 2. 计算 PIS 总到达时间
            % 时间 = (光子历史飞行距离 + 散射点到Rx的虚拟距离) / c
            t_pis_arrival = (current_dist_traveled + dist_to_rx) / param.c_water;
            
            % 3. 放入时间桶
            bin_idx = floor((t_pis_arrival - t_min) / dt) + 1;
            if bin_idx >= 1 && bin_idx <= N_bins
                h_time(bin_idx) = h_time(bin_idx) + energy_pis;
            end
        end
        
        if order == param.n_max, break; end
        
        weight = weight * param.albedo;
        if weight < 1e-9, break; end
        
        if strcmp(param.phase_func, 'HG')
            dir = sample_dir_HG(dir, param.g_HG);
        else
            dir = sample_dir_HG(dir, 0.9);% 还没写别的相函数，这里也得改，以后不用反变换法求角度了
        end
        
        d_step = -log(rand()) / param.coef_c;
        
        % 飞向下一个点
        [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
            Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, ...
            Screen_Chain, k_wave, x_axis, dx, N_grid);
        
        % 更新历史距离
        current_dist_traveled = current_dist_traveled + step_len;
    end
    
    if mod(p, N_packets/10) == 0, fprintf('   进度: %.0f%%\n', p/N_packets*100); end
end
toc;

%% ================= 4. 结果处理: 时延拓展与总损耗 =================
% h_time 现在包含了所有光子的能量分布 (Impulse Response)
% 归一化 (Power per packet)
h_time_norm = h_time / N_packets;

% 1. 计算总接收功率 (线性求和)
P_rx_total = sum(h_time_norm);

% 2. 转换为 dB (总路径衰减)
Path_Loss_dB = 10 * log10(P_rx_total);

% 3. 计算时延拓展 (Delay Spread) - 使用线性 h(t)
% 时间轴
T = param.T_bins; 
% 平均时延 (Mean Delay)
if P_rx_total > 0
    tau_mean = sum(T .* h_time_norm) / P_rx_total;
    % RMS 时延 (RMS Delay Spread)
    tau_rms = sqrt( sum( ((T - tau_mean).^2) .* h_time_norm ) / P_rx_total );
else
    tau_mean = 0; tau_rms = 0;
end

fprintf('\n=== 仿真结果 ===\n');
fprintf('总接收功率: %.4e (线性)\n', P_rx_total);
fprintf('总路径衰减: %.2f dB\n', Path_Loss_dB);
fprintf('平均时延:   %.4f ns\n', tau_mean * 1e9);
fprintf('RMS时延拓展: %.4f ns\n', tau_rms * 1e9);

% 4. 绘图 (虽然计算用线性，但绘图通常看 dB 或者线性归一化)
figure;

plot(T, h_time_norm, 'b-', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Received Power (Linear)');
title(['Channel Impulse Response (Linear) - RMS: ' num2str(tau_rms*1e9, '%.2f') 'ns']);
grid on; xlim([t_min, t_max]);
% figure;
% plot(T*1e9, 10*log10(h_time_norm + 1e-30), 'r-', 'LineWidth', 1.5);
% xlabel('Time (ns)'); ylabel('Power (dB)');
% title('Channel Impulse Response (Log Scale)');
% grid on; xlim([t_min*1e9, t_min*1e9 + 10]);

%% ================= 5. 核心追踪函数 =================
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, ...
    Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, ...
    Screen_Chain, k_wave, x_axis, dx, N_grid)

    hit_flag = false;
    total_len = 0;
    remaining_dist = dist_limit;
    
    while remaining_dist > 1e-6
        min_dist = remaining_dist;
        event_type = 'none';
        target_idx = -1;
        
        % A. 检查相位屏
        for i = 1:length(Screen_Chain)
            denom = dot(dir, Screen_Chain(i).Normal);
            if abs(denom) > 1e-6
                t = dot(Screen_Chain(i).Center - pos, Screen_Chain(i).Normal) / denom;
                if t > 1e-6 && t < min_dist
                    min_dist = t;
                    event_type = 'screen';
                    target_idx = i;
                end
            end
        end
        
        % B. 检查接收机
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
        
        % 移动
        pos = pos + dir * min_dist;
        remaining_dist = remaining_dist - min_dist;
        total_len = total_len + min_dist; % 累加路径长度
        
        % 事件
        if strcmp(event_type, 'rx')
            dist_to_center = norm(pos - Rx_Pos);
            angle_inc = acos(dot(-dir, Rx_Normal));
            if dist_to_center <= Rx_Aperture/2 && angle_inc <= Rx_FOV
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
%%
function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    lambda_nm = lambda * 1e9; 
    k_wave = 2 * pi / lambda; 
    
    a1 = 1.779e-4;  a2 = -1.05e-6;  a3 = 1.6e-8; 
    a4 = -2.02e-6;  a5 = 1.155e-2;  a6 = -4.23e-3;
    A = a2*S + 2*a3*T*S + 2*a4*T + a6/lambda_nm; 
    B = a1 + a2*T + a3*T^2 + a5/lambda_nm; 
    
    rho = 1025; mu = 1e-3; cp = 4182; sigma_T = 0.6;      
    Pr = mu * cp / sigma_T;  
    Sc = mu^2 / (5.954e-15 * (T + 273.15) * rho); 
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1);
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    alpha_c = 2.6e-4; beta_c = 7.6e-4;  
    R_rho = alpha_c * abs(H_ratio) / beta_c;
    if R_rho >= 1, d_r = R_rho + sqrt(R_rho)*sqrt(R_rho-1);
    elseif R_rho >= 0.5, d_r = 1.85*R_rho - 0.85;
    else, d_r = 0.15*R_rho; end
    
    chi_S = chi_T * d_r / (H_ratio^2);
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    
    Phi_Hill = @(K, chi_M, c_M) 0.033 * chi_M * epsilon^(-1/3) .* ...
        (K.^2 + 1e-25).^(-11/6) .* ... 
        exp(-176.90 * (K*eta).^2 .* c_M^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + ...
                       B^2 * Phi_Hill(K, chi_S, c_S) + ...
                       2*A*B * Phi_Hill(K, chi_TS, c_TS));
end
%%
function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2 * pi / D;
    dx = D / N;
    kx = (-N/2 : N/2-1) * dk;
    [KX, KY] = meshgrid(kx, kx);
    K_grid = sqrt(KX.^2 + KY.^2);
    K_grid(N/2+1, N/2+1) = 1e-10; 
    Phi_n_val = Phi_n_func(K_grid);
    Phi_n_val(N/2+1, N/2+1) = 0; 
    F_phi = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_val;
    noise = (randn(N) + 1i * randn(N)) / sqrt(2);
    C_nm = noise .* sqrt(F_phi) * dk;
    phase_high = real(ifft2(ifftshift(C_nm))) * N^2;
    phase_low = zeros(N, N);
    [xx, yy] = meshgrid((-N/2 : N/2-1) * dx);
    n_sub = 3; 
    for p = 1:n_sub
        dk_p = dk / (3^p); 
        for m = -1:1
            for n = -1:1
                if (m==0 && n==0), continue; end
                kx_p = m * dk_p; ky_p = n * dk_p;
                k_p = sqrt(kx_p^2 + ky_p^2);
                Phi_n_p = Phi_n_func(k_p);
                F_phi_p = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_p;
                amp = sqrt(F_phi_p) * dk_p;
                r_c = (randn(1) + 1i * randn(1)) / sqrt(2);
                phase_low = phase_low + real( r_c * amp * exp(1i * (kx_p * xx + ky_p * yy)) );
            end
        end
    end
    phase_screen = phase_high + phase_low;
end
%%
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
%%
function p = pdf_HG(cos_theta, g)
    p = (1 - g^2) ./ (4*pi * (1 + g^2 - 2*g*cos_theta).^1.5);
end
%%
function p = pdf_Rayleigh(cos_theta, gamma)
    p = 3 * (1 + 3*gamma + (1-gamma)*cos_theta^2) / (16 * pi * (1 + 2*gamma));
end
%%
function p = pdf_Mie(cos_theta, g, f)
    val = (1 - g^2)/2 * ( 1./(1 + g^2 - 2*g*cos_theta).^1.5 + f * 0.5 * (3*cos_theta^2 - 1) / (1 + g^2)^1.5 );
    p = val / (2*pi);
end
%% 当采用HG函数，用反变换法求得天顶角,HG函数好的地方在于有解析反解形式，不像别的得用数值反解方法，慢的一批，但我待会得换成TTHG，这个拟合程度更高
function dir_new = sample_dir_HG(dir, g)
    rand_val = rand();
    if abs(g) < 1e-6, cos_t = 2*rand_val - 1;
    else, temp = (1 - g^2) / (1 - g + 2*g*rand_val); cos_t = (1 + g^2 - temp^2) / (2*g); end
    cos_t = max(min(cos_t, 1), -1); 
    theta = acos(cos_t); 
    phi = 2*pi*rand();
    dir_new = rotate_direction(dir, theta, phi);
end