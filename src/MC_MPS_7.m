%% 水下/大气光通信混合仿真: MC-MPS (Final Version)
% 这个我是想换成反变换法求散射角

clc; clear; close all;

%% ================= 1. 参数初始化 =================
param = struct();

% --- [关键] 模式选择 ---
% 可选: 'Haltrin', 'TTHG', 'Ding', 'HG', 'Petzold'
param.phase_func = 'Haltrin';  
param.n_max = 3;          

% --- 光源与几何 ---
lambda_nm = 532; 
lambda = lambda_nm * 1e-9;
k_wave = 2*pi/lambda;

w0 = 0.01;                  % 束腰半径 (m)
div_angle = 1e-3;           % 发散角 (rad)

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

% --- 介质参数 (通用) ---
param.c_water = 2.25e8;     % 水中光速 (m/s)
param.coef_a = 0.114;       % 吸收系数 (参考 Clear Ocean a=0.114)
param.coef_b = 0.1514 - 0.114; % 默认散射系数 (会被特定模型覆盖)

% ================= 相函数参数配置 (基于 scatter_test3) =================

if strcmp(param.phase_func, 'Haltrin')
    % Haltrin 经验公式参数 (Clear Ocean)
    param.c_e = 0.1514;
    param.a_e = 0.114;
    param.b_e = param.c_e - param.a_e;
    param.coef_b = param.b_e; % 更新仿真用的散射系数
    
    % 计算 Haltrin 系数
    albedo_e = param.b_e / param.c_e;
    param.q_e = 2.598 + 17.748*sqrt(param.b_e) - 16.722*param.b_e + 5.932*param.b_e*sqrt(param.b_e);
    param.k1 = 1.188 - 0.688*albedo_e;
    param.k2 = 0.1 * (3.07 - 1.90*albedo_e);
    param.k3 = 0.01 * (4.58 - 3.02*albedo_e);
    param.k4 = 0.001 * (3.24 - 2.25*albedo_e);
    param.k5 = 0.0001 * (0.84 - 0.61*albedo_e);
    
    % [重要] 预计算归一化因子 b_calc
    % 因为 Haltrin 公式给出的是 VSF，我们需要积分得到散射系数 b，
    % 然后用 P = 4pi * VSF / b 来归一化相函数
    fprintf('正在初始化 Haltrin 模型 (积分归一化)...\n');
    th_integ = linspace(0, pi, 2000); % 积分网格
    th_deg = th_integ * 180 / pi;
    term = 1 + (-1)^1*param.k1*th_deg.^0.5 + (-1)^2*param.k2*th_deg.^1.0 + ...
           (-1)^3*param.k3*th_deg.^1.5 + (-1)^4*param.k4*th_deg.^2.0 + ...
           (-1)^5*param.k5*th_deg.^2.5;
    vsf_val = exp(param.q_e * term);
    % b = 2pi * int(VSF * sin(theta))
    param.b_haltrin_norm = 2 * pi * trapz(th_integ, vsf_val .* sin(th_integ));
    fprintf('Haltrin 归一化完成, 计算出的 b = %.4f (经验值 b_e = %.4f)\n', param.b_haltrin_norm, param.b_e);

elseif strcmp(param.phase_func, 'TTHG')
    % TTHG 参数 (0-180度 Clear Ocean 拟合)
    param.coef_b = 0.037; % 示例值，可修改
    param.alpha_TTHG = 0.4435;
    param.g1_TTHG = 0.9900;
    param.g2_TTHG = 0.8238;

elseif strcmp(param.phase_func, 'Ding')
    % Ding 参数 (0-180度 Clear Ocean 拟合)
    param.coef_b = 0.037;
    param.w_mie_ding = 0.8758;
    param.gamma_ding = 0.0000;
    param.g_ding = 0.9803;
    param.f_ding = 4.6721; % 修正因子
    
elseif strcmp(param.phase_func, 'HG')
    % Single HG 参数
    param.coef_b = 0.037;
    param.g_HG = 0.9807; 

elseif strcmp(param.phase_func, 'Petzold')
    % Petzold 实测数据插值
    param.coef_b = 0.037;
    if exist('petzold_ocean.mat', 'file')
        load('petzold_ocean.mat', 'petzold_ocean');
        theta_raw = petzold_ocean(:, 1);
        vsf_raw = petzold_ocean(:, 4); % Clear Ocean
        
        [theta_raw, sort_idx] = sort(theta_raw);
        vsf_raw = vsf_raw(sort_idx);
        
        if theta_raw(1) > 0, theta_raw = [0; theta_raw]; vsf_raw = [vsf_raw(1); vsf_raw]; end
        if theta_raw(end) < pi, theta_raw = [theta_raw; pi]; vsf_raw = [vsf_raw(end); vsf_raw]; end
        
        b_calc = 2 * pi * trapz(theta_raw, vsf_raw .* sin(theta_raw));
        param.petzold_theta = theta_raw;
        param.petzold_p = vsf_raw / b_calc; % 1/sr (未乘4pi，注意后续处理统一乘4pi)
        % 修正：为了配合统一接口，这里让 P 的积分为 1/4pi? 
        % 不，通常 P_HG 积分为 1。如果 vsf/b_calc，则积分为 1/4pi。
        % 为了代码统一，我们让 param.petzold_p 存储 4pi*VSF/b，即标准的归一化相函数
        param.petzold_p = 4 * pi * vsf_raw / b_calc; 
    else
        error('未找到 petzold_ocean.mat');
    end
end

% 计算总衰减
param.coef_c = param.coef_a + param.coef_b; 
param.albedo = param.coef_b / param.coef_c; 

% --- 湍流与相位屏 (保持不变) ---
T_avg = 20; S_avg = 35; epsilon = 1e-6; chi_T = 1e-8; eta = 1e-3; H_ratio = -20;
N_screens = 20; D_screen = 2.0; N_grid = 256; delta_z_screen = Link_Dist / N_screens; 

% --- 时间轴设置 ---
dt = 1e-9; t_min = dt; t_max= 1e-6;                       
param.T_bins = t_min : dt : t_max;
N_bins = length(param.T_bins);
h_time = zeros(1, N_bins); 

% --- 仿真控制 ---
N_packets = 1e5; % 光子数

%% ================= 2. 预计算: 相位屏链 (保持不变) =================
fprintf('1. 生成相位屏链...\n');
% ... (此处省略相位屏生成代码，与原版一致，为节省篇幅) ...
% 如果需要完整运行，请保留原有的相位屏生成代码块
% 这里为了演示核心逻辑，直接跳过生成步骤，实际运行时请补全
Screen_Chain = []; % 占位，如需启用湍流请取消注释下方代码
% [代码块略，请保持原文件中的相位屏生成部分]

%% ================= 3. 主仿真循环 =================
fprintf('3. 开始仿真 (Mode: %s)...\n', param.phase_func);
tic;

% 预计算辅助向量
n_vec = Link_Dir;
if abs(n_vec(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
u_vec = cross(up_temp, n_vec); u_vec = u_vec / norm(u_vec);
v_vec = cross(n_vec, u_vec);   v_vec = v_vec / norm(v_vec);

for p = 1:N_packets
    % --- 发射初始化 ---
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos_local = r0*cos(phi0)*u_vec + r0*sin(phi0)*v_vec;
    pos_init = Tx_Pos + pos_local;
    
    theta_div_physics = lambda / (pi * w0);
    U = theta_div_physics * sqrt(-0.5 * log(rand()));
    psi_ini = 2 * pi * rand;
    dir_init = rotate_direction(mu_T, U, psi_ini);
    
    weight_init = 1.0;
    
    % ================= [分支 A] 直射 =================
    % 简化逻辑：如果没有相位屏，直接算距离
    % 如果有相位屏，请调用 ray_march_generic
    path_len_ballistic = Link_Dist; % 简化
    hit_flag = true; % 假设对准
    
    if hit_flag
        attenuation = exp(-param.coef_c * path_len_ballistic);
        energy_val = weight_init * attenuation;
        t_arrival = path_len_ballistic / param.c_water;
        bin_idx = floor((t_arrival - t_min) / dt) + 1;
        if bin_idx >= 1 && bin_idx <= N_bins
            h_time(bin_idx) = h_time(bin_idx) + energy_val;
        end
    end
    
    % ================= [分支 B] 散射 (重点修改部分) =================
    pos = pos_init; 
    dir = dir_init;
    weight = weight_init;
    current_dist_traveled = 0;
    
    d_step = -log(rand()) / param.coef_c;
    
    % 简化移动逻辑 (不含相位屏)
    pos_new = pos + dir * d_step;
    step_len = d_step;
    current_dist_traveled = current_dist_traveled + step_len;
    
    if dot(pos_new - Tx_Pos, Link_Dir) >= Link_Dist
        continue; 
    end
    pos = pos_new;
    
    % --- 高阶散射循环 ---
    for order = 1 : param.n_max
        vec_to_rx = Rx_Pos - pos;
        dist_to_rx = norm(vec_to_rx);
        dir_to_rx = vec_to_rx / dist_to_rx;
        
        % --- 1. PIS 估计 (计算 P_phase) ---
        if acos(dot(-dir_to_rx, Rx_Normal)) <= Rx_FOV
            cos_theta_s = dot(dir, dir_to_rx);
            % 获取相函数值 (单位 1/sr, 归一化为 4pi积分)
            p_phase_val = get_phase_val(cos_theta_s, param);
            
            % PIS 公式中的 p_phase 是概率密度 PDF(Omega)，积分为1
            % 我们的 get_phase_val 返回的是 P (积分为4pi)
            % 所以需要除以 4pi 才能变成 standard PDF (1/sr)
            p_pdf_val = p_phase_val / (4*pi); 
            
            omega = Rx_Area / (dist_to_rx^2) * abs(dot(dir_to_rx, Rx_Normal));
            prob_survival = exp(-param.coef_c * dist_to_rx);
            
            energy_pis = weight * param.albedo * min(1, p_pdf_val * omega) * prob_survival;
            
            t_pis_arrival = (current_dist_traveled + dist_to_rx) / param.c_water;
            bin_idx = floor((t_pis_arrival - t_min) / dt) + 1;
            if bin_idx >= 1 && bin_idx <= N_bins
                h_time(bin_idx) = h_time(bin_idx) + energy_pis;
            end
        end
        
        if order == param.n_max, break; end
        
        weight = weight * param.albedo;
        if weight < 1e-9, break; end
        
        % --- 2. 散射采样 (Importance Sampling) ---
        theta_i = pi * rand; % [0, pi]
        phi = 2 * pi * rand; % [0, 2pi]
        
        % 计算目标分布的 P 值 (1/sr, 4pi归一化)
        p_val_target = get_phase_val(cos(theta_i), param);
        
        % 权重修正公式: 
        % Target PDF(theta) = P(theta)/4pi * 2pi * sin(theta) = 0.5 * P * sin
        % Proposal PDF(theta) = 1/pi
        % Weight = Target / Proposal = 0.5 * pi * P * sin(theta)
        % 等等，之前的推导是:
        % W = P(theta)[1/sr] * 2pi^2 * sin(theta) ? 
        % 让我们复查一下：
        % Target (Angle Density) f(theta) = P(theta_norm_4pi) / 4pi * 2pi * sin(theta) = 0.5 * P * sin(theta)
        % Proposal (Angle Density) q(theta) = 1/pi
        % Ratio = (0.5 * P * sin) / (1/pi) = 0.5 * pi * P * sin(theta)
        
        % [注意] 之前的代码用了 2*pi^2，那是假设 P 是归一化为 1 的？
        % 不，之前的代码 pdf_HG 返回的是 (1-g^2)/(4pi... )，这积分为 1。
        % 如果 get_phase_val 返回的是积分为 4pi 的 P (如 Haltrin)，那公式得变。
        
        % 统一口径：get_phase_val 返回的 P 都是 积分为 4pi 的形式 (类似 scatter_test3 的图2)
        % 那么对应的 PDF(1/sr) = P / 4pi
        % 对应的 PDF(theta) = (P/4pi) * 2pi * sin(theta) = 0.5 * P * sin(theta)
        % 权重因子 = (0.5 * P * sin) / (1/pi) = 0.5 * pi * P * sin(theta)
        
        weight_factor = 0.5 * pi * p_val_target * sin(theta_i);
        
        weight = weight * weight_factor;
        dir = rotate_direction(dir, theta_i, phi);
        
        d_step = -log(rand()) / param.coef_c;
        pos = pos + dir * d_step;
        current_dist_traveled = current_dist_traveled + d_step;
    end
end
toc;

% 结果绘图 (略，保持原样)
P_rx_total = sum(h_time);
fprintf('Total Power: %.4e\n', P_rx_total);
figure; plot(param.T_bins, h_time); title(['IRF - ' param.phase_func]);


%% ================= 4. 核心辅助函数 =================

% --- [核心] 统一相函数接口 ---
% 输出: P(theta) [sr^-1], 归一化为全空间积分 = 4pi
function p_val = get_phase_val(cos_theta, param)
    if strcmp(param.phase_func, 'Haltrin')
        % 将 cos(theta) 转为角度制 (degrees)
        theta_rad = acos(cos_theta);
        theta_deg = theta_rad * 180 / pi;
        
        % 代入经验公式
        t = theta_deg;
        term = 1 + (-1)^1*param.k1*t.^0.5 + (-1)^2*param.k2*t.^1.0 + ...
               (-1)^3*param.k3*t.^1.5 + (-1)^4*param.k4*t.^2.0 + ...
               (-1)^5*param.k5*t.^2.5;
        vsf = exp(param.q_e * term);
        
        % 转换为 P (P = 4pi * VSF / b)
        p_val = 4 * pi * vsf / param.b_haltrin_norm;
        
    elseif strcmp(param.phase_func, 'TTHG')
        p1 = pdf_HG_norm_4pi(cos_theta, param.g1_TTHG);
        p2 = pdf_HG_norm_4pi(cos_theta, param.g2_TTHG);
        p_val = param.alpha_TTHG * p1 + (1 - param.alpha_TTHG) * p2;
        
    elseif strcmp(param.phase_func, 'Ding')
        p_ray = pdf_Ray_norm_4pi(cos_theta, param.gamma_ding);
        p_ghg = pdf_GHG_norm_4pi(cos_theta, param.g_ding, param.f_ding);
        p_val = (1 - param.w_mie_ding) * p_ray + param.w_mie_ding * p_ghg;
        
    elseif strcmp(param.phase_func, 'HG')
        p_val = pdf_HG_norm_4pi(cos_theta, param.g_HG);
        
    elseif strcmp(param.phase_func, 'Petzold')
        theta_rad = acos(cos_theta);
        % 线性插值 (petzold_p 已经是 4pi 归一化的了)
        p_val = interp1(param.petzold_theta, param.petzold_p, theta_rad, 'linear', 0);
    else
        p_val = 1; % 默认各向同性
    end
end

% --- 基础相函数 (归一化为 4pi) ---
function p = pdf_HG_norm_4pi(cos_theta, g)
    % 标准 HG 公式 (积分为 4pi)
    p = (1 - g^2) ./ (1 + g^2 - 2*g*cos_theta).^1.5;
end

function p = pdf_Ray_norm_4pi(cos_theta, gamma)
    % 广义瑞利 (积分为 4pi)
    % 原 pdf_Rayleigh 返回的是 1/sr (积分为1)，分子系数是 3/16pi
    % 乘以 4pi 后 -> 分子系数 3/4
    p = (3/4) * (1 + 3*gamma + (1-gamma)*cos_theta^2) / (1 + 2*gamma);
end

function p = pdf_GHG_norm_4pi(cos_theta, g, f)
    % 广义 HG (积分为 4pi)
    % 原 pdf_Mie 返回的是 1/sr (积分为1)，除以了 4pi
    % 这里直接用无 4pi 的版本
    term1 = (1 - g^2) ./ (1 + g^2 - 2*g*cos_theta).^1.5;
    term2 = f * 0.5 * (1 - g^2) * (3*cos_theta^2 - 1) ./ (1 + g^2)^1.5;
    p = term1 + term2;
end

% --- 几何旋转函数 ---
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