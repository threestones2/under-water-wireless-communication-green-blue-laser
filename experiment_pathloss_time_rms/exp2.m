%% UWOC Channel Characteristics vs. Link Distance (Integrated Comparison - Corrected)
clc; clear; close all;

% ================= 1. 全局参数初始化 =================
N_packets = 5e4;           % 建议出正式图时设为 1e5，调试可设为 1e4
n_max = 10;                % 最大散射阶数
lambda = 514e-9;           % 光波长 (m)
k_wave = 2 * pi / lambda;
w0 = 0.1;                  % 束腰半径 (m)
div_angle = 3 * pi / 180;  % 发散角 (rad)
Rx_Aperture = 0.2;         % 接收机孔径 (m)
Rx_FOV = 10 * pi / 180;    % 接收机视场角 (rad)
Rx_Area = pi * (Rx_Aperture / 2)^2;
c_water = 2.25e8;          % 水中光速 (m/s)

% Clear Ocean 参数 (Petzold)
param.coef_a = 0.114; 
param.coef_b = 0.0374; 
param.coef_c = 0.1514;
% %Petzold近岸海水
% param.coef_c = 0.398;
% param.coef_a = 0.179;
% param.coef_b=0.219;

% Petzold港口海水
% param.coef_c = 2.190;
% param.coef_a = 0.366;
% param.coef_b=1.824;
param.albedo = param.coef_b / param.coef_c;
param = calc_haltrin_params(param);

% 自适应计算水体固有发散角 (用于湍流波束惩罚)
theta_peak = 1e-6; 
p_peak = pdf_Empirical(cos(theta_peak), param);
opts = optimset('Display', 'off');
try
    param.theta_water = fzero(@(th) pdf_Empirical(cos(th), param) - p_peak*exp(-2), [1e-6, 0.1], opts);
catch
    param.theta_water = 0.02; % 异常回退基准值
end

% 距离设置
distances = 10:10:50; 
num_L = length(distances);

% 结果存储矩阵
PL_WCIMC = zeros(1, num_L); Tau_WCIMC = zeros(1, num_L);
PL_HGPIS = zeros(1, num_L); Tau_HGPIS = zeros(1, num_L);
PL_MCIPIS = zeros(1, num_L); Tau_MCIPIS = zeros(1, num_L);
PL_MCS   = zeros(1, num_L); Tau_MCS   = zeros(1, num_L);

% 湍流等效参数预估 (弱湍流设定)
Cn2_eq = 1e-11; % 简化 OTOPS 的等效 Cn2 以避免逐距离重新积分

fprintf('=== 开始传输距离对比仿真实验 ===\n');

for idx = 1:num_L
    L = distances(idx);
    fprintf('\n正在仿真传输距离: %d m (%d/%d)\n', L, idx, num_L);
    
    Tx_Pos = [0, 0, 0];
    Rx_Pos = [0, L, 0];
    Link_Dir = [0, 1, 0];
    Rx_Normal = [0, -1, 0];
    
    z_R = (pi * w0^2) / lambda; 
    W_L_diff = w0 * sqrt(1 + (L / z_R)^2);
    geom_factor = (2 * Rx_Area) / (pi * W_L_diff^2);
    P_ballistic_base = exp(-param.coef_c * L) * min(1, geom_factor);
    
    % 时间轴动态划分
    t_min = L / c_water; 
    dt = 1e-10; 
    t_max = t_min + 1e-6; 
    T_bins = t_min : dt : t_max;
    N_bins = length(T_bins);
    
    % LUT 初始化 (MCS使用)
    th_axis = linspace(0, pi, 5000); pdf_vals = zeros(size(th_axis));
    for i = 1:length(th_axis), pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i)); end
    cdf_vals = cumtrapz(th_axis, pdf_vals); cdf_vals = cdf_vals / cdf_vals(end);
    
    % HG 自适应 g 初始化
    P_max_val = pdf_Empirical(1.0, param);
    C_val = 4 * pi * P_max_val;
    g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);
    
    % 相干长度初始化 (WCI-MC使用)
    rho0_Link = (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5);
    
    % 初始化各算法时间统计直方图
    h_wcimc = zeros(1, N_bins); h_hgpis = zeros(1, N_bins);
    h_mcipis = zeros(1, N_bins); h_mcs = zeros(1, N_bins);
    
    % ================= [算法基准修正] 直射能量补偿 =================
    bin_b = floor((L/c_water - t_min)/dt) + 1;
    
    % WCI-MC 的直射部分附加湍流衰减
    rho0_ballistic = (0.545 * k_wave^2 * Cn2_eq * L)^(-3/5);
    W_turb_LT = 2 * L / (k_wave * rho0_ballistic);
    corr_factor = max(0, 1 - 0.37 * (rho0_ballistic / (2*w0))^(1/3));
    W_spot_sq = W_L_diff^2 + (W_turb_LT * corr_factor)^2;
    WCI_geom = (2 * Rx_Area) / (pi * W_spot_sq);
    P_ballistic_WCI = exp(-param.coef_c * L) * min(1, WCI_geom);
    
    if bin_b >= 1 && bin_b <= N_bins
        h_wcimc(bin_b)  = h_wcimc(bin_b)  + P_ballistic_WCI * N_packets;
        h_hgpis(bin_b)  = h_hgpis(bin_b)  + P_ballistic_base * N_packets;
        h_mcipis(bin_b) = h_mcipis(bin_b) + P_ballistic_base * N_packets;
        h_mcs(bin_b)    = h_mcs(bin_b)    + P_ballistic_base * N_packets;
    end
    
    for p = 1:N_packets
        if mod(p, 10000) == 0, fprintf('  光子进度: %d / %d\n', p, N_packets); end
        
        % [核心修正] 采用高斯光束发散角统一采样初始出射方向
        U_init = div_angle * sqrt(-0.5 * log(rand()));
        psi_init = 2 * pi * rand();
        dir_global_init = rotate_direction(Link_Dir, U_init, psi_init);
        
        % ================= [算法1] WCI-MC (含湍流与波束拓展) =================
        pos = Tx_Pos; dir = dir_global_init; weight = 1.0; dist_travel = 0;
        
        for ord = 1:n_max
            d_s = -log(rand()) / param.coef_c;
            % 简化直射追踪（因已在外部统计质心宏观衰减，游走阶段聚焦散射特性）
            pos = pos + dir * d_s; dist_travel = dist_travel + d_s;
            if dot(pos - Tx_Pos, Link_Dir) >= L, break; end
            
            % PIS 与 波前退化惩罚
            vec2rx = Rx_Pos - pos; d2rx = norm(vec2rx); d2rx_dir = vec2rx / d2rx;
            if acos(dot(-d2rx_dir, Rx_Normal)) <= Rx_FOV/2
                cos_th_s = dot(dir, d2rx_dir); p_ph = pdf_Empirical(cos_th_s, param);
                omega = Rx_Area / (d2rx^2) * abs(dot(d2rx_dir, Rx_Normal));
                base_w = weight * param.albedo * min(1, p_ph * omega) * exp(-param.coef_c * d2rx);
                
                if base_w > 1e-15
                    % 湍流波前畸变惩罚
                    rho_0 = rho0_Link * (L / d2rx)^(3/5);
                    th_turb_LT = 2 / (k_wave * rho_0);
                    cf = max(0, 1 - 0.37 * (rho_0 / Rx_Aperture)^(1/3));
                    th_turb_ST = th_turb_LT * cf;
                    penalty = (param.theta_water^2) / (param.theta_water^2 + th_turb_ST^2);
                    
                    % 偏轴指向惩罚近似计算 (由于非相干散射光主要受发散角控制)
                    r_wander = d2rx * acos(dot(dir, d2rx_dir)); % 近似偏移量
                    W_tot_sq = (d2rx^2) * (param.theta_water^2 + th_turb_ST^2);
                    point_loss = exp(-2 * r_wander^2 / W_tot_sq);
                    
                    final_w = base_w * penalty * point_loss;
                    t_arr = (dist_travel + d2rx) / c_water;
                    b_idx = floor((t_arr - t_min)/dt) + 1;
                    if b_idx >= 1 && b_idx <= N_bins, h_wcimc(b_idx) = h_wcimc(b_idx) + final_w; end
                end
            end
            
            weight = weight * param.albedo; if weight < 1e-6, break; end
            % 真实相函数采样 (作为后续游走的基础)
            th_s = interp1(cdf_vals, th_axis, rand(), 'linear', 'extrap');
            dir = rotate_direction(dir, th_s, 2*pi*rand());
        end
        
        % ================= [算法2] Adaptive HG-PIS (无湍流) =================
        pos = Tx_Pos; dir = dir_global_init; weight = 1.0; dist_travel = 0;
        
        for ord = 1:n_max
            d_s = -log(rand()) / param.coef_c;
            pos = pos + dir * d_s; dist_travel = dist_travel + d_s;
            if dot(pos - Tx_Pos, Link_Dir) >= L, break; end
            
            vec2rx = Rx_Pos - pos; d2rx = norm(vec2rx); d2rx_dir = vec2rx / d2rx;
            if acos(dot(-d2rx_dir, Rx_Normal)) <= Rx_FOV/2
                cos_th_s = dot(dir, d2rx_dir); p_ph = pdf_Empirical(cos_th_s, param);
                omega = Rx_Area / (d2rx^2) * abs(dot(d2rx_dir, Rx_Normal));
                base_w = weight * param.albedo * min(1, p_ph * omega) * exp(-param.coef_c * d2rx);
                t_arr = (dist_travel + d2rx) / c_water;
                b_idx = floor((t_arr - t_min)/dt) + 1;
                if b_idx >= 1 && b_idx <= N_bins, h_hgpis(b_idx) = h_hgpis(b_idx) + base_w; end
            end
            
            weight = weight * param.albedo; if weight < 1e-6, break; end
            xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi);
            cos_th_i = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1); th_i = acos(cos_th_i);
            
            q_HG = (1 - g_prop^2) / (4 * pi * (1 + g_prop^2 - 2 * g_prop * cos_th_i)^1.5);
            p_H = pdf_Empirical(cos_th_i, param);
            weight = weight * min(1e3, p_H / q_HG);
            dir = rotate_direction(dir, th_i, 2*pi*rand());
        end
        
        % ================= [算法3] MCI-PIS (均匀采样) =================
        pos = Tx_Pos; dir = dir_global_init; weight = 1.0; dist_travel = 0;
        
        for ord = 1:n_max
            d_s = -log(rand()) / param.coef_c;
            pos = pos + dir * d_s; dist_travel = dist_travel + d_s;
            if dot(pos - Tx_Pos, Link_Dir) >= L, break; end
            
            th_i = pi * rand(); phi_i = 2 * pi * rand();
            p_val = pdf_Empirical(cos(th_i), param);
            w_fac = min(1.0, p_val * 2 * pi^2 * sin(th_i));
            weight = weight * param.albedo * w_fac;
            dir = rotate_direction(dir, th_i, phi_i);
            
            vec2rx = Rx_Pos - pos; d2rx = norm(vec2rx); d2rx_dir = vec2rx / d2rx;
            if acos(dot(-d2rx_dir, Rx_Normal)) <= Rx_FOV/2
                cos_th_s = dot(dir, d2rx_dir); p_ph = pdf_Empirical(cos_th_s, param);
                omega = Rx_Area / (d2rx^2) * abs(dot(d2rx_dir, Rx_Normal));
                base_w = weight * min(1, p_ph * omega) * exp(-param.coef_c * d2rx);
                t_arr = (dist_travel + d2rx) / c_water;
                b_idx = floor((t_arr - t_min)/dt) + 1;
                if b_idx >= 1 && b_idx <= N_bins, h_mcipis(b_idx) = h_mcipis(b_idx) + base_w; end
            end
        end
        
        % ================= [算法4] MCS (传统LUT) =================
        pos = Tx_Pos; dir = dir_global_init; weight = 1.0; dist_travel = 0;
        
        for ord = 1:n_max
            d_s = -log(rand()) / param.coef_c;
            pos = pos + dir * d_s; dist_travel = dist_travel + d_s;
            weight = weight * param.albedo; if weight < 1e-6, break; end
            if dot(pos - Tx_Pos, Link_Dir) >= L, break; end
            
            vec2rx = Rx_Pos - pos; d2rx = norm(vec2rx); d2rx_dir = vec2rx / d2rx;
            if acos(dot(-d2rx_dir, Rx_Normal)) <= Rx_FOV/2
                cos_th_s = dot(dir, d2rx_dir); p_ph = pdf_Empirical(cos_th_s, param);
                omega = Rx_Area / (d2rx^2) * abs(dot(d2rx_dir, Rx_Normal));
                base_w = weight * min(1, p_ph * omega) * exp(-param.coef_c * d2rx);
                t_arr = (dist_travel + d2rx) / c_water;
                b_idx = floor((t_arr - t_min)/dt) + 1;
                if b_idx >= 1 && b_idx <= N_bins, h_mcs(b_idx) = h_mcs(b_idx) + base_w; end
            end
            
            th_s = interp1(cdf_vals, th_axis, rand(), 'linear', 'extrap');
            dir = rotate_direction(dir, th_s, 2*pi*rand());
        end
    end
    
    % 数据整理与时延计算
    H_arrays = {h_wcimc/N_packets, h_hgpis/N_packets, h_mcipis/N_packets, h_mcs/N_packets};
    PL_arr = zeros(1,4); Tau_arr = zeros(1,4);
    
    for i = 1:4
        P_tot = sum(H_arrays{i});
        PL_arr(i) = 10 * log10(max(P_tot, 1e-300));
        if P_tot > 0
            t_mean = sum(T_bins .* H_arrays{i}) / P_tot;
            Tau_arr(i) = sqrt(sum(((T_bins - t_mean).^2) .* H_arrays{i}) / P_tot) * 1e9; % 转化为 ns
        end
    end
    
    PL_WCIMC(idx) = PL_arr(1); Tau_WCIMC(idx) = Tau_arr(1);
    PL_HGPIS(idx) = PL_arr(2); Tau_HGPIS(idx) = Tau_arr(2);
    PL_MCIPIS(idx) = PL_arr(3); Tau_MCIPIS(idx) = Tau_arr(3);
    PL_MCS(idx)   = PL_arr(4); Tau_MCS(idx)   = Tau_arr(4);
end

%% ================= 3. 结果绘图 =================
figure('Name', 'Path Loss vs Distance', 'Color', 'w', 'Position', [100, 100, 700, 500]);
plot(distances, PL_WCIMC, '-ro', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
plot(distances, PL_HGPIS, '--bs', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(distances, PL_MCIPIS, '-.g^', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(distances, PL_MCS, ':kd', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Transmission Distance L (m)', 'FontSize', 12);
ylabel('Path Loss (dB)', 'FontSize', 12);
title('Path Loss vs. Transmission Distance', 'FontSize', 14);
legend('WCI-MC (Proposed, w/ Turb)', 'Adaptive HG-PIS (w/o Turb)', 'MCI-PIS (Uniform)', 'MCS (LUT)', 'Location', 'southwest');

figure('Name', 'RMS Delay Spread vs Distance', 'Color', 'w', 'Position', [150, 150, 700, 500]);
plot(distances, Tau_WCIMC, '-ro', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
plot(distances, Tau_HGPIS, '--bs', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(distances, Tau_MCIPIS, '-.g^', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(distances, Tau_MCS, ':kd', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Transmission Distance L (m)', 'FontSize', 12);
ylabel('RMS Delay Spread \tau_{rms} (ns)', 'FontSize', 12);
title('RMS Delay Spread vs. Transmission Distance', 'FontSize', 14);
legend('WCI-MC (Proposed, w/ Turb)', 'Adaptive HG-PIS (w/o Turb)', 'MCI-PIS (Uniform)', 'MCS (LUT)', 'Location', 'northwest');


%% ================= 辅助函数区 =================
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo;
    p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); 
    p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); 
    p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000
        t = max(th(i)*180/pi, 1e-6);
        val(i) = exp(p.q_e*(1 - p.k1*t^0.5 + p.k2*t^1.0 - p.k3*t^1.5 + p.k4*t^2.0 - p.k5*t^2.5)); 
    end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end

function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(max(min(cos_theta, 1), -1)) * 180 / pi, 1e-6);
    term = 1 - param.k1*t_deg^0.5 + param.k2*t_deg^1.0 - param.k3*t_deg^1.5 + param.k4*t_deg^2.0 - param.k5*t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end

function nd = rotate_direction(d, t, p)
    denom = sqrt(1 - d(3)^2);
    if denom < 1e-10, nd = [sin(t)*cos(p), sin(t)*sin(p), sign(d(3))*cos(t)];
    else, nd = [sin(t)/denom*(d(1)*d(3)*cos(p)-d(2)*sin(p))+d(1)*cos(t), sin(t)/denom*(d(2)*d(3)*cos(p)+d(1)*sin(p))+d(2)*cos(t), -sin(t)*cos(p)*denom+d(3)*cos(t)]; end
    nd = nd / norm(nd);
end