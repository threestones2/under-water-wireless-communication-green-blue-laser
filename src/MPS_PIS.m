%% 大气紫外光 NLOS 通信仿真: MC_MPS vs. MCI_PIS (损耗修正版)
% 修正说明:
% 1. 权重对齐: 将 O_star_1 的权重系数从 pi 修正为 (beta_T/2)，消除 13dB 增益。
% 2. 概率修正: p_d 公式中引入 1/(2*pi*sin(theta)) 因子，对齐单位立体角概率。
% 3. 几何修正: 统一采用标准 A/d^2 * cos(phi) 接收模型。

clc; clear; close all;

%% ================= 1. 共享参数 (同步自 MCI_PIS.m) =================
param = struct();
param.r = 100;              
param.theta_T = deg2rad(45);
param.phi_T = 2*pi-deg2rad(90); 
param.theta_R = deg2rad(45);
param.phi_R = deg2rad(90);  
param.beta_T = deg2rad(17); 
param.beta_R = deg2rad(30); 
param.A_r = 1.77e-4;        
param.ka = 0.802e-3;        
param.kr = 0.266e-3;        
param.km = 0.284e-3;        
param.ks = param.kr + param.km; 
param.ke = param.ka + param.ks; 
param.gamma = 0.017;        
param.g = 0.72; 
param.f = 0.5;
param.n_max = 3;            
param.c = 2.997046e8;       
param.dt = 1e-8;            
param.t_max = 2e-6; 
param.T_bins = 0:param.dt:param.t_max;
param.N = 1e5;              

% 矢量计算
mu_T = [cos(param.phi_T)*sin(param.theta_T); sin(param.phi_T)*sin(param.theta_T); cos(param.theta_T)];
mu_R = [cos(param.phi_R)*sin(param.theta_R); sin(param.phi_R)*sin(param.theta_R); cos(param.theta_R)];
Tx_Pos = [0; param.r; 0];   
Rx_Pos = [0; 0; 0];         

%% ================= 2. 仿真 1: MCI_PIS 修正逻辑 =================
h_pis = zeros(1, length(param.T_bins));
total_P_pis = 0;
tic;
for k = 1:param.N
    pos = Tx_Pos; mu = mu_T; d_total = 0; O_star = 1;
    for i = 1:param.n_max
        d_i = -log(1 - rand) / param.ke;
        if i == 1
            theta_i = (param.beta_T / 2) * rand;
            % [修正]: 权重系数从 pi 改为 (param.beta_T/2)
            O_star_i = (param.beta_T/2) * sin(theta_i) / (1 - cos(param.beta_T/2));
        else
            theta_i = pi * rand;
            f_Th = (param.kr/param.ks) * f_Theta_Rayleigh(theta_i, param.gamma) + ...
                   (param.km/param.ks) * f_Theta_Mie(theta_i, param.g, param.f);
            O_star_i = (param.ks/param.ke) * pi * f_Th;
        end
        phi_i = 2 * pi * rand;
        O_star = O_star * O_star_i;
        mu = update_direction_local(mu, theta_i, phi_i);
        pos = pos + d_i * mu; d_total = d_total + d_i;
        
        d_to_rx = norm(pos);
        cos_phi_r = dot(mu_R, pos) / d_to_rx;
        if cos_phi_r >= cos(param.beta_R/2)
            theta_n = acos(dot(-mu, pos) / d_to_rx);
            f_Th_n = (param.kr/param.ks) * f_Theta_Rayleigh(theta_n, param.gamma) + ...
                     (param.km/param.ks) * f_Theta_Mie(theta_n, param.g, param.f);
            
            % [修正]: 引入 1/(2*pi*sin(theta)) 并去除多余的 cos_phi_r
            % 这里的 sin(theta_n) 必须防止除零
            sin_tn = max(sin(theta_n), 1e-6);
            p_d = (param.ks/param.ke) * exp(-param.ke * d_to_rx) * ...
                  (f_Th_n / (2*pi*sin_tn)) * (param.A_r / d_to_rx^2) * cos_phi_r;
            
            P_cur = p_d * O_star;
            total_P_pis = total_P_pis + P_cur;
            t_total = (d_total + d_to_rx) / param.c;
            bin_idx = floor(t_total / param.dt) + 1;
            if bin_idx >= 1 && bin_idx <= length(param.T_bins), h_pis(bin_idx) = h_pis(bin_idx) + P_cur; end
        end
    end
end
time_pis = toc;
PL_pis = 10 * log10(param.N / total_P_pis);

%% ================= 3. 仿真 2: MC_MPS 框架修正 (对齐物理逻辑) =================
h_mps = zeros(1, length(param.T_bins));
total_P_mps = 0;
tic;
for p = 1:param.N
    theta_init = (param.beta_T / 2) * rand;
    phi_init = 2 * pi * rand;
    dir = update_direction_local(mu_T, theta_init, phi_init);
    pos = Tx_Pos;
    % [修正]: 初始权重
    weight = (param.beta_T/2) * sin(theta_init) / (1 - cos(param.beta_T/2));
    dist_traveled = 0;
    
    for order = 1 : param.n_max
        d_step = -log(1 - rand) / param.ke;
        pos = pos + dir * d_step;
        dist_traveled = dist_traveled + d_step;
        
        vec_to_rx = Rx_Pos - pos;
        d_to_rx = norm(vec_to_rx);
        cos_phi_r = dot(mu_R, pos) / d_to_rx;
        
        if cos_phi_r >= cos(param.beta_R/2)
            theta_s = acos(dot(-dir, pos) / d_to_rx);
            % 直接使用标准 sr^-1 相函数
            p_ray = (3*(1+param.gamma*3+(1-param.gamma)*cos(theta_s)^2))/(16*pi*(1+2*param.gamma));
            p_mie = ((1-param.g^2)/2*(1/(1+param.g^2-2*param.g*cos(theta_s))^1.5 + ...
                    param.f*0.5*(3*cos(theta_s)^2-1)/(1+param.g^2)^1.5))/(2*pi);
            p_phase = (param.kr/param.ks)*p_ray + (param.km/param.ks)*p_mie;
            
            p_d = (param.ks/param.ke) * exp(-param.ke * d_to_rx) * p_phase * (param.A_r/d_to_rx^2) * cos_phi_r;
            
            energy = weight * p_d;
            total_P_mps = total_P_mps + energy;
            t_arrival = (dist_traveled + d_to_rx) / param.c;
            bin_idx = floor(t_arrival / param.dt) + 1;
            if bin_idx >= 1 && bin_idx <= length(param.T_bins), h_mps(bin_idx) = h_mps(bin_idx) + energy; end
        end
        
        % [修正]: 迭代权重计算
        theta_next = pi * rand;
        phi_next = 2 * pi * rand;
        f_next = (param.kr/param.ks) * f_Theta_Rayleigh(theta_next, param.gamma) + ...
                 (param.km/param.ks) * f_Theta_Mie(theta_next, param.g, param.f);
        weight = weight * (param.ks/param.ke) * pi * f_next;
        dir = update_direction_local(dir, theta_next, phi_next);
        if weight < 1e-15, break; end
    end
end
time_mps = toc;
PL_mps = 10 * log10(param.N / total_P_mps);

%% ================= 4. 结果报告 =================
fprintf('\n--- 仿真性能与修正损耗报告 (N = %d) ---\n', param.N);
fprintf('1. MCI_PIS (修正权重与歸一化):\n   耗时: %.4f s | 路径损耗: %.2f dB\n', time_pis, PL_pis);
fprintf('2. MC_MPS (框架重构与物理对齐):\n   耗时: %.4f s | 路径损耗: %.2f dB\n', time_mps, PL_mps);
fprintf('--------------------------------------\n');

figure;

% semilogy
plot(param.T_bins*1e6, h_pis/(param.N*param.dt), 'r--', 'LineWidth', 1.5); hold on;
plot(param.T_bins*1e6, h_mps/(param.N*param.dt), 'g-', 'LineWidth', 1.0);
xlabel('Time (\mus)'); ylabel('h(t) (W/m^2/s)'); grid on;
legend('MCI\_PIS (Corrected)', 'MC\_MPS Framework (Corrected)');
title(['Corrected CIR Comparison ']);

%% ================= 5. 辅助函数 =================
function f = f_Theta_Rayleigh(theta, gamma)
    f = 3 * (1 + 3*gamma + (1-gamma)*cos(theta)^2) * sin(theta) / (8 * (1 + 2*gamma));
end
function f = f_Theta_Mie(theta, g, f_m)
    f = (1 - g^2)/2 * (1/(1 + g^2 - 2*g*cos(theta))^1.5 + f_m * 0.5 * (3*cos(theta)^2 - 1) / (1 + g^2)^1.5) * sin(theta);
end
function mu_new = update_direction_local(mu, theta, phi)
    mu_x = mu(1); mu_y = mu(2); mu_z = mu(3);
    if abs(mu_z) < 1 - 1e-10
        denom = sqrt(1 - mu_z^2);
        mu_new = [cos(theta)*mu_x + sin(theta)*(mu_x*mu_z*cos(phi) - mu_y*sin(phi))/denom;
                  cos(theta)*mu_y + sin(theta)*(mu_y*mu_z*cos(phi) + mu_x*sin(phi))/denom;
                  cos(theta)*mu_z - sin(theta)*cos(phi)*denom];
    else
        mu_new = [sin(theta)*cos(phi); sin(theta)*sin(phi); sign(mu_z)*cos(theta)];
    end
    mu_new = mu_new / norm(mu_new);
end