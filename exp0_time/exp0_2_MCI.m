%% 2. MCI-PIS (Pure Directed Reception + Uniform PIS)
clc; clear; close all;

water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; coef_b_arr = [0.0374, 0.219, 1.824]; coef_c_arr = [0.1514, 0.398, 2.190];
N_packets = 1e5; n_max = 200; 

% --- 实际硬件系统几何与光学参数 (对齐 exp3) ---
r = 10;                     % 通信距离 20 m
lambda = 514e-9;             
w0 = 0.002;                  
div_angle = 0.1 * pi / 180;  
theta_half_div_physics = div_angle / 2;
rx_radius = 0.01;           
beta_R = 5 * pi / 180;      % 接收机全视场角 5 度
A_r = pi*(rx_radius)^2; 

mu_T = [0, 0, 1]; mu_R = [0, 0, -1]; rx_pos = [0, 0, r];
if abs(mu_T(3)) < 0.9, up_temp_Tx = [0, 0, 1]; else, up_temp_Tx = [1, 0, 0]; end
u_vec_Tx = cross(up_temp_Tx, mu_T); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(mu_T, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

fprintf('--- 启动 MCI-PIS 仿真 (均匀 PIS 基准) ---\n');

for wt_idx = 1:3
    ka = coef_a_arr(wt_idx); ks = coef_b_arr(wt_idx); ke = coef_c_arr(wt_idx);
    param.coef_b = ks; param.albedo = ks/ke;
    param = calc_haltrin_params(param);
    
    P_rx_total_accum = 0;
    tic;
    rng(123456, 'twister'); 
    for k = 1:N_packets
        r0 = w0 * sqrt(-0.5*log(rand())); 
        phi0 = 2*pi*rand();
        pos_local = r0*cos(phi0)*u_vec_Tx + r0*sin(phi0)*v_vec_Tx;
        pos_init = [0, 0, 0] + pos_local;
        
        U_init = theta_half_div_physics * sqrt(-0.5 * log(rand()));
        psi_init = 2 * pi * rand();
        dir_init = rotate_direction(mu_T, U_init, psi_init);
        
        P_packet_rx = 0;
        
        % 2. 0阶直射路径 (对齐 exp3 的纯几何截断)
        pos_centroid = [0, 0, 0];
        dir_centroid = mu_T;
        
        cos_theta = dot(dir_centroid, mu_T);
        if cos_theta > 0
            dist_to_plane = r / cos_theta;
            pos_end = pos_centroid + dir_centroid * dist_to_plane;
            r_wander = norm(pos_end - rx_pos);
            cos_rx_tilt = abs(dot(dir_centroid, mu_R));
            
            if acos(cos_rx_tilt) <= beta_R / 2
                if r_wander <= rx_radius
                    w_ballistic = exp(-ke * dist_to_plane);
                    P_packet_rx = P_packet_rx + w_ballistic;
                end
            end
        end
        
        % 3. 1~n阶散射演化 (均匀采样)
        pos = pos_init; 
        mu = dir_init; 
        O_star = 1.0;
        
        for i = 1:n_max
            d_i = -log(rand) / ks;
            
            theta_i = pi * rand();
            phi_i = 2 * pi * rand;
            p_val = pdf_Empirical(cos(theta_i), param);
            
            weight_factor = min(1.0, p_val * 2 * pi^2 * sin(theta_i)); 
            O_star_i = exp(-ka * d_i) * weight_factor;
            O_star = O_star * O_star_i;
            
            mu = rotate_direction(mu, theta_i, phi_i);
            pos = pos + d_i * mu;
            if dot(pos, mu_T) >= r, break; end
            if O_star < 1e-9
                if rand() > 0.1, break; else, O_star = O_star * 10; end
            end
            
            vec_to_rx = rx_pos - pos; d_to_rx = norm(vec_to_rx); cos_phi_r = dot(-mu_R, vec_to_rx) / d_to_rx;
            if cos_phi_r >= cos(beta_R/2)
                theta_n = acos(dot(mu, vec_to_rx) / d_to_rx);
                f_Theta_n = pdf_Empirical(cos(theta_n), param);
                Omega_r = A_r / d_to_rx^2 * cos_phi_r;
                p_d = exp(-ke * d_to_rx) * cos_phi_r * min(1, f_Theta_n * Omega_r);
                P_packet_rx = P_packet_rx + p_d * O_star;
            end
        end
        
        P_rx_total_accum = P_rx_total_accum + P_packet_rx;
    end
    t_run = toc;
    
    P_rx_mean = P_rx_total_accum / N_packets;
    fprintf('%-15s | Time: %6.2f s | Path Loss: %6.2f dB\n', water_types{wt_idx}, t_run, -10*log10(max(P_rx_mean, 1e-300)));
end

% 辅助函数
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo;
    p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e);
    p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000, val(i) = exp(p.q_e*(1-p.k1*max(th(i)*180/pi,1e-6)^0.5 + p.k2*max(th(i)*180/pi,1e-6)^1.0 - p.k3*max(th(i)*180/pi,1e-6)^1.5 + p.k4*max(th(i)*180/pi,1e-6)^2.0 - p.k5*max(th(i)*180/pi,1e-6)^2.5)); end
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