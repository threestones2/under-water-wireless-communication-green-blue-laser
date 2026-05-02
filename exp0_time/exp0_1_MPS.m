%% 1. WCI-MC (Wavefront-Coupled + Adaptive HG-PIS Accelerated)
clc; clear; close all;

water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr = [0.1514, 0.398, 2.190];
N_packets = 1e5; n_max = 200; % 保持仿真光子数不变，对齐散射阶数

% --- 实际硬件系统几何与光学参数 (对齐 exp3) ---
lambda = 514e-9;             
w0 = 0.002;                  % 初始束腰半径 5 mm
div_angle = 0.1 * pi / 180;  % 发射源全束散角 0.1 度
theta_half_div_physics = div_angle / 2;
Rx_Aperture = 0.01;          % 接收端孔径直径 5 cm
Rx_FOV = 5 * pi / 180;       % 接收机全视场角 5 度

% --- 空间链路布局 ---
Tx_Pos = [0, 0, 0]; 
Rx_Pos = [0, 10, 0];         % 通信距离 20 m
Link_Dist = norm(Rx_Pos - Tx_Pos); 
Link_Dir = (Rx_Pos - Tx_Pos) / Link_Dist; 
Rx_Normal = -Link_Dir;
Rx_Area = pi*(Rx_Aperture/2)^2;

% 生成发射端局部发射平面的正交基
if abs(Link_Dir(3)) < 0.9, up_temp_Tx = [0, 0, 1]; else, up_temp_Tx = [1, 0, 0]; end
u_vec_Tx = cross(up_temp_Tx, Link_Dir); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(Link_Dir, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

fprintf('--- 启动 WCI-MC 仿真 (自适应最优 HG-PIS 加速版) ---\n');

for wt_idx = 1:3
    param.coef_a = coef_a_arr(wt_idx); param.coef_b = coef_b_arr(wt_idx); param.coef_c = coef_c_arr(wt_idx);
    param.albedo = param.coef_b / param.coef_c;
    param = calc_haltrin_params(param);
    
    P_max = pdf_Empirical(1.0, param); 
    C_val = 4 * pi * P_max;
    g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);
    
    P_rx_total_accum = 0;
    tic;
    rng(123456, 'twister');
    for p = 1:N_packets
        r0 = w0 * sqrt(-0.5*log(rand())); 
        phi0 = 2*pi*rand();
        pos_local = r0*cos(phi0)*u_vec_Tx + r0*sin(phi0)*v_vec_Tx;
        pos_init = Tx_Pos + pos_local;
        
        U_init = theta_half_div_physics * sqrt(-0.5 * log(rand()));
        psi_init = 2 * pi * rand();
        dir_init = rotate_direction(Link_Dir, U_init, psi_init);
        
        weight_init = 1.0;
        P_packet_rx = 0;
        
        % 2. 0阶直射路径 (对齐 exp3 的纯几何截断)
        pos_centroid = Tx_Pos;
        dir_centroid = Link_Dir;
        
        cos_theta = dot(dir_centroid, Link_Dir);
        if cos_theta > 0
            dist_to_plane = Link_Dist / cos_theta;
            pos_end = pos_centroid + dir_centroid * dist_to_plane;
            r_wander = norm(pos_end - Rx_Pos);
            cos_rx_tilt = abs(dot(dir_centroid, Rx_Normal));
            
            if acos(cos_rx_tilt) <= Rx_FOV / 2
                if r_wander <= Rx_Aperture / 2
                    w_ballistic = exp(-param.coef_c * dist_to_plane);
                    P_packet_rx = P_packet_rx + w_ballistic;
                end
            end
        end
        
        % 3. 1~n阶散射路径演化
        pos = pos_init; 
        dir = dir_init; 
        weight = weight_init;
        
        for order = 1 : n_max
            d_step = -log(rand()) / param.coef_b;
            pos = pos + dir * d_step;
            weight = weight * exp(-param.coef_a * d_step);
            
            if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist, break; end
            % 引入无偏轮盘赌
            if weight < 1e-9
                if rand() > 0.1, break; else, weight = weight * 10; end
            end
            
            vec_to_rx = Rx_Pos - pos; dist_to_rx = norm(vec_to_rx); dir_to_rx = vec_to_rx / dist_to_rx;
            if acos(dot(-dir_to_rx, Rx_Normal)) <= Rx_FOV/2
                cos_theta_s = dot(dir, dir_to_rx);
                p_phase = pdf_Empirical(cos_theta_s, param);
                omega = Rx_Area / (dist_to_rx^2) * abs(dot(dir_to_rx, Rx_Normal));
                prob_survival = exp(-param.coef_c * dist_to_rx);
                P_packet_rx = P_packet_rx + weight * min(1, p_phase * omega) * prob_survival;
            end
            
            xi = rand();
            term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi);
            cos_th_i = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1);
            th_i = acos(cos_th_i);
            
            q_HG = (1 - g_prop^2) / (4 * pi * (1 + g_prop^2 - 2 * g_prop * cos_th_i)^1.5);
            p_Haltrin = pdf_Empirical(cos_th_i, param);
            weight_factor = min(1e3, p_Haltrin / q_HG);
            weight = weight * weight_factor;
            
            dir = rotate_direction(dir, th_i, 2*pi*rand());
        end
        
        P_rx_total_accum = P_rx_total_accum + P_packet_rx;
    end
    t_run = toc;
    
    P_rx_mean = P_rx_total_accum / N_packets;
    fprintf('%-15s | Time: %6.2f s | Path Loss: %6.2f dB | Adaptive g: %.4f\n', water_types{wt_idx}, t_run, -10*log10(max(P_rx_mean, 1e-300)), g_prop);
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