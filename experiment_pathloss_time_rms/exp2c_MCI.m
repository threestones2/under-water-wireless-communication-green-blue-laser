%% 3. MCI-PIS vs Distance (Uniform PIS, Without Turbulence)
clc; clear; close all;

% ================= 参数初始化 =================
dist_array = 10:10:50;     % 传输距离数组 (m)
N_packets = 1e5;           
n_max = 10;                
lambda = 514e-9; w0 = 0.1; div_angle = 3 * pi / 180;
Rx_Aperture = 0.2; Rx_FOV = 10 * pi / 180; Rx_Area = pi * (Rx_Aperture / 2)^2;
c_water = 2.25e8; z_R = (pi * w0^2) / lambda;

param.coef_a = 0.114; param.coef_b = 0.0374; param.coef_c = 0.1514;      
param.albedo = param.coef_b / param.coef_c; param = calc_haltrin_params(param);

fprintf('=== 启动 MCI-PIS 距离演化仿真 ===\n');

for d_idx = 1:length(dist_array)
    Link_Dist = dist_array(d_idx);
    Tx_Pos = [0, 0, 0]; Rx_Pos = [0, Link_Dist, 0]; Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
    
    t_min = Link_Dist / c_water; dt = 1e-10; t_max = t_min + 1e-6;
    T_bins = t_min : dt : t_max; N_bins = length(T_bins); h_time = zeros(1, N_bins);
    W_spot = w0 * sqrt(1 + (Link_Dist / z_R)^2);
    
    tic;
    for p = 1:N_packets
        U_init = div_angle * sqrt(-0.5 * log(rand())); psi_init = 2 * pi * rand();
        dir_init = rotate_direction(Link_Dir, U_init, psi_init);
        pos = Tx_Pos; dir = dir_init; weight = 1.0; dist_travel = 0;
        
        % --- 0阶直射路径微观评估 ---
        cos_theta = dot(dir, Link_Dir);
        if cos_theta > 0
            dist_to_plane = Link_Dist / cos_theta; pos_end = Tx_Pos + dir * dist_to_plane; r_wander = norm(pos_end - Rx_Pos);
            cos_rx_tilt = abs(dot(dir, Rx_Normal)); 
            
            geom_loss = (2 * Rx_Area * cos_rx_tilt) / (pi * W_spot^2); 
            if acos(cos_rx_tilt) <= Rx_FOV / 2
                pointing_loss = exp(-2 * r_wander^2 / W_spot^2);
                w_ballistic = exp(-param.coef_c * dist_to_plane) * min(1, geom_loss * pointing_loss);
                b_idx = floor((dist_to_plane/c_water - t_min)/dt) + 1;
                if b_idx >= 1 && b_idx <= N_bins, h_time(b_idx) = h_time(b_idx) + w_ballistic; end
            end
        end
        
        % --- 1~n阶散射路径 ---
        for order = 1 : n_max
            d_step = -log(rand()) / param.coef_c;
            pos = pos + dir * d_step; dist_travel = dist_travel + d_step;
            if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist, break; end
            
            th_i = pi * rand(); phi_i = 2 * pi * rand();
            weight = weight * param.albedo * min(1.0, pdf_Empirical(cos(th_i), param) * 2 * pi^2 * sin(th_i));
            dir = rotate_direction(dir, th_i, phi_i);
            
            vec2rx = Rx_Pos - pos; dist2rx = norm(vec2rx); dir2rx = vec2rx / dist2rx;
            if acos(dot(-dir2rx, Rx_Normal)) <= Rx_FOV/2
                w_rx = weight * min(1, pdf_Empirical(dot(dir, dir2rx), param) * (Rx_Area / dist2rx^2 * abs(dot(dir2rx, Rx_Normal)))) * exp(-param.coef_c * dist2rx);
                b_idx = floor(((dist_travel + dist2rx)/c_water - t_min)/dt) + 1;
                if b_idx >= 1 && b_idx <= N_bins, h_time(b_idx) = h_time(b_idx) + w_rx; end
            end
        end
    end
    t_run = toc;
    
    h_norm = h_time / N_packets; P_rx = sum(h_norm); tau_rms = 0;
    if P_rx > 0, t_mean = sum(T_bins .* h_norm) / P_rx; tau_rms = sqrt( sum(((T_bins - t_mean).^2) .* h_norm) / P_rx ); end
    fprintf('L = %2d m | Time: %5.2fs | Path Loss: %7.2f dB | RMS Delay: %.4f ns\n', Link_Dist, t_run, 10*log10(max(P_rx, 1e-300)), tau_rms*1e9);
end

% 辅助函数
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000, t=max(th(i)*180/pi,1e-6); val(i)=exp(p.q_e*(1-p.k1*t^0.5+p.k2*t^1.0-p.k3*t^1.5+p.k4*t^2.0-p.k5*t^2.5)); end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(max(min(cos_theta, 1), -1)) * 180 / pi, 1e-6); term = 1 - param.k1*t_deg^0.5 + param.k2*t_deg^1.0 - param.k3*t_deg^1.5 + param.k4*t_deg^2.0 - param.k5*t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end
function nd = rotate_direction(d, t, p)
    denom = sqrt(1 - d(3)^2);
    if denom < 1e-10, nd = [sin(t)*cos(p), sin(t)*sin(p), sign(d(3))*cos(t)]; else, nd = [sin(t)/denom*(d(1)*d(3)*cos(p)-d(2)*sin(p))+d(1)*cos(t), sin(t)/denom*(d(2)*d(3)*cos(p)+d(1)*sin(p))+d(2)*cos(t), -sin(t)*cos(p)*denom+d(3)*cos(t)]; end
    nd = nd / norm(nd);
end