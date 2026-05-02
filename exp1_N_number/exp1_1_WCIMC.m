%% Exp2: WCI-MC (Wavefront-Coupled + Adaptive HG-PIS, With Turbulence)
clc; clear; close all;

% ================= 物理参数初始化 =================
N_arr = [1e3, 1e4, 1e5, 1e6, 1e7];
num_N = length(N_arr);

lambda = 514e-9; k_wave = 2 * pi / lambda;
w0 = 0.002; 
div_angle = 0.1 * pi / 180; 
theta_half_div_physics = div_angle / 2;
Rx_Aperture = 0.05; 
Rx_FOV = 40 * pi / 180; 
Rx_Area = pi * (Rx_Aperture / 2)^2;
L = 20; % 传输距离 20m

% Coastal Ocean (Petzold 514nm)
param.coef_a = 0.179; param.coef_b = 0.219; param.coef_c = 0.398;
param.albedo = param.coef_b / param.coef_c;
param = calc_haltrin_params(param);

% 弱湍流参数
Cn2 = 1e-15; 

Tx_Pos = [0, 0, 0]; Rx_Pos = [0, L, 0];
Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
if abs(Link_Dir(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
u_vec_Tx = cross(up_temp, Link_Dir); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(Link_Dir, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

W_L_geo = w0 + L * tan(theta_half_div_physics);
rho0_Link = (0.545 * k_wave^2 * Cn2 * L)^(-3/5);

P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max;
g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);

PL_arr = zeros(1, num_N); Time_arr = zeros(1, num_N);

fprintf('--- 运行 WCI-MC 收敛性测试 ---\n');
for idx = 1:num_N
    N_packets = N_arr(idx); P_rx_accum = 0; tic;
    
    for p = 1:N_packets
        r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
        pos_init = Tx_Pos + r0*cos(phi0)*u_vec_Tx + r0*sin(phi0)*v_vec_Tx;
        U_init = theta_half_div_physics * sqrt(-0.5 * log(rand()));
        dir_init = rotate_direction(Link_Dir, U_init, 2*pi*rand());
        
        pos = pos_init; dir = dir_init; weight = 1.0; P_packet = 0;
        
        % 0阶直射
        cos_th = dot(dir_init, Link_Dir);
        if cos_th > 0
            d_plane = L / cos_th; pos_end = Tx_Pos + dir_init * d_plane;
            cos_rx_tilt = abs(dot(dir_init, Rx_Normal));
            if acos(cos_rx_tilt) <= Rx_FOV / 2
                r_wander = norm(pos_end - Rx_Pos);
                r_eff = (Rx_Aperture / 2) * sqrt(cos_rx_tilt);
                
                % 湍流惩罚
                rho0_ballistic = rho0_Link * (L / d_plane)^(3/5);
                W_turb_LT = 2 * d_plane / (k_wave * rho0_ballistic);
                cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*w0))^(1/3));
                W_spot = sqrt(W_L_geo^2 + (W_turb_LT * cf)^2);
                
                a_p = 2 * r_wander / W_spot; b_p = 2 * r_eff / W_spot;
                recv_frac = 1 - marcumq(a_p, b_p);
                P_packet = P_packet + exp(-param.coef_c * d_plane) * recv_frac;
            end
        end
        
        % 多重散射 (1~10阶)
        for ord = 1:10
            d_step = -log(rand()) / param.coef_b; % Method B
            pos = pos + dir * d_step; weight = weight * exp(-param.coef_a * d_step);
            if dot(pos - Tx_Pos, Link_Dir) >= L || weight < 1e-9, break; end
            
            vec2rx = Rx_Pos - pos; d2rx = norm(vec2rx); dir2rx = vec2rx / d2rx;
            if acos(dot(-dir2rx, Rx_Normal)) <= Rx_FOV/2
                p_ph = pdf_Empirical(dot(dir, dir2rx), param);
                omega = Rx_Area / (d2rx^2) * abs(dot(dir2rx, Rx_Normal));
                P_packet = P_packet + weight * min(1, p_ph * omega) * exp(-param.coef_c * d2rx);
            end
            
            xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi);
            cos_th_i = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1); th_i = acos(cos_th_i);
            weight = weight * min(1e3, pdf_Empirical(cos_th_i, param) / ((1 - g_prop^2) / (4 * pi * (1 + g_prop^2 - 2 * g_prop * cos_th_i)^1.5)));
            dir = rotate_direction(dir, th_i, 2*pi*rand());
        end
        P_rx_accum = P_rx_accum + P_packet;
    end
    
    Time_arr(idx) = toc; PL_arr(idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
    fprintf('N = 10^%d | Time: %.2f s | Path Loss: %.2f dB\n', log10(N_packets), Time_arr(idx), -PL_arr(idx));
end
save('data_exp2_WCIMC.mat', 'N_arr', 'PL_arr', 'Time_arr');

% --- Helper Functions ---
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