%% Exp1: Adaptive HG-PIS (Without Turbulence, Zero-Allocation Scalarized)
clc; clear; close all;

N_arr = [1e3, 1e4, 1e5, 1e6, 1e7]; num_N = length(N_arr);
lambda = 514e-9; 
w0 = 0.002; 
div_angle = 0.1 * pi / 180; 
theta_half_div_physics = div_angle / 2;
Rx_Aperture = 0.05; 
Rx_FOV = 5 * pi / 180; 
Rx_Area = pi * (Rx_Aperture / 2)^2; 
L = 20;

param.coef_a = 0.179; 
param.coef_b = 0.219; 
param.coef_c = 0.398; 
param.albedo = param.coef_b / param.coef_c;
param = calc_haltrin_params(param);

Tx_Pos = [0, 0, 0]; Rx_Pos = [0, L, 0]; Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
if abs(Link_Dir(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
u_vec_Tx = cross(up_temp, Link_Dir); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx); 
v_vec_Tx = cross(Link_Dir, u_vec_Tx); v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

% 预计算判定常数与解包参考基矢 (消除循环内解包开销)
cos_FOV_half = cos(Rx_FOV / 2);
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3);
Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3);
Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
Ux = u_vec_Tx(1); Uy = u_vec_Tx(2); Uz = u_vec_Tx(3);
Vx = v_vec_Tx(1); Vy = v_vec_Tx(2); Vz = v_vec_Tx(3);

P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max;
g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);

PL_arr = zeros(1, num_N); Time_arr = zeros(1, num_N);

fprintf('--- 运行 Adaptive HG-PIS 收敛性测试 (零分配标量加速) ---\n');
for idx = 1:num_N
    N_packets = N_arr(idx); P_rx_accum = 0; tic;
    rng(123456, 'twister'); 
    
    for p = 1:N_packets
        r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
        cp0 = cos(phi0); sp0 = sin(phi0);
        
        % 标量化初始位置
        p1 = Tx + r0*cp0*Ux + r0*sp0*Vx;
        p2 = Ty + r0*cp0*Uy + r0*sp0*Vy;
        p3 = Tz + r0*cp0*Uz + r0*sp0*Vz;
        
        U_init = theta_half_div_physics * sqrt(-0.5 * log(rand()));
        dir = rotate_direction_fast(Link_Dir, cos(U_init), sin(U_init), 2*pi*rand());
        d1 = dir(1); d2 = dir(2); d3 = dir(3);
        
        weight = 1.0; P_packet = 0;
        
        % 0阶直射路径 (纯几何截断)
        cos_th = d1*Lx + d2*Ly + d3*Lz;
        if cos_th > 0
            d_plane = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / cos_th; 
            pos_end_1 = p1 + d1 * d_plane; 
            pos_end_2 = p2 + d2 * d_plane; 
            pos_end_3 = p3 + d3 * d_plane; 
            
            cos_rx_tilt = abs(d1*Nx + d2*Ny + d3*Nz);
            if cos_rx_tilt >= cos_FOV_half
                r_hit_sq = (pos_end_1 - Rx)^2 + (pos_end_2 - Ry)^2 + (pos_end_3 - Rz)^2;
                if r_hit_sq <= Rx_Aperture_half_sq
                    P_packet = P_packet + exp(-param.coef_c * d_plane);
                end
            end
        end
        
        % 多重散射分支
        for ord = 1:10
            d_step = -log(rand()) / param.coef_b; 
            p1 = p1 + d1 * d_step; p2 = p2 + d2 * d_step; p3 = p3 + d3 * d_step; 
            weight = weight * exp(-param.coef_a * d_step);
            
            if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L, break; end 
            if weight < 1e-9
                if rand() > 0.1, break; else, weight = weight * 10; end
            end
            
            vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; 
            d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2;
            d2rx = sqrt(d2rx_sq); 
            dir2rx_1 = vec2rx_1 / d2rx; dir2rx_2 = vec2rx_2 / d2rx; dir2rx_3 = vec2rx_3 / d2rx;
            
            cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
            if cos_inc >= cos_FOV_half
                omega = Rx_Area / d2rx_sq * cos_inc;
                cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                P_packet = P_packet + weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * d2rx);
            end
            
            xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); 
            cos_th_i = (1 + g_prop^2 - term^2) / (2 * g_prop);
            if cos_th_i > 1, cos_th_i = 1; elseif cos_th_i < -1, cos_th_i = -1; end
            
            denom_core = 1 + g_prop^2 - 2 * g_prop * cos_th_i;
            q_HG = (1 - g_prop^2) / (4 * pi * denom_core * sqrt(denom_core));
            p_val = pdf_Empirical(cos_th_i, param);
            
            weight = weight * (p_val / q_HG); % 解除 1e3 限制
            
            dir = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1 - cos_th_i^2), 2*pi*rand());
            d1 = dir(1); d2 = dir(2); d3 = dir(3);
        end
        P_rx_accum = P_rx_accum + P_packet;
    end
    Time_arr(idx) = toc; PL_arr(idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
    fprintf('N = 10^%d | Time: %.2f s | Path Loss: %.2f dB\n', log10(N_packets), Time_arr(idx), -PL_arr(idx));
end
save('data_exp1_HGPIS.mat', 'N_arr', 'PL_arr', 'Time_arr');

% ================= 辅助函数 =================
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000
        t = max(th(i)*180/pi, 1e-6); sq_t = sqrt(t); t_sq = t*t;
        term = 1 - p.k1*sq_t + p.k2*t - p.k3*t*sq_t + p.k4*t_sq - p.k5*t_sq*sq_t;
        val(i) = exp(p.q_e * term); 
    end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6); sq_t = sqrt(t_deg); t_sq = t_deg * t_deg;
    term = 1 - param.k1*sq_t + param.k2*t_deg - param.k3*t_deg*sq_t + param.k4*t_sq - param.k5*t_sq*sq_t;
    p = exp(param.q_e * term) / param.b_emp_norm;
end
function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2); sp = sin(psi_s); cp = cos(psi_s);
    if denom < 1e-10, nd = [st*cp, st*sp, sign(d(3))*ct]; else
        nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, -st*cp*denom + d(3)*ct]; 
    end
end