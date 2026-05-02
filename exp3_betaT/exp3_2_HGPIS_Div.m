%% Exp3.3.1: HG-PIS Algorithm for Different Divergence Angles (No Turbulence)
clc; clear; close all;

N_packets = 1e5; n_max = 200; 
dist_axis = 5:5:60; num_dist = length(dist_axis);

w0 = 0.002; Rx_Aperture = 0.05; Rx_FOV = 180 * pi / 180; Rx_Area = pi * (Rx_Aperture / 2)^2;
cos_FOV_half = cos(Rx_FOV / 2); Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

% [修复点] 显式声明独立标量 coef_c 以供后续保存
coef_c = 0.1514;
param.coef_a = 0.114; param.coef_b = 0.0374; param.coef_c = coef_c;
param.albedo = param.coef_b / param.coef_c;
param = calc_haltrin_params(param);
P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max; 
g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);

div_angles_deg = [0.1, 1.0, 5.0]; num_div = length(div_angles_deg);
PL_Cell_HG = cell(1, num_div);

fprintf('--- 运行 Exp3.3.1: Adaptive HG-PIS (束散角遍历) ---\n');
for i = 1:num_div
    theta_half_div = (div_angles_deg(i) * pi / 180) / 2;
    fprintf('--- 当前束散角: %.2f deg ---\n', div_angles_deg(i));
    
    PL_arr_temp = zeros(1, num_dist);
    for d_idx = 1:num_dist
        L = dist_axis(d_idx);
        Tx = 0; Ty = 0; Tz = 0; Rx = 0; Ry = L; Rz = 0;
        Lx = 0; Ly = 1; Lz = 0; Nx = 0; Ny = -1; Nz = 0;
        Ux = 1; Uy = 0; Uz = 0; Vx = 0; Vy = 0; Vz = 1;
        
        P_rx_accum = 0; tic; rng(123456, 'twister'); 
        
        for p = 1:N_packets
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand(); cp0 = cos(phi0); sp0 = sin(phi0);
            p1 = Tx + r0*cp0*Ux + r0*sp0*Vx; p2 = Ty + r0*cp0*Uy + r0*sp0*Vy; p3 = Tz + r0*cp0*Uz + r0*sp0*Vz;
            
            U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
            dir = rotate_direction_fast([Lx,Ly,Lz], cos(U_init), sin(U_init), 2*pi*rand());
            d1 = dir(1); d2 = dir(2); d3 = dir(3);
            weight = 1.0; P_packet = 0;
            
% ================= 直射路径处理 (引入几何拓展) =================
            cos_th = d1*Lx + d2*Ly + d3*Lz;
            if cos_th > 0
                % 计算射线到接收平面的距离及交点
                d_plane = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / cos_th; 
                pos_end_1 = p1 + d1 * d_plane; 
                pos_end_2 = p2 + d2 * d_plane; 
                pos_end_3 = p3 + d3 * d_plane; 
                
                cos_rx_tilt = abs(d1*Nx + d2*Ny + d3*Nz);
                if cos_rx_tilt >= cos_FOV_half
                    % 计算光子落点到接收机中心的径向距离
                    r_wp = sqrt((pos_end_1 - Rx)^2 + (pos_end_2 - Ry)^2 + (pos_end_3 - Rz)^2);
                    
                    % 考虑倾斜入射时的有效接收半径
                    r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);
                    
                    % 计算纯几何拓展下的光斑半径 (无湍流影响)
                    W_geo = w0 + d_plane * tan(theta_half_div);
                    W_spot = W_geo; 
                    
                    % 利用 Marcum Q 函数计算空间展宽后的孔径截获概率
                    if W_spot > 1e-6
                        recv_frac = 1 - marcumq(2 * r_wp / W_spot, 2 * r_eff / W_spot);
                    else
                        recv_frac = double(r_wp <= r_eff);
                    end
                    
                    % 累加直射路径衰减与接收权重
                    P_packet = P_packet + exp(-param.coef_c * d_plane) * recv_frac;
                end
            end
            % ===============================================================
            
            for ord = 1:n_max
                d_step = -log(rand()) / param.coef_b; 
                p1 = p1 + d1 * d_step; p2 = p2 + d2 * d_step; p3 = p3 + d3 * d_step; 
                weight = weight * exp(-param.coef_a * d_step);
                if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L, break; end 
                if weight < 1e-9, if rand() > 0.1, break; else, weight = weight * 10; end; end
                
                vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; 
                d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; d2rx = sqrt(d2rx_sq); 
                dir2rx_1 = vec2rx_1 / d2rx; dir2rx_2 = vec2rx_2 / d2rx; dir2rx_3 = vec2rx_3 / d2rx;
                cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                
                if cos_inc >= cos_FOV_half
                    omega = Rx_Area / d2rx_sq * cos_inc; cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                    P_packet = P_packet + weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * d2rx);
                end
                
                xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); 
                cos_th_i = (1 + g_prop^2 - term^2) / (2 * g_prop);
                if cos_th_i > 1, cos_th_i = 1; elseif cos_th_i < -1, cos_th_i = -1; end
                denom_core = 1 + g_prop^2 - 2 * g_prop * cos_th_i;
                q_HG = (1 - g_prop^2) / (4 * pi * denom_core * sqrt(denom_core));
                p_val = pdf_Empirical(cos_th_i, param);
                weight = weight * (p_val / q_HG);
                
                dir = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1 - cos_th_i^2), 2*pi*rand());
                d1 = dir(1); d2 = dir(2); d3 = dir(3);
            end
            P_rx_accum = P_rx_accum + P_packet;
        end
        PL_arr_temp(d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
        fprintf('    L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', L, toc, -PL_arr_temp(d_idx));
    end
    PL_Cell_HG{i} = PL_arr_temp;
end
save('data_exp3_div_HGPIS.mat', 'dist_axis', 'div_angles_deg', 'PL_Cell_HG', 'coef_c');
fprintf('已保存数据至 data_exp3_div_HGPIS.mat\n');

% --- 辅助函数 ---
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000, t = max(th(i)*180/pi, 1e-6); sq_t = sqrt(t); t_sq = t*t; term = 1 - p.k1*sq_t + p.k2*t - p.k3*t*sq_t + p.k4*t_sq - p.k5*t_sq*sq_t; val(i) = exp(p.q_e * term); end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6); sq_t = sqrt(t_deg); t_sq = t_deg * t_deg;
    term = 1 - param.k1*sq_t + param.k2*t_deg - param.k3*t_deg*sq_t + param.k4*t_sq - param.k5*t_sq*sq_t; p = exp(param.q_e * term) / param.b_emp_norm;
end
function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2); sp = sin(psi_s); cp = cos(psi_s);
    if denom < 1e-10, nd = [st*cp, st*sp, sign(d(3))*ct]; else, nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, -st*cp*denom + d(3)*ct]; end
end