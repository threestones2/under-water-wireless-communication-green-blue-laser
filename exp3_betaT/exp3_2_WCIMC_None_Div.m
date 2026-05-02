%% Exp3.3.1: WCI-MC Algorithm for Different Divergence Angles (No Turbulence Baseline)
% 理论更新: 修正MC框架下的无湍流基准。移除单光子的几何光斑拓展(已由初始采样暗含)，
% 将接收波包展宽严格约束为由有效接收孔径决定的衍射极限。
clc; clear; close all;

% 强制切换至本脚本所在的物理目录
cd(fileparts(mfilename('fullpath')));

N_packets = 1e5; n_max = 200; 
dist_axis = 5:5:60; num_dist = length(dist_axis);

% 光学及物理常数补充
lambda = 514e-9; k_wave = 2 * pi / lambda;

w0 = 0.002; Rx_Aperture = 0.01; Rx_FOV = 5 * pi / 180; Rx_Area = pi * (Rx_Aperture / 2)^2;
cos_FOV_half = cos(Rx_FOV / 2); r_rx = Rx_Aperture / 2;

coef_c = 0.1514;
param.coef_a = 0.114; param.coef_b = 0.0374; param.coef_c = coef_c;
param.albedo = param.coef_b / param.coef_c;
param = calc_haltrin_params(param);

P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max; 
g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);

div_angles_deg = [0.05, 0.1, 1]; num_div = length(div_angles_deg);
PL_Cell_WCI_None = cell(1, num_div);

fprintf('--- 运行 Exp3.3.1: WCI-MC (无湍流基准 - 纯MC追踪 + 衍射极限) ---\n');
for i = 1:num_div
    theta_half_div = (div_angles_deg(i) * pi / 180) / 2;
    fprintf('--- 当前束散角: %.2f deg ---\n', div_angles_deg(i));
    
    PL_arr_temp = zeros(1, num_dist);
    for d_idx = 1:num_dist
        L = dist_axis(d_idx);
        rng(123456+L, 'twister'); 
        Tx_Pos = [0, 0, 0]; Rx_Pos = [0, L, 0]; Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
        Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3);
        Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
        Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3);
        Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
        Ux = 1; Uy = 0; Uz = 0; Vx = 0; Vy = 0; Vz = 1;
        
        P_accum = 0; tic; 
        rng(123456, 'twister'); 
        
        for p = 1:N_packets
            % 初始发射参数 (已包含几何束散角的统计采样)
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand(); cp0 = cos(phi0); sp0 = sin(phi0);
            p1 = Tx + r0*cp0*Ux + r0*sp0*Vx; p2 = Ty + r0*cp0*Uy + r0*sp0*Vy; p3 = Tz + r0*cp0*Uz + r0*sp0*Vz;
            
            U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
            dir = rotate_direction_fast([Lx,Ly,Lz], cos(U_init), sin(U_init), 2*pi*rand());
            d1 = dir(1); d2 = dir(2); d3 = dir(3);
            weight = 1.0; P_packet = 0;
            
            % ================= 直射路径处理 =================
            cos_th = d1*Lx + d2*Ly + d3*Lz;
            if cos_th > 0
                d_plane = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / cos_th;
                pos_end_1 = p1 + d1 * d_plane; 
                pos_end_2 = p2 + d2 * d_plane; 
                pos_end_3 = p3 + d3 * d_plane; 
                
                cos_rx_tilt = abs(d1*Nx + d2*Ny + d3*Nz);
                if cos_rx_tilt >= cos_FOV_half
                    vec_x = pos_end_1 - Rx; vec_y = pos_end_2 - Ry; vec_z = pos_end_3 - Rz;
                    cx = vec_y*d3 - vec_z*d2; cy = vec_z*d1 - vec_x*d3; cz = vec_x*d2 - vec_y*d1;
                    r_wander_perp = sqrt(cx^2 + cy^2 + cz^2);
                    
                    r_eff = r_rx * sqrt(cos_rx_tilt);
                    
                    % 修正: 使用基于接收机有效孔径的衍射极限展宽，取代宏观几何束散展宽
                    W_spot_diff = d_plane / (k_wave * r_eff);
                    
                    if W_spot_diff > 1e-8
                        recv_frac = 1 - marcumq(2 * r_wander_perp / W_spot_diff, 2 * r_eff / W_spot_diff);
                    else
                        recv_frac = double(r_wander_perp <= r_eff);
                    end
                    P_packet = P_packet + exp(-param.coef_c * d_plane) * recv_frac;
                end
            end
            
            % ================= 多重散射演化 =================
            for ord = 1:n_max
                d_step = -log(rand()) / param.coef_b; 
                p1 = p1 + d1 * d_step; p2 = p2 + d2 * d_step; p3 = p3 + d3 * d_step; 
                weight = weight * exp(-param.coef_a * d_step);
                
                if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L || weight < 1e-15, break; end 
                
                vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; 
                d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; d2rx = sqrt(d2rx_sq); 
                dir2rx_1 = vec2rx_1 / d2rx; dir2rx_2 = vec2rx_2 / d2rx; dir2rx_3 = vec2rx_3 / d2rx;
                cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                
                if cos_inc < cos_FOV_half
                    break; 
                end
                
                if cos_inc >= cos_FOV_half
                    omega = Rx_Area / d2rx_sq * cos_inc; cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                    base_w = weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * d2rx);
                    if base_w > 1e-15
                        d_plane_v = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / (dir2rx_1*Lx + dir2rx_2*Ly + dir2rx_3*Lz);
                        p1_v = p1 + dir2rx_1 * d_plane_v; p2_v = p2 + dir2rx_2 * d_plane_v; p3_v = p3 + dir2rx_3 * d_plane_v;
                        vec_x_v = p1_v - Rx; vec_y_v = p2_v - Ry; vec_z_v = p3_v - Rz;
                        r_wp = sqrt(vec_x_v^2 + vec_y_v^2 + vec_z_v^2);
                        r_eff_v = r_rx * sqrt(cos_inc);
                        
                        % 修正: 散射支路同样适用基础的孔径衍射极限评估
                        W_ST_diff = 2 * d_plane_v / (k_wave * r_eff_v);
                        
                        if W_ST_diff > 1e-8
                            point_loss = 1 - marcumq(2 * r_wp / W_ST_diff, 2 * r_eff_v / W_ST_diff);
                        else
                            point_loss = double(r_wp <= r_eff_v);
                        end
                        
                        P_packet = P_packet + base_w * point_loss;
                    end
                end
                
                xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); 
                cos_th_i = (1 + g_prop^2 - term^2) / (2 * g_prop);
                if cos_th_i > 1, cos_th_i = 1; elseif cos_th_i < -1, cos_th_i = -1; end
                denom_core = 1 + g_prop^2 - 2 * g_prop * cos_th_i;
                q_HG = (1 - g_prop^2) / (4 * pi * denom_core * sqrt(denom_core));
                p_val = pdf_Empirical(cos_th_i, param);
                
                weight = weight * max(0.5, min(2, p_val / q_HG));
                
                dir = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1 - cos_th_i^2), 2*pi*rand());
                d1 = dir(1); d2 = dir(2); d3 = dir(3);
            end
            P_accum = P_accum + P_packet;
        end
        
        PL_arr_temp(d_idx) = 10 * log10(max(P_accum / N_packets, 1e-300));
        fprintf('    L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', L, toc, -PL_arr_temp(d_idx));
    end
    PL_Cell_WCI_None{i} = PL_arr_temp;
end
save('data_exp3_div_WCIMC_None.mat', 'dist_axis', 'div_angles_deg', 'PL_Cell_WCI_None', 'coef_c');
fprintf('已保存无湍流解析基准数据至 data_exp3_div_WCIMC_None.mat\n');

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