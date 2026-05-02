%% Exp3.3.1: Standard MCS Algorithm for Different Divergence Angles
% 理论更新: 加入FOV超出截断机制；直射路径对齐基于衍射极限的孔径平滑接收，散射路径保留标准立体角局部估计。
clc; clear; close all;

% 强制切换至本脚本所在的物理目录
cd(fileparts(mfilename('fullpath')));

N_packets = 1e5; n_max = 200; 
dist_axis = 5:5:60; num_dist = length(dist_axis);

% 物理常数补充
lambda = 514e-9; k_wave = 2 * pi / lambda;

w0 = 0.002; Rx_Aperture = 0.01; Rx_FOV = 5 * pi / 180; Rx_Area = pi * (Rx_Aperture / 2)^2;
cos_FOV_half = cos(Rx_FOV / 2); r_rx = Rx_Aperture / 2;

% 显式声明独立标量 coef_c 以供后续保存
coef_c = 0.1514;
param.coef_a = 0.114; param.coef_b = 0.0374; param.coef_c = coef_c;
param.albedo = param.coef_b / param.coef_c;
param = calc_haltrin_params(param);

% 预构建 O(1) LUT
th_axis = linspace(0, pi, 50000); pdf_vals = zeros(size(th_axis));
for i = 1:length(th_axis), pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i)); end
cdf_vals = cumtrapz(th_axis, pdf_vals); cdf_vals = cdf_vals / cdf_vals(end);
[cdf_uniq, idx_uniq] = unique(cdf_vals); LUT_SIZE = 100000; P_grid = linspace(0, 1, LUT_SIZE);
th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap');
cos_th_LUT = cos(th_LUT); sin_th_LUT = sin(th_LUT);

div_angles_deg = [0.05, 0.1, 1.0]; num_div = length(div_angles_deg);
PL_Cell_MCS = cell(1, num_div);

fprintf('--- 运行 Exp3.3.1: Standard MCS (束散角遍历 - FOV截断与衍射对齐) ---\n');
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
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand(); 
            p1 = Tx + r0*cos(phi0)*Ux + r0*sin(phi0)*Vx; p2 = Ty + r0*cos(phi0)*Uy + r0*sin(phi0)*Vy; p3 = Tz + r0*cos(phi0)*Uz + r0*sin(phi0)*Vz;
            
            U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
            dir = rotate_direction_fast([Lx,Ly,Lz], cos(U_init), sin(U_init), 2*pi*rand());
            d1 = dir(1); d2 = dir(2); d3 = dir(3);
            weight = 1.0; P_packet = 0;
            
            % ================= 直射路径处理 =================
            cos_th = d1*Lx + d2*Ly + d3*Lz;
            if cos_th > 0
                d_plane = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / cos_th; 
                pos_end_1 = p1 + d1 * d_plane; pos_end_2 = p2 + d2 * d_plane; pos_end_3 = p3 + d3 * d_plane; 
                cos_rx_tilt = abs(d1*Nx + d2*Ny + d3*Nz);
                if cos_rx_tilt >= cos_FOV_half
                    vec_x = pos_end_1 - Rx; vec_y = pos_end_2 - Ry; vec_z = pos_end_3 - Rz;
                    % 计算光子落点到接收机中心的垂直偏移距离
                    cx = vec_y*d3 - vec_z*d2; cy = vec_z*d1 - vec_x*d3; cz = vec_x*d2 - vec_y*d1;
                    r_wander_perp = sqrt(cx^2 + cy^2 + cz^2);
                    
                    r_eff = r_rx * sqrt(cos_rx_tilt);
                    
                    % 采用与无湍流基准一致的衍射极限展宽
                    W_spot_diff = 2 * d_plane / (k_wave * r_eff);
                    
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
                
                % 权重低于截断阈值或超出传输轴线则终止
                if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L || weight < 1e-15, break; end 
                
                vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; 
                d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; d2rx = sqrt(d2rx_sq); 
                dir2rx_1 = vec2rx_1 / d2rx; dir2rx_2 = vec2rx_2 / d2rx; dir2rx_3 = vec2rx_3 / d2rx;
                cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                
                % 新增: FOV 隔离截断
                if cos_inc < cos_FOV_half
                    break; 
                end
                
                if cos_inc >= cos_FOV_half
                    cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                    % 保留标准 MCS 原本的立体角局部估计法
                    P_packet = P_packet + weight * min(1, pdf_Empirical(cos_scatter, param) * (Rx_Area / d2rx_sq * cos_inc)) * exp(-param.coef_c * d2rx);
                end
                
                % 采用 LUT 纯随机采样散射角（非重要性采样，无需设定阈值比）
                u_rand = rand(); idx_lut = floor(u_rand * (LUT_SIZE - 1)) + 1;
                ct_s = cos_th_LUT(idx_lut); st_s = sin_th_LUT(idx_lut);
                dir = rotate_direction_fast([d1, d2, d3], ct_s, st_s, 2*pi*rand());
                d1 = dir(1); d2 = dir(2); d3 = dir(3);
            end
            P_rx_accum = P_rx_accum + P_packet;
        end
        PL_arr_temp(d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
        fprintf('    L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', L, toc, -PL_arr_temp(d_idx));
    end
    PL_Cell_MCS{i} = PL_arr_temp;
end
save('data_exp3_div_MCS.mat', 'dist_axis', 'div_angles_deg', 'PL_Cell_MCS', 'coef_c');
fprintf('已保存数据至 data_exp3_div_MCS.mat\n');

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
    nd = nd / norm(nd);
end