%% Exp2: MCS (Traditional LUT Baseline, Extreme O(1) LUT Acceleration)
% 优化声明: 
% 1. 保留三维广义矢量运算，完全兼容非视距(NLOS)与收发未对准场景。
% 2. 引入 O(1) 复杂度等概率反函数查找表(Inverse CDF LUT)，避免差值算法的耗时。
% 3. 解除矩阵运算中的临时内存分配 (如 .^2 和 sum)，恢复 JIT 极速编译。
% 4. 算力优化: 引入 FOV 空间截断剔除策略，切断游走出视场锥的无效光子追踪。
clc; clear; close all;

% % ================= 物理参数初始化 =================
dist_cell = {5:10:85, 5:10:55, 5:5:20}; 
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr = [0.1514, 0.398, 2.190];
num_W = length(water_types); 

N_packets = 1e5; 
n_max = 200; 
lambda = 514e-9; 

% --- 核心物理边界 ---
w0 = 0.002;                                
div_angle = 10 * pi / 180;                
theta_half_div = div_angle / 2; 
Rx_Aperture = 0.05;                        
Rx_FOV = 40 * pi / 180;                     
Rx_Area = pi * (Rx_Aperture / 2)^2;

% 预计算判定常数
cos_FOV_half = cos(Rx_FOV / 2);
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

PL_Cell_MCS = cell(1, num_W);

fprintf('--- 运行 Exp2: MCS (传统 LUT + O(1)极速标量版 + FOV空间截断) ---\n');
for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx); param.coef_b = coef_b_arr(w_idx); param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c; param = calc_haltrin_params(param);
    
    th_axis = linspace(0, pi, 50000); pdf_vals = zeros(size(th_axis));
    for i = 1:length(th_axis), pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i)); end
    cdf_vals = cumtrapz(th_axis, pdf_vals); cdf_vals = cdf_vals / cdf_vals(end);
    
    % =======================================================
    % 构建 O(1) 复杂度的等概率反函数哈希查找表
    % =======================================================
    [cdf_uniq, idx_uniq] = unique(cdf_vals);
    LUT_SIZE = 100000;
    P_grid = linspace(0, 1, LUT_SIZE);
    th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap');
    % 将三角函数运算前置，彻底移出千万级主循环
    cos_th_LUT = cos(th_LUT);
    sin_th_LUT = sin(th_LUT);
    
    dist_arr = dist_cell{w_idx};
    num_D = length(dist_arr);
    PL_arr_temp = zeros(1, num_D);
    
    for d_idx = 1:num_D
        L = dist_arr(d_idx);
        Tx_Pos = [0, 0, 0]; 
        Rx_Pos = [0, L, 0]; 
        Link_Dir = [0, 1, 0]; 
        Rx_Normal = [0, -1, 0];
        u_vec_Tx = [1, 0, 0]; 
        v_vec_Tx = [0, 0, 1]; 
        
        % 提取通用矢量的标量分量，避免结构体或数组索引开销
        Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3);
        Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
        Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3);
        Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
        
        P_rx_accum = 0; tic;
        rng(123456, 'twister'); 
        
        for p = 1:N_packets
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
            pos = [Tx + r0*cos(phi0)*u_vec_Tx(1) + r0*sin(phi0)*v_vec_Tx(1), ...
                   Ty + r0*cos(phi0)*u_vec_Tx(2) + r0*sin(phi0)*v_vec_Tx(2), ...
                   Tz + r0*cos(phi0)*u_vec_Tx(3) + r0*sin(phi0)*v_vec_Tx(3)];
            
            U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
            ct_init = cos(U_init); st_init = sin(U_init);
            dir = rotate_direction_fast(Link_Dir, ct_init, st_init, 2*pi*rand());
            
            weight = 1.0; P_packet = 0;
            
            % =======================================================
            % 0阶直射路径 (泛化三维矢量零分配解算)
            % =======================================================
            cos_th = dir(1)*Lx + dir(2)*Ly + dir(3)*Lz;
            if cos_th > 0
                d_plane = ((Rx - pos(1))*Lx + (Ry - pos(2))*Ly + (Rz - pos(3))*Lz) / cos_th; 
                pos_end_1 = pos(1) + dir(1) * d_plane; 
                pos_end_2 = pos(2) + dir(2) * d_plane; 
                pos_end_3 = pos(3) + dir(3) * d_plane; 
                
                cos_rx_tilt = abs(dir(1)*Nx + dir(2)*Ny + dir(3)*Nz);
                
                if cos_rx_tilt >= cos_FOV_half
                    r_hit_sq = (pos_end_1 - Rx)^2 + (pos_end_2 - Ry)^2 + (pos_end_3 - Rz)^2;
                    if r_hit_sq <= Rx_Aperture_half_sq
                        P_packet = P_packet + exp(-param.coef_c * d_plane);
                    end
                end
            end
            
            for ord = 1:n_max
                d_step = -log(rand()) / param.coef_b; 
                pos(1) = pos(1) + dir(1) * d_step; 
                pos(2) = pos(2) + dir(2) * d_step; 
                pos(3) = pos(3) + dir(3) * d_step; 
                
                weight = weight * exp(-param.coef_a * d_step);
                
                z_progression = (pos(1) - Tx)*Lx + (pos(2) - Ty)*Ly + (pos(3) - Tz)*Lz;
                if z_progression >= L || weight < 1e-20
                    break; 
                end 
                
                % =======================================================
                % Local Estimation 的无临时数组标量展开
                % =======================================================
                vec2rx_1 = Rx - pos(1); 
                vec2rx_2 = Ry - pos(2); 
                vec2rx_3 = Rz - pos(3); 
                
                d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2;
                d2rx = sqrt(d2rx_sq); 
                
                dir2rx_1 = vec2rx_1 / d2rx;
                dir2rx_2 = vec2rx_2 / d2rx;
                dir2rx_3 = vec2rx_3 / d2rx;
                
                cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                
                if cos_inc < cos_FOV_half
                    break; 
                end
                
                cos_scatter = dir(1)*dir2rx_1 + dir(2)*dir2rx_2 + dir(3)*dir2rx_3;
                P_packet = P_packet + weight * min(1, pdf_Empirical(cos_scatter, param) * (Rx_Area / d2rx_sq * cos_inc)) * exp(-param.coef_c * d2rx);
                
                % =======================================================
                % O(1) 极速哈希表角度抽样
                % =======================================================
                u_rand = rand();
                idx_lut = floor(u_rand * (LUT_SIZE - 1)) + 1;
                ct_s = cos_th_LUT(idx_lut);
                st_s = sin_th_LUT(idx_lut);
                
                dir = rotate_direction_fast(dir, ct_s, st_s, 2*pi*rand());
            end
            P_rx_accum = P_rx_accum + P_packet;
        end
        PL_arr_temp(d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
        fprintf('  %s | L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', water_types{w_idx}, L, toc, -PL_arr_temp(d_idx));
    end
    PL_Cell_MCS{w_idx} = PL_arr_temp;
end
save('data_exp2_MCS.mat', 'dist_cell', 'PL_Cell_MCS');

% ================= 核心辅助函数区域 =================

function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000
        t = max(th(i)*180/pi, 1e-6); 
        sq_t = sqrt(t); t_sq = t*t;
        term = 1 - p.k1*sq_t + p.k2*t - p.k3*t*sq_t + p.k4*t_sq - p.k5*t_sq*sq_t;
        val(i) = exp(p.q_e * term); 
    end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); 
    param = p;
end

function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6); 
    sq_t = sqrt(t_deg); t_sq = t_deg * t_deg;
    term = 1 - param.k1*sq_t + param.k2*t_deg - param.k3*t_deg*sq_t + param.k4*t_sq - param.k5*t_sq*sq_t;
    p = exp(param.q_e * term) / param.b_emp_norm;
end

function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2);
    sp = sin(psi_s); cp = cos(psi_s);
    if denom < 1e-10
        nd = [st*cp, st*sp, sign(d(3))*ct]; 
    else
        nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, ...
              st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, ...
             -st*cp*denom + d(3)*ct]; 
    end
    nd = nd / norm(nd);
end