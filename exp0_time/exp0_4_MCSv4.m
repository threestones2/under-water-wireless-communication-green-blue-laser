%% 3. MCS (Local Estimation + O(1) Fast Hash LUT Accelerated)
clc; clear; close all;
% 强制将 MATLAB 的当前工作文件夹切换至本脚本所在的物理目录
cd(fileparts(mfilename('fullpath')));
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr = [0.1514, 0.398, 2.190];
num_photons = 1e5; max_scattering = 200; 

% --- 实际硬件系统几何与光学参数 (对齐 exp3) ---
r = 10;                      
lambda = 514e-9;             
w0 = 0.002;                  
div_angle = 0.1 * pi / 180;  
theta_half_div_physics = div_angle / 2;
rx_radius = 0.01;           
phi2 = 5 * pi / 180;         
rx_area = pi*(rx_radius)^2; 

rx_pos = [0, 0, r]; miu0 = [0, 0, 1]; miur = [0, 0, -1];

if abs(miu0(3)) < 0.9, up_temp_Tx = [0, 0, 1]; else, up_temp_Tx = [1, 0, 0]; end
u_vec_Tx = cross(up_temp_Tx, miu0); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(miu0, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

fprintf('--- 启动 MCS 仿真 (传统物理基准 - O(1) 哈希查表法极限加速) ---\n');

for wt_idx = 1:3
    k_a = coef_a_arr(wt_idx); k_s = coef_b_arr(wt_idx); k_e = coef_c_arr(wt_idx);
    param.coef_b = k_s; param.albedo = k_s/k_e;
    param = calc_haltrin_params(param);
    
    % =====================================================================
    % 【核心优化】：构建 O(1) 等概率反函数哈希查找表
    % =====================================================================
    th_axis = linspace(0, pi, 50000); pdf_vals = zeros(size(th_axis));
    for i = 1:length(th_axis)
        pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i));
    end
    cdf_vals = cumtrapz(th_axis, pdf_vals); 
    cdf_vals = cdf_vals / cdf_vals(end);
    
    [cdf_uniq, idx_uniq] = unique(cdf_vals);
    LUT_SIZE = 100000;
    P_grid = linspace(0, 1, LUT_SIZE);
    
    % 在循环外部进行一次性的逆变换插值
    th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap');
    
    % 彻底消灭循环内的三角函数调用，将计算前置
    cos_th_LUT = cos(th_LUT);
    sin_th_LUT = sin(th_LUT);
    % =====================================================================
    
    P_rx_total_accum = 0;
    tic;
    for i = 1:num_photons
        r0 = w0 * sqrt(-0.5*log(rand())); 
        phi0 = 2*pi*rand();
        pos_local = r0*cos(phi0)*u_vec_Tx + r0*sin(phi0)*v_vec_Tx;
        pos_init = [0, 0, 0] + pos_local;
        
        U_init = theta_half_div_physics * sqrt(-0.5 * log(rand()));
        psi_init = 2 * pi * rand();
        
        % 初始发射也采用高速版的旋转函数
        dir_init = rotate_direction_fast(miu0, cos(U_init), sin(U_init), psi_init);
        
        P_packet_rx = 0;
        
        % 2. 0阶直射路径 (对齐 exp3 的纯几何截断)
        pos_centroid = [0, 0, 0];
        dir_centroid = miu0;
        
        cos_theta = dot(dir_centroid, miu0);
        if cos_theta > 0
            dist_to_plane = r / cos_theta;
            pos_end = pos_centroid + dir_centroid * dist_to_plane;
            r_wander = norm(pos_end - rx_pos);
            cos_rx_tilt = abs(dot(dir_centroid, miur));
            
            if acos(cos_rx_tilt) <= phi2 / 2
                if r_wander <= rx_radius
                    w_ballistic = exp(-k_e * dist_to_plane);
                    P_packet_rx = P_packet_rx + w_ballistic;
                end
            end
        end
        
        % 3. 1~n阶散射演化
        pos = pos_init; 
        direction = dir_init; 
        survival_prob = 1.0;
        
        for scnt = 1:max_scattering
            delta_s = -log(rand) / k_s;
            pos = pos + direction * delta_s;
            survival_prob = survival_prob * exp(-k_a * delta_s);
            
            if survival_prob < 1e-9
                if rand() > 0.1, break; else, survival_prob = survival_prob * 10; end
            end 
            if dot(pos, miu0) >= r, break; end
            
            vec_rx = rx_pos - pos; dist_rx = norm(vec_rx);
            if dot(vec_rx/dist_rx, -miur) >= cos(phi2/2)
                cos_theta_s = dot(vec_rx/dist_rx, direction);
                P_phase = pdf_Empirical(cos_theta_s, param);
                solid_angle = rx_area * abs(dot(vec_rx/dist_rx, -miur)) / (dist_rx^2);
                prob_reception = min(1.0, P_phase * solid_angle);
                P_packet_rx = P_packet_rx + survival_prob * prob_reception * exp(-k_e * dist_rx);
            end
            
            % =========================================================================
            % 【核心优化】：O(1) 极速哈希表抽样
            % 生成 [1, LUT_SIZE] 的整数索引，直接提取预计算的正余弦值，彻底绕开插值
            % =========================================================================
            u_rand = rand();
            idx_lut = floor(u_rand * (LUT_SIZE - 1)) + 1;
            ct_s = cos_th_LUT(idx_lut);
            st_s = sin_th_LUT(idx_lut);
            
            direction = rotate_direction_fast(direction, ct_s, st_s, 2 * pi * rand());
            % =========================================================================
        end
        
        P_rx_total_accum = P_rx_total_accum + P_packet_rx;
    end
    t_run = toc;
    
    P_rx_mean = P_rx_total_accum / num_photons;
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
    cos_theta = max(min(cos_theta, 1), -1); 
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6);
    term = 1 + (-1)^1 * param.k1 * t_deg.^0.5 + ...
           (-1)^2 * param.k2 * t_deg.^1.0 + ...
           (-1)^3 * param.k3 * t_deg.^1.5 + ...
           (-1)^4 * param.k4 * t_deg.^2.0 + ...
           (-1)^5 * param.k5 * t_deg.^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end

% 改用接收正余弦值的极速旋转函数
function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2);
    sp = sin(psi_s); cp = cos(psi_s);
    if denom < 1e-10, nd = [st*cp, st*sp, sign(d(3))*ct];
    else, nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, ...
                st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, ...
                -st*cp*denom + d(3)*ct]; end
    nd = nd / norm(nd);
end