%% Exp2: Adaptive HG-PIS (Without Turbulence Baseline, Zero-Allocation Scalarized)
% 优化声明: 
% 1. 保留三维广义矢量运算，完全兼容非视距(NLOS)与收发未对准场景。
% 2. 引入不等式余弦域等价判定、分数次幂降维与解析正余弦映射。
% 3. 全局实施无分配标量展开 (Zero-Allocation Scalarization)，消灭所有 dot, norm, sum 及隐式数组创建。
clc; clear; close all;

% ================= 物理参数初始化 =================
dist_cell = {5:10:85, 5:10:55, 5:5:20}; 
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr = [0.1514, 0.398, 2.190];

% ================= 物理参数初始化 =================
% dist_cell = {5:10:65, 5:10:55, 5:5:20}; 
% water_types = {'Jerlov IB', 'Jerlov II', 'Jerlov 5C'};
% coef_a_arr = [0.0389, 0.0389, 0.1376]; 
% coef_b_arr = [0.0579, 0.4020, 1.4928]; 
% coef_c_arr = [0.0968, 0.4409, 1.6304];

num_W = length(water_types); 
N_packets = 1e5; 
n_max = 200; 

% --- 核心物理与几何边界 ---
w0 = 0.02;                                
div_angle = 0.1 * pi / 180;                
theta_half_div = div_angle / 2; 
Rx_Aperture = 0.01;                        
Rx_FOV = 5 * pi / 180;                     
Rx_Area = pi * (Rx_Aperture / 2)^2;

% [底层优化] 预计算判定常数
cos_FOV_half = cos(Rx_FOV / 2);
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

PL_Cell_HG = cell(1, num_W);

fprintf('--- 运行 Exp2: Adaptive HG-PIS (无湍流, 零分配标量加速版) ---\n');
for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx); param.coef_b = coef_b_arr(w_idx); param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c; 
    param = calc_haltrin_params(param);
    
    P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max; 
    g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);
    
    dist_arr = dist_cell{w_idx};
    num_D = length(dist_arr);
    PL_arr_temp = zeros(1, num_D);
    
    for d_idx = 1:num_D
        L = dist_arr(d_idx);
        Tx_Pos = [0, 0, 0]; Rx_Pos = [0, L, 0]; Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
        u_vec_Tx = [1, 0, 0]; v_vec_Tx = [0, 0, 1]; 
        
        % [底层优化] 预提取结构标量，避免循环内的数组解包
        Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3);
        Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
        Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3);
        Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
        Ux = u_vec_Tx(1); Uy = u_vec_Tx(2); Uz = u_vec_Tx(3);
        Vx = v_vec_Tx(1); Vy = v_vec_Tx(2); Vz = v_vec_Tx(3);
        
        P_rx_accum = 0; tic;
        rng(123456, 'twister'); 
        
        for p = 1:N_packets
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
            cp0 = cos(phi0); sp0 = sin(phi0);
            
            % 标量化初始位置
            p1 = Tx + r0*cp0*Ux + r0*sp0*Vx;
            p2 = Ty + r0*cp0*Uy + r0*sp0*Vy;
            p3 = Tz + r0*cp0*Uz + r0*sp0*Vz;
            
            U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
            dir = rotate_direction_fast(Link_Dir, cos(U_init), sin(U_init), 2*pi*rand());
            d1 = dir(1); d2 = dir(2); d3 = dir(3);
            
            weight = 1.0; P_packet = 0;
            
            % =======================================================
            % 0阶直射路径 (零分配标量代数展开)
            % =======================================================
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
            
            % =======================================================
            % 多重散射分支 (HG-PIS 零分配标量化)
            % =======================================================
            for ord = 1:n_max
                d_step = -log(rand()) / param.coef_b; 
                p1 = p1 + d1 * d_step; 
                p2 = p2 + d2 * d_step; 
                p3 = p3 + d3 * d_step; 
                weight = weight * exp(-param.coef_a * d_step);
                
                if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L
                    break; 
                end 
                
                if weight < 1e-9
                    if rand() > 0.1, 
                        break; 
                    else, 
                        weight = weight * 10; 
                    end
                end
                
                % --- 虚拟强制接收 (Local Estimation) ---
                vec2rx_1 = Rx - p1; 
                vec2rx_2 = Ry - p2; 
                vec2rx_3 = Rz - p3; 
                
                d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2;
                d2rx = sqrt(d2rx_sq); 
                dir2rx_1 = vec2rx_1 / d2rx;
                dir2rx_2 = vec2rx_2 / d2rx;
                dir2rx_3 = vec2rx_3 / d2rx;
                
                cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                if cos_inc >= cos_FOV_half
                    omega = Rx_Area / d2rx_sq * cos_inc;
                    cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                    P_packet = P_packet + weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * d2rx);
                end
                
                % --- 自适应 HG-PIS ---
                xi = rand(); 
                term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); 
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
        fprintf('  %s | L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', water_types{w_idx}, L, toc, -PL_arr_temp(d_idx));
    end
    PL_Cell_HG{w_idx} = PL_arr_temp;
end
save('data_exp2_HGPIS.mat', 'dist_cell', 'PL_Cell_HG');

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
    sq_t = sqrt(t_deg);
    t_sq = t_deg * t_deg;
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
end