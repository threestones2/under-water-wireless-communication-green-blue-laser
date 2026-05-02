%% Exp2: WCI-MC_None (No Turbulence Baseline via Unified Wavelet Degeneration)
% 架构更新: 
% 1. 理论统一: 废除0阶宏观解析 W_geo，将直射路径拉回 MC 循环。宏观光斑发散由光子角采样自然涌现。
% 2. 极限退化: 无湍流态下 W_ST -> 0，波前补偿完全退化为硬孔径 (Hard Aperture) O(1) 几何求交判定。
% 3. 方差对齐: 确保本基准模型与湍流 WCI-MC 模型的随机数序列与采样方差完全一致。
% 4. 算力优化: 引入 FOV 空间截断剔除策略，彻底切断游走出视场外的无效光子追踪。

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= 物理参数初始化 =================
dist_cell = {5:10:85, 5:10:55, 5:5:20}; 
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr = [0.1514, 0.398, 2.190];
num_W = length(water_types);

N_packets = 1e5; n_max = 200; 
lambda = 514e-9; k_wave = 2 * pi / lambda;

% --- 核心边界 ---
w0 = 0.002;                                
div_angle = 0.1 * pi / 180;                
theta_half_div = div_angle / 2; 
Rx_Aperture = 0.05;                        
Rx_FOV = 20 * pi / 180;                     
Rx_Area = pi * (Rx_Aperture / 2)^2;
cos_FOV_half = cos(Rx_FOV / 2);
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

PL_Cell_WCI_None = cell(1, num_W);

fprintf('--- 运行 Exp2: WCI-MC (无湍流物理基准 | 统一波束微元极限退化) ---\n');

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
        
        Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3);
        Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
        Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3);
        Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
        Ux = u_vec_Tx(1); Uy = u_vec_Tx(2); Uz = u_vec_Tx(3);
        Vx = v_vec_Tx(1); Vy = v_vec_Tx(2); Vz = v_vec_Tx(3);
        
        P_rx_accum = 0; 
        tic; 
        rng(123456, 'twister'); 
        
        for p = 1:N_packets
            % =================================================================
            % [初始抽样]: 发射端空间和角度抽样，自然涌现 W_geo
            % =================================================================
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
            p1 = Tx + r0*cos(phi0)*Ux + r0*sin(phi0)*Vx;
            p2 = Ty + r0*cos(phi0)*Uy + r0*sin(phi0)*Vy;
            p3 = Tz + r0*cos(phi0)*Uz + r0*sin(phi0)*Vz;
            
            U_init = theta_half_div * sqrt(-0.5 * log(rand()));
            dir = rotate_direction_fast(Link_Dir, cos(U_init), sin(U_init), 2*pi*rand());
            d1 = dir(1); d2 = dir(2); d3 = dir(3);
            
            weight = 1.0; 
            P_packet = 0;
            
            % =================================================================
            % [0阶直射路径]: 无湍流环境极限退化，使用解析几何求交验证硬孔径
            % =================================================================
            cos_th = d1*Lx + d2*Ly + d3*Lz;
            if cos_th > 1e-10
                path_len_ballistic = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / cos_th;
                if path_len_ballistic > 0
                    p1_end = p1 + d1 * path_len_ballistic;
                    p2_end = p2 + d2 * path_len_ballistic;
                    p3_end = p3 + d3 * path_len_ballistic;
                    
                    cos_rx_tilt = abs(d1*Nx + d2*Ny + d3*Nz);
                    if cos_rx_tilt >= cos_FOV_half
                        vec_x = p1_end - Rx; vec_y = p2_end - Ry; vec_z = p3_end - Rz;
                        cx = vec_y*d3 - vec_z*d2; cy = vec_z*d1 - vec_x*d3; cz = vec_x*d2 - vec_y*d1;
                        r_wander_perp = sqrt(cx^2 + cy^2 + cz^2);
                        r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);
                        
                        % 微元波束展宽 W_ST -> 0，直接判定刚体碰撞边界
                        if r_wander_perp <= r_eff
                            P_packet = P_packet + exp(-param.coef_c * path_len_ballistic);
                        end
                    end
                end
            end
            
            % =================================================================
            % [多重散射路径]: PIS 演化，无湍流下纯吸收与散射主导
            % =================================================================
            for ord = 1:n_max
                % 游走一步
                d_step = -log(rand()) / param.coef_b; 
                p1 = p1 + d1 * d_step; p2 = p2 + d2 * d_step; p3 = p3 + d3 * d_step;
                weight = weight * exp(-param.coef_a * d_step);
                
                % 越界判断
                if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L || weight < 1e-15
                    break; 
                end 
                
                % 强制接收估计 (Virtual Ray Local Estimation)
                vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3;
                d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; dist2rx = sqrt(d2rx_sq); 
                dir2rx_1 = vec2rx_1/dist2rx; dir2rx_2 = vec2rx_2/dist2rx; dir2rx_3 = vec2rx_3/dist2rx;
                
                cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                
                % =========================================================
                % 【核心修复：FOV 空间截断剔除策略】
                % 如果当前散射点超出了接收机 FOV 向后投影的几何视场锥
                % 立即舍弃该光子，终止追踪，切断无效算力消耗
                if cos_inc < cos_FOV_half
                    break; 
                end
                % =========================================================
                
                % 如果在 FOV 内，则继续接收计算
                omega = Rx_Area / d2rx_sq * cos_inc;
                cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                % 累加到达接收面的高阶散射能量 (无湍流截断惩罚)
                P_packet = P_packet + weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * dist2rx);
                
                
                % HG-PIS 散射角更新
                xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); 
                cos_th_i = max(min((1 + g_prop^2 - term^2) / (2 * g_prop), 1), -1); 
                denom_core = 1 + g_prop^2 - 2 * g_prop * cos_th_i;
                q_HG = (1 - g_prop^2) / (4 * pi * denom_core * sqrt(denom_core));
                p_val = pdf_Empirical(cos_th_i, param);                   
                weight = weight * max(0.5,min(2, p_val / q_HG));
                %weight = weight *  p_val / q_HG;
                
                dir_new = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1 - cos_th_i^2), 2*pi*rand());
                d1 = dir_new(1); d2 = dir_new(2); d3 = dir_new(3);
            end
            P_rx_accum = P_rx_accum + P_packet;
        end
        
        PL_arr_temp(d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
        fprintf('  %s | L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', water_types{w_idx}, L, toc, -PL_arr_temp(d_idx));
    end
    PL_Cell_WCI_None{w_idx} = PL_arr_temp;
end

save('data_exp2_WCIMC_None.mat', 'dist_cell', 'PL_Cell_WCI_None', 'water_types');
fprintf('\n统一理论无湍流基准数据已成功保存至 data_exp2_WCIMC_None.mat！\n');

% === 辅助函数区域 (保持不变) ===
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000, t = max(th(i)*180/pi, 1e-6); val(i) = exp(p.q_e * (1 - p.k1*sqrt(t) + p.k2*t - p.k3*t^1.5 + p.k4*t^2 - p.k5*t^2.5)); end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6); 
    p = exp(param.q_e * (1 - param.k1*sqrt(t_deg) + param.k2*t_deg - param.k3*t_deg^1.5 + param.k4*t_deg^2 - param.k5*t_deg^2.5)) / param.b_emp_norm; 
end
function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2); sp = sin(psi_s); cp = cos(psi_s); 
    if denom < 1e-10
        nd = [st*cp, st*sp, sign(d(3))*ct]; 
    else
        nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, -st*cp*denom + d(3)*ct]; 
    end
end