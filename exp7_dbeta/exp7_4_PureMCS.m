%% Exp7: Pure Physical MCS (Analog MC, No Local Estimation)
% 实验目的: 验证完全基于真实物理碰撞与硬孔径几何相交的纯粹蒙特卡洛模型
% 机制说明: 剥离 PIS，光子必须在三维空间中物理撞击接收机才会被记录

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= 物理与仿真参数初始化 =================
dist_cell = {10:10:60, 10:10:60, 5:5:20}; 
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr =[0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr =[0.1514, 0.398, 2.190];
num_W = length(water_types);

off_axis_angles = [0, 0.3, 1];
num_A = length(off_axis_angles);

N_packets = 1e7; n_max = 200; 

w0 = 0.002;                                
div_angle = 0.1 * pi / 180;                
theta_half_div = div_angle / 2; 
Rx_Aperture = 0.01;                        
Rx_FOV = 5 * pi / 180;                     

cos_FOV_half = cos(Rx_FOV / 2);
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

PL_Cell_PureMCS = cell(1, num_W);

fprintf('--- 运行 Exp7: Pure Analog MCS (纯物理无偏转硬孔径接收) ---\n');

for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx); param.coef_b = coef_b_arr(w_idx); param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c; 
    param = calc_haltrin_params(param);
    
    % 构建 Haltrin 相函数的高精度 LUT
    th_axis = logspace(-6, log10(pi), 50000); th_axis = [0, th_axis]; 
    pdf_vals = zeros(size(th_axis));
    for i = 1:length(th_axis), pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i)); end
    cdf_vals = cumtrapz(th_axis, pdf_vals); cdf_vals = cdf_vals / cdf_vals(end);
    
    [cdf_uniq, idx_uniq] = unique(cdf_vals);
    LUT_SIZE = 100000; P_grid = linspace(0, 1, LUT_SIZE);
    th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap');
    cos_th_LUT = cos(th_LUT); sin_th_LUT = sin(th_LUT);
    
    dist_arr = dist_cell{w_idx};
    num_D = length(dist_arr);
    PL_matrix = zeros(num_A, num_D);
    
    for a_idx = 1:num_A
        theta_err = off_axis_angles(a_idx) * pi / 180;
        % 发射端存在物理偏轴对准误差
        Nominal_Dir = [sin(theta_err), cos(theta_err), 0]; 
        
        for d_idx = 1:num_D
            L = dist_arr(d_idx);
            
            P_rx_accum = 0; tic;
            rng(123456, 'twister'); 
            
            for p = 1:N_packets
                % 初始发射参数
                r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
                pos = [r0*cos(phi0)*cos(theta_err), r0*cos(phi0)*sin(theta_err), r0*sin(phi0)];
                
                U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
                dir = rotate_direction_fast(Nominal_Dir, cos(U_init), sin(U_init), 2*pi*rand());
                
                weight = 1.0; 
                P_packet = 0;
                
                % 统一物理演化循环 (直射与多重散射)
                for ord = 1:n_max
                    d_step = -log(rand()) / param.coef_b; 
                    
                    % ---------------------------------------------------------
                    % 核心逻辑：物理硬孔径求交判定 (无任何强制偏转)
                    % ---------------------------------------------------------
                    if dir(2) > 0 % 光子朝向接收平面 (y = L) 运动
                        t_hit = (L - pos(2)) / dir(2);
                        if t_hit > 0 && t_hit <= d_step
                            % 光子的物理轨迹穿透了接收平面
                            x_hit = pos(1) + t_hit * dir(1);
                            z_hit = pos(3) + t_hit * dir(3);
                            r_hit_sq = x_hit^2 + z_hit^2;
                            
                            % 判定是否击中有限尺寸的接收机并满足 FOV
                            if r_hit_sq <= Rx_Aperture_half_sq && dir(2) >= cos_FOV_half
                                P_packet = weight * exp(-param.coef_a * t_hit);
                            end
                            % 无论是否击中靶心，光子均已越过接收面，终止追踪
                            break; 
                        end
                    end
                    
                    % ---------------------------------------------------------
                    % 若未穿透接收面，执行物理空间步进与真实散射
                    % ---------------------------------------------------------
                    pos = pos + dir * d_step; 
                    weight = weight * exp(-param.coef_a * d_step);
                    
                    % 边界清理：背向飞出过远或能量耗尽
                    if pos(2) < -10 || weight < 1e-15
                        break; 
                    end 
                    
                    % LUT O(1) 物理散射角抽样
                    u_rand = rand();
                    idx_lut = floor(u_rand * (LUT_SIZE - 1)) + 1;
                    ct_s = cos_th_LUT(idx_lut);
                    st_s = sin_th_LUT(idx_lut);
                    
                    dir = rotate_direction_fast(dir, ct_s, st_s, 2*pi*rand());
                end
                P_rx_accum = P_rx_accum + P_packet;
            end
            % 记录结果，若无任何光子命中，输出极小值对应的 300 dB
            PL_matrix(a_idx, d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
            fprintf('  %s | Off-axis: %.1f | L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', ...
                water_types{w_idx}, off_axis_angles(a_idx), L, toc, -PL_matrix(a_idx, d_idx));
        end
    end
    PL_Cell_PureMCS{w_idx} = PL_matrix;
end

save('data_exp7_MCS.mat', 'dist_cell', 'PL_Cell_MCS');
fprintf('Exp7 MCS 偏轴数据已保存！\n');

% === 辅助函数区域 ===
function param = calc_haltrin_params(p)
    b_e=p.coef_b; al_e=p.albedo; p.q_e=2.598+17.748*sqrt(b_e)-16.722*b_e+5.932*b_e*sqrt(b_e); 
    p.k1=1.188-0.688*al_e; p.k2=0.1*(3.07-1.90*al_e); p.k3=0.01*(4.58-3.02*al_e); 
    p.k4=0.001*(3.24-2.25*al_e); p.k5=0.0001*(0.84-0.61*al_e); param = p;
end

function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6); 
    sq_t = sqrt(t_deg); t_sq = t_deg * t_deg;
    term = 1 - param.k1*sq_t + param.k2*t_deg - param.k3*t_deg*sq_t + param.k4*t_sq - param.k5*t_sq*sq_t;
    p = exp(param.q_e * term); % 归一化在LUT生成时处理
end

function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2); sp = sin(psi_s); cp = cos(psi_s);
    if denom < 1e-10
        nd = [st*cp, st*sp, sign(d(3))*ct]; 
    else
        nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, ...
              st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, -st*cp*denom + d(3)*ct]; 
    end
    nd = nd / norm(nd);
end