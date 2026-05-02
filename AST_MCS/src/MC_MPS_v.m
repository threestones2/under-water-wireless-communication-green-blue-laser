% %% 水下/大气光通信混合仿真: MC-MPS (Turbulence Comparison)
% %  Comparison: With Turbulence vs. Without Turbulence
% %  Unified Petzold Clear Ocean parameters
% clc; clear; close all;
% 
% %% ================= 1. 参数初始化 =================
% param = struct();
% % --- 模式选择 ---
% param.phase_func = 'Empirical';  
% param.n_max = 10;          
% % --- 波长 ---
% lambda_nm = 514; 
% lambda = lambda_nm * 1e-9;
% k_wave = 2*pi/lambda;
% w0 = 0.01;                  % 束腰半径 (m)
% % div_angle = 1e-3;           % 发散角 (rad)，采用物理衍射极限theta_div_physics = lambda / (pi * w0) =1.63 *10^(-5)
% % --- 3D 空间布局 ---
% Tx_Pos = [0, 0, 0];         
% Rx_Pos = [0, 100, 0];       % 传输距离
% Link_Vec = Rx_Pos - Tx_Pos;
% Link_Dist = norm(Link_Vec);
% Link_Dir = Link_Vec / Link_Dist;
% mu_T = Link_Dir;            
% Rx_Normal = -Link_Dir;
% Rx_Aperture = 0.02;          % 孔径 2cm
% Rx_FOV = 3 * pi/180;        
% Rx_Area = pi*(Rx_Aperture/2)^2;
% 
% % --- 介质参数 (Petzold Clear Ocean) ---
% param.c_water = 2.25e8;     
% param.coef_c = 0.1514;      
% param.coef_a = 0.114;       
% param.coef_b = 0.0374; 
% param.albedo = param.coef_b / param.coef_c;
% 
% % --- 相函数参数 (Empirical - Petzold) ---
% param.c_e = param.coef_c; param.a_e = param.coef_a; b_e = param.coef_b; albedo_e = param.albedo;
% param.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
% param.k1 = 1.188 - 0.688*albedo_e; param.k2 = 0.1 * (3.07 - 1.90*albedo_e);
% param.k3 = 0.01 * (4.58 - 3.02*albedo_e); param.k4 = 0.001 * (3.24 - 2.25*albedo_e);
% param.k5 = 0.0001 * (0.84 - 0.61*albedo_e);
% 
% % 预计算归一化系数
% th_test = linspace(0, pi, 2000); val_test = zeros(size(th_test));
% for i=1:length(th_test)
%     t_deg = th_test(i) * 180 / pi; if t_deg < 1e-6, t_deg = 1e-6; end
%     term = 1 + (-1)^1*param.k1*t_deg^0.5 + (-1)^2*param.k2*t_deg^1.0 + ...
%            (-1)^3*param.k3*t_deg^1.5 + (-1)^4*param.k4*t_deg^2.0 + (-1)^5*param.k5*t_deg^2.5;
%     val_test(i) = exp(param.q_e * term);
% end
% param.b_emp_norm = 2 * pi * trapz(th_test, val_test .* sin(th_test));
% 
% % --- 湍流与相位屏 ---
% T_avg = 20; S_avg = 35; epsilon = 1e-6; chi_T = 1e-8; eta = 1e-3; H_ratio = -20;
% N_screens = 20;              
% D_screen = 2.0;              
% N_grid = 256;                
% delta_z_screen = Link_Dist / N_screens; 
% 
% % --- 时间轴 ---
% dt = 1e-10; t_min = dt; t_max= 1e-6;                        
% param.T_bins = t_min : dt : t_max;
% N_bins = length(param.T_bins);
% 
% % --- 仿真控制 ---
% N_packets = 1e5; % 可根据需要调整
% scenarios = {'With Turbulence', 'Without Turbulence'};
% results = struct(); % 存储结果
% 
% %% ================= 2. 对比仿真循环 =================
% dx = D_screen / N_grid;
% x_axis = (-N_grid/2 : N_grid/2-1) * dx;
% [Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
% 
% % 预先计算几何方向向量 (Screen Orientation)
% n_vec = Link_Dir;
% if abs(n_vec(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
% u_vec = cross(up_temp, n_vec); u_vec = u_vec / norm(u_vec);
% v_vec = cross(n_vec, u_vec);   v_vec = v_vec / norm(v_vec);
% 
% for s_idx = 1:2
%     scenario_name = scenarios{s_idx};
%     fprintf('=== Running Scenario: %s ===\n', scenario_name);
% 
%     % --- Step 2.1: 生成/重置相位屏 ---
%     Screen_Chain = repmat(struct('Center', [0,0,0], 'Normal', [0,0,1], ...
%                              'u_vec', [], 'v_vec', [], 'grad_x', [], 'grad_y', []), 1, N_screens);
% 
%     for i = 1:N_screens
%         Center_Pos = Tx_Pos + i * delta_z_screen * Link_Dir;
%         Screen_Chain(i).Center = Center_Pos;
%         Screen_Chain(i).Normal = Link_Dir;
%         Screen_Chain(i).u_vec = u_vec; 
%         Screen_Chain(i).v_vec = v_vec;
% 
%         if strcmp(scenario_name, 'With Turbulence')
%             % 生成湍流相位屏
%             phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
%             [gx, gy] = gradient(phi, dx);
%             Screen_Chain(i).grad_x = gx; 
%             Screen_Chain(i).grad_y = gy;
%         else
%             % 无湍流：梯度设为全零 (Flat Phase)
%             Screen_Chain(i).grad_x = zeros(N_grid, N_grid); 
%             Screen_Chain(i).grad_y = zeros(N_grid, N_grid);
%         end
%     end
% 
%     % --- Step 2.2: MC-MPS 传输 ---
%     h_time = 1e-14*ones(1, N_bins); 
% 
%     tic;
%     for p = 1:N_packets
%         % 光子初始化，后续改成光包，这部分初始化应该删去，但是其实影响不大，这部分初始化方向几乎都是轴向，用的准直高斯光束
%         r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
%         pos_local = r0*cos(phi0)*u_vec + r0*sin(phi0)*v_vec;
%         pos_init = Tx_Pos + pos_local;
%         theta_div_physics = lambda / (pi * w0); 
%         U = theta_div_physics * sqrt(-0.5 * log(rand()));
%         psi_ini = 2 * pi * rand;
%         dir_init = rotate_direction(mu_T, U, psi_ini);
%         weight_init = 1.0;
% 
%         % [Ballistic Branch - Modified for Beam Spreading]
%         pos = pos_init; dir = dir_init;
% 
%         % 1. 强制追踪到接收平面（使用一个巨大的虚拟孔径，确保只检测平面相交）
%         Huge_Aperture = 1e5; % 足够大，确保一定"命中"平面
%         [pos_end, dir_end, plane_hit, path_len_ballistic] = ray_march_generic(pos, dir, 1e9, ...
%             Rx_Pos, Huge_Aperture, pi, Rx_Normal, true, ...
%             Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
% 
%         % 2. 只有当光线确实向前传输并到达接收平面时才计算
%         if plane_hit
%             % 计算光束中心偏离接收机中心的距离 (Beam Wander)
%             r_wander = norm(pos_end - Rx_Pos);
% 
%             % 3. 计算光束在该处的有效束宽 (Beam Spreading)
%             % 包含：初始束腰 + 自由衍射 + 湍流扩展
%             if strcmp(scenario_name, 'With Turbulence')
%                  % 这里需要估算等效Cn2，根据OTOPS参数预估一个典型值
%                  % 或者根据相位屏梯度统计值动态计算（较复杂），此处用参数化方法
%                  % 假设弱湍流环境 (可根据实际情况调整)
%                  Cn2_eff = 1e-12; 
%                  W_spot = calc_beam_spot_size(w0, lambda, path_len_ballistic, Cn2_eff);
%             else
%                  Cn2_none = 0;
%                  W_spot = calc_beam_spot_size(w0, lambda, path_len_ballistic, Cn2_none);
%             end
% 
%             % 4. 计算接收能量 (高斯光束在圆孔径上的积分近似)
%             % 假设孔径半径 << 光斑半径 (在强湍流或远距离下通常成立)
%             % 接收功率 P_rx = P_total * (2 * A_rx / (pi * W^2)) * exp(-2 * r^2 / W^2)
% 
%             % 几何衰减因子 (Energy Dilution)
%             geometric_loss = (2 * Rx_Area) / (pi * W_spot^2); 
% 
%             % 偏轴衰减因子 (Off-axis Loss)
%             pointing_loss = exp(-2 * r_wander^2 / W_spot^2);
% 
%             % 水体吸收衰减
%             attenuation = exp(-param.coef_c * path_len_ballistic);
% 
%             % 5. 最终能量权重 (限制最大值为1，防止近场奇异点)
%             final_weight = weight_init * attenuation * min(1, geometric_loss * pointing_loss);
% 
%             t_arrival = path_len_ballistic / param.c_water;
%             bin_idx = floor((t_arrival - t_min) / dt) + 1;
%             if bin_idx >= 1 && bin_idx <= N_bins
%                 h_time(bin_idx) = h_time(bin_idx) + final_weight;
%             end
%         end
% 
%         % [Scattering Branch]
%         pos = pos_init; dir = dir_init; weight = weight_init;
%         current_dist_traveled = 0;
%         d_step = -log(rand()) / param.coef_c;
% 
%         [pos_new, dir_new, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
%             Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, ...
%             Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
% 
%         current_dist_traveled = current_dist_traveled + step_len;
%         if dot(pos_new - Tx_Pos, Link_Dir) >= Link_Dist, continue; end
%         pos = pos_new; dir = dir_new;
% 
%         for order = 1 : param.n_max
%             vec_to_rx = Rx_Pos - pos; dist_to_rx = norm(vec_to_rx); dir_to_rx = vec_to_rx / dist_to_rx;
% 
%            % PIS 估计 (混合修正版：虚拟射线质心追踪 + 波动解析展宽)
%             is_in_FOV_theta = acos(dot(-dir_to_rx, Rx_Normal));
%             if is_in_FOV_theta <= Rx_FOV/2
%                 cos_theta_s = dot(dir, dir_to_rx);
%                 p_phase = pdf_Empirical(cos_theta_s, param);
% 
%                 omega = Rx_Area / (dist_to_rx^2) * abs(dot(dir_to_rx, Rx_Normal));
%                 prob_survival = exp(-param.coef_c * dist_to_rx);
%                 base_energy_pis = weight * param.albedo * min(1, p_phase * omega) * prob_survival;
% 
%                 if base_energy_pis > 1e-15 % 剔除极小权重以加速
%                      if strcmp(scenario_name, 'With Turbulence')
%                         % 1. 追踪光束质心轨迹 (必须使用虚拟大孔径，确保计算交点)
%                         Huge_Aperture = 1e5;
%                         [pos_end, ~, v_hit_flag, v_len] = ray_march_generic(pos, dir_to_rx, dist_to_rx + 1e-3, ...
%                             Rx_Pos, Huge_Aperture, pi, Rx_Normal, true, ...
%                             Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
% 
%                         % 2. 质心命中虚拟平面，计算波前破碎导致的能量稀释与偏轴惩罚
%                         if v_hit_flag
%                             % 计算散射光线的质心在接收面上的漂移距离
%                             r_wander = norm(pos_end - Rx_Pos);
% 
%                             % 设定等效湍流强度 (根据 OTOPS 或典型值设定)
%                             Cn2_eq = 1e-12; 
% 
%                             % 计算从当前散射点到接收机距离上的球面波相干长度
%                             rho_0 = (0.545 * k_wave^2 * Cn2_eq * v_len)^(-3/5); 
% 
%                             % 计算长期湍流发散角
%                             theta_turb_LT = 2 / (k_wave * rho_0); 
% 
%                             % 转换为短期湍流发散角 (特征去倾斜尺度采用 Rx_Aperture)
%                             correction_factor = max(0, 1 - 0.37 * (rho_0 / Rx_Aperture)^(1/3));
%                             theta_turb_ST = theta_turb_LT * correction_factor;
% 
%                             % 水体固有相函数的等效发散角 (Petzold 前向散射区约 0.02 rad)
%                             theta_water = 0.02; 
% 
%                             % 中心能量密度稀释因子 (Dilution Factor)
%                             penalty_factor = (theta_water^2) / (theta_water^2 + theta_turb_ST^2);
% 
%                             % 光斑在接收面上的等效平方半径 (水体散射展宽 + 湍流短期展宽)
%                             W_total_sq = (v_len^2) * (theta_water^2 + theta_turb_ST^2);
% 
%                             % 偏轴指向损耗 (Off-axis Pointing Loss)
%                             pointing_loss = exp(-2 * r_wander^2 / W_total_sq);
% 
%                             % 施加双重惩罚计算最终接收能量
%                             final_energy_pis = base_energy_pis * penalty_factor * pointing_loss;
% 
%                             % 使用累积真实路径长度计算到达时间
%                             t_pis_arrival = (current_dist_traveled + v_len) / param.c_water;
%                             bin_idx = floor((t_pis_arrival - t_min) / dt) + 1;
%                             if bin_idx >= 1 && bin_idx <= N_bins
%                                 h_time(bin_idx) = h_time(bin_idx) + final_energy_pis;
%                             end
%                         end
%                     else
%                         % 无湍流情况：纯理想几何直线，无惩罚
%                         t_pis_arrival = (current_dist_traveled + dist_to_rx) / param.c_water;
%                         bin_idx = floor((t_pis_arrival - t_min) / dt) + 1;
%                         if bin_idx >= 1 && bin_idx <= N_bins
%                             h_time(bin_idx) = h_time(bin_idx) + base_energy_pis;
%                         end
%                     end
%                 end
%             end
% 
%             if order == param.n_max, break; end
%             weight = weight * param.albedo;
%             if weight < 1e-9, break; end
% 
%             % 散射偏转
%             theta_i = pi * rand; phi_scat = 2 * pi * rand; cos_t = cos(theta_i);
%             p_val = pdf_Empirical(cos_t, param);
%             weight_factor = 2 * pi^2 * p_val * sin(theta_i);
%             weight = weight * weight_factor;
%             dir = rotate_direction(dir, theta_i, phi_scat);
% 
%             d_step = -log(rand()) / param.coef_c;
%             [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
%                 Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, ...
%                 Screen_Chain, k_wave, x_axis, dx, N_grid, Tx_Pos, Link_Dir, delta_z_screen);
%             current_dist_traveled = current_dist_traveled + step_len;
%         end
%     end
%     t_run = toc;
%     fprintf('   Completed in %.2f s\n', t_run);
% 
%     % --- Step 2.3: 存储结果 ---
%     h_time_norm = h_time / N_packets;
%     P_rx_total = sum(h_time_norm);
% 
%     results(s_idx).h_time = h_time_norm;
%     results(s_idx).P_rx = P_rx_total;
%     results(s_idx).Path_Loss_dB = 10 * log10(P_rx_total);
%     if P_rx_total > 0
%         tau_mean = sum(param.T_bins .* h_time_norm) / P_rx_total;
%         tau_rms = sqrt( sum( ((param.T_bins - tau_mean).^2) .* h_time_norm ) / P_rx_total );
%         results(s_idx).tau_rms = tau_rms;
%     else
%         results(s_idx).tau_rms = 0;
%     end
% end

%% ================= 3. 对比绘图与分析 =================
fprintf('\n=== 仿真对比结果 (Tx-Rx: %.1fm) ===\n', Link_Dist);
fprintf('%-20s | Path Loss (dB) | RMS Delay (ns)\n', 'Scenario');
fprintf('---------------------------------------------------\n');
for i = 1:2
    fprintf('%-20s | %14.4f | %14.4f\n', scenarios{i}, results(i).Path_Loss_dB, results(i).tau_rms*1e9);
end

% --- Plot 1: CIR 对比 (Zoomed) ---
figure('Name', 'CIR Comparison', 'Color', 'w');
% plot(param.T_bins*1e9, 10*log10(results(1).h_time), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With Turbulence');
semilogy(param.T_bins*1e9, 10*log10(results(1).h_time), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With Turbulence');
hold on;
plot(param.T_bins*1e9, 10*log10(results(2).h_time), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Without Turbulence');
grid on;
legend('Location', 'best');
xlabel('Time (ns)');
ylabel('Received Power (dB)');
title(['CIR Comparison: Turbulence vs. No Turbulence (d=' num2str(Link_Dist) 'm)']);
xlim([400, 600]); % 聚焦在直射光到达附近
ylim([-160, max(10*log10([results.h_time]))+5]); 

% --- Plot 2: 归一化能量分布对比 (线性) ---
figure('Name', 'Linear Power Comparison', 'Color', 'w');
area(param.T_bins*1e9, results(2).h_time, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'No Turb');
hold on;
plot(param.T_bins*1e9, results(1).h_time, 'r-', 'LineWidth', 1.5, 'DisplayName', 'With Turb');
grid on;
legend;
xlabel('Time (ns)'); ylabel('Normalized Power (W)');
title('Linear Impulse Response Comparison');
xlim([430, 460]); % 高精度聚焦波峰

%% ================= 4. 辅助函数 =================

% --- [优化后] O(1) 复杂度 Ray Marching ---
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, ...
    Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, ...
    Screen_Chain, k_wave, x_axis, dx, N_grid, ...
    Tx_Pos, Link_Dir, delta_z_screen) % <--- 新增几何参数

    hit_flag = false;
    total_len = 0;
    remaining_dist = dist_limit;
    
    N_screens = length(Screen_Chain);
    
    while remaining_dist > 1e-6
        min_dist = remaining_dist;
        event_type = 'none';
        target_idx = -1;
        
        % --- [优化核心]: 直接几何计算目标相位屏，替代循环 ---
        % 计算光子在链路方向上的投影距离
        dist_projected = dot(pos - Tx_Pos, Link_Dir);
        % 计算当前所在的区间索引 (floor向下取整)
        current_layer_idx = floor(dist_projected / delta_z_screen);
        
        target_screen_idx = -1;
        cos_theta = dot(dir, Link_Dir);
        
        % 根据飞行方向决定检测哪一个屏
        if cos_theta > 1e-6       % 向前飞: 检测下一屏
            target_screen_idx = current_layer_idx + 1;
        elseif cos_theta < -1e-6  % 向后飞: 检测当前屏(即身后的屏)
            target_screen_idx = current_layer_idx;
        end
        
        % 仅当目标索引有效时进行相交检测
        if target_screen_idx >= 1 && target_screen_idx <= N_screens
            scr = Screen_Chain(target_screen_idx);
            denom = dot(dir, scr.Normal);
            if abs(denom) > 1e-6
                t = dot(scr.Center - pos, scr.Normal) / denom;
                if t > 1e-6 && t < min_dist
                    min_dist = t;
                    event_type = 'screen';
                    target_idx = target_screen_idx;
                end
            end
        end
        % --------------------------------------------------
        
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx) > 1e-6
                t_rx = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
                if t_rx > 1e-6 && t_rx <= min_dist
                    min_dist = t_rx;
                    event_type = 'rx';
                end
            end
        end
        
        pos = pos + dir * min_dist;
        remaining_dist = remaining_dist - min_dist;
        total_len = total_len + min_dist; 
        
        if strcmp(event_type, 'rx')
            dist_to_center = norm(pos - Rx_Pos);
            angle_inc = acos(dot(-dir, Rx_Normal));
            if dist_to_center <= Rx_Aperture/2 && angle_inc <= Rx_FOV/2
                hit_flag = true;
                return;
            end
            
        elseif strcmp(event_type, 'screen')
            scr = Screen_Chain(target_idx);
            vec_on_plane = pos - scr.Center;
            loc_u = dot(vec_on_plane, scr.u_vec);
            loc_v = dot(vec_on_plane, scr.v_vec);
            idx_x = mod(round((loc_u - x_axis(1))/dx), N_grid) + 1;
            idx_y = mod(round((loc_v - x_axis(1))/dx), N_grid) + 1;
            g_u = scr.grad_x(idx_y, idx_x); 
            g_v = scr.grad_y(idx_y, idx_x);
            delta_vec = (g_u * scr.u_vec + g_v * scr.v_vec) / k_wave;
            dir = dir + delta_vec; 
            dir = dir / norm(dir);
        end
    end
end

function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    % GET_OTOPS_SPECTRUM_HANDLE 计算 OTOPS 折射率谱的函数句柄
    % 基于 Yao et al., JOSA A, 2020.
    
    %% 1. 基础参数
    lambda_nm = lambda * 1e9; 
    k_wave = 2 * pi / lambda; 
    
    %% 2. 计算 OTOPS 系数 A 和 B
    % 注意：此处变量命名沿用 Wen et al. (2023) 的下标习惯
    % a5 对应 Yao 论文中的 a6 (盐度项), a6 对应 a7 (温度项)
    a1 = 1.779e-4;  a2 = -1.05e-6;  a3 = 1.6e-8; 
    a4 = -2.02e-6;  a5 = 1.155e-2;  a6 = -4.23e-3;
    
    A = a2*S + 2*a3*T*S + 2*a4*T + a6/lambda_nm; 
    B = a1 + a2*T + a3*T^2 + a5/lambda_nm; 
    
    %% 3. [修正] 流体热力学参数 (基于 Yao 2020 Appendix A)
    % 原始代码使用了固定值，这里改为动态计算以支持任意 T, S
    
    % 辅助变量
    T_k = T + 273.15; % 开尔文温度
    s_frac = S * 1e-3; % 盐度比例 (ppt -> ratio)
    
    % (1) 比热容 cp (Eq. A1-A2) [J/(kg K)]
    a11 = 5.328 - 9.76e-2*S + 4.04e-4*S^2;
    a12 = -6.913e-3 + 7.351e-4*S - 3.15e-6*S^2;
    a13 = 9.6e-6 - 1.927e-6*S + 8.23e-9*S^2;
    a14 = 2.5e-9 + 1.666e-9*S - 7.125e-12*S^2;
    cp = 1000 * (a11 + a12*T + a13*T^2 + a14*T^3); 
    
    % (2) 密度 rho (Eq. A8-A10) [kg/m^3]
    rho_T = 9.9992293295e2 + 2.0341179217e-2*T - 6.1624591598e-3*T^2 + ...
            2.2614664708e-5*T^3 - 4.6570659168e-8*T^4;
    rho_S = s_frac * (8.0200240891e2 - 2.0005183488*T + 1.6771024982e-2*T^2 - ...
            3.0600536746e-5*T^3 - 1.6132224742e-5*T*S);
    rho = rho_T + rho_S;
    
    % (3) 动力粘度 mu (Eq. A5-A7) [N s / m^2]
    mu_0 = (0.15700386464*(T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5;
    a21 = 1.5409136040 + 1.9981117208e-2*T - 9.5203865864e-5*T^2;
    a22 = 7.9739318223 - 7.561456881e-2*T + 4.7237011074e-4*T^2;
    mu = mu_0 * (1 + a21*s_frac + a22*s_frac^2);
    
    % (4) 热导率 sigma_T (Eq. A3-A4) [W/(m K)]
    T_b = 1.00024 * T;
    S_b = S / 1.00472;
    term1 = log10(240 + 0.0002*S_b);
    term2 = 0.434 * (2.3 - (343.5 + 0.037*S_b)/(T_b + 273.15));
    term3 = (1 - (T_b + 273.15)/(647.3 + 0.03*S_b))^(1/3);
    log_sigma = term1 - 3 + term2 * term3; % 修正公式结构
    % 注意：Yao论文Eq A3 写法比较复杂，这里简化处理或直接使用近似值
    % 考虑到公式实现的复杂性和排版歧义，如果对精度要求不极高，
    % 也可以保留 sigma_T = 0.6 或使用海水标准值，但 mu, rho, cp 建议用上面的公式
    sigma_T = 10^log_sigma; 
    % 备用方案：若上述公式报错，可回退到 sigma_T = 0.6;
    
    %% 4. 计算 Pr 和 Sc
    Pr = mu * cp / sigma_T;  
    Sc = mu^2 / (5.954e-15 * T_k * rho); 
    
    % Hill 参数 (Eq. 23)
    c_T = 0.072^(4/3) * Pr^(-1); 
    c_S = 0.072^(4/3) * Sc^(-1);
    c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    
    %% 5. Eddy Diffusivity Ratio (d_r)
    % 注：严格来说 alpha_c 和 beta_c 也应通过 TEOS-10 计算
    % 但此处保留您的近似常数以维持代码独立性
    alpha_c = 2.6e-4; beta_c = 7.6e-4;  
    R_rho = alpha_c * abs(H_ratio) / beta_c;
    
    if R_rho >= 1
        d_r = R_rho + sqrt(R_rho)*sqrt(R_rho-1);
    elseif R_rho >= 0.5
        d_r = 1.85*R_rho - 0.85;
    else
        d_r = 0.15*R_rho; 
    end
    
    chi_S = chi_T * d_r / (H_ratio^2);
    chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    
    %% 6. [修正] 定义谱函数 (Hill Spectrum)
    % Eq. 20: 系数应为 (0.72 / 4pi)
    coeff_Hill = 0.72 / (4 * pi); % ≈ 0.0573
    
    Phi_Hill = @(K, chi_M, c_M) coeff_Hill * chi_M * epsilon^(-1/3) .* ...
        (K.^2 + 1e-25).^(-11/6) .* ... 
        exp(-176.90 * (K*eta).^2 .* c_M^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_M^(0.02) - 18.18*(K*eta).^(0.55).*c_M^(0.04));
    
    % 返回总谱句柄
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + ...
                       B^2 * Phi_Hill(K, chi_S, c_S) + ...
                       2*A*B * Phi_Hill(K, chi_TS, c_TS));
end

function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    % GEN_SCREEN_FROM_SPECTRUM 根据给定的谱函数生成相位屏
    %
    % 输入:
    %   Phi_n_func - 折射率谱函数句柄 @(K)
    %   D          - 相位屏边长 (m)
    %   N          - 网格点数
    %   k_wave     - 光波数 (rad/m)
    %   delta_z    - 该相位屏代表的传输距离/厚度 (m)
    
    dk = 2 * pi / D;
    dx = D / N;
    
    %% 1. 高频部分 (FFT)
    kx = (-N/2 : N/2-1) * dk;
    [KX, KY] = meshgrid(kx, kx);
    K_grid = sqrt(KX.^2 + KY.^2);
    K_grid(N/2+1, N/2+1) = 1e-10; % 避开奇异点
    
    % 调用传入的谱函数计算 Phi_n
    Phi_n_val = Phi_n_func(K_grid);
    Phi_n_val(N/2+1, N/2+1) = 0; 
    
    % 转换为相位谱 F_phi = 2 * pi^2 * k^2 * dz * Phi_n
    F_phi = 2 * pi * k_wave^2 * delta_z * Phi_n_val;
    
    % 生成复高斯噪声
    noise = (randn(N) + 1i * randn(N));%不应该除以根号2，因为后面取了实部
    
    % 频域滤波
    C_nm = noise .* sqrt(F_phi) * dk;
    
    % 逆变换 (注意 ifft2 的缩放)
    phase_high = real(ifftshift(ifft2(ifftshift(C_nm))) )*N^2;
    
    %% 2. 次谐波补偿 (Subharmonics)
    phase_low = zeros(N, N);
    [xx, yy] = meshgrid((-N/2 : N/2-1) * dx);
    n_sub = 3; % 补偿级数
    
    for p = 1:n_sub
        dk_p = dk / (3^p); % 频率步长指数衰减
        for m = -1:1
            for n = -1:1
                if (m==0 && n==0), continue; end
                
                kx_p = m * dk_p;
                ky_p = n * dk_p;
                k_p = sqrt(kx_p^2 + ky_p^2);
                
                % 计算低频点的谱值
                Phi_n_p = Phi_n_func(k_p);
                F_phi_p = 2 * pi^2 * k_wave^2 * delta_z * Phi_n_p;
                
                amp = sqrt(F_phi_p) * dk_p;
                r_c = (randn(1) + 1i * randn(1));
                
                phase_low = phase_low + real( r_c * amp * exp(1i * (kx_p * xx + ky_p * yy)) );
            end
        end
    end
    
    %% 3. 合成
    phase_screen = phase_high + phase_low;
end

function new_dir = rotate_direction(dir, theta_s, psi_s)
    mu_x = dir(1); mu_y = dir(2); mu_z = dir(3);
    denom = sqrt(1 - mu_z^2);
    if denom < 1e-10
        if mu_z > 0
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), cos(theta_s)];
        else
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), -cos(theta_s)];
        end
    else
        sin_theta = sin(theta_s); cos_theta = cos(theta_s);
        cos_psi = cos(psi_s); sin_psi = sin(psi_s);
        new_dir_x = sin_theta/denom * (mu_x*mu_z*cos_psi - mu_y*sin_psi) + mu_x*cos_theta;
        new_dir_y = sin_theta/denom * (mu_y*mu_z*cos_psi + mu_x*sin_psi) + mu_y*cos_theta;
        new_dir_z = -sin_theta*cos_psi*denom + mu_z*cos_theta;
        new_dir = [new_dir_x, new_dir_y, new_dir_z];
    end
    new_dir = new_dir / norm(new_dir);
end

% --- PDF 函数 (1/sr) ---
function p = pdf_HG(cos_theta, g)
    p = (1 - g^2) ./ (4*pi * (1 + g^2 - 2*g*cos_theta).^1.5);
end

function p = pdf_Rayleigh(cos_theta, gamma)
    p = 3 * (1 + 3*gamma + (1-gamma)*cos_theta^2) / (16 * pi * (1 + 2*gamma));
end

function p = pdf_Mie(cos_theta, g, f)
    val = (1 - g^2)/2 * ( 1./(1 + g^2 - 2*g*cos_theta).^1.5 + f * 0.5 * (3*cos_theta^2 - 1) / (1 + g^2)^1.5 );
    p = val / (2*pi);
end

function p = pdf_TTHG(cos_theta, alpha, g1, g2)
    p1 = pdf_HG(cos_theta, g1);
    p2 = pdf_HG(cos_theta, g2);
    p = alpha * p1 + (1 - alpha) * p2;
end

function p = pdf_Empirical(cos_theta, param)
    cos_theta = max(min(cos_theta, 1), -1);
    theta_rad = acos(cos_theta);
    t_deg = theta_rad * 180 / pi;
    if t_deg < 1e-6, t_deg = 1e-6; end
    
    term = 1 + (-1)^1*param.k1*t_deg^0.5 + ...
           (-1)^2*param.k2*t_deg^1.0 + ...
           (-1)^3*param.k3*t_deg^1.5 + ...
           (-1)^4*param.k4*t_deg^2.0 + ...
           (-1)^5*param.k5*t_deg^2.5;
       
    VSF = exp(param.q_e * term);
    p = VSF / param.b_emp_norm;
end

function W_L = calc_beam_spot_size(w0, lambda, L, Cn2)
    % 计算高斯光束在距离 L 处的短期有效光斑半径 (Short-term spot size)
    % 基于 Yura (1973) Eq.(16), (20), (21)
    
    k = 2 * pi / lambda;
    D = 2 * w0; % 初始光束直径
    
    % 1. 自由空间衍射束宽 (Gaussian Beam Propagation)
    z_R = (pi * w0^2) / lambda;
    W_diff = w0 * sqrt(1 + (L / z_R)^2);
    
    % 2. 湍流引起的短期束宽扩展 (Short-term Beam Spread)
    if Cn2 > 1e-15
        % 长期的球面波相干长度 rho_0
        rho_0 = (0.545 * k^2 * Cn2 * L)^(-3/5);
        
        % 长期的湍流附加展宽 (Yura 1973, Eq. 21)
        W_turb_LT = 2 * L / (k * rho_0);
        
        % 短期的湍流附加展宽 (Yura 1973, Eq. 20)
        % 修正因子: 剔除了 Beam Wander 的贡献
        if rho_0 < D
            % 使用 Yura 近似公式，当相干长度小于初始孔径时短期修正尤为重要
            correction_factor = 1 - 0.37 * (rho_0 / D)^(1/3);
            % 防止极端情况下因子变为负数 (确保物理意义)
            correction_factor = max(correction_factor, 0); 
            W_turb_ST = W_turb_LT * correction_factor;
        else
            % 若相干长度大于孔径，湍流较弱，长短期差异极小
            W_turb_ST = W_turb_LT;
        end
        
        % 短期总有效束宽 (Quadratic Sum)
        W_L = sqrt(W_diff^2 + W_turb_ST^2);
    else
        W_L = W_diff;
    end
end
