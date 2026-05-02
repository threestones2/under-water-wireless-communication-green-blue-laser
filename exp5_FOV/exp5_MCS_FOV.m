%% Exp8_2: MCS FOV Simulation (Traditional LUT Baseline)
% 实验目的: 分析无湍流时，接收视场角 (FOV) 对传统MCS纯散射能量收集能力的影响
clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= 物理与仿真参数初始化 =================
dist_cell = {10:10:60, 10:10:60, 5:5:20}; 
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr =[0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr =[0.1514, 0.398, 2.190];
num_W = length(water_types);

% 设定接收视场角遍历数组 (单位: 度)
fov_angles_deg = [5, 10, 40, 90];
num_F = length(fov_angles_deg);

N_packets = 1e5; n_max = 200; 

w0 = 0.002;                                
div_angle = 0.1 * pi / 180;                
theta_half_div = div_angle / 2; 
Rx_Aperture = 0.01;                        
Rx_Area = pi * (Rx_Aperture / 2)^2;
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

PL_Cell_MCS_FOV = cell(1, num_W);

fprintf('--- 运行 Exp8_2: MCS FOV分析 (O(1) LUT 极速标量版) ---\n');

% ================= 蒙特卡洛散射过程 =================
for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx); param.coef_b = coef_b_arr(w_idx); param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c; 
    param = calc_haltrin_params(param);
    
    th_axis = linspace(0, pi, 50000); pdf_vals = zeros(size(th_axis));
    for i = 1:length(th_axis)
        pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i)); 
    end
    cdf_vals = cumtrapz(th_axis, pdf_vals); cdf_vals = cdf_vals / cdf_vals(end);
    [cdf_uniq, idx_uniq] = unique(cdf_vals); LUT_SIZE = 100000;
    P_grid = linspace(0, 1, LUT_SIZE);
    th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap');
    cos_th_LUT = cos(th_LUT); sin_th_LUT = sin(th_LUT);
    
    dist_arr = dist_cell{w_idx}; num_D = length(dist_arr);
    PL_matrix = zeros(num_F, num_D);
    
    for f_idx = 1:num_F
        % 动态更新 FOV 相关判定参数
        Rx_FOV = fov_angles_deg(f_idx) * pi / 180;
        cos_FOV_half = cos(Rx_FOV / 2);
        
        % 固定无偏轴误差
        theta_err = 0;
        Tx_Dir = [sin(theta_err), cos(theta_err), 0];
        u_vec_Tx = [cos(theta_err), -sin(theta_err), 0]; 
        v_vec_Tx = [0, 0, 1];
        
        fprintf('\n[MCS] | %s | 接收 FOV: %.1f deg\n', water_types{w_idx}, fov_angles_deg(f_idx));
        
        for d_idx = 1:num_D
            L = dist_arr(d_idx);
            Tx_Pos =[0, 0, 0]; Rx_Pos = [0, L, 0]; Rx_Normal =[0, -1, 0];
            
            Lx = Tx_Dir(1); Ly = Tx_Dir(2); Lz = Tx_Dir(3);
            Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
            Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3); 
            Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
            
            P_rx_accum = 0; tic; 
            rng(123456, 'twister'); 
            
            for p = 1:N_packets
                r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand(); cp0 = cos(phi0); sp0 = sin(phi0);
                pos =[Tx + r0*cp0*u_vec_Tx(1) + r0*sp0*v_vec_Tx(1), ...
                       Ty + r0*cp0*u_vec_Tx(2) + r0*sp0*v_vec_Tx(2), ...
                       Tz + r0*cp0*u_vec_Tx(3) + r0*sp0*v_vec_Tx(3)];
                
                U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
                dir = rotate_direction_fast(Tx_Dir, cos(U_init), sin(U_init), 2*pi*rand());
                
                weight = 1.0; P_packet = 0;
                
                % --- 0阶直射路径分析 ---
                denom_rx = dir(1)*Nx + dir(2)*Ny + dir(3)*Nz;
                if abs(denom_rx) > 1e-10
                    d_plane = ((Rx - pos(1))*Nx + (Ry - pos(2))*Ny + (Rz - pos(3))*Nz) / denom_rx; 
                    if d_plane > 0
                        pos_end_1 = pos(1) + dir(1) * d_plane; 
                        pos_end_2 = pos(2) + dir(2) * d_plane; 
                        pos_end_3 = pos(3) + dir(3) * d_plane; 
                        
                        cos_rx_tilt = abs(dir(1)*Nx + dir(2)*Ny + dir(3)*Nz);
                        if cos_rx_tilt >= cos_FOV_half
                            if (pos_end_1-Rx)^2 + (pos_end_2-Ry)^2 + (pos_end_3-Rz)^2 <= Rx_Aperture_half_sq
                                P_packet = P_packet + exp(-param.coef_c * d_plane);
                            end
                        end
                    end
                end
                
                % --- 散射游走与能量收集 ---
                for ord = 1:n_max
                    d_step = -log(rand()) / param.coef_b; 
                    pos(1) = pos(1) + dir(1)*d_step; 
                    pos(2) = pos(2) + dir(2)*d_step; 
                    pos(3) = pos(3) + dir(3)*d_step; 
                    weight = weight * exp(-param.coef_a * d_step);
                    
                    if (pos(1)-Tx)*Lx + (pos(2)-Ty)*Ly + (pos(3)-Tz)*Lz >= L || weight < 1e-15
                        break; 
                    end 
                    
                    vec2rx_1 = Rx - pos(1); vec2rx_2 = Ry - pos(2); vec2rx_3 = Rz - pos(3); 
                    d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; d2rx = sqrt(d2rx_sq); 
                    dir2rx_1 = vec2rx_1/d2rx; dir2rx_2 = vec2rx_2/d2rx; dir2rx_3 = vec2rx_3/d2rx;
                    
                    cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);

                    if cos_inc >= cos_FOV_half
                        cos_scatter = dir(1)*dir2rx_1 + dir(2)*dir2rx_2 + dir(3)*dir2rx_3;
                        P_packet = P_packet + weight * min(1, pdf_Empirical(cos_scatter, param) * (Rx_Area/d2rx_sq * cos_inc)) * exp(-param.coef_c * d2rx);
                    end
                    
                    u_rand = rand(); idx_lut = floor(u_rand * (LUT_SIZE - 1)) + 1;
                    dir = rotate_direction_fast(dir, cos_th_LUT(idx_lut), sin_th_LUT(idx_lut), 2*pi*rand());
                end
                P_rx_accum = P_rx_accum + P_packet;
            end
            PL_matrix(f_idx, d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
            fprintf('    L=%2dm | Time: %5.2fs | Path Loss: %6.2f dB\n', L, toc, -PL_matrix(f_idx, d_idx));
        end
    end
    PL_Cell_MCS_FOV{w_idx} = PL_matrix;
end
save('data_exp8_MCS_FOV.mat', 'dist_cell', 'PL_Cell_MCS_FOV', 'water_types', 'fov_angles_deg');
fprintf('Exp8_2 MCS FOV 数据已保存！\n');

% ================= 辅助函数区域 =================
function param = calc_haltrin_params(p)
    b_e=p.coef_b; al_e=p.albedo; p.q_e=2.598+17.748*sqrt(b_e)-16.722*b_e+5.932*b_e*sqrt(b_e); p.k1=1.188-0.688*al_e; p.k2=0.1*(3.07-1.90*al_e); p.k3=0.01*(4.58-3.02*al_e); p.k4=0.001*(3.24-2.25*al_e); p.k5=0.0001*(0.84-0.61*al_e);
    th=linspace(0,pi,2000); val=zeros(size(th)); for i=1:2000, t=max(th(i)*180/pi,1e-6); val(i)=exp(p.q_e*(1-p.k1*t^0.5+p.k2*t^1.0-p.k3*t^1.5+p.k4*t^2.0-p.k5*t^2.5)); end; p.b_emp_norm=2*pi*trapz(th,val.*sin(th)); param=p;
end
function p = pdf_Empirical(cos_theta, param), t_deg=max(acos(max(min(cos_theta,1),-1))*180/pi,1e-6); p=exp(param.q_e*(1-param.k1*t_deg^0.5+param.k2*t_deg^1.0-param.k3*t_deg^1.5+param.k4*t_deg^2.0-param.k5*t_deg^2.5))/param.b_emp_norm; end
function nd = rotate_direction_fast(d, ct, st, psi_s), denom=sqrt(1-d(3)^2); sp=sin(psi_s); cp=cos(psi_s); if denom<1e-10, nd=[st*cp,st*sp,sign(d(3))*ct]; else, nd=[st/denom*(d(1)*d(3)*cp-d(2)*sp)+d(1)*ct, st/denom*(d(2)*d(3)*cp+d(1)*sp)+d(2)*ct, -st*cp*denom+d(3)*ct]; end; nd=nd/norm(nd); end