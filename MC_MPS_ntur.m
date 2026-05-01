%% 2. 水下光通信仿真: MC-MPS (无湍流纯散射版本)
clc; clear; close all;

%% ================= 参数初始化 =================
param = struct();
param.n_max = 10;          

lambda = 514 * 1e-9;
w0 = 0.1;                  
div_angle = 10*pi/180;     

Tx_Pos = [0, 0, 0];         
Rx_Pos = [0, 100, 0];       
Link_Vec = Rx_Pos - Tx_Pos;
Link_Dist = norm(Link_Vec);
Link_Dir = Link_Vec / Link_Dist;

mu_T = Link_Dir;            
Rx_Normal = -Link_Dir;

if abs(mu_T(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
u_vec_Tx = cross(up_temp, mu_T); u_vec_Tx = u_vec_Tx / norm(u_vec_Tx);
v_vec_Tx = cross(mu_T, u_vec_Tx);   v_vec_Tx = v_vec_Tx / norm(v_vec_Tx);

Rx_Aperture = 0.2;          
Rx_FOV = 20 * pi/180;        
Rx_Area = pi*(Rx_Aperture/2)^2;

param.c_water = 2.25e8;     
param.coef_c = 0.1514;      
param.coef_a = 0.114;       
param.coef_b = 0.0374; 
param.albedo = param.coef_b / param.coef_c;

param.q_e = 2.598 + 17.748*sqrt(param.coef_b) - 16.722*param.coef_b + 5.932*param.coef_b*sqrt(param.coef_b);
param.k1 = 1.188 - 0.688*param.albedo; param.k2 = 0.1 * (3.07 - 1.90*param.albedo);
param.k3 = 0.01 * (4.58 - 3.02*param.albedo); param.k4 = 0.001 * (3.24 - 2.25*param.albedo);
param.k5 = 0.0001 * (0.84 - 0.61*param.albedo);

th_test = linspace(0, pi, 2000); val_test = zeros(size(th_test));
for i=1:length(th_test)
    t_deg = max(th_test(i) * 180 / pi, 1e-6);
    term = 1 + (-1)^1*param.k1*t_deg^0.5 + (-1)^2*param.k2*t_deg^1.0 + ...
           (-1)^3*param.k3*t_deg^1.5 + (-1)^4*param.k4*t_deg^2.0 + (-1)^5*param.k5*t_deg^2.5;
    val_test(i) = exp(param.q_e * term);
end
param.b_emp_norm = 2 * pi * trapz(th_test, val_test .* sin(th_test));

dt = 1e-10; t_min = dt; t_max= 1e-6; param.T_bins = t_min : dt : t_max;
N_bins = length(param.T_bins);
N_packets = 1e5;

%% ================= 仿真核心 =================
h_time = 1e-12*ones(1, N_bins); 
tic;
for p = 1:N_packets
    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
    pos_init = Tx_Pos + r0*cos(phi0)*u_vec_Tx + r0*sin(phi0)*v_vec_Tx;
    dir_init = rotate_direction(mu_T, div_angle * sqrt(-0.5 * log(rand())), 2 * pi * rand);
    
    % [直射光]
    [pos_end, dir_end, plane_hit, path_len] = ray_march_simple(pos_init, dir_init, 1e9, Rx_Pos, 1e5, pi, Rx_Normal, true);
    if plane_hit
        W_spot = w0 * sqrt(1 + (path_len / ((pi * w0^2) / lambda))^2);
        cos_tilt = abs(dot(dir_end, Rx_Normal));
        geom_loss = (2 * Rx_Area * cos_tilt) / (pi * W_spot^2); 
        if acos(cos_tilt) > Rx_FOV / 2, geom_loss = 0; end
        
        weight = exp(-param.coef_c * path_len) * min(1, geom_loss * exp(-2 * norm(pos_end - Rx_Pos)^2 / W_spot^2));
        idx = floor((path_len/param.c_water - t_min) / dt) + 1;
        if idx >= 1 && idx <= N_bins, h_time(idx) = h_time(idx) + weight; end
    end

    % [散射光]
    pos = pos_init; dir = dir_init; weight = 1.0; dist_travel = 0;
    d_step = -log(rand()) / param.coef_c;
    [pos, dir, ~, s_len] = ray_march_simple(pos, dir, d_step, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false);
    dist_travel = dist_travel + s_len;
    
    for order = 1 : param.n_max
        if dot(pos - Tx_Pos, Link_Dir) >= Link_Dist, break; end
        v_rx = Rx_Pos - pos; d_rx = norm(v_rx); d_dir = v_rx / d_rx;
        
        if acos(dot(-d_dir, Rx_Normal)) <= Rx_FOV/2
            base_eng = weight * param.albedo * min(1, pdf_Empirical(dot(dir, d_dir), param) * (Rx_Area/d_rx^2*abs(dot(d_dir, Rx_Normal)))) * exp(-param.coef_c * d_rx);
            if base_eng > 1e-15
                idx = floor(((dist_travel + d_rx)/param.c_water - t_min) / dt) + 1;
                if idx >= 1 && idx <= N_bins, h_time(idx) = h_time(idx) + base_eng; end
            end
        end
        
        weight = weight * param.albedo; if weight < 1e-9, break; end
        th_i = pi * rand; p_scat = 2 * pi * rand;
        weight = weight * (2 * pi^2 * pdf_Empirical(cos(th_i), param) * sin(th_i));
        dir = rotate_direction(dir, th_i, p_scat);
        
        [pos, dir, ~, s_len] = ray_march_simple(pos, dir, -log(rand())/param.coef_c, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false);
        dist_travel = dist_travel + s_len;
    end
end
fprintf('Completed in %.2f s\n', toc);

% 绘图
figure('Color', 'w'); plot(param.T_bins*1e9, 10*log10(h_time/N_packets), 'b--'); grid on;
xlabel('Time (ns)'); ylabel('Received Power (dB)'); title('CIR - Without Turbulence');

%% ================= 辅助函数 =================
function [pos, dir, hit_flag, total_len] = ray_march_simple(pos, dir, dist_limit, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check)
    hit_flag = false; total_len = dist_limit;
    if enable_hit_check
        denom_rx = dot(dir, Rx_Normal);
        if abs(denom_rx) > 1e-6
            t_rx = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
            if t_rx > 1e-6 && t_rx <= dist_limit
                total_len = t_rx; hit_flag = true;
                if norm(pos + dir*t_rx - Rx_Pos) > Rx_Aperture/2 || acos(dot(-dir, Rx_Normal)) > Rx_FOV/2
                    hit_flag = false;
                end
            end
        end
    end
    pos = pos + dir * total_len;
end
function new_dir = rotate_direction(dir, theta, psi)
    mz = dir(3); denom = sqrt(1 - mz^2);
    if denom < 1e-10, new_dir = [sin(theta)*cos(psi), sin(theta)*sin(psi), sign(mz)*cos(theta)];
    else, new_dir = [sin(theta)/denom*(dir(1)*mz*cos(psi) - dir(2)*sin(psi))+dir(1)*cos(theta), ...
                     sin(theta)/denom*(dir(2)*mz*cos(psi) + dir(1)*sin(psi))+dir(2)*cos(theta), -sin(theta)*cos(psi)*denom+mz*cos(theta)];
    end
    new_dir = new_dir / norm(new_dir);
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(max(min(cos_theta, 1), -1)) * 180 / pi, 1e-6);
    term = 1 + (-1)^1*param.k1*t_deg^0.5 + (-1)^2*param.k2*t_deg^1.0 + ...
           (-1)^3*param.k3*t_deg^1.5 + (-1)^4*param.k4*t_deg^2.0 + (-1)^5*param.k5*t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end