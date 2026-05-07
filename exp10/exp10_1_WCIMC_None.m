%% exp10_1: WCI-MC (No Turbulence) — 光子数扫描
%  无湍流 WCI-MC：弹道光子 + 散射虚拟光子，无 Marcum Q 统计模型
%  对比不同光子数下的路径损耗收敛行为

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= 物理参数 (以 MCS_Physical 为准) =================
dist_cell = {[15, 55], [15, 35], [10, 15]};
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366];
coef_b_arr = [0.0374, 0.219, 1.824];
coef_c_arr = [0.1514, 0.398, 2.190];
num_W = length(water_types);

N_packets_arr = [1e3, 1e4, 1e5, 1e6, 1e7];
num_N = length(N_packets_arr);
n_max = 200;

lambda = 514e-9; k_wave = 2 * pi / lambda;
w0 = 0.002;
div_angle = 0.1 * pi / 180;
theta_half_div = div_angle / 2;
Rx_Aperture = 0.01;
Rx_FOV = 20 * pi / 180;
Rx_Area = pi * (Rx_Aperture / 2)^2;

cos_FOV_half = cos(Rx_FOV / 2);
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

% ================= 无湍流：哑元相位屏 =================
Grad_X_3D = zeros(1, 1, 1); Grad_Y_3D = zeros(1, 1, 1);
Screen_Z_1D = [1e10]; dx_s = 1; N_grid_s = 1;
x_axis = [0]; delta_z_screen = 1e10;

PL_Cell = cell(1, num_W);   % {water_type}(distance, N_idx)
Time_Cell = cell(1, num_W);

fprintf('--- exp10_1: WCI-MC (无湍流) 光子数扫描 ---\n');

for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx);
    param.coef_b = coef_b_arr(w_idx);
    param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c;
    param = calc_haltrin_params(param);

    P_max = pdf_Empirical(1.0, param);
    C_val = 4 * pi * P_max;
    g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);

    dist_arr = dist_cell{w_idx};
    num_D = length(dist_arr);
    PL_matrix = zeros(num_D, num_N);
    Time_matrix = zeros(num_D, num_N);

    for d_idx = 1:num_D
        L = dist_arr(d_idx);
        Tx_Pos = [0, 0, 0]; Rx_Pos = [0, L, 0];
        Link_Dir = [0, 1, 0]; Rx_Normal = [0, -1, 0];
        u_vec_Tx = [1, 0, 0]; v_vec_Tx = [0, 0, 1];

        Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3);
        Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
        Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3);
        Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
        Ux = u_vec_Tx(1); Uy = u_vec_Tx(2); Uz = u_vec_Tx(3);
        Vx = v_vec_Tx(1); Vy = v_vec_Tx(2); Vz = v_vec_Tx(3);

        for n_idx = 1:num_N
            N_packets = N_packets_arr(n_idx);
            P_rx_accum = 0; tic;
            rng(123456, 'twister');

            for p = 1:N_packets
                r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
                cp0 = cos(phi0); sp0 = sin(phi0);

                p1 = Tx + r0*cp0*Ux + r0*sp0*Vx;
                p2 = Ty + r0*cp0*Uy + r0*sp0*Vy;
                p3 = Tz + r0*cp0*Uz + r0*sp0*Vz;

                U_init = theta_half_div * sqrt(-0.5 * log(rand()));
                dir = rotate_direction_fast(Link_Dir, cos(U_init), sin(U_init), 2*pi*rand());
                d1 = dir(1); d2 = dir(2); d3 = dir(3);

                weight = 1.0; P_packet = 0;

                % --- 弹道光子 (无湍流：纯几何截断) ---
                [p1_end, p2_end, p3_end, d1_end, d2_end, d3_end, plane_hit, path_len_ballistic] = ...
                    ray_march_flat_scalar(p1, p2, p3, d1, d2, d3, 1e9, ...
                    Rx, Ry, Rz, 1e9, -1, Nx, Ny, Nz, true, ...
                    Grad_X_3D, Grad_Y_3D, Screen_Z_1D, ...
                    Ux, Uy, Uz, Vx, Vy, Vz, k_wave, ...
                    x_axis(1), dx_s, N_grid_s, ...
                    Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen);

                if plane_hit
                    cos_rx_tilt = abs(d1_end*Nx + d2_end*Ny + d3_end*Nz);
                    if cos_rx_tilt >= cos_FOV_half
                        vec_x = p1_end - Rx; vec_y = p2_end - Ry; vec_z = p3_end - Rz;
                        cx = vec_y*d3_end - vec_z*d2_end;
                        cy = vec_z*d1_end - vec_x*d3_end;
                        cz = vec_x*d2_end - vec_y*d1_end;
                        r_wander_perp = sqrt(cx^2 + cy^2 + cz^2);
                        r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);
                        recv_frac = double(r_wander_perp <= r_eff);
                        P_packet = P_packet + exp(-param.coef_c * path_len_ballistic) * recv_frac;
                    end
                end

                % --- 散射光子 ---
                for ord = 1:n_max
                    d_step = -log(rand()) / param.coef_b;
                    [p1, p2, p3, d1, d2, d3, ~, step_len] = ...
                        ray_march_flat_scalar(p1, p2, p3, d1, d2, d3, d_step, ...
                        Rx, Ry, Rz, Rx_Aperture_half_sq, cos_FOV_half, ...
                        Nx, Ny, Nz, false, ...
                        Grad_X_3D, Grad_Y_3D, Screen_Z_1D, ...
                        Ux, Uy, Uz, Vx, Vy, Vz, k_wave, ...
                        x_axis(1), dx_s, N_grid_s, ...
                        Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen);

                    weight = weight * exp(-param.coef_a * step_len);

                    if (p1-Tx)*Lx + (p2-Ty)*Ly + (p3-Tz)*Lz >= L || weight < 1e-4
                        break;
                    end

                    vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3;
                    d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2;
                    dist2rx = sqrt(d2rx_sq);
                    dir2rx_1 = vec2rx_1 / dist2rx;
                    dir2rx_2 = vec2rx_2 / dist2rx;
                    dir2rx_3 = vec2rx_3 / dist2rx;

                    cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                    if cos_inc < cos_FOV_half, break; end

                    omega = Rx_Area / d2rx_sq * cos_inc;
                    cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                    base_w = weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * dist2rx);

                    if base_w > 1e-16
                        [p1_v, p2_v, p3_v, d1_v, d2_v, d3_v, v_hit, ~] = ...
                            ray_march_flat_scalar(p1, p2, p3, dir2rx_1, dir2rx_2, dir2rx_3, dist2rx+1e-1, ...
                            Rx, Ry, Rz, 1e10, -1, Nx, Ny, Nz, true, ...
                            Grad_X_3D, Grad_Y_3D, Screen_Z_1D, ...
                            Ux, Uy, Uz, Vx, Vy, Vz, k_wave, ...
                            x_axis(1), dx_s, N_grid_s, ...
                            Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen);

                        if v_hit
                            cos_inc_v = abs(d1_v*Nx + d2_v*Ny + d3_v*Nz);
                            r_eff_v = (Rx_Aperture/2) * sqrt(cos_inc_v);
                            vec_x = p1_v - Rx; vec_y = p2_v - Ry; vec_z = p3_v - Rz;
                            cx = vec_y*d3_v - vec_z*d2_v;
                            cy = vec_z*d1_v - vec_x*d3_v;
                            cz = vec_x*d2_v - vec_y*d1_v;
                            r_wp = sqrt(cx^2 + cy^2 + cz^2);
                            point_loss = double(r_wp <= r_eff_v);
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

                    dir_new = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1-cos_th_i^2), 2*pi*rand());
                    d1 = dir_new(1); d2 = dir_new(2); d3 = dir_new(3);
                end
                P_rx_accum = P_rx_accum + P_packet;
            end

            PL_matrix(d_idx, n_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
            Time_matrix(d_idx, n_idx) = toc;
            fprintf('  %s | L=%2dm | N=%7d | Time=%5.2fs | PL=%6.2f dB\n', ...
                water_types{w_idx}, L, N_packets, Time_matrix(d_idx, n_idx), -PL_matrix(d_idx, n_idx));
        end
    end
    PL_Cell{w_idx} = PL_matrix;
    Time_Cell{w_idx} = Time_matrix;
end

save('data_exp10_WCIMC_None.mat', 'dist_cell', 'PL_Cell', 'Time_Cell', 'N_packets_arr', 'water_types');
fprintf('已保存至 data_exp10_WCIMC_None.mat\n');

% ================= 辅助函数 =================
function [p1, p2, p3, d1, d2, d3, hit_flag, total_len] = ray_march_flat_scalar(...
        p1, p2, p3, d1, d2, d3, dist_limit, Rx, Ry, Rz, ...
        Rx_Aperture_half_sq, cos_FOV_half, Nx, Ny, Nz, enable_hit_check, ...
        Grad_X_3D, Grad_Y_3D, Screen_Z_1D, ...
        Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x0, dx, N_grid, ...
        Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen)
    hit_flag = false; total_len = 0; rem_dist = dist_limit;
    N_screens = length(Screen_Z_1D);
    while rem_dist > 1e-9
        dir_z = d1*Lx + d2*Ly + d3*Lz;
        z_pos = (p1-Tx)*Lx + (p2-Ty)*Ly + (p3-Tz)*Lz;
        t_scr = inf; target_idx = -1;
        if dir_z > 1e-10
            target_idx = floor((z_pos + 1e-9) / delta_z_screen) + 1;
            if target_idx < 1, target_idx = 1; end
            if target_idx <= N_screens
                t_scr = (Screen_Z_1D(target_idx) - z_pos) / dir_z;
            end
        elseif dir_z < -1e-10
            target_idx = ceil((z_pos - 1e-9) / delta_z_screen) - 1;
            if target_idx > N_screens, target_idx = N_screens; end
            if target_idx >= 1
                t_scr = (Screen_Z_1D(target_idx) - z_pos) / dir_z;
            end
        end
        t_rx = inf;
        if enable_hit_check
            denom_rx = d1*Nx + d2*Ny + d3*Nz;
            if abs(denom_rx) > 1e-10
                t_temp = ((Rx-p1)*Nx + (Ry-p2)*Ny + (Rz-p3)*Nz) / denom_rx;
                if t_temp > -1e-7, t_rx = max(0, t_temp); end
            end
        end
        [min_dist, event_idx] = min([rem_dist, t_rx, t_scr]);
        p1 = p1 + d1 * min_dist;
        p2 = p2 + d2 * min_dist;
        p3 = p3 + d3 * min_dist;
        rem_dist = rem_dist - min_dist;
        total_len = total_len + min_dist;
        if event_idx == 2
            if (p1-Rx)^2 + (p2-Ry)^2 + (p3-Rz)^2 <= Rx_Aperture_half_sq && ...
               -(d1*Nx + d2*Ny + d3*Nz) >= cos_FOV_half
                hit_flag = true;
            end
            break;
        elseif event_idx == 3
            loc_u = (p1-Tx)*Ux + (p2-Ty)*Uy + (p3-Tz)*Uz;
            loc_v = (p1-Tx)*Vx + (p2-Ty)*Vy + (p3-Tz)*Vz;
            idx_x = mod(round((loc_u - x0)/dx), N_grid) + 1;
            idx_y = mod(round((loc_v - x0)/dx), N_grid) + 1;
            gx = Grad_X_3D(idx_y, idx_x, target_idx);
            gy = Grad_Y_3D(idx_y, idx_x, target_idx);
            d1 = d1 + (gx*Ux + gy*Vx) / k_wave;
            d2 = d2 + (gx*Uy + gy*Vy) / k_wave;
            d3 = d3 + (gx*Uz + gy*Vz) / k_wave;
            d_norm = sqrt(d1^2 + d2^2 + d3^2);
            d1 = d1/d_norm; d2 = d2/d_norm; d3 = d3/d_norm;
        elseif event_idx == 1
            break;
        end
    end
end

function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo;
    p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e);
    p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e);
    p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i = 1:2000
        t = max(th(i)*180/pi, 1e-6); sq_t = sqrt(t); t_sq = t*t;
        term = 1 - p.k1*sq_t + p.k2*t - p.k3*t*sq_t + p.k4*t_sq - p.k5*t_sq*sq_t;
        val(i) = exp(p.q_e * term);
    end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end

function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6);
    sq_t = sqrt(t_deg); t_sq = t_deg * t_deg;
    term = 1 - param.k1*sq_t + param.k2*t_deg - param.k3*t_deg*sq_t ...
           + param.k4*t_sq - param.k5*t_sq*sq_t;
    p = exp(param.q_e * term) / param.b_emp_norm;
end

function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2); sp = sin(psi_s); cp = cos(psi_s);
    if denom < 1e-10
        nd = [st*cp, st*sp, sign(d(3))*ct];
    else
        nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, ...
              st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, ...
              -st*cp*denom + d(3)*ct];
    end
end
