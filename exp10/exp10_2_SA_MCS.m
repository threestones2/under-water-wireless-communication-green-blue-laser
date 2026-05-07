%% exp10_2: SA-MCS (Semi-Analytical MCS) — 光子数扫描
%  半解析 MCS：纯几何弹道 + 每散射节点局部估计贡献 + LUT 散射角抽样
%  无 ray_march 开销，纯解析几何，面向嵌入式高效率信道估计

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

lambda = 514e-9;
w0 = 0.002;
div_angle = 0.1 * pi / 180;
theta_half_div = div_angle / 2;
Rx_Aperture = 0.01;
Rx_FOV = 20 * pi / 180;
Rx_Area = pi * (Rx_Aperture / 2)^2;

cos_FOV_half = cos(Rx_FOV / 2);
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

PL_Cell = cell(1, num_W);
Time_Cell = cell(1, num_W);

fprintf('--- exp10_2: SA-MCS (半解析 Local Estimation) 光子数扫描 ---\n');

for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx);
    param.coef_b = coef_b_arr(w_idx);
    param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c;
    param = calc_haltrin_params(param);

    % 构建散射相函数 O(1) LUT
    th_axis = linspace(0, pi, 50000); pdf_vals = zeros(size(th_axis));
    for i = 1:length(th_axis)
        pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i));
    end
    cdf_vals = cumtrapz(th_axis, pdf_vals);
    cdf_vals = cdf_vals / cdf_vals(end);
    [cdf_uniq, idx_uniq] = unique(cdf_vals);
    LUT_SIZE = 100000;
    P_grid = linspace(0, 1, LUT_SIZE);
    th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap');
    cos_th_LUT = cos(th_LUT);
    sin_th_LUT = sin(th_LUT);

    dist_arr = dist_cell{w_idx};
    num_D = length(dist_arr);
    PL_matrix = zeros(num_D, num_N);
    Time_matrix = zeros(num_D, num_N);

    for d_idx = 1:num_D
        L = dist_arr(d_idx);
        rx_pos = [0, L, 0];
        link_dir = [0, 1, 0];
        rx_normal = [0, -1, 0];
        u_vec_Tx = [1, 0, 0];
        v_vec_Tx = [0, 0, 1];

        for n_idx = 1:num_N
            N_packets = N_packets_arr(n_idx);
            P_rx_accum = 0; tic;
            rng(123456, 'twister');

            for p = 1:N_packets
                % --- 光源采样 ---
                r0 = w0 * sqrt(-0.5*log(rand()));
                phi0 = 2*pi*rand();
                pos = [r0*cos(phi0), 0, r0*sin(phi0)];

                U_init = theta_half_div * sqrt(-0.5 * log(rand()));
                dir = rotate_direction_fast(link_dir, cos(U_init), sin(U_init), 2*pi*rand());

                weight = 1.0; P_packet = 0;

                % --- 弹道路径 (纯几何，无 ray_march) ---
                cos_th_bal = dir(2);  % dot(dir, link_dir), link_dir=[0,1,0]
                if cos_th_bal > 0
                    d_plane = (L - pos(2)) / cos_th_bal;
                    pos_end = pos + dir * d_plane;
                    r_hit_sq = pos_end(1)^2 + pos_end(3)^2;
                    cos_rx_tilt = abs(dir(2));  % abs(dot(dir, rx_normal)), rx_normal=[0,-1,0]

                    if cos_rx_tilt >= cos_FOV_half && r_hit_sq <= Rx_Aperture_half_sq
                        P_packet = P_packet + exp(-param.coef_c * d_plane);
                    end
                end

                % --- 散射游走 ---
                for ord = 1:n_max
                    d_step = -log(rand()) / param.coef_b;
                    pos = pos + dir * d_step;
                    weight = weight * exp(-param.coef_a * d_step);

                    if pos(2) >= L || weight < 1e-4
                        break;
                    end

                    % Local Estimation: 散射节点向接收器的半解析贡献
                    vec2rx = rx_pos - pos;
                    dist2rx = norm(vec2rx);
                    dir2rx = vec2rx / dist2rx;

                    cos_inc = dir2rx(2);  % dot(dir2rx, -rx_normal), -rx_normal=[0,1,0]
                    if cos_inc >= cos_FOV_half
                        cos_scatter = dot(dir, dir2rx);
                        solid_angle = Rx_Area / (dist2rx^2) * cos_inc;
                        base_w = weight * min(1.0, pdf_Empirical(cos_scatter, param) * solid_angle) * exp(-param.coef_c * dist2rx);
                        P_packet = P_packet + base_w;
                    end

                    % O(1) LUT 散射角抽样
                    idx_lut = floor(rand() * (LUT_SIZE - 1)) + 1;
                    ct_s = cos_th_LUT(idx_lut);
                    st_s = sin_th_LUT(idx_lut);
                    dir = rotate_direction_fast(dir, ct_s, st_s, 2*pi*rand());
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

save('data_exp10_SA_MCS.mat', 'dist_cell', 'PL_Cell', 'Time_Cell', 'N_packets_arr', 'water_types');
fprintf('已保存至 data_exp10_SA_MCS.mat\n');

% ================= 辅助函数 =================
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
