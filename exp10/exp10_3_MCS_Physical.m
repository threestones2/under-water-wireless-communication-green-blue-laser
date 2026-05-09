%% exp10_3: MCS Physical (无湍流, 1e7 光子参考真值)
%  纯物理接收 MCS：光子必须几何命中接收孔径，无任何半解析近似
%  作为参考真值与 WCI-MC / SA-MCS 比较收敛精度

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= 物理参数 (以 MCS_Physical 为准) =================
dist_cell = {[15, 55], [15, 35], [10, 15]};
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366];
coef_b_arr = [0.0374, 0.219, 1.824];
coef_c_arr = [0.1514, 0.398, 2.190];
num_W = length(water_types);

N_packets = 1e8;
n_max = 200;

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

fprintf('--- exp10_3: MCS Physical (无湍流, 参考真值, N=%d) ---\n', N_packets);

for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx);
    param.coef_b = coef_b_arr(w_idx);
    param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c;
    param = calc_haltrin_params(param);

    % 构建 LUT (散射角抽样)
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
    PL_arr = zeros(1, num_D);
    Time_arr = zeros(1, num_D);

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

            % --- 弹道检查 (阶数 0，命中孔径才跳过散射) ---
            denom_bal = d1*Nx + d2*Ny + d3*Nz;
            t_bal = -1;
            if abs(denom_bal) > 1e-10
                t_bal = ((Rx-p1)*Nx + (Ry-p2)*Ny + (Rz-p3)*Nz) / denom_bal;
                if t_bal > 0
                    p_end_1 = p1 + d1 * t_bal;
                    p_end_2 = p2 + d2 * t_bal;
                    p_end_3 = p3 + d3 * t_bal;
                    cos_rx_tilt = abs(denom_bal);
                    if cos_rx_tilt >= cos_FOV_half
                        vec_x = p_end_1 - Rx; vec_y = p_end_2 - Ry; vec_z = p_end_3 - Rz;
                        cx = vec_y*d3 - vec_z*d2;
                        cy = vec_z*d1 - vec_x*d3;
                        cz = vec_x*d2 - vec_y*d1;
                        r_wander_perp = sqrt(cx^2 + cy^2 + cz^2);
                        r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);
                        if r_wander_perp <= r_eff
                            P_rx_accum = P_rx_accum + exp(-param.coef_c * t_bal);
                            continue;  % 弹道命中，跳过散射
                        end
                    end
                end
            end

            % --- 散射游走 (阶数 >= 1) ---
            for ord = 1:n_max
                d_step = -log(rand()) / param.coef_b;

                p1_prev = p1; p2_prev = p2; p3_prev = p3;

                p1 = p1 + d1 * d_step;
                p2 = p2 + d2 * d_step;
                p3 = p3 + d3 * d_step;

                if (p1-Tx)*Lx + (p2-Ty)*Ly + (p3-Tz)*Lz >= L
                    if ord > 1
                        denom = d1*Nx + d2*Ny + d3*Nz;
                        if abs(denom) > 1e-10
                            t_rx = ((Rx-p1_prev)*Nx + (Ry-p2_prev)*Ny + (Rz-p3_prev)*Nz) / denom;
                            if t_rx > 0 && t_rx <= d_step
                                px = p1_prev + d1 * t_rx;
                                py = p2_prev + d2 * t_rx;
                                pz = p3_prev + d3 * t_rx;
                                r2 = (px-Rx)^2 + (py-Ry)^2 + (pz-Rz)^2;
                                if r2 <= Rx_Aperture_half_sq && -(denom) >= cos_FOV_half
                                    P_packet = P_packet + weight * exp(-param.coef_c * t_rx);
                                end
                            end
                        end
                    end
                    break;
                end

                weight = weight * exp(-param.coef_a * d_step);

                if weight < 1e-8
                    break;
                end

                % O(1) LUT 散射角抽样
                idx_lut = floor(rand() * (LUT_SIZE - 1)) + 1;
                ct_s = cos_th_LUT(idx_lut);
                st_s = sin_th_LUT(idx_lut);
                dir_new = rotate_direction_fast([d1, d2, d3], ct_s, st_s, 2*pi*rand());
                d1 = dir_new(1); d2 = dir_new(2); d3 = dir_new(3);
            end
            P_rx_accum = P_rx_accum + P_packet;
        end

        PL_arr(d_idx) = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
        Time_arr(d_idx) = toc;
        fprintf('  %s | L=%2dm | N=%d | Time=%5.2fs | PL=%6.2f dB\n', ...
            water_types{w_idx}, L, N_packets, Time_arr(d_idx), -PL_arr(d_idx));
    end
    PL_Cell{w_idx} = PL_arr;
    Time_Cell{w_idx} = Time_arr;
end

save('data_exp10_MCS_Physical.mat', 'dist_cell', 'PL_Cell', 'Time_Cell', 'N_packets', 'water_types');
fprintf('已保存参考真值至 data_exp10_MCS_Physical.mat\n');

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
