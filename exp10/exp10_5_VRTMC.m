%% exp10_5: VRT-MC (With Turbulence) — 光子数扫描
%  有湍流 VRT-MC：弹道光子 + 散射虚拟光子 + 直射路径湍流解析惩罚
%  执行逻辑参考 WCIMC (内联几何)，HG-PIS 固定 g=0.97
%  湍流参数同步自 AST_MCS/src/MC_MPS_v6.m
%  计时从脚本起始开始 (含全部准备工作: Haltrin + OTOPS 谱计算)

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));
t_global = tic;  % 全局计时起点

% ================= 物理参数 =================
dist_cell = {[15, 55], [15, 35], [10, 15]};
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
% 对应 Jerlov 水质：Clear Ocean → IB, Coastal → II, Turbid Harbor → 7C
coef_a_arr = [0.0512, 0.0512, 0.2224];
coef_b_arr = [0.0557, 0.3789, 2.4868];
coef_c_arr = [0.1069, 0.4301, 2.7092];
num_W = length(water_types);

N_packets_arr = [1e3, 1e4, 1e5];
num_N = length(N_packets_arr);
n_max = 200;

lambda_nm = 514;
lambda = lambda_nm * 1e-9;
k_wave = 2 * pi / lambda;

w0 = 0.002;
div_angle = 0.1 * pi / 180;
theta_half_div = div_angle / 2;
Rx_Aperture = 0.01;
Rx_FOV = 20 * pi / 180;
Rx_Area = pi * (Rx_Aperture / 2)^2;

cos_FOV_half = cos(Rx_FOV / 2);

% VRT-MC: HG-PIS 固定 g=0.97
g_prop = 0.97;

% ================= 湍流参数 (同步自 AST_MCS/src/MC_MPS_v6.m) =================
T_avg = 20;
S_avg = 35;
H_ratio = -1;
epsilon = 1e-10;
chi_T = 1e-5;

fprintf('--- exp10_5: VRT-MC (有湍流) 光子数扫描 ---\n');
fprintf('Turbulence: T=%.0f C S=%.0f ppt H_ratio=%d epsilon=%.0e chi_T=%.0e\n', ...
    T_avg, S_avg, H_ratio, epsilon, chi_T);

% ================= 预计算 Haltrin 参数 (全部水质, 计入总时间) =================
param_cell = cell(1, num_W);
for w_idx = 1:num_W
    p = struct();
    p.coef_a = coef_a_arr(w_idx);
    p.coef_b = coef_b_arr(w_idx);
    p.coef_c = coef_c_arr(w_idx);
    p.albedo = p.coef_b / p.coef_c;
    param_cell{w_idx} = calc_haltrin_params(p);
end

% ================= 预计算 OTOPS 湍流相干参数 (全部距离, 计入总时间) =================
all_dists = unique([dist_cell{:}]);
[Phi_func, ~, eta_physical] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, H_ratio);

fprintf('预计算 OTOPS 湍流相干参数...\n');
fprintf('  内尺度 eta = %.4f mm\n', eta_physical * 1000);

rho0_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:length(all_dists)
    [rho0_map(all_dists(i)), ~] = calc_turb_coherence_params(Phi_func, k_wave, all_dists(i), eta_physical);
    fprintf('  L=%2dm: rho0=%.4fm\n', all_dists(i), rho0_map(all_dists(i)));
end
fprintf('Prep done at t = %.2f s\n\n', toc(t_global));

PL_Cell = cell(1, num_W);
Time_Cell = cell(1, num_W);

for w_idx = 1:num_W
    param = param_cell{w_idx};
    dist_arr = dist_cell{w_idx};
    num_D = length(dist_arr);
    PL_matrix = zeros(num_D, num_N);
    Time_matrix = zeros(num_D, num_N);

    for d_idx = 1:num_D
        L = dist_arr(d_idx);
        rho0_Link = rho0_map(L);

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
            P_rx_accum = 0;
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

                % --- 弹道光子 (纯几何截断 + 湍流解析惩罚) ---
                denom_bal = d1*Nx + d2*Ny + d3*Nz;
                if abs(denom_bal) > 1e-10
                    t_bal = ((Rx-p1)*Nx + (Ry-p2)*Ny + (Rz-p3)*Nz) / denom_bal;
                    if t_bal > 0
                        p_end_1 = p1 + d1 * t_bal;
                        p_end_2 = p2 + d2 * t_bal;
                        p_end_3 = p3 + d3 * t_bal;
                        cos_rx_tilt = abs(d1*Nx + d2*Ny + d3*Nz);
                        if cos_rx_tilt >= cos_FOV_half
                            vec_x = p_end_1 - Rx; vec_y = p_end_2 - Ry; vec_z = p_end_3 - Rz;
                            cx = vec_y*d3 - vec_z*d2;
                            cy = vec_z*d1 - vec_x*d3;
                            cz = vec_x*d2 - vec_y*d1;
                            r_wander_perp = sqrt(cx^2 + cy^2 + cz^2);
                            r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);

                            % 几何光束宽度
                            W_geo = sqrt(w0^2 + (t_bal * tan(theta_half_div))^2);

                            % 方案A: 湍流解析惩罚 (Likelihood Ratio)
                            rho0_ballistic = rho0_Link * (L / t_bal)^(3/5);
                            W_turb_LT = 2 * t_bal / (k_wave * rho0_ballistic);
                            cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*w0))^(1/3));
                            W_ST = W_turb_LT * cf;

                            W_eff = sqrt(W_geo^2 + W_ST^2);
                            exponent = -2 * (r_wander_perp^2) * (1/(W_eff^2) - 1/(W_geo^2));
                            weight_factor = (W_geo^2 / W_eff^2) * exp(exponent);

                            if r_wander_perp <= r_eff
                                P_packet = P_packet + exp(-param.coef_c * t_bal) * weight_factor;
                            end
                        end
                    end
                end

                % --- 散射光子 (直线步进, 虚拟光子局部估计, 散射分量无湍流) ---
                for ord = 1:n_max
                    d_step = -log(rand()) / param.coef_b;
                    p1 = p1 + d1 * d_step;
                    p2 = p2 + d2 * d_step;
                    p3 = p3 + d3 * d_step;

                    weight = weight * exp(-param.coef_a * d_step);

                    if (p1-Tx)*Lx + (p2-Ty)*Ly + (p3-Tz)*Lz >= L || weight < 1e-8
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

                    P_packet = P_packet + base_w;

                    % HG-PIS 散射角采样 (固定 g=0.97)
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
            Time_matrix(d_idx, n_idx) = toc(t_global);  % 累计时间, 含全部准备
            fprintf('  %s | L=%2dm | N=%7d | Time=%5.2fs | PL=%6.2f dB\n', ...
                water_types{w_idx}, L, N_packets, Time_matrix(d_idx, n_idx), -PL_matrix(d_idx, n_idx));
        end
    end
    PL_Cell{w_idx} = PL_matrix;
    Time_Cell{w_idx} = Time_matrix;
end

save('data_exp10_VRTMC.mat', 'dist_cell', 'PL_Cell', 'Time_Cell', 'N_packets_arr', 'water_types');
fprintf('已保存至 data_exp10_VRTMC.mat\n');

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

function [Phi_n_func, k_wave, eta_physical] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, H_ratio)
    lambda_nm = lambda * 1e9; k_wave = 2 * pi / lambda;
    A = -1.05e-6*S + 2*1.6e-8*T*S - 2*2.02e-6*T - 4.23e-3/lambda_nm;
    B = 1.779e-4 - 1.05e-6*T + 1.6e-8*T^2 + 1.155e-2/lambda_nm;
    s_f = S*1e-3; T_k = T+273.15;
    cp = 1000 * ((5.328 - 9.76e-2 * S + 4.04e-4 * S^2) + (-6.913e-3 + 7.351e-4 * S - 3.15e-6 * S^2) * T + (9.6e-6 - 1.927e-6 * S + 8.23e-9 * S^2) * T^2 + (2.5e-9 + 1.666e-9 * S - 7.125e-12 * S^2) * T^3);
    rho = (9.9992293295e2 + 2.0341179217e-2 * T - 6.1624591598e-3 * T^2 + 2.2614664708e-5 * T^3 - 4.6570659168e-8 * T^4) + s_f * (8.0200240891e2 - 2.0005183488 * T + 1.6771024982e-2 * T^2 - 3.0600536746e-5 * T^3 - 1.6132224742e-5 * T * S);
    mu = ((0.15700386464 * (T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5) * (1 + (1.5409136040 + 1.9981117208e-2 * T - 9.5203865864e-5 * T^2) * s_f + (7.9739318223 - 7.561456881e-2 * T + 4.7237011074e-4 * T^2) * s_f^2);
    sigma_T = 10^(log10(240 + 0.0002 * (S / 1.00472)) - 3 + 0.434 * (2.3 - (343.5 + 0.037 * (S / 1.00472)) / (1.00024 * T + 273.15)) * (1 - (1.00024 * T + 273.15) / (647.3 + 0.03 * (S / 1.00472)))^(1/3));
    Pr = mu * cp / sigma_T; Sc = mu^2 / (5.954e-15 * T_k * rho);
    c_T = 0.072^(4/3) / Pr; c_S = 0.072^(4/3) / Sc; c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    R_rho = 2.6e-4 * abs(H_ratio) / 7.6e-4;
    if R_rho >= 1, d_r = R_rho + sqrt(R_rho)*sqrt(R_rho - 1); elseif R_rho >= 0.5, d_r = 1.85*R_rho - 0.85; else, d_r = 0.15*R_rho; end
    chi_S = chi_T * d_r / (H_ratio^2); chi_TS = chi_T * (1 + d_r) / (2 * H_ratio);

    nu_kinematic = mu / rho; eta_physical = (nu_kinematic^3 / epsilon)^(0.25);
    Phi_Hill = @(K, chi_M, c_M) (0.72 / (4 * pi)) * chi_M * epsilon^(-1/3) .* (K.^2 + 1e-25).^(-11/6) .* exp(-176.90 * (K * eta_physical).^2 .* c_M^(0.96)) .* (1 + 21.61 * (K * eta_physical).^(0.61) .* c_M^(0.02) - 18.18 * (K * eta_physical).^(0.55) .* c_M^(0.04));
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + B^2 * Phi_Hill(K, chi_S, c_S) + 2 * A * B * Phi_Hill(K, chi_TS, c_TS));
end

function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L, eta_physical)
    kappa_eval = 0.01 / eta_physical;
    Cn2_eq = Phi_n_func(kappa_eval) / (0.033 * kappa_eval^(-11/3));
    try
        opts = optimset('Display', 'off');
        rho0_exact = fzero(@(rho) 8*pi^2*k_wave^2*L*integral2(@(K, xi) K.*Phi_n_func(K).*(1-besselj(0, K.*rho.*xi)), 1e-2, 1e5, 0, 1, 'Method', 'iterated', 'RelTol', 1e-3) - 2, (0.545*k_wave^2*Cn2_eq*L)^(-3/5), opts);
    catch
        rho0_exact = (0.545*k_wave^2*Cn2_eq*L)^(-3/5);
    end
end
