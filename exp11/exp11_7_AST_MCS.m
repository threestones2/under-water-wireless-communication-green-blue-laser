%% exp11_7: AST-MCS (MC_MPS_v6 算法) — Jerlov IB 距离扫描
%  AST-MCS: 纯几何步进 + 直射路径湍流解析惩罚 (方案A)
%  对比: 有湍流 vs 无湍流
%  湍流参数同步自 AST_MCS/src/MC_MPS_v6.m
%  光学/几何参数与 exp11 其他脚本一致以公平对比

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= Jerlov IB 水质参数 =================
water_type = 'Jerlov IB';
coef_a = 0.0512;
coef_b = 0.0557;
coef_c = 0.1069;

% ================= 距离扫描 =================
dist_arr = [10, 20, 30, 40, 50, 60, 70, 80];
num_D = length(dist_arr);

% ================= 光学 / 几何参数 (与 exp11 一致) =================
N_packets = 1e5;
n_max = 200;

lambda_nm = 514;
lambda = lambda_nm * 1e-9;
k_wave = 2 * pi / lambda;

w0 = 0.002;
div_angle = 0.1 * pi / 180;
theta_half_div = div_angle / 2;
Rx_Aperture = 0.05;
Rx_Radius = Rx_Aperture / 2;
Rx_FOV = 20 * pi / 180;
Rx_Area = pi * Rx_Radius^2;

c_water = 2.237e8;

% --- 发射端坐标系 (Link_Dir = [0,1,0]) ---
Tx_Pos_ref = [0, 0, 0];
Link_Dir_ref = [0, 1, 0];
Rx_Normal_ref = [0, -1, 0];
u_vec_Tx_ref = [1, 0, 0];
v_vec_Tx_ref = [0, 0, 1];

% ================= 湍流参数 (同步自 AST_MCS/src/MC_MPS_v6.m) =================
T_avg = 20;
S_avg = 35;
H_ratio = -1;
epsilon = 1e-10;
chi_T = 1e-5;

% ================= 计算 Haltrin 参数 =================
param = struct();
param.coef_a = coef_a;
param.coef_b = coef_b;
param.coef_c = coef_c;
param.albedo = coef_b / coef_c;
param = calc_haltrin_params(param);

% HG-PIS 提议分布: 固定 g=0.97 (MC_MPS_v6 方式)
param.g_prop = 0.97;

fprintf('=== exp11_7: AST-MCS (MC_MPS_v6 算法) Jerlov IB 距离扫描 ===\n');
fprintf('Water: %s | a=%.4f b=%.4f c=%.4f\n', water_type, coef_a, coef_b, coef_c);
fprintf('N=%d | w0=%.1fmm | FOV=%.0f deg | Aperture=%.0fmm\n', ...
    N_packets, w0*1e3, Rx_FOV*180/pi, Rx_Aperture*1e3);
fprintf('Turbulence params: T=%.0f C S=%.0f ppt H_ratio=%d epsilon=%.0e chi_T=%.0e\n', ...
    T_avg, S_avg, H_ratio, epsilon, chi_T);
fprintf('Algorithm: pure geometric stepping + analytical turbulence penalty (no phase screens)\n\n');

% ================= 预计算 OTOPS 湍流相干参数 (每个距离) =================
[Phi_func, ~, eta_physical] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, H_ratio);

fprintf('预计算 OTOPS 湍流相干参数...\n');
fprintf('  内尺度 eta = %.4f mm\n', eta_physical * 1000);

rho0_arr = zeros(1, num_D);
Cn2_arr = zeros(1, num_D);
for d_idx = 1:num_D
    [rho0_arr(d_idx), Cn2_arr(d_idx)] = calc_turb_coherence_params(Phi_func, k_wave, dist_arr(d_idx), eta_physical);
    fprintf('  L=%2dm: rho0=%.4fm, Cn2=%.2e m^(-2/3)\n', dist_arr(d_idx), rho0_arr(d_idx), Cn2_arr(d_idx));
end
fprintf('\n');

% ================= 主仿真循环 =================
scenarios = {'None', 'Turb'};
data_files = {'data_exp11_AST_MCS_None.mat', 'data_exp11_AST_MCS_Turb.mat'};

for s_idx = 1:2
    scenario_name = scenarios{s_idx};
    fprintf('=== Running Scenario: AST-MCS (%s Turbulence) ===\n', scenario_name);

    PL_arr = zeros(1, num_D);
    Time_arr = zeros(1, num_D);

    for d_idx = 1:num_D
        L = dist_arr(d_idx);
        rho0_Link = rho0_arr(d_idx);

        Tx_Pos = [0, 0, 0];
        Rx_Pos = [0, L, 0];
        Link_Dir = Link_Dir_ref;
        Rx_Normal = Rx_Normal_ref;
        u_vec_Tx = u_vec_Tx_ref;
        v_vec_Tx = v_vec_Tx_ref;

        Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3);
        Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3);
        Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3);
        Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3);
        Ux = u_vec_Tx(1); Uy = u_vec_Tx(2); Uz = u_vec_Tx(3);
        Vx = v_vec_Tx(1); Vy = v_vec_Tx(2); Vz = v_vec_Tx(3);

        P_rx_accum = 0;
        rng(123456, 'twister');
        tic;

        for p = 1:N_packets
            % --- 初始光子采样 (高斯光束) ---
            r0 = w0 * sqrt(-0.5 * log(rand()));
            phi0 = 2 * pi * rand();
            pos_local = r0 * cos(phi0) * u_vec_Tx + r0 * sin(phi0) * v_vec_Tx;
            pos_init = Tx_Pos + pos_local;

            U_init = theta_half_div * sqrt(-0.5 * log(rand()));
            dir_init = rotate_direction(Link_Dir, U_init, 2 * pi * rand());

            weight_init = 1.0;
            Huge_Aperture = 1e5;

            % ==== 0阶直射路径 (纯几何直线穿透) ====
            [pos_end, dir_end, plane_hit, path_len_ballistic] = ...
                ray_march_flat(pos_init, dir_init, 1e9, Rx_Pos, Huge_Aperture, pi, Rx_Normal, true);

            if plane_hit
                cos_rx_tilt = abs(dot(dir_end, Rx_Normal));
                if acos(cos_rx_tilt) <= Rx_FOV / 2
                    r_wander_perp = norm(cross(pos_end - Rx_Pos, dir_end));
                    r_eff = Rx_Radius * sqrt(cos_rx_tilt);

                    % 方案A: 几何截断判定 + PDF重要性采样权重补偿
                    W_geo = (w0.^2 + (path_len_ballistic * tan(theta_half_div)).^2).^0.5;

                    if strcmp(scenario_name, 'Turb')
                        % 1. 计算湍流引起的短期展宽 W_ST
                        rho0_ballistic = rho0_Link * (L / path_len_ballistic)^(3/5);
                        W_turb_LT = 2 * path_len_ballistic / (k_wave * rho0_ballistic);
                        cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*w0))^(1/3));
                        W_ST = W_turb_LT * cf;

                        % 2. 计算包含湍流的有效光斑大小 W_eff
                        W_eff = sqrt(W_geo^2 + W_ST^2);

                        % 3. 核心：计算重要性采样权重因子 (Likelihood Ratio)
                        exponent = -2 * (r_wander_perp^2) * (1/(W_eff^2) - 1/(W_geo^2));
                        weight_factor = (W_geo^2 / W_eff^2) * exp(exponent);
                    else
                        weight_factor = 1.0;
                    end

                    is_received = double(r_wander_perp <= r_eff);
                    P_rx_accum = P_rx_accum + weight_init * exp(-param.coef_c * path_len_ballistic) * is_received * weight_factor;
                end
            end

            % ==== 多重散射分支 PIS 追踪 (纯几何步进) ====
            pos = pos_init;
            dir = dir_init;
            current_dist = 0;
            weight_pis_likelihood = weight_init;

            % 初始物理游走
            d_step = -log(rand()) / param.coef_b;
            [pos, dir, ~, step_len] = ray_march_flat(pos, dir, d_step, Rx_Pos, Huge_Aperture, Rx_FOV, Rx_Normal, false);
            current_dist = current_dist + step_len;

            if dot(pos - Tx_Pos, Link_Dir) >= L
                continue;
            end

            for order = 1:n_max
                % ---- 虚拟强制接收 (Local Estimation) ----
                vec_rx = Rx_Pos - pos;
                dist_rx = norm(vec_rx);
                dir_rx = vec_rx / dist_rx;

                if acos(dot(-dir_rx, Rx_Normal)) <= Rx_FOV/2
                    p_phase = pdf_Empirical(dot(dir, dir_rx), param);
                    omega = Rx_Area / (dist_rx^2) * abs(dot(dir_rx, Rx_Normal));

                    current_physical_weight = weight_pis_likelihood * exp(-param.coef_a * current_dist);
                    base_energy = current_physical_weight * min(1, p_phase * omega) * exp(-param.coef_c * dist_rx);

                    if base_energy > 1e-15
                        [pos_v_end, ~, v_hit, v_len] = ray_march_flat(pos, dir_rx, dist_rx + 1e-1, Rx_Pos, Huge_Aperture, pi, Rx_Normal, true);

                        if v_hit
                            cos_theta_inc = abs(dot(dir_rx, Rx_Normal));
                            r_eff_v = Rx_Radius * sqrt(cos_theta_inc);
                            r_wander_perp_v = norm(cross(pos_v_end - Rx_Pos, dir_rx));
                            f_spread = double(r_wander_perp_v <= r_eff_v);
                            P_rx_accum = P_rx_accum + base_energy * f_spread;
                        end
                    end
                end

                % ---- 轮盘赌 ----
                current_physical_weight = weight_pis_likelihood * exp(-param.coef_a * current_dist);
                if current_physical_weight < 1e-9
                    if rand() > 0.1
                        break;
                    else
                        weight_pis_likelihood = weight_pis_likelihood / 0.1;
                    end
                end

                % ---- HG-PIS 生成提议散射角 (固定 g=0.97) ----
                xi_rand = rand();
                term = (1 - param.g_prop^2) / (1 - param.g_prop + 2 * param.g_prop * xi_rand);
                cos_t = max(min((1 + param.g_prop^2 - term^2) / (2 * param.g_prop), 1), -1);
                theta_i = acos(cos_t);

                q_HG = (1 - param.g_prop^2) / (4 * pi * (1 + param.g_prop^2 - 2*param.g_prop*cos_t)^1.5);
                p_val = pdf_Empirical(cos_t, param);

                weight_factor = min(1e3, p_val / q_HG);
                weight_pis_likelihood = weight_pis_likelihood * weight_factor;

                dir = rotate_direction(dir, theta_i, 2 * pi * rand());
                d_step = -log(rand()) / param.coef_b;

                [pos, dir, ~, step_len] = ray_march_flat(pos, dir, d_step, Rx_Pos, Huge_Aperture, Rx_FOV, Rx_Normal, false);
                current_dist = current_dist + step_len;

                if dot(pos - Tx_Pos, Link_Dir) >= L
                    break;
                end
            end
        end

        t_run = toc;
        PL_val = 10 * log10(max(P_rx_accum / N_packets, 1e-300));
        PL_arr(d_idx) = PL_val;
        Time_arr(d_idx) = t_run;
        fprintf('  L=%2dm | PL=%7.2f dB | Time=%6.1fs\n', L, -PL_val, t_run);
    end

    % 保存数据
    save(data_files{s_idx}, 'dist_arr', 'PL_arr', 'Time_arr', 'N_packets', 'water_type');
    fprintf('已保存至 %s\n\n', data_files{s_idx});
end

fprintf('=== exp11_7 完成 ===\n');

% ================= 核心辅助函数 =================

function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo;
    p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e);
    p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e);
    p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i = 1:length(th)
        t = max(th(i)*180/pi, 1e-6);
        term = 1 - p.k1*sqrt(t) + p.k2*t - p.k3*t.*sqrt(t) + p.k4*t.^2 - p.k5*t.^2.*sqrt(t);
        val(i) = exp(p.q_e * term);
    end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th));
    param = p;
end

function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(max(min(cos_theta, 1), -1)) * 180 / pi, 1e-6);
    term = 1 - param.k1*t_deg^0.5 + param.k2*t_deg^1.0 - param.k3*t_deg^1.5 + param.k4*t_deg^2.0 - param.k5*t_deg^2.5;
    p = exp(param.q_e * term) / param.b_emp_norm;
end

function [pos, dir, hit_flag, total_len] = ray_march_flat(pos, dir, dist_limit, Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check)
    hit_flag = false;
    total_len = 0;
    rem_dist = dist_limit;

    while rem_dist > 1e-9
        t_rx = inf;
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx) > 1e-10
                t_temp = dot(Rx_Pos - pos, Rx_Normal) / denom_rx;
                if t_temp > -1e-7
                    t_rx = max(0, t_temp);
                end
            end
        end

        [min_dist, event_idx] = min([rem_dist, t_rx]);

        pos = pos + dir * min_dist;
        rem_dist = rem_dist - min_dist;
        total_len = total_len + min_dist;

        if event_idx == 2
            if norm(pos - Rx_Pos) <= Rx_Aperture/2 && acos(dot(-dir, Rx_Normal)) <= Rx_FOV/2
                hit_flag = true;
            end
            break;
        elseif event_idx == 1
            break;
        end
    end
end

function nd = rotate_direction(d, t, p)
    denom = sqrt(1 - d(3)^2);
    if denom < 1e-10
        nd = [sin(t)*cos(p), sin(t)*sin(p), sign(d(3))*cos(t)];
    else
        nd = [sin(t)/denom*(d(1)*d(3)*cos(p)-d(2)*sin(p))+d(1)*cos(t), ...
              sin(t)/denom*(d(2)*d(3)*cos(p)+d(1)*sin(p))+d(2)*cos(t), ...
              -sin(t)*cos(p)*denom+d(3)*cos(t)];
    end
    nd = nd / norm(nd);
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
