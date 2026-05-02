%% Exp2_PDF: Turbulence Intensity PDF Simulation (20m)
clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= 用户选择开关 =================
WATER_FLAG = 1; % 1: 'Clear Ocean'
TURB_FLAG = 1;  % 1: 'Paper_Fig2a'

% ================= 参数初始化 =================
L = 50;                              % 固定传输距离 (m)
N_packets = 1000;                    
N_realizations = 200;                

water_types_all = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_all = [0.114, 0.179, 0.366]; 
coef_b_all = [0.0374, 0.219, 1.824]; 
coef_c_all = [0.1514, 0.398, 2.190];

water_types = {water_types_all{WATER_FLAG}};
coef_a_arr = coef_a_all(WATER_FLAG); 
coef_b_arr = coef_b_all(WATER_FLAG); 
coef_c_arr = coef_c_all(WATER_FLAG);
num_W = 1; n_max = 200; 

lambda = 514e-9; k_wave = 2 * pi / lambda;
w0 = 0.02;                          
theta_half_div = 0.1 * pi / 180 / 2; 
Rx_Aperture = 0.005;                  
Rx_FOV = 5 * pi / 180;             
Rx_Area = pi * (Rx_Aperture / 2)^2;
cos_FOV_half = cos(Rx_FOV / 2);      
Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

turb_params_all(1).T_avg = 20; turb_params_all(1).S_avg = 35; turb_params_all(1).H_ratio = -0.3;     
turb_params_all(1).epsilon = 1e-4; turb_params_all(1).chi_T = 1e-7; turb_params_all(1).name = 'Paper_Fig2a';

turb_params_all(2).T_avg = 20; turb_params_all(2).S_avg = 35; turb_params_all(2).H_ratio = -20; 
turb_params_all(2).epsilon = 1e-9; turb_params_all(2).chi_T = 1e-7; turb_params_all(2).name = 'Weak';

turb_params = turb_params_all(TURB_FLAG);
num_T = 1;

N_screens = 10; 
D_screen = 1; N_grid = 2^8; dx = D_screen / N_grid; x_axis = (-N_grid/2 : N_grid/2-1) * dx;

Intensity_Samples = zeros(num_T, num_W, N_realizations);

fprintf('--- 运行光强 PDF 统计仿真---\n');
for w = 1:num_W
    param.coef_a = coef_a_arr(w); param.coef_b = coef_b_arr(w); param.coef_c = coef_c_arr(w);
    param.albedo = param.coef_b / param.coef_c;
    param = calc_haltrin_params(param);
    P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max;
    g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);
    Tx = 0; Ty = 0; Tz = 0; Rx = 0; Ry = L; Rz = 0;
    Lx = 0; Ly = 1; Lz = 0; Nx = 0; Ny = -1; Nz = 0;
    Ux = 1; Uy = 0; Uz = 0; Vx = 0; Vy = 0; Vz = 1;
    delta_z_screen = L / N_screens;
    Screen_Z_1D = (1:N_screens) * delta_z_screen;
    
    for t_idx = 1:num_T
        tp = turb_params(t_idx);
        [Phi_func, ~, eta_physical] = get_OTOPS_spectrum_handle(lambda, tp.T_avg, tp.S_avg, tp.epsilon, tp.chi_T, tp.H_ratio);
        [rho0_Link, ~] = calc_turb_coherence_params(Phi_func, k_wave, L, eta_physical);
        tic;
        for real_idx = 1:N_realizations
            rng(real_idx * 1000 + WATER_FLAG * 10 + TURB_FLAG + L, 'twister'); 
            Grad_X_3D = zeros(N_grid, N_grid, N_screens);
            Grad_Y_3D = zeros(N_grid, N_grid, N_screens);
            for scr = 1:N_screens
                phi_screen = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
                [gx, gy] = gradient(phi_screen, dx);
                Grad_X_3D(:, :, scr) = gx; Grad_Y_3D(:, :, scr) = gy;
            end
            P_rx_accum = 0; 
            for p = 1:N_packets
                r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand(); cp0 = cos(phi0); sp0 = sin(phi0);
                p1 = Tx + r0*cp0*Ux + r0*sp0*Vx; p2 = Ty + r0*cp0*Uy + r0*sp0*Vy; p3 = Tz + r0*cp0*Uz + r0*sp0*Vz;
                U_init = theta_half_div * sqrt(-0.5 * log(rand()));
                dir = rotate_direction_fast([Lx,Ly,Lz], cos(U_init), sin(U_init), 2*pi*rand());
                d1 = dir(1); d2 = dir(2); d3 = dir(3); weight = 1.0; P_packet = 0;
                [p1_end, p2_end, p3_end, d1_end, d2_end, d3_end, plane_hit, path_len_ballistic] = ...
                    ray_march_flat_scalar_v3(p1, p2, p3, d1, d2, d3, 1e9, Rx, Ry, Rz, 1e10, -1, Nx, Ny, Nz, true, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x_axis(1), dx, N_grid, Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen, N_screens);
                if plane_hit
                    cos_rx_tilt = abs(d1_end*Nx + d2_end*Ny + d3_end*Nz);
                    if cos_rx_tilt >= cos_FOV_half
                        vec_x = p1_end - Rx; vec_y = p2_end - Ry; vec_z = p3_end - Rz;
                        cx = vec_y*d3_end - vec_z*d2_end; cy = vec_z*d1_end - vec_x*d3_end; cz = vec_x*d2_end - vec_y*d1_end;
                        r_wander_perp = sqrt(cx^2 + cy^2 + cz^2); r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);
                        rho0_ballistic = rho0_Link * (L / path_len_ballistic)^(3/5);
                        W_turb_LT = 2 * path_len_ballistic / (k_wave * rho0_ballistic);
                        cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*r_eff))^(1/3));
                        W_spot = W_turb_LT * cf;
                        if W_spot > 1e-6, recv_frac = 1 - marcumq(2*r_wander_perp/W_spot, 2*r_eff/W_spot); else, recv_frac = double(r_wander_perp <= r_eff); end
                        P_packet = P_packet + exp(-param.coef_c * path_len_ballistic) * recv_frac;
                    end
                end
                for ord = 1:n_max
                    d_step = -log(rand()) / param.coef_b; 
                    [p1, p2, p3, d1, d2, d3, ~, step_len] = ray_march_flat_scalar_v3(p1, p2, p3, d1, d2, d3, d_step, Rx, Ry, Rz, Rx_Aperture_half_sq, cos_FOV_half, Nx, Ny, Nz, false, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x_axis(1), dx, N_grid, Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen, N_screens);
                    weight = weight * exp(-param.coef_a * step_len);
                    if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L || weight < 1e-15, break; end 
                    vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; dist2rx = sqrt(d2rx_sq); 
                    dir2rx_1 = vec2rx_1 / dist2rx; dir2rx_2 = vec2rx_2 / dist2rx; dir2rx_3 = vec2rx_3 / dist2rx; cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                    if cos_inc < cos_FOV_half, break; end
                    omega = Rx_Area / d2rx_sq * cos_inc; cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                    base_w = weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * dist2rx);
                    if base_w > 1e-15
                        [p1_v, p2_v, p3_v, d1_v, d2_v, d3_v, v_hit, v_len] = ray_march_flat_scalar_v3(p1, p2, p3, dir2rx_1, dir2rx_2, dir2rx_3, dist2rx+1e-1, Rx, Ry, Rz, 1e10, -1, Nx, Ny, Nz, true, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x_axis(1), dx, N_grid, Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen, N_screens);
                        if v_hit
                            cos_inc_v = abs(d1_v*Nx + d2_v*Ny + d3_v*Nz); r_eff_v = (Rx_Aperture/2) * sqrt(cos_inc_v);
                            vec_x = p1_v - Rx; vec_y = p2_v - Ry; vec_z = p3_v - Rz; r_wp = sqrt((vec_y*d3_v - vec_z*d2_v)^2 + (vec_z*d1_v - vec_x*d3_v)^2 + (vec_x*d2_v - vec_y*d1_v)^2);
                            rho_0_sph = rho0_Link * (L / v_len)^(3/5); W_ST = (2 * v_len / (k_wave * rho_0_sph)) * max(0, 1 - 0.37*(rho_0_sph/(2*r_eff_v))^(1/3));
                            if W_ST > 1e-6, point_loss = 1 - marcumq(2*r_wp/W_ST, 2*r_eff_v/W_ST); else, point_loss = double(r_wp <= r_eff_v); end
                            P_packet = P_packet + base_w * point_loss;
                        end
                    end
                    xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); cos_th_i = (1 + g_prop^2 - term^2) / (2 * g_prop);
                    if cos_th_i > 1, cos_th_i = 1; elseif cos_th_i < -1, cos_th_i = -1; end
                    denom_core = 1 + g_prop^2 - 2 * g_prop * cos_th_i; q_HG = (1 - g_prop^2) / (4 * pi * denom_core * sqrt(denom_core));
                    weight = weight * max(0.5,min(2, pdf_Empirical(cos_th_i, param) / q_HG));
                    dir_new = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1 - cos_th_i^2), 2*pi*rand()); d1 = dir_new(1); d2 = dir_new(2); d3 = dir_new(3);
                end
                P_rx_accum = P_rx_accum + P_packet;
            end
            Intensity_Samples(1, 1, real_idx) = P_rx_accum / N_packets;
            if mod(real_idx, 50) == 0, fprintf('    完成 %d/%d 次采样...\n', real_idx, N_realizations); end
        end
        fprintf('  -> 本项耗时: %.2f s\n', toc);
    end
end

% save('data_exp2_PDF_20m.mat', 'Intensity_Samples', 'turb_params', 'water_types', 'L');

save('data_exp2_PDF_50m.mat', 'Intensity_Samples', 'turb_params', 'water_types', 'L');

%% ================= 基础物理计算函数支持 =================
function [p1, p2, p3, d1, d2, d3, hit_flag, total_len] = ray_march_flat_scalar_v3(p1, p2, p3, d1, d2, d3, dist_limit, Rx, Ry, Rz, Rx_Aperture_half_sq, cos_FOV_half, Nx, Ny, Nz, enable_hit_check, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x0, dx, N_grid, Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen, N_screens)
    hit_flag = false; total_len = 0; rem_dist = dist_limit;
    while rem_dist > 1e-9
        dir_z = d1*Lx + d2*Ly + d3*Lz; z_pos = (p1-Tx)*Lx + (p2-Ty)*Ly + (p3-Tz)*Lz; t_scr = inf; target_idx = -1;
        if dir_z > 1e-10
            target_idx = floor((z_pos + 1e-9) / delta_z_screen) + 1; if target_idx < 1, target_idx = 1; end
            if target_idx <= N_screens, t_scr = (Screen_Z_1D(target_idx) - z_pos) / dir_z; end
        elseif dir_z < -1e-10
            target_idx = ceil((z_pos - 1e-9) / delta_z_screen) - 1; if target_idx > N_screens, target_idx = N_screens; end
            if target_idx >= 1, t_scr = (Screen_Z_1D(target_idx) - z_pos) / dir_z; end
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
        p1 = p1 + d1 * min_dist; p2 = p2 + d2 * min_dist; p3 = p3 + d3 * min_dist;
        rem_dist = rem_dist - min_dist; total_len = total_len + min_dist;
        if event_idx == 2 
            if (p1-Rx)^2 + (p2-Ry)^2 + (p3-Rz)^2 <= Rx_Aperture_half_sq && -(d1*Nx + d2*Ny + d3*Nz) >= cos_FOV_half, hit_flag = true; end
            break; 
        elseif event_idx == 3 
            loc_u = (p1-Tx)*Ux + (p2-Ty)*Uy + (p3-Tz)*Uz; loc_v = (p1-Tx)*Vx + (p2-Ty)*Vy + (p3-Tz)*Vz;
            idx_x = mod(round((loc_u - x0)/dx), N_grid) + 1; idx_y = mod(round((loc_v - x0)/dx), N_grid) + 1;
            gx = Grad_X_3D(idx_y, idx_x, target_idx); gy = Grad_Y_3D(idx_y, idx_x, target_idx);
            d1 = d1 + (gx*Ux + gy*Vx) / k_wave; d2 = d2 + (gx*Uy + gy*Vy) / k_wave; d3 = d3 + (gx*Uz + gy*Vz) / k_wave;
            d_norm = sqrt(d1^2 + d2^2 + d3^2); d1 = d1/d_norm; d2 = d2/d_norm; d3 = d3/d_norm;
        elseif event_idx == 1
            break;
        end
    end
end
function param = calc_haltrin_params(p)
    b_e = p.coef_b; al_e = p.albedo; p.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
    p.k1 = 1.188 - 0.688*al_e; p.k2 = 0.1*(3.07 - 1.90*al_e); p.k3 = 0.01*(4.58 - 3.02*al_e); p.k4 = 0.001*(3.24 - 2.25*al_e); p.k5 = 0.0001*(0.84 - 0.61*al_e);
    th = linspace(0, pi, 2000); val = zeros(size(th));
    for i=1:2000, t = max(th(i)*180/pi, 1e-6); sq_t = sqrt(t); t_sq = t*t; term = 1 - p.k1*sq_t + p.k2*t - p.k3*t*sq_t + p.k4*t_sq - p.k5*t_sq*sq_t; val(i) = exp(p.q_e * term); end
    p.b_emp_norm = 2*pi*trapz(th, val.*sin(th)); param = p;
end
function p = pdf_Empirical(cos_theta, param)
    t_deg = max(acos(cos_theta) * 180 / pi, 1e-6); sq_t = sqrt(t_deg); t_sq = t_deg * t_deg;
    term = 1 - param.k1*sq_t + param.k2*t_deg - param.k3*t_deg*sq_t + param.k4*t_sq - param.k5*t_sq*sq_t; p = exp(param.q_e * term) / param.b_emp_norm;
end
function nd = rotate_direction_fast(d, ct, st, psi_s)
    denom = sqrt(1 - d(3)^2); sp = sin(psi_s); cp = cos(psi_s);
    if denom < 1e-10, nd = [st*cp, st*sp, sign(d(3))*ct]; else, nd = [st/denom*(d(1)*d(3)*cp - d(2)*sp) + d(1)*ct, st/denom*(d(2)*d(3)*cp + d(1)*sp) + d(2)*ct, -st*cp*denom + d(3)*ct]; end
end
function [Phi_n_func, k_wave, eta_physical] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, H_ratio)
    lambda_nm = lambda * 1e9; k_wave = 2 * pi / lambda; 
    A = -1.05e-6*S + 2*1.6e-8*T*S - 2*2.02e-6*T - 4.23e-3/lambda_nm; B = 1.779e-4 - 1.05e-6*T + 1.6e-8*T^2 + 1.155e-2/lambda_nm; s_f = S*1e-3; T_k = T+273.15;
    cp = 1000 * ((5.328 - 9.76e-2 * S + 4.04e-4 * S^2) + (-6.913e-3 + 7.351e-4 * S - 3.15e-6 * S^2) * T + (9.6e-6 - 1.927e-6 * S + 8.23e-9 * S^2) * T^2 + (2.5e-9 + 1.666e-9 * S - 7.125e-12 * S^2) * T^3); 
    rho = (9.9992293295e2 + 2.0341179217e-2 * T - 6.1624591598e-3 * T^2 + 2.2614664708e-5 * T^3 - 4.6570659168e-8 * T^4) + s_f * (8.0200240891e2 - 2.0005183488 * T + 1.6771024982e-2 * T^2 - 3.0600536746e-5 * T^3 - 1.6132224742e-5 * T * S);
    mu = ((0.15700386464 * (T + 64.992620050)^2 - 91.296496657)^(-1) + 4.2844324477e-5) * (1 + (1.5409136040 + 1.9981117208e-2 * T - 9.5203865864e-5 * T^2) * s_f + (7.9739318223 - 7.561456881e-2 * T + 4.7237011074e-4 * T^2) * s_f^2);
    sigma_T = 10^(log10(240 + 0.0002 * (S / 1.00472)) - 3 + 0.434 * (2.3 - (343.5 + 0.037 * (S / 1.00472)) / (1.00024 * T + 273.15)) * (1 - (1.00024 * T + 273.15) / (647.3 + 0.03 * (S / 1.00472)))^(1/3)); 
    Pr = mu * cp / sigma_T; Sc = mu^2 / (5.954e-15 * T_k * rho); 
    c_T = 0.072^(4/3) / Pr; c_S = 0.072^(4/3) / Sc; c_TS = 0.072^(4/3) * (Pr + Sc) / (2 * Pr * Sc);
    R_rho = 2.6e-4 * abs(H_ratio) / 7.6e-4;
    if R_rho >= 1, d_r = R_rho + sqrt(R_rho)*sqrt(R_rho - 1); elseif R_rho >= 0.5, d_r = 1.85*R_rho - 0.85; else, d_r = 0.15*R_rho; end
    chi_S = chi_T * d_r / (H_ratio^2); chi_TS = chi_T * (1 + d_r) / (2 * H_ratio); 
    nu_kinematic = mu / rho; 
    % eta_physical = (nu_kinematic^3 / epsilon)^(0.25); 
    eta_physical = 5e-3;
    Phi_Hill = @(K, chi_M, c_M) (0.72 / (4 * pi)) * chi_M * epsilon^(-1/3) .* (K.^2 + 1e-25).^(-11/6) .* exp(-176.90 * (K * eta_physical).^2 .* c_M^(0.96)) .* (1 + 21.61 * (K * eta_physical).^(0.61) .* c_M^(0.02) - 18.18 * (K * eta_physical).^(0.55) .* c_M^(0.04));
    Phi_n_func = @(K) (A^2 * Phi_Hill(K, chi_T, c_T) + B^2 * Phi_Hill(K, chi_S, c_S) + 2 * A * B * Phi_Hill(K, chi_TS, c_TS));
end
function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk = 2*pi/D; dx = D/N; kx = (-N/2 : N/2-1)*dk; [KX, KY] = meshgrid(kx, kx); K_grid = sqrt(KX.^2 + KY.^2); K_grid(N/2+1, N/2+1) = 1e-10;
    C_nm = (randn(N)+1i*randn(N)).*sqrt(2*pi*k_wave^2*delta_z*Phi_n_func(K_grid))*dk; phase_high = real(ifftshift(ifft2(ifftshift(C_nm)))) * N^2;
    phase_low = zeros(N, N); [xx, yy] = meshgrid((-N/2 : N/2-1)*dx);
    for p = 1:3, dk_p = dk/(3^p);
        for m = -1:1, for n = -1:1, if m==0&&n==0, continue; end
            k_p = sqrt((m*dk_p)^2+(n*dk_p)^2); phase_low = phase_low + real((randn(1)+1i*randn(1))*sqrt(2*pi^2*k_wave^2*delta_z*Phi_n_func(k_p))*dk_p * exp(1i*(m*dk_p*xx+n*dk_p*yy)));
        end; end
    end
    phase_screen = phase_high + phase_low;
end
function [rho0_exact, Cn2_eq] = calc_turb_coherence_params(Phi_n_func, k_wave, L, eta_physical)
    kappa_eval = 0.01 / eta_physical; Cn2_eq = Phi_n_func(kappa_eval) / (0.033 * kappa_eval^(-11/3));
    try, opts = optimset('Display', 'off'); rho0_exact = fzero(@(rho) 8*pi^2*k_wave^2*L*integral2(@(K, xi) K.*Phi_n_func(K).*(1-besselj(0, K.*rho.*xi)), 1e-2, 1e5, 0, 1, 'Method', 'iterated', 'RelTol', 1e-3) - 2, (0.545*k_wave^2*Cn2_eq*L)^(-3/5), opts); catch, rho0_exact = (0.545*k_wave^2*Cn2_eq*L)^(-3/5); end
end