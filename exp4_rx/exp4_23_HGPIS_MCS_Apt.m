%% Exp4: HG-PIS and MCS Algorithms for Different Receiver Apertures
clc; clear; close all;

N_packets = 1e5; n_max = 200; dist_axis = 5:5:60; num_dist = length(dist_axis);
w0 = 0.002; div_angle = 0.1 * pi / 180; theta_half_div = div_angle / 2; 
Rx_FOV = 5 * pi / 180; cos_FOV_half = cos(Rx_FOV / 2);

coef_c = 0.1514; param.coef_a = 0.114; param.coef_b = 0.0374; param.coef_c = coef_c; param.albedo = param.coef_b / coef_c;
param = calc_haltrin_params(param); P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max; g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);

% O(1) LUT for MCS
th_axis = linspace(0, pi, 50000); pdf_vals = zeros(size(th_axis));
for i = 1:length(th_axis), pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i)); end
cdf_vals = cumtrapz(th_axis, pdf_vals); cdf_vals = cdf_vals / cdf_vals(end);
[cdf_uniq, idx_uniq] = unique(cdf_vals); LUT_SIZE = 100000; P_grid = linspace(0, 1, LUT_SIZE);
th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap'); cos_th_LUT = cos(th_LUT); sin_th_LUT = sin(th_LUT);

apertures_m = [0.01, 0.05, 0.10]; num_apt = length(apertures_m);
PL_Cell_HG = cell(1, num_apt); PL_Cell_MCS = cell(1, num_apt);

for i = 1:num_apt
    Rx_Aperture = apertures_m(i); r_rx = Rx_Aperture / 2; Rx_Area = pi * r_rx^2; Rx_Aperture_half_sq = r_rx^2;
    fprintf('--- 孔径: %.2f m ---\n', Rx_Aperture);
    PL_HG_temp = zeros(1, num_dist); PL_MCS_temp = zeros(1, num_dist);
    
    for d_idx = 1:num_dist
        L = dist_axis(d_idx); Tx = 0; Ty = 0; Tz = 0; Rx = 0; Ry = L; Rz = 0;
        Lx = 0; Ly = 1; Lz = 0; Nx = 0; Ny = -1; Nz = 0; Ux = 1; Uy = 0; Uz = 0; Vx = 0; Vy = 0; Vz = 1;
        
        P_hg_accum = 0; P_mcs_accum = 0; tic; rng(123456, 'twister'); 
        for p = 1:N_packets
            r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand(); p1 = Tx + r0*cos(phi0)*Ux + r0*sin(phi0)*Vx; p2 = Ty + r0*cos(phi0)*Uy + r0*sin(phi0)*Vy; p3 = Tz + r0*cos(phi0)*Uz + r0*sin(phi0)*Vz;
            p1_mcs = p1; p2_mcs = p2; p3_mcs = p3;
            U_init = theta_half_div * sqrt(-0.5 * log(rand())); dir = rotate_direction_fast([Lx,Ly,Lz], cos(U_init), sin(U_init), 2*pi*rand());
            d1 = dir(1); d2 = dir(2); d3 = dir(3); d1_mcs = d1; d2_mcs = d2; d3_mcs = d3;
            w_hg = 1.0; w_mcs = 1.0; P_packet_hg = 0; P_packet_mcs = 0;
            
            % Zero-order hit
            cos_th = d1*Lx + d2*Ly + d3*Lz;
            if cos_th > 0
                d_plane = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / cos_th; pos_end_1 = p1 + d1 * d_plane; pos_end_2 = p2 + d2 * d_plane; pos_end_3 = p3 + d3 * d_plane; 
                if abs(d1*Nx + d2*Ny + d3*Nz) >= cos_FOV_half && (pos_end_1 - Rx)^2 + (pos_end_2 - Ry)^2 + (pos_end_3 - Rz)^2 <= Rx_Aperture_half_sq
                    P_packet_hg = exp(-param.coef_c * d_plane); P_packet_mcs = P_packet_hg;
                end
            end
            
            % HG-PIS Scatter
            for ord = 1:n_max
                d_step = -log(rand()) / param.coef_b; p1 = p1 + d1 * d_step; p2 = p2 + d2 * d_step; p3 = p3 + d3 * d_step; w_hg = w_hg * exp(-param.coef_a * d_step);
                if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L, break; end; if w_hg < 1e-9, if rand() > 0.1, break; else, w_hg = w_hg * 10; end; end
                vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; d2rx = sqrt(d2rx_sq); dir2rx_1 = vec2rx_1 / d2rx; dir2rx_2 = vec2rx_2 / d2rx; dir2rx_3 = vec2rx_3 / d2rx; cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                if cos_inc >= cos_FOV_half, P_packet_hg = P_packet_hg + w_hg * min(1, pdf_Empirical(d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3, param) * (Rx_Area / d2rx_sq * cos_inc)) * exp(-param.coef_c * d2rx); end
                xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); cos_th_i = (1 + g_prop^2 - term^2) / (2 * g_prop); if cos_th_i > 1, cos_th_i = 1; elseif cos_th_i < -1, cos_th_i = -1; end
                denom_core = 1 + g_prop^2 - 2 * g_prop * cos_th_i; q_HG = (1 - g_prop^2) / (4 * pi * denom_core * sqrt(denom_core)); w_hg = w_hg * (pdf_Empirical(cos_th_i, param) / q_HG);
                dir_new = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1 - cos_th_i^2), 2*pi*rand()); d1 = dir_new(1); d2 = dir_new(2); d3 = dir_new(3);
            end
            P_hg_accum = P_hg_accum + P_packet_hg;
            
            % MCS Scatter
            rng(123456 + p*100, 'twister'); % ensure independence
            for ord = 1:n_max
                d_step = -log(rand()) / param.coef_b; p1_mcs = p1_mcs + d1_mcs * d_step; p2_mcs = p2_mcs + d2_mcs * d_step; p3_mcs = p3_mcs + d3_mcs * d_step; w_mcs = w_mcs * exp(-param.coef_a * d_step);
                if (p1_mcs - Tx)*Lx + (p2_mcs - Ty)*Ly + (p3_mcs - Tz)*Lz >= L || w_mcs < 1e-15, break; end 
                vec2rx_1 = Rx - p1_mcs; vec2rx_2 = Ry - p2_mcs; vec2rx_3 = Rz - p3_mcs; d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; d2rx = sqrt(d2rx_sq); dir2rx_1 = vec2rx_1 / d2rx; dir2rx_2 = vec2rx_2 / d2rx; dir2rx_3 = vec2rx_3 / d2rx; cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                if cos_inc >= cos_FOV_half, P_packet_mcs = P_packet_mcs + w_mcs * min(1, pdf_Empirical(d1_mcs*dir2rx_1 + d2_mcs*dir2rx_2 + d3_mcs*dir2rx_3, param) * (Rx_Area / d2rx_sq * cos_inc)) * exp(-param.coef_c * d2rx); end
                u_rand = rand(); idx_lut = floor(u_rand * (LUT_SIZE - 1)) + 1; ct_s = cos_th_LUT(idx_lut); st_s = sin_th_LUT(idx_lut);
                dir_new = rotate_direction_fast([d1_mcs, d2_mcs, d3_mcs], ct_s, st_s, 2*pi*rand()); d1_mcs = dir_new(1); d2_mcs = dir_new(2); d3_mcs = dir_new(3);
            end
            P_mcs_accum = P_mcs_accum + P_packet_mcs;
        end
        PL_HG_temp(d_idx) = 10 * log10(max(P_hg_accum / N_packets, 1e-300));
        PL_MCS_temp(d_idx) = 10 * log10(max(P_mcs_accum / N_packets, 1e-300));
        fprintf('    L=%2dm | Time: %5.2fs | HG: %6.2f dB | MCS: %6.2f dB\n', L, toc, -PL_HG_temp(d_idx), -PL_MCS_temp(d_idx));
    end
    PL_Cell_HG{i} = PL_HG_temp; PL_Cell_MCS{i} = PL_MCS_temp;
end
save('data_exp4_apt_HGPIS.mat', 'dist_axis', 'apertures_m', 'PL_Cell_HG', 'coef_c');
save('data_exp4_apt_MCS.mat', 'dist_axis', 'apertures_m', 'PL_Cell_MCS', 'coef_c');
% === Helper Functions ===
function param=calc_haltrin_params(p), b_e=p.coef_b; al_e=p.albedo; p.q_e=2.598+17.748*sqrt(b_e)-16.722*b_e+5.932*b_e*sqrt(b_e); p.k1=1.188-0.688*al_e; p.k2=0.1*(3.07-1.90*al_e); p.k3=0.01*(4.58-3.02*al_e); p.k4=0.001*(3.24-2.25*al_e); p.k5=0.0001*(0.84-0.61*al_e); th=linspace(0,pi,2000); val=zeros(size(th)); for i=1:2000, t=max(th(i)*180/pi,1e-6); sq_t=sqrt(t); t_sq=t*t; term=1-p.k1*sq_t+p.k2*t-p.k3*t*sq_t+p.k4*t_sq-p.k5*t_sq*sq_t; val(i)=exp(p.q_e*term); end; p.b_emp_norm=2*pi*trapz(th,val.*sin(th)); param=p; end
function p=pdf_Empirical(cos_theta,param), t_deg=max(acos(cos_theta)*180/pi,1e-6); sq_t=sqrt(t_deg); t_sq=t_deg*t_deg; term=1-param.k1*sq_t+param.k2*t_deg-param.k3*t_deg*sq_t+param.k4*t_sq-param.k5*t_sq*sq_t; p=exp(param.q_e*term)/param.b_emp_norm; end
function nd=rotate_direction_fast(d,ct,st,psi_s), denom=sqrt(1-d(3)^2); sp=sin(psi_s); cp=cos(psi_s); if denom<1e-10, nd=[st*cp,st*sp,sign(d(3))*ct]; else, nd=[st/denom*(d(1)*d(3)*cp-d(2)*sp)+d(1)*ct, st/denom*(d(2)*d(3)*cp+d(1)*sp)+d(2)*ct, -st*cp*denom+d(3)*ct]; end; nd=nd/norm(nd); end