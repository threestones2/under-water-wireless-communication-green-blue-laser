%% Exp8: SA-MCS CIR Simulation (支持 Haltrin 与 HG 相函数切换)
clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

N_packets = 1e7; 
n_max = 200; 
c_water = 2.237e8;

% --- 核心工程参数设定 ---
w0 = 0.002;                  
div_angle = 1 * pi / 180;  
theta_half_div = div_angle / 2; 

Rx_Aperture = 0.01;          
r_rx = Rx_Aperture / 2; 
Rx_Area = pi * r_rx^2; 
Rx_Aperture_half_sq = r_rx^2;
Rx_FOV = 90* pi / 180;     
cos_FOV_half = cos(Rx_FOV / 2);

% ================= 相函数模式切换 =================
use_HG_phase = false;         % true: 使用 HG 相函数; false: 使用 Haltrin 查表
g_HG = 0.924;                % HG 相函数的非对称因子
% =================================================

water_params = struct();
water_params(1).name = 'Clear Ocean (50m)';   water_params(1).a = 0.114; water_params(1).b = 0.0374; water_params(1).c = 0.1514; water_params(1).L = 50;
water_params(2).name = 'Coastal Ocean (20m)'; water_params(2).a = 0.179; water_params(2).b = 0.219;  water_params(2).c = 0.398;  water_params(2).L = 20;
num_water = length(water_params);

dt = 0.01e-9;         
delta_t_max = 10e-9;  

CIR_MCS_Scattered = cell(1, num_water);
E_Ballistic_Total = cell(1, num_water); 
Time_Axis = cell(1, num_water);

for w = 1:num_water
    wp = water_params(w);
    L = wp.L;
    
    param.coef_a = wp.a; param.coef_b = wp.b; param.coef_c = wp.c; param.albedo = wp.b / wp.c;
    
    if ~use_HG_phase
        param = calc_haltrin_params(param);
        th_axis = linspace(0, pi, 50000); pdf_vals = zeros(size(th_axis));
        for i = 1:length(th_axis), pdf_vals(i) = pdf_Empirical(cos(th_axis(i)), param) * 2 * pi * sin(th_axis(i)); end
        cdf_vals = cumtrapz(th_axis, pdf_vals); cdf_vals = cdf_vals / cdf_vals(end);
        [cdf_uniq, idx_uniq] = unique(cdf_vals); LUT_SIZE = 100000; P_grid = linspace(0, 1, LUT_SIZE);
        th_LUT = interp1(cdf_uniq, th_axis(idx_uniq), P_grid, 'linear', 'extrap'); 
        cos_th_LUT = cos(th_LUT); sin_th_LUT = sin(th_LUT);
    end

    t_LOS = L / c_water; t_min = t_LOS; t_max = t_LOS + delta_t_max;
    T_bins = t_min : dt : t_max; N_bins = length(T_bins); Time_Axis{w} = T_bins - t_LOS; 
    
    h_time_MCS = zeros(1, N_bins); E_ballistic_accum = 0;         
    

    if use_HG_phase
        phase_str = sprintf('HG (g=%.3f)', g_HG);
    else
        phase_str = 'Haltrin';
    end
    
    fprintf('\n======================================================\n');
    fprintf('>>> SA-MCS CIR 仿真 | 相函数: %s | %s <<<\n', phase_str, wp.name);
    % ---------------------------------
    
    Tx = 0; Ty = 0; Tz = 0; Rx = 0; Ry = L; Rz = 0;
    Lx = 0; Ly = 1; Lz = 0; Nx = 0; Ny = -1; Nz = 0; 
    
    tic; rng(123456, 'twister'); 
    for p = 1:N_packets
        r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
        p1 = Tx + r0*cos(phi0); p2 = Ty; p3 = Tz + r0*sin(phi0);
        U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
        dir = rotate_direction_fast([Lx,Ly,Lz], cos(U_init), sin(U_init), 2*pi*rand()); 
        d1 = dir(1); d2 = dir(2); d3 = dir(3); 
        
        w_mcs = 1.0; current_dist = 0;
        
        cos_th = d1*Lx + d2*Ly + d3*Lz;
        if cos_th > 0
            d_plane = ((Rx - p1)*Lx + (Ry - p2)*Ly + (Rz - p3)*Lz) / cos_th; 
            pos_e1 = p1 + d1 * d_plane; pos_e2 = p2 + d2 * d_plane; pos_e3 = p3 + d3 * d_plane; 
            if abs(d1*Nx + d2*Ny + d3*Nz) >= cos_FOV_half && (pos_e1 - Rx)^2 + (pos_e2 - Ry)^2 + (pos_e3 - Rz)^2 <= Rx_Aperture_half_sq
                E_ballistic_accum = E_ballistic_accum + exp(-wp.c * d_plane); 
            end
        end
        
        for ord = 1:n_max
            d_step = -log(rand()) / wp.b; 
            p1 = p1 + d1 * d_step; p2 = p2 + d2 * d_step; p3 = p3 + d3 * d_step; 
            w_mcs = w_mcs * exp(-wp.a * d_step); 
            current_dist = current_dist + d_step;
            
            if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L || w_mcs < 1e-15, break; end 
            
            vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; 
            d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; dist2rx = sqrt(d2rx_sq); 
            dir2rx_1 = vec2rx_1 / dist2rx; dir2rx_2 = vec2rx_2 / dist2rx; dir2rx_3 = vec2rx_3 / dist2rx; 
            cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
            dist_total = current_dist + dist2rx;
            
            if cos_inc >= cos_FOV_half && dist_total > (L + 1e-9)
                cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                
                % [核心替换] 解析求取 PIS 的相函数权重
                if use_HG_phase
                    denom_core = 1 + g_HG^2 - 2 * g_HG * cos_scatter;
                    pdf_val = (1 - g_HG^2) / (4 * pi * denom_core * sqrt(denom_core));
                else
                    pdf_val = pdf_Empirical(cos_scatter, param);
                end
                
                omega = Rx_Area / d2rx_sq * cos_inc;
                base_w = w_mcs * min(1, pdf_val * omega) * exp(-wp.c * dist2rx); 
                
                idx = floor((dist_total / c_water - t_min)/dt) + 1;
                if idx >= 1 && idx <= N_bins, h_time_MCS(idx) = h_time_MCS(idx) + base_w; end
            end
            
            if w_mcs < 1e-9
                if rand() > 0.1, break; else, w_mcs = w_mcs * 10; end
            end
            
            % [核心替换] 解析求取物理散射角
            if use_HG_phase
                xi = rand();
                term = (1 - g_HG^2) / (1 - g_HG + 2 * g_HG * xi);
                cos_th_s = (1 + g_HG^2 - term^2) / (2 * g_HG);
                if cos_th_s > 1, cos_th_s = 1; elseif cos_th_s < -1, cos_th_s = -1; end
                sin_th_s = sqrt(1 - cos_th_s^2);
            else
                idx_lut = min(LUT_SIZE, floor(rand() * LUT_SIZE) + 1); 
                cos_th_s = cos_th_LUT(idx_lut); sin_th_s = sin_th_LUT(idx_lut);
            end
            
            dir_new = rotate_direction_fast([d1, d2, d3], cos_th_s, sin_th_s, 2*pi*rand()); 
            d1 = dir_new(1); d2 = dir_new(2); d3 = dir_new(3);
        end
    end
    CIR_MCS_Scattered{w} = h_time_MCS / N_packets;
    E_Ballistic_Total{w} = E_ballistic_accum / N_packets;
    fprintf('  计算完成, 耗时: %5.2f s\n', toc);
end

save('data_exp8_CIR_MCS.mat', 'water_params', 'Time_Axis', 'CIR_MCS_Scattered', 'E_Ballistic_Total', 'dt');
disp('SA-MCS CIR 数据已保存！');

% --- 辅助函数区域 ---
function param=calc_haltrin_params(p), b_e=p.coef_b; al_e=p.albedo; p.q_e=2.598+17.748*sqrt(b_e)-16.722*b_e+5.932*b_e*sqrt(b_e); p.k1=1.188-0.688*al_e; p.k2=0.1*(3.07-1.90*al_e); p.k3=0.01*(4.58-3.02*al_e); p.k4=0.001*(3.24-2.25*al_e); p.k5=0.0001*(0.84-0.61*al_e); th=linspace(0,pi,2000); val=zeros(size(th)); for i=1:2000, t=max(th(i)*180/pi,1e-6); sq_t=sqrt(t); t_sq=t*t; term=1-p.k1*sq_t+p.k2*t-p.k3*t*sq_t+p.k4*t_sq-p.k5*t_sq*sq_t; val(i)=exp(p.q_e*term); end; p.b_emp_norm=2*pi*trapz(th,val.*sin(th)); param=p; end
function p=pdf_Empirical(cos_theta,param), t_deg=max(acos(cos_theta)*180/pi,1e-6); sq_t=sqrt(t_deg); t_sq=t_deg*t_deg; term=1-param.k1*sq_t+param.k2*t_deg-param.k3*t_deg*sq_t+param.k4*t_sq-param.k5*t_sq*sq_t; p=exp(param.q_e*term)/param.b_emp_norm; end
function nd=rotate_direction_fast(d,ct,st,psi_s), denom=sqrt(1-d(3)^2); sp=sin(psi_s); cp=cos(psi_s); if denom<1e-10, nd=[st*cp,st*sp,sign(d(3))*ct]; else, nd=[st/denom*(d(1)*d(3)*cp-d(2)*sp)+d(1)*ct, st/denom*(d(2)*d(3)*cp+d(1)*sp)+d(2)*ct, -st*cp*denom+d(3)*ct]; end; nd=nd/norm(nd); end