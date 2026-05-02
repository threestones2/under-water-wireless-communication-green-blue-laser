%% Exp7: WCI-MC Pointing Error Simulation (Strong Turbulence - Ensemble Average)
% 实验目的: 分析收发偏轴误差对链路衰减的影响 (包含强湍流影响的系综平均无偏估计)
clc; clear; close all;

% ================= 物理与仿真参数初始化 =================
dist_cell = {10:10:60, 10:10:60, 5:5:20};
water_types = {'Clear Ocean', 'Coastal Ocean', 'Turbid Harbor'};
coef_a_arr = [0.114, 0.179, 0.366]; 
coef_b_arr = [0.0374, 0.219, 1.824]; 
coef_c_arr = [0.1514, 0.398, 2.190];
num_W = length(water_types);

off_axis_angles = [0.1, 1, 10];
num_A = length(off_axis_angles);

% --- 系综平均参数设置 ---
N_total_packets = 1e5;                   % 总光子数保持10^5不变
N_ensemble = 20;                         % 系综实现次数（生成20组不同随机相位的湍流屏）
N_per_ens = N_total_packets / N_ensemble; % 每个系综分配的光子数 (5000)
n_max = 200; 

lambda = 514e-9; k_wave = 2 * pi / lambda;

w0 = 0.002; div_angle = 0.1 * pi / 180; theta_half_div = div_angle / 2; 
Rx_Aperture = 0.01; Rx_FOV = 5 * pi / 180; Rx_Area = pi * (Rx_Aperture / 2)^2;
cos_FOV_half = cos(Rx_FOV / 2); Rx_Aperture_half_sq = (Rx_Aperture / 2)^2;

% ================= 强湍流环境设置 =================
turb_params.T_avg = 20; turb_params.S_avg = 35; turb_params.H_ratio = -1; 
turb_params.epsilon = 1e-10; turb_params.chi_T = 1e-5; turb_params.name = 'Strong';

N_screens = 20; D_screen = 1; N_grid = 2^8; dx = D_screen / N_grid; x_axis = (-N_grid/2 : N_grid/2-1) * dx;

PL_Cell_WCI_Strong = cell(1, num_W);

fprintf('--- 运行 Exp7: WCI-MC (Strong Turbulence - Ensemble Average) ---\n');
[Phi_func, ~, eta_physical] = get_OTOPS_spectrum_handle(lambda, turb_params.T_avg, turb_params.S_avg, turb_params.epsilon, turb_params.chi_T, turb_params.H_ratio);

for w_idx = 1:num_W
    param.coef_a = coef_a_arr(w_idx); param.coef_b = coef_b_arr(w_idx); param.coef_c = coef_c_arr(w_idx);
    param.albedo = param.coef_b / param.coef_c; param = calc_haltrin_params(param);
    P_max = pdf_Empirical(1.0, param); C_val = 4 * pi * P_max;
    g_prop = 1 + (1 - sqrt(8*C_val + 1)) / (2*C_val);
    
    dist_arr = dist_cell{w_idx}; num_D = length(dist_arr);
    PL_matrix = zeros(num_A, num_D);
    
    for a_idx = 1:num_A
        theta_err = off_axis_angles(a_idx) * pi / 180;
        
        % 统一采用严格的坐标旋转
        Link_Dir = [sin(theta_err), cos(theta_err), 0];
        u_vec_Tx = [cos(theta_err), -sin(theta_err), 0]; 
        v_vec_Tx = [0, 0, 1];
        fprintf('>>> 水质: %s | 偏轴角度: %.1f deg <<<\n', water_types{w_idx}, off_axis_angles(a_idx));
        
        for d_idx = 1:num_D
            L = dist_arr(d_idx);
            Tx_Pos = [0, 0, 0]; Rx_Pos = [0, L, 0]; Rx_Normal = [0, -1, 0];
            [rho0_Link, ~] = calc_turb_coherence_params(Phi_func, k_wave, L, eta_physical); 
            delta_z_screen = L / N_screens;
            
            Lx = Link_Dir(1); Ly = Link_Dir(2); Lz = Link_Dir(3); 
            Nx = Rx_Normal(1); Ny = Rx_Normal(2); Nz = Rx_Normal(3); 
            Tx = Tx_Pos(1); Ty = Tx_Pos(2); Tz = Tx_Pos(3); 
            Rx = Rx_Pos(1); Ry = Rx_Pos(2); Rz = Rx_Pos(3); 
            Ux = u_vec_Tx(1); Uy = u_vec_Tx(2); Uz = u_vec_Tx(3); 
            Vx = v_vec_Tx(1); Vy = v_vec_Tx(2); Vz = v_vec_Tx(3);
            
            P_rx_accum_total = 0; tic; 
            
            % --- 系综循环开始 ---
            for ens_idx = 1:N_ensemble
                % 构建复合种子，确保各距离、偏轴、系综间的流形独立且可复现
                rng_seed = 123456 + w_idx*10000 + a_idx*1000 + d_idx*100 + ens_idx;
                rng(rng_seed, 'twister');
                
                % 预生成当前系综的相位屏及其梯度分布
                Grad_X_3D = zeros(N_grid, N_grid, N_screens); 
                Grad_Y_3D = zeros(N_grid, N_grid, N_screens); 
                Screen_Z_1D = zeros(1, N_screens);
                for i = 1:N_screens
                    phi_screen = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
                    [gx, gy] = gradient(phi_screen, dx); 
                    Grad_X_3D(:,:,i) = gx; 
                    Grad_Y_3D(:,:,i) = gy; 
                    Screen_Z_1D(i) = i * delta_z_screen; 
                end
                
                P_rx_accum_ens = 0; 
                
                % 当前系综的光子传输
                for p = 1:N_per_ens
                    r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand(); cp0 = cos(phi0); sp0 = sin(phi0);
                    p1 = Tx + r0*cp0*Ux + r0*sp0*Vx; p2 = Ty + r0*cp0*Uy + r0*sp0*Vy; p3 = Tz + r0*cp0*Uz + r0*sp0*Vz;
                    U_init = theta_half_div * sqrt(-0.5 * log(rand())); 
                    dir = rotate_direction_fast(Link_Dir, cos(U_init), sin(U_init), 2*pi*rand());
                    d1 = dir(1); d2 = dir(2); d3 = dir(3); weight = 1.0; P_packet = 0;
                    
                    % --- 1. 直射光 (Ballistic) 分支 ---
                    [p1_end, p2_end, p3_end, d1_end, d2_end, d3_end, plane_hit, path_len_ballistic] = ...
                        ray_march_flat_scalar(p1, p2, p3, d1, d2, d3, 1e9, Rx, Ry, Rz, 1e10, -1, Nx, Ny, Nz, true, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x_axis(1), dx, N_grid, Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen);
                    
                    if plane_hit
                        cos_rx_tilt = abs(d1_end*Nx + d2_end*Ny + d3_end*Nz);
                        if cos_rx_tilt >= cos_FOV_half
                            vec_x = p1_end - Rx; vec_y = p2_end - Ry; vec_z = p3_end - Rz;
                            cx = vec_y*d3_end - vec_z*d2_end; cy = vec_z*d1_end - vec_x*d3_end; cz = vec_x*d2_end - vec_y*d1_end;
                            r_wander_perp = sqrt(cx^2 + cy^2 + cz^2); 
                            r_eff = (Rx_Aperture/2) * sqrt(cos_rx_tilt);
                            
                            rho0_ballistic = rho0_Link * (L / path_len_ballistic)^(3/5);
                            W_turb_LT = 2 * path_len_ballistic / (k_wave * rho0_ballistic);
                            cf = max(0, 1 - 0.37 * (rho0_ballistic / (2*r_eff))^(1/3));
                            W_spot = W_turb_LT * cf; 
                            
                            if W_spot > 1e-6 
                                recv_frac = 1 - marcumq(2*r_wander_perp/W_spot, 2*r_eff/W_spot); 
                            else 
                                recv_frac = double(r_wander_perp <= r_eff); 
                            end
                            P_packet = P_packet + exp(-param.coef_c * path_len_ballistic) * recv_frac;
                        end
                    end
                    
                    % --- 2. 散射光 (Scattered) 游走分支 ---
                    for ord = 1:n_max
                        d_step = -log(rand()) / param.coef_b; 
                        [p1, p2, p3, d1, d2, d3, ~, step_len] = ray_march_flat_scalar(p1, p2, p3, d1, d2, d3, d_step, Rx, Ry, Rz, Rx_Aperture_half_sq, cos_FOV_half, Nx, Ny, Nz, false, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x_axis(1), dx, N_grid, Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen);
                        
                        weight = weight * exp(-param.coef_a * step_len); 
                        if (p1 - Tx)*Lx + (p2 - Ty)*Ly + (p3 - Tz)*Lz >= L || weight < 1e-15, break; end 
                        
                        vec2rx_1 = Rx - p1; vec2rx_2 = Ry - p2; vec2rx_3 = Rz - p3; 
                        d2rx_sq = vec2rx_1^2 + vec2rx_2^2 + vec2rx_3^2; dist2rx = sqrt(d2rx_sq); 
                        dir2rx_1 = vec2rx_1 / dist2rx; dir2rx_2 = vec2rx_2 / dist2rx; dir2rx_3 = vec2rx_3 / dist2rx;
                        cos_inc = -(dir2rx_1*Nx + dir2rx_2*Ny + dir2rx_3*Nz);
                        
                        if cos_inc >= cos_FOV_half
                            omega = Rx_Area / d2rx_sq * cos_inc; cos_scatter = d1*dir2rx_1 + d2*dir2rx_2 + d3*dir2rx_3;
                            base_w = weight * min(1, pdf_Empirical(cos_scatter, param) * omega) * exp(-param.coef_c * dist2rx);
                            
                            if base_w > 1e-16
                                [p1_v, p2_v, p3_v, d1_v, d2_v, d3_v, v_hit, v_len] = ray_march_flat_scalar(p1, p2, p3, dir2rx_1, dir2rx_2, dir2rx_3, dist2rx+1e-1, Rx, Ry, Rz, 1e10, -1, Nx, Ny, Nz, true, Grad_X_3D, Grad_Y_3D, Screen_Z_1D, Ux, Uy, Uz, Vx, Vy, Vz, k_wave, x_axis(1), dx, N_grid, Tx, Ty, Tz, Lx, Ly, Lz, delta_z_screen);
                                
                                if v_hit
                                    cos_inc_v = abs(d1_v*Nx + d2_v*Ny + d3_v*Nz); r_eff_v = (Rx_Aperture/2) * sqrt(cos_inc_v);
                                    vec_x = p1_v - Rx; vec_y = p2_v - Ry; vec_z = p3_v - Rz; cx = vec_y*d3_v - vec_z*d2_v; cy = vec_z*d1_v - vec_x*d3_v; cz = vec_x*d2_v - vec_y*d1_v; r_wp = sqrt(cx^2 + cy^2 + cz^2);
                                    
                                    point_loss = double(r_wp <= r_eff_v); 
                                    P_packet = P_packet + base_w * point_loss;
                                end
                            end
                        end
                        
                        xi = rand(); term = (1 - g_prop^2) / (1 - g_prop + 2 * g_prop * xi); cos_th_i = (1 + g_prop^2 - term^2) / (2 * g_prop); 
                        if cos_th_i > 1, cos_th_i = 1; elseif cos_th_i < -1, cos_th_i = -1; end
                        denom_core = 1 + g_prop^2 - 2 * g_prop * cos_th_i; q_HG = (1 - g_prop^2) / (4 * pi * denom_core * sqrt(denom_core)); p_val = pdf_Empirical(cos_th_i, param); 
                        
                        weight = weight * max(0.5,min(2, p_val / q_HG)); 
                        dir_new = rotate_direction_fast([d1, d2, d3], cos_th_i, sqrt(1 - cos_th_i^2), 2*pi*rand()); d1 = dir_new(1); d2 = dir_new(2); d3 = dir_new(3);
                    end
                    P_rx_accum_ens = P_rx_accum_ens + P_packet;
                end
                
                % 累加系综能量
                P_rx_accum_total = P_rx_accum_total + P_rx_accum_ens;
            end
            % --- 系综循环结束 ---
            
            % 使用总包数归一化，计算系综平均后的路径损耗
            PL_matrix(a_idx, d_idx) = 10 * log10(max(P_rx_accum_total / N_total_packets, 1e-300));
            fprintf('    L=%2dm | Time: %5.2fs | 系综平均 Path Loss: %6.2f dB\n', L, toc, -PL_matrix(a_idx, d_idx));
        end
    end
    PL_Cell_WCI_Strong{w_idx} = PL_matrix;
end
save('data_exp7_WCIMC_Strong_Ensemble.mat', 'dist_cell', 'PL_Cell_WCI_Strong', 'water_types', 'off_axis_angles');

% ================= 辅助函数区域 =================
function [p1,p2,p3,d1,d2,d3,hit_flag,total_len]=ray_march_flat_scalar(p1,p2,p3,d1,d2,d3,dist_limit,Rx,Ry,Rz,Rx_Aperture_half_sq,cos_FOV_half,Nx,Ny,Nz,enable_hit_check,Grad_X_3D,Grad_Y_3D,Screen_Z_1D,Ux,Uy,Uz,Vx,Vy,Vz,k_wave,x0,dx,N_grid,Tx,Ty,Tz,Lx,Ly,Lz,delta_z_screen)
hit_flag=false; total_len=0; rem_dist=dist_limit; N_screens=length(Screen_Z_1D);
while rem_dist>1e-9
    dir_z=d1*Lx+d2*Ly+d3*Lz; z_pos=(p1-Tx)*Lx+(p2-Ty)*Ly+(p3-Tz)*Lz; t_scr=inf; target_idx=-1;
    if dir_z>1e-10, target_idx=floor((z_pos+1e-9)/delta_z_screen)+1; if target_idx<1, target_idx=1; end; if target_idx<=N_screens, t_scr=(Screen_Z_1D(target_idx)-z_pos)/dir_z; end
    elseif dir_z<-1e-10, target_idx=ceil((z_pos-1e-9)/delta_z_screen)-1; if target_idx>N_screens, target_idx=N_screens; end; if target_idx>=1, t_scr=(Screen_Z_1D(target_idx)-z_pos)/dir_z; end; end
    t_rx=inf; 
    if enable_hit_check
        denom_rx=d1*Nx+d2*Ny+d3*Nz; 
        if abs(denom_rx)>1e-10
            t_temp=((Rx-p1)*Nx+(Ry-p2)*Ny+(Rz-p3)*Nz)/denom_rx; 
            if t_temp>-1e-7, t_rx=max(0,t_temp); end
        end
    end
    [min_dist,event_idx]=min([rem_dist,t_rx,t_scr]); 
    p1=p1+d1*min_dist; p2=p2+d2*min_dist; p3=p3+d3*min_dist; 
    rem_dist=rem_dist-min_dist; total_len=total_len+min_dist;
    if event_idx==2
        if (p1-Rx)^2+(p2-Ry)^2+(p3-Rz)^2<=Rx_Aperture_half_sq && -(d1*Nx+d2*Ny+d3*Nz)>=cos_FOV_half
            hit_flag=true; 
        end
        break; 
    elseif event_idx==3
        loc_u=(p1-Tx)*Ux+(p2-Ty)*Uy+(p3-Tz)*Uz; loc_v=(p1-Tx)*Vx+(p2-Ty)*Vy+(p3-Tz)*Vz; 
        idx_x=mod(round((loc_u-x0)/dx),N_grid)+1; idx_y=mod(round((loc_v-x0)/dx),N_grid)+1; 
        gx=Grad_X_3D(idx_y,idx_x,target_idx); gy=Grad_Y_3D(idx_y,idx_x,target_idx); 
        d1=d1+(gx*Ux+gy*Vx)/k_wave; d2=d2+(gx*Uy+gy*Vy)/k_wave; d3=d3+(gx*Uz+gy*Vz)/k_wave; 
        d_norm=sqrt(d1^2+d2^2+d3^2); d1=d1/d_norm; d2=d2/d_norm; d3=d3/d_norm;
    elseif event_idx==1
        break; 
    end
end
end

function param=calc_haltrin_params(p), b_e=p.coef_b; al_e=p.albedo; p.q_e=2.598+17.748*sqrt(b_e)-16.722*b_e+5.932*b_e*sqrt(b_e); p.k1=1.188-0.688*al_e; p.k2=0.1*(3.07-1.90*al_e); p.k3=0.01*(4.58-3.02*al_e); p.k4=0.001*(3.24-2.25*al_e); p.k5=0.0001*(0.84-0.61*al_e); th=linspace(0,pi,2000); val=zeros(size(th)); for i=1:2000, t=max(th(i)*180/pi,1e-6); sq_t=sqrt(t); t_sq=t*t; term=1-p.k1*sq_t+p.k2*t-p.k3*t*sq_t+p.k4*t_sq-p.k5*t_sq*sq_t; val(i)=exp(p.q_e*term); end; p.b_emp_norm=2*pi*trapz(th,val.*sin(th)); param=p; end
function p=pdf_Empirical(cos_theta,param), t_deg=max(acos(cos_theta)*180/pi,1e-6); sq_t=sqrt(t_deg); t_sq=t_deg*t_deg; term=1-param.k1*sq_t+param.k2*t_deg-param.k3*t_deg*sq_t+param.k4*t_sq-param.k5*t_sq*sq_t; p=exp(param.q_e*term)/param.b_emp_norm; end
function nd=rotate_direction_fast(d,ct,st,psi_s), denom=sqrt(1-d(3)^2); sp=sin(psi_s); cp=cos(psi_s); if denom<1e-10, nd=[st*cp,st*sp,sign(d(3))*ct]; else, nd=[st/denom*(d(1)*d(3)*cp-d(2)*sp)+d(1)*ct, st/denom*(d(2)*d(3)*cp+d(1)*sp)+d(2)*ct, -st*cp*denom+d(3)*ct]; end; end
function [Phi_n,k_wave,eta]=get_OTOPS_spectrum_handle(lam,T,S,eps,chi_T,H_r), lam_nm=lam*1e9; k_wave=2*pi/lam; A=-1.05e-6*S+2*1.6e-8*T*S-2*2.02e-6*T-4.23e-3/lam_nm; B=1.779e-4-1.05e-6*T+1.6e-8*T^2+1.155e-2/lam_nm; s_f=S*1e-3; T_k=T+273.15; cp=1000*((5.328-9.76e-2*S+4.04e-4*S^2)+(-6.913e-3+7.351e-4*S-3.15e-6*S^2)*T+(9.6e-6-1.927e-6*S+8.23e-9*S^2)*T^2+(2.5e-9+1.666e-9*S-7.125e-12*S^2)*T^3); rho=(9.9992293295e2+2.0341179217e-2*T-6.1624591598e-3*T^2+2.2614664708e-5*T^3-4.6570659168e-8*T^4)+s_f*(8.0200240891e2-2.0005183488*T+1.6771024982e-2*T^2-3.0600536746e-5*T^3-1.6132224742e-5*T*S); mu=((0.15700386464*(T+64.992620050)^2-91.296496657)^(-1)+4.2844324477e-5)*(1+(1.5409136040+1.9981117208e-2*T-9.5203865864e-5*T^2)*s_f+(7.9739318223-7.561456881e-2*T+4.7237011074e-4*T^2)*s_f^2); sigma_T=10^(log10(240+0.0002*(S/1.00472))-3+0.434*(2.3-(343.5+0.037*(S/1.00472))/(1.00024*T+273.15))*(1-(1.00024*T+273.15)/(647.3+0.03*(S/1.00472)))^(1/3)); Pr=mu*cp/sigma_T; Sc=mu^2/(5.954e-15*T_k*rho); c_T=0.072^(4/3)/Pr; c_S=0.072^(4/3)/Sc; c_TS=0.072^(4/3)*(Pr+Sc)/(2*Pr*Sc); R_rho=2.6e-4*abs(H_r)/7.6e-4; if R_rho>=1, d_r=R_rho+sqrt(R_rho)*sqrt(R_rho-1); elseif R_rho>=0.5, d_r=1.85*R_rho-0.85; else, d_r=0.15*R_rho; end; chi_S=chi_T*d_r/(H_r^2); chi_TS=chi_T*(1+d_r)/(2*H_r); eta=(mu/rho)^0.75/eps^0.25; Phi_Hill=@(K,chi_M,c_M) (0.72/(4*pi))*chi_M*eps^(-1/3).*(K.^2+1e-25).^(-11/6).*exp(-176.90*(K*eta).^2.*c_M^(0.96)).*(1+21.61*(K*eta).^(0.61).*c_M^(0.02)-18.18*(K*eta).^(0.55).*c_M^(0.04)); Phi_n=@(K) (A^2*Phi_Hill(K,chi_T,c_T)+B^2*Phi_Hill(K,chi_S,c_S)+2*A*B*Phi_Hill(K,chi_TS,c_TS)); end
function screen=gen_screen_from_spectrum(Phi_n,D,N,k_wave,dz), dk=2*pi/D; dx=D/N; kx=(-N/2:N/2-1)*dk; [KX,KY]=meshgrid(kx,kx); K_grid=sqrt(KX.^2+KY.^2); K_grid(N/2+1,N/2+1)=1e-10; C_nm=(randn(N)+1i*randn(N)).*sqrt(2*pi*k_wave^2*dz*Phi_n(K_grid))*dk; phase_high=real(ifftshift(ifft2(ifftshift(C_nm))))*N^2; phase_low=zeros(N,N); [xx,yy]=meshgrid((-N/2:N/2-1)*dx); for p=1:3, dk_p=dk/(3^p); for m=-1:1, for n=-1:1, if m==0&&n==0, continue; end; k_p=sqrt((m*dk_p)^2+(n*dk_p)^2); phase_low=phase_low+real((randn(1)+1i*randn(1))*sqrt(2*pi^2*k_wave^2*dz*Phi_n(k_p))*dk_p*exp(1i*(m*dk_p*xx+n*dk_p*yy))); end; end; end; screen=phase_high+phase_low; end
function [rho0,Cn2_eq]=calc_turb_coherence_params(Phi_n,k_wave,L,eta), kappa_eval=0.01/eta; Cn2_eq=Phi_n(kappa_eval)/(0.033*kappa_eval^(-11/3)); try opts=optimset('Display','off'); rho0=fzero(@(rho) 8*pi^2*k_wave^2*L*integral2(@(K,xi) K.*Phi_n(K).*(1-besselj(0,K.*rho.*xi)),1e-2,1e5,0,1,'Method','iterated','RelTol',1e-3)-2,(0.545*k_wave^2*Cn2_eq*L)^(-3/5),opts); catch, rho0=(0.545*k_wave^2*Cn2_eq*L)^(-3/5); end; end