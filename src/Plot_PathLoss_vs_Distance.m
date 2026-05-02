%% Path Loss & Mean Delay vs. Distance (Unified Petzold Clear Ocean)
% 功能: 
% 1. 扫描不同传输距离下的路径衰减 (Path Loss)
% 2. 扫描不同传输距离下的平均到达时延 (Mean Delay / Simulation Time)
% 3. 对比四种相函数 (HG, TTHG, Mix, Empirical) 的差异
% 
% 解释: "仿真时间"在此处解读为物理层的"平均传播时延"(Mean Propagation Delay)，
% 即光信号穿过水下信道所需的平均飞行时间。

clc; clear; close all;

%% ================= 1. 扫描参数设置 =================
% --- 距离扫描范围 (米) ---
Dist_List = 10 : 10 : 100;  % 距离范围 [10m, 100m]

% --- 待对比的相函数列表 ---
Func_Names = {'HG', 'TTHG', 'Mix', 'Empirical'};

% --- 仿真光子数 ---
% 建议: 调试用 1e4, 论文出图建议 1e5 或更高
N_packets = 2e4; 

% --- 结果存储 ---
% 行: 模式, 列: 距离点
Result_Loss  = zeros(length(Func_Names), length(Dist_List));
Result_Delay = zeros(length(Func_Names), length(Dist_List));

%% ================= 2. 主循环 =================
fprintf('=== 开始 Path Loss & Delay 仿真 ===\n');
fprintf('光学参数: Petzold Clear Ocean (c=0.1514, a=0.114)\n');
fprintf('扫描距离: %d m ~ %d m\n', min(Dist_List), max(Dist_List));

total_steps = length(Func_Names) * length(Dist_List);
current_step = 0;
t_start_all = tic;

for f_idx = 1 : length(Func_Names)
    p_func = Func_Names{f_idx};
    
    for d_idx = 1 : length(Dist_List)
        dist = Dist_List(d_idx);
        current_step = current_step + 1;
        
        fprintf('[%d/%d] 仿真中: Mode=%-9s | Dist=%3d m ... ', ...
            current_step, total_steps, p_func, dist);
        
        % === 调用核心仿真函数 (返回结构体) ===
        res = run_simulation_core(dist, p_func, N_packets);
        
        % 记录数据
        Result_Loss(f_idx, d_idx)  = res.path_loss_db;
        Result_Delay(f_idx, d_idx) = res.mean_delay;
        
        fprintf('Loss=%.2f dB, Delay=%.2f ns\n', ...
            res.path_loss_db, res.mean_delay*1e9);
    end
    fprintf('------------------------------------------------\n');
end

toc(t_start_all);

%% ================= 3. 绘图结果 =================
colors = lines(length(Func_Names));
markers = {'o-', 's-', '^-', 'd-'};
lw = 1.5; % 线宽
ms = 6;   % 标记大小

% --- 图 1: 路径衰减 vs 距离 ---
figure('Name', 'Path Loss & Delay Analysis', 'Color', 'w', 'Position', [100, 100, 1000, 450]);

subplot(1, 2, 1);
for f_idx = 1 : length(Func_Names)
    plot(Dist_List, Result_Loss(f_idx, :), markers{f_idx}, ...
        'Color', colors(f_idx, :), 'LineWidth', lw, 'MarkerSize', ms, ...
        'DisplayName', Func_Names{f_idx});
    hold on;
end
grid on;
xlabel('Transmission Distance (m)', 'FontWeight', 'bold');
ylabel('Path Loss (dB)', 'FontWeight', 'bold');
title('Fig 1: Path Loss vs. Distance');
legend('Location', 'SouthEast', 'FontSize', 9);
xlim([min(Dist_List)-5, max(Dist_List)+5]);

% --- 图 2: 平均时延 vs 距离 ---
subplot(1, 2, 2);
for f_idx = 1 : length(Func_Names)
    % 将时延转换为纳秒 (ns) 显示，或者保留秒(s)
    delay_ns = Result_Delay(f_idx, :) * 1e9; 
    
    plot(Dist_List, delay_ns, markers{f_idx}, ...
        'Color', colors(f_idx, :), 'LineWidth', lw, 'MarkerSize', ms, ...
        'DisplayName', Func_Names{f_idx});
    hold on;
end

% 添加理论直射时间参考线 (t = L / c_water)
c_water = 2.25e8;
ref_time_ns = (Dist_List / c_water) * 1e9;
plot(Dist_List, ref_time_ns, 'k--', 'LineWidth', 1.0, 'DisplayName', 'Ballistic Limit (L/c)');

grid on;
xlabel('Transmission Distance (m)', 'FontWeight', 'bold');
ylabel('Mean Propagation Delay (ns)', 'FontWeight', 'bold');
title('Fig 2: Mean Delay vs. Distance');
legend('Location', 'NorthWest', 'FontSize', 9);
xlim([min(Dist_List)-5, max(Dist_List)+5]);

%% ================= 4. 核心仿真函数 (封装) =================
function res = run_simulation_core(link_dist, phase_func_mode, N_packets)
    % 结果结构体初始化
    res = struct('path_loss_db', 0, 'mean_delay', 0);

    % --- 1. 基础参数 ---
    param = struct();
    param.phase_func = phase_func_mode;  
    param.n_max = 8;  % 最大散射阶数
    
    lambda_nm = 514; 
    lambda = lambda_nm * 1e-9;
    k_wave = 2*pi/lambda;
    w0 = 0.01;                  
    
    Tx_Pos = [0, 0, 0];         
    Rx_Pos = [0, link_dist, 0]; % 动态接收机位置
    Link_Vec = Rx_Pos - Tx_Pos;
    Link_Dir = Link_Vec / link_dist;
    
    Rx_Normal = -Link_Dir;      
    Rx_Aperture = 0.2;          
    Rx_FOV = 30 * pi/180;       
    Rx_Area = pi*(Rx_Aperture/2)^2;
    
    % --- 2. 介质参数 (Petzold Clear Ocean 统一值) ---
    param.c_water = 2.25e8;     
    param.coef_c = 0.1514;                
    param.coef_a = 0.114;                 
    param.coef_b = param.coef_c - param.coef_a; 
    param.albedo = param.coef_b / param.coef_c; 
    
    % --- 3. 相函数参数 (统一自 scatter_test3.m) ---
    if strcmp(param.phase_func, 'HG')
        param.g_HG = 0.9707; 
    elseif strcmp(param.phase_func, 'Mix')
        w_mie = 0.8440; 
        param.coef_kr = param.coef_b * (1 - w_mie);
        param.coef_km = param.coef_b * w_mie;
        param.gamma = 0.017; 
        param.g_mie = 0.9814; 
        param.f_mie = 13.0101;      
    elseif strcmp(param.phase_func, 'TTHG')
        param.alpha_TTHG = 0.4437;
        param.g1_TTHG = 0.9900;
        param.g2_TTHG = 0.8232;
    elseif strcmp(param.phase_func, 'Empirical')
        param.c_e = param.coef_c; 
        param.a_e = param.coef_a;
        b_e = param.coef_b; albedo_e = param.albedo;
        param.q_e = 2.598 + 17.748*sqrt(b_e) - 16.722*b_e + 5.932*b_e*sqrt(b_e);
        param.k1 = 1.188 - 0.688*albedo_e;
        param.k2 = 0.1 * (3.07 - 1.90*albedo_e);
        param.k3 = 0.01 * (4.58 - 3.02*albedo_e);
        param.k4 = 0.001 * (3.24 - 2.25*albedo_e);
        param.k5 = 0.0001 * (0.84 - 0.61*albedo_e);
        
        % 积分计算归一化因子
        th_test = linspace(0, pi, 1000);
        val_test = zeros(size(th_test));
        for i=1:length(th_test)
            t_deg = th_test(i) * 180 / pi; if t_deg<1e-6, t_deg=1e-6; end
            term = 1 + (-1)^1*param.k1*t_deg^0.5 + (-1)^2*param.k2*t_deg^1.0 + ...
                   (-1)^3*param.k3*t_deg^1.5 + (-1)^4*param.k4*t_deg^2.0 + (-1)^5*param.k5*t_deg^2.5;
            val_test(i) = exp(param.q_e * term);
        end
        param.b_emp_norm = 2 * pi * trapz(th_test, val_test .* sin(th_test));
    end
    
    % --- 4. 湍流与相位屏 ---
    % 简化设置以加快扫描速度
    T_avg=20; S_avg=35; epsilon=1e-6; chi_T=1e-8; eta=1e-3; H_ratio=-20;
    N_screens = 10; 
    D_screen = 2.0; N_grid = 128; 
    delta_z_screen = link_dist / N_screens;
    
    [Phi_func, ~] = get_OTOPS_spectrum_handle(lambda, T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
    dx = D_screen / N_grid;
    x_axis = (-N_grid/2 : N_grid/2-1) * dx;
    
    Screen_Chain = repmat(struct('Center',[],'Normal',[],'grad_x',[],'grad_y',[],'u_vec',[],'v_vec',[]), 1, N_screens);
    if abs(Link_Dir(3)) < 0.9, up_temp = [0, 0, 1]; else, up_temp = [1, 0, 0]; end
    u_vec = cross(up_temp, Link_Dir); u_vec = u_vec / norm(u_vec);
    v_vec = cross(Link_Dir, u_vec);   v_vec = v_vec / norm(v_vec);
    
    for i = 1:N_screens
        phi = gen_screen_from_spectrum(Phi_func, D_screen, N_grid, k_wave, delta_z_screen);
        [gx, gy] = gradient(phi, dx);
        Screen_Chain(i).Center = Tx_Pos + i * delta_z_screen * Link_Dir;
        Screen_Chain(i).Normal = Link_Dir;
        Screen_Chain(i).grad_x = gx; Screen_Chain(i).grad_y = gy;
        Screen_Chain(i).u_vec = u_vec; Screen_Chain(i).v_vec = v_vec;
    end
    
    % --- 5. MC 追踪 ---
    sum_power = 0;
    sum_power_time = 0; % 用于计算加权平均时间 (Power * Time)
    
    for p = 1:N_packets
        % 发射
        r0 = w0 * sqrt(-0.5*log(rand())); phi0 = 2*pi*rand();
        pos_local = r0*cos(phi0)*u_vec + r0*sin(phi0)*v_vec;
        pos_init = Tx_Pos + pos_local;
        
        theta_div = lambda / (pi * w0);
        U = theta_div * sqrt(-0.5 * log(rand())); psi_ini = 2*pi*rand;
        dir_init = rotate_direction(Link_Dir, U, psi_ini);
        
        % 分支A: 直射 (Ballistic)
        [~, ~, hit_A, len_A] = ray_march_generic(pos_init, dir_init, 1e9, ...
            Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, true, Screen_Chain, k_wave, x_axis, dx, N_grid);
        
        if hit_A
            e_val = exp(-param.coef_c * len_A);
            t_val = len_A / param.c_water;
            
            sum_power = sum_power + e_val;
            sum_power_time = sum_power_time + e_val * t_val;
        end
        
        % 分支B: 散射 (Scattered)
        pos = pos_init; dir = dir_init; weight = 1.0;
        current_len = 0; % 累计飞行距离
        
        d_step = -log(rand()) / param.coef_c;
        [pos_new, dir_new, ~, step_len] = ray_march_generic(pos, dir, d_step, ...
            Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, Screen_Chain, k_wave, x_axis, dx, N_grid);
        
        current_len = current_len + step_len;
        if dot(pos_new - Tx_Pos, Link_Dir) >= link_dist, continue; end
        
        pos = pos_new; dir = dir_new;
        
        for order = 1 : param.n_max
            vec_to_rx = Rx_Pos - pos; dist_to_rx = norm(vec_to_rx); dir_to_rx = vec_to_rx/dist_to_rx;
            
            % PIS (Point Intersection Estimation)
            if acos(dot(-dir_to_rx, Rx_Normal)) <= Rx_FOV
                cos_ts = dot(dir, dir_to_rx);
                % 相函数计算
                if strcmp(param.phase_func,'HG'), p_val=pdf_HG(cos_ts,param.g_HG);
                elseif strcmp(param.phase_func,'TTHG'), p_val=pdf_TTHG(cos_ts,param.alpha_TTHG,param.g1_TTHG,param.g2_TTHG);
                elseif strcmp(param.phase_func,'Mix')
                    p_ray=pdf_Rayleigh(cos_ts,param.gamma); p_mie=pdf_Mie(cos_ts,param.g_mie,param.f_mie);
                    p_val=(param.coef_kr*p_ray+param.coef_km*p_mie)/param.coef_b;
                elseif strcmp(param.phase_func,'Empirical'), p_val=pdf_Empirical(cos_ts, param);
                end
                
                omega = Rx_Area/(dist_to_rx^2) * abs(dot(dir_to_rx, Rx_Normal));
                E_pis = weight * param.albedo * min(1, p_val*omega) * exp(-param.coef_c * dist_to_rx);
                
                % 累计能量与时间矩
                t_pis = (current_len + dist_to_rx) / param.c_water;
                sum_power = sum_power + E_pis;
                sum_power_time = sum_power_time + E_pis * t_pis;
            end
            
            % 散射游走
            weight = weight * param.albedo;
            if weight < 1e-9, break; end
            
            % 采样新方向
            theta_i = pi*rand; phi_i = 2*pi*rand; cos_ti = cos(theta_i);
            
            if strcmp(param.phase_func,'HG'), p_samp=pdf_HG(cos_ti,param.g_HG);
            elseif strcmp(param.phase_func,'TTHG'), p_samp=pdf_TTHG(cos_ti,param.alpha_TTHG,param.g1_TTHG,param.g2_TTHG);
            elseif strcmp(param.phase_func,'Mix')
                 p_ray=pdf_Rayleigh(cos_ti,param.gamma); p_mie=pdf_Mie(cos_ti,param.g_mie,param.f_mie);
                 p_samp=(param.coef_kr*p_ray+param.coef_km*p_mie)/param.coef_b;
            elseif strcmp(param.phase_func,'Empirical'), p_samp=pdf_Empirical(cos_ti, param);
            end
            
            weight = weight * (2 * pi^2 * p_samp * sin(theta_i));
            dir = rotate_direction(dir, theta_i, phi_i);
            d_next = -log(rand())/param.coef_c;
            
            [pos, dir, ~, step_len] = ray_march_generic(pos, dir, d_next, ...
                Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, false, Screen_Chain, k_wave, x_axis, dx, N_grid);
            
            current_len = current_len + step_len;
        end
    end
    
    % --- 6. 统计计算 ---
    if sum_power > 0
        P_rx_norm = sum_power / N_packets;
        res.path_loss_db = 10 * log10(P_rx_norm + 1e-30);
        res.mean_delay = sum_power_time / sum_power; % 加权平均时间
    else
        res.path_loss_db = -300; % 无信号
        res.mean_delay = link_dist / param.c_water; % 默认给直射时间
    end
end

%% ================= 5. 辅助通用函数 =================
function [pos, dir, hit_flag, total_len] = ray_march_generic(pos, dir, dist_limit, ...
    Rx_Pos, Rx_Aperture, Rx_FOV, Rx_Normal, enable_hit_check, ...
    Screen_Chain, k_wave, x_axis, dx, N_grid)
    hit_flag = false; total_len = 0; remaining_dist = dist_limit;
    while remaining_dist > 1e-6
        min_dist = remaining_dist; event_type = 'none'; target_idx = -1;
        for i = 1:length(Screen_Chain)
            denom = dot(dir, Screen_Chain(i).Normal);
            if abs(denom)>1e-6
                t = dot(Screen_Chain(i).Center - pos, Screen_Chain(i).Normal)/denom;
                if t>1e-6 && t<min_dist, min_dist=t; event_type='screen'; target_idx=i; end
            end
        end
        if enable_hit_check
            denom_rx = dot(dir, Rx_Normal);
            if abs(denom_rx)>1e-6
                t_rx = dot(Rx_Pos - pos, Rx_Normal)/denom_rx;
                if t_rx>1e-6 && t_rx<=min_dist, min_dist=t_rx; event_type='rx'; end
            end
        end
        pos = pos + dir*min_dist; remaining_dist = remaining_dist - min_dist; total_len = total_len + min_dist;
        if strcmp(event_type, 'rx')
            if norm(pos-Rx_Pos)<=Rx_Aperture/2 && acos(dot(-dir, Rx_Normal))<=Rx_FOV, hit_flag=true; return; end
        elseif strcmp(event_type, 'screen')
            scr=Screen_Chain(target_idx); vec=pos-scr.Center;
            lx=dot(vec,scr.u_vec); ly=dot(vec,scr.v_vec);
            ix=mod(round((lx-x_axis(1))/dx), N_grid)+1; iy=mod(round((ly-x_axis(1))/dx), N_grid)+1;
            d_vec = (scr.grad_x(iy,ix)*scr.u_vec + scr.grad_y(iy,ix)*scr.v_vec)/k_wave;
            dir = dir + d_vec; dir = dir/norm(dir);
        end
    end
end

function new_dir = rotate_direction(dir, theta_s, psi_s)
    mu_x=dir(1); mu_y=dir(2); mu_z=dir(3); denom=sqrt(1-mu_z^2);
    if denom<1e-10
        if mu_z>0, new_dir=[sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), cos(theta_s)];
        else, new_dir=[sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), -cos(theta_s)]; end
    else
        st=sin(theta_s); ct=cos(theta_s); cp=cos(psi_s); sp=sin(psi_s);
        nx=st/denom*(mu_x*mu_z*cp-mu_y*sp)+mu_x*ct;
        ny=st/denom*(mu_y*mu_z*cp+mu_x*sp)+mu_y*ct;
        nz=-st*cp*denom+mu_z*ct;
        new_dir=[nx,ny,nz];
    end
    new_dir=new_dir/norm(new_dir);
end

% --- 相函数 PDF 定义 ---
function p = pdf_HG(c, g), p = (1-g^2)./(4*pi*(1+g^2-2*g*c).^1.5); end
function p = pdf_Rayleigh(c, gam), p = 3*(1+3*gam+(1-gam)*c.^2)/(16*pi*(1+2*gam)); end
function p = pdf_Mie(c, g, f), val=(1-g^2)/2*(1./(1+g^2-2*g*c).^1.5+f*0.5*(3*c.^2-1)/(1+g^2)^1.5); p=val/(2*pi); end
function p = pdf_TTHG(c, a, g1, g2), p = a*pdf_HG(c,g1) + (1-a)*pdf_HG(c,g2); end
function p = pdf_Empirical(c, pm)
    c=max(min(c,1),-1); t_deg=acos(c)*180/pi; if t_deg<1e-6,t_deg=1e-6;end
    term=1+(-1)*pm.k1*t_deg^0.5+pm.k2*t_deg-pm.k3*t_deg^1.5+pm.k4*t_deg^2-pm.k5*t_deg^2.5;
    p = exp(pm.q_e*term)/pm.b_emp_norm;
end

% --- 湍流谱 ---
function [Phi_n_func, k_wave] = get_OTOPS_spectrum_handle(lambda, T, S, epsilon, chi_T, eta, H_ratio)
    lambda_nm = lambda * 1e9; k_wave = 2 * pi / lambda; 
    a1=1.779e-4; a2=-1.05e-6; a3=1.6e-8; a4=-2.02e-6; a5=1.155e-2; a6=-4.23e-3;
    A=a2*S+2*a3*T*S+2*a4*T+a6/lambda_nm; B=a1+a2*T+a3*T^2+a5/lambda_nm; 
    Phi_Hill=@(K) 0.033*chi_T*epsilon^(-1/3)*(K.^2+1e-25).^(-11/6); % 简化用于快速绘图
    Phi_n_func=@(K) (A+B)^2 * Phi_Hill(K); 
end
function phase_screen = gen_screen_from_spectrum(Phi_n_func, D, N, k_wave, delta_z)
    dk=2*pi/D; kx=(-N/2:N/2-1)*dk; [KX,KY]=meshgrid(kx,kx); K_grid=sqrt(KX.^2+KY.^2); 
    K_grid(N/2+1,N/2+1)=1e-10; Phi_val=Phi_n_func(K_grid); Phi_val(N/2+1,N/2+1)=0;
    F_phi=2*pi^2*k_wave^2*delta_z*Phi_val; 
    phase_screen=real(ifft2(ifftshift((randn(N)+1i*randn(N)).*sqrt(F_phi)*dk)))*N^2;
end