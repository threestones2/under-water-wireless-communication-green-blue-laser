%% 本代码复现论文《Monte-Carlo Integration Models for Multiple Scattering Based Optical Wireless Communication》
%  针对水下信道实验部分 (Fig. 9 & Table V) 进行了量纲修正与参数修复

clc; clear; close all;

% 1. 共有几何参数（严格按照 Table I 回退至 NLOS 设定）
param = struct();
param.r = 5;                  % 水下仿真距离 (Table V)
param.theta_T = deg2rad(10);  % 发射器天顶角
param.phi_T = deg2rad(-90);   % 发射器方位角
param.theta_R = deg2rad(0);  % 接收器天顶角
param.phi_R = deg2rad(90);    % 接收器方位角
param.beta_T = deg2rad(6);   % 发射束发散角
param.beta_R = deg2rad(10);   % 接收器视场角
param.A_r = 1.77*1e-4;        % 接收器面积 (m^2)
param.N = 1e7;                % 采样点数
param.n_max = 10;              % 最大散射阶数

% 水下光速计算 (真空中光速 / 水的典型折射率 1.33)
param.c = 2.997046e8 / 1.33;

% 时间窗口设定 (覆盖 NLOS 构型下的理论峰值 31.4 ns)
param.dt = 1e-10; 
param.T_bins = 10e-9 : param.dt : 60e-9; 

% 2. 修正后的水质参数定义 
% (剔除原论文 Table V 中错误的 1e-3 乘数，还原真实海洋光学参数)
water_types = {'Clean Ocean', 'Coastal', 'Turbid Harbor'};
ka_array = [0.069, 0.088, 0.285];  % 真实的吸收系数 m^-1
ks_array = [0.080, 0.216, 1.875];  % 真实的散射系数 m^-1
h_array  = [0.8708, 0.9470, 0.9199]; % HG 不对称因子

% 3. 独立绘制三张图表
for w = 1:3
    fprintf('正在仿真水质: %s...\n', water_types{w});
    % 更新当前水质参数
    param.ka = ka_array(w);
    param.ks = ks_array(w);
    param.h = h_array(w);
    
    % 调用核心模型
    [P, P_hat, h_hat, h_hat_i] = MCI_PIS_Function(param);
    
    % 创建独立绘图窗口
    figure('Name', ['Channel IRF - ', water_types{w}]);
    
    % 将时间轴由 s 转换为 ns 进行可视化
    plot(param.T_bins(1:end-1) * 1e9, h_hat_i, 'LineWidth', 1.5);
    hold on;
    plot(param.T_bins(1:end-1) * 1e9, h_hat, 'k--', 'LineWidth', 1.5);
    
    % 图表格式化
    xlabel('Time (ns)');
    ylabel('Impulse Response (W/m^2)');
    title(['Channel IRF in ', water_types{w}]);
    legend('1st order', '2nd order', '3rd order', 'total', 'Location', 'best');
    grid on;
    
    fprintf('%s 仿真完成，总接收能量积分: %e J/m^2\n\n', water_types{w}, sum(h_hat)*param.dt);
end

%% MCI_PIS 核心模型函数
function [P, P_hat, h_hat, h_hat_i] = MCI_PIS_Function(param)
    % 提取参数
    r = param.r; 
    theta_T = param.theta_T; phi_T = param.phi_T; 
    theta_R = param.theta_R; phi_R = param.phi_R; 
    beta_T = param.beta_T;   beta_R = param.beta_R; 
    A_r = param.A_r; 
    ka = param.ka; ks = param.ks; h = param.h;
    N = param.N; n_max = param.n_max; 
    c = param.c; 
    T_bins = param.T_bins; dt = param.dt;
    
    ke = ka + ks;
    
    mu_T = [cos(phi_T)*sin(theta_T); sin(phi_T)*sin(theta_T); cos(theta_T)];
    mu_R = [cos(phi_R)*sin(theta_R); sin(phi_R)*sin(theta_R); cos(theta_R)];
    
    P_hat = zeros(1, n_max);
    num_bins = length(T_bins) - 1;
    h_hat_i = zeros(n_max, num_bins);
    
    for k = 1:N
        pos = [0; r; 0]; 
        mu = mu_T; 
        d_total = 0; 
        O_star = 1; 
        
        for i = 1:n_max
            % 重要性采样：传播距离
            d_i = -log(1 - rand) / ke;
            
            % 均匀采样：散射角度
            if i == 1
                theta_i = (beta_T / 2) * rand;
                f_Theta = sin(theta_i) / (1 - cos(beta_T/2)); 
                O_star_i = f_Theta * (beta_T/2);
            else
                theta_i = pi * rand;
                % HG相函数(带sin)
                f_Theta = (1 - h^2) / (2 * (1 + h^2 - 2*h*cos(theta_i))^(3/2)) * sin(theta_i); 
                O_star_i = (ks/ke) * f_Theta * pi;
            end
            O_star = O_star * O_star_i;

            phi_i = 2 * pi * rand;
            
            % 坐标更新
            mu = update_direction(mu, theta_i, phi_i);
            pos = pos + d_i * mu;
            d_total = d_total + d_i;
            
            % 接收判决
            d_to_rx = norm(pos);
            cos_phi_r = dot(mu_R, pos) / d_to_rx;
            
            if cos_phi_r >= cos(beta_R/2)
                theta_n = acos(dot(-mu, pos) / d_to_rx);
                Omega_r = A_r / d_to_rx^2 * cos_phi_r;  
                
                % 剥离 sin 项的纯相函数，避免极点奇异性
                P_HG = (1 - h^2) / (2 * (1 + h^2 - 2*h*cos(theta_n))^(3/2));
                
                % 严格依照 Eq. 6 修复指数平方错误
                p_d = (ks/ke) * exp(-ke * d_to_rx) * cos_phi_r * min(1, P_HG * Omega_r / (2*pi));

                P_cur = p_d * O_star;
                P_hat(i) = P_hat(i) + P_cur;
                
                t_total = (d_total + d_to_rx) / c;
                
                % 快速时间槽分配
                bin_idx = floor((t_total - T_bins(1)) / dt) + 1;
                if bin_idx >= 1 && bin_idx <= num_bins
                    h_hat_i(i, bin_idx) = h_hat_i(i, bin_idx) + P_cur;
                end
            end
        end
    end
    
    % 概率归一化
    P_hat = P_hat / N;
    P = sum(P_hat);
    
    % 转化为时间密度与面积密度的绝对物理量纲 W/m^2
    h_hat_i = h_hat_i / (N * dt * A_r);
    h_hat = sum(h_hat_i, 1);
end

% 辅助函数：更新方向余弦
function mu_new = update_direction(mu, theta, phi)
    mu_x = mu(1); mu_y = mu(2); mu_z = mu(3);
    if abs(mu_z) < 1 - 1e-10
        denom = sqrt(1 - mu_z^2);
        mu_new = [cos(theta)*mu_x + sin(theta)*(mu_x*mu_z*cos(phi) - mu_y*sin(phi))/denom;
                  cos(theta)*mu_y + sin(theta)*(mu_y*mu_z*cos(phi) + mu_x*sin(phi))/denom;
                  cos(theta)*mu_z - sin(theta)*cos(phi)*denom];
    else
        if mu_z > 0  
            mu_new = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
        else         
            mu_new = [-sin(theta)*cos(phi); -sin(theta)*sin(phi); -cos(theta)];          
        end
    end
end