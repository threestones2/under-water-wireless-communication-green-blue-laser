%% 本代码复现论文《Monte-Carlo Integration Models for Multiple  Scattering Based Optical Wireless Communication》
%  修改版：MCI_PIS_4
%  修改内容：将 i>1 阶散射的采样方式，从“均匀采样”改为 MCS 中的“反变换法重要性采样”,收敛速度慢的一批
%  这相当于实现了 MCI-IS (Importance Sampling) 的核心逻辑

function [P, P_hat, h_hat, h_hat_i] = MSI_PIS_Function(param)
    % 输入参数:
    % (参数列表保持不变，请参考原文件)
    
    % 参数设置
    r = param.r; 
    theta_T = param.theta_T; phi_T = param.phi_T; 
    theta_R = param.theta_R; phi_R = param.phi_R; 
    beta_T = param.beta_T; beta_R = param.beta_R; 
    A_r = param.A_r; 
    ka = param.ka; kr = param.kr; km = param.km; 
    gamma = param.gamma; g = param.g; f = param.f; 
    N = param.N; n_max = param.n_max; 
    c = param.c; t_max = param.t_max; T_bins = param.T_bins; 
    r_R = sqrt(A_r / pi); 
    
    % 计算散射系数和消光系数
    ks = kr + km;
    ke = ka + ks;
    
    % 预计算散射权重，用于混合相位函数
    w_ray = kr / ks;
    w_mie = km / ks;

    % 计算方向余弦
    mu_T = [cos(phi_T)*sin(theta_T); sin(phi_T)*sin(theta_T); cos(theta_T)];
    mu_R = [cos(phi_R)*sin(theta_R); sin(phi_R)*sin(theta_R); cos(theta_R)];
    
    % 初始化
    P_hat = zeros(1, n_max);
    h_hat_i = zeros(n_max, length(T_bins));
    
    % 主循环：对每个采样点
    for k = 1:N
        % 初始化光子位置和方向
        pos = [0; r; 0]; 
        mu = mu_T; 

        d_total = 0; 
        O_star = 1; % 目标函数值
        
        % 对每个散射阶
        for i = 1:n_max
            % 1. 采样传播距离（保持指数分布 - 重要性采样）
            d_i = -log(1 - rand) / ke;
            
            % 2. 采样散射天顶角 (theta_i) 与 计算权重 (O_star_i)
            if i == 1
                % --- 第1阶：保持原有的光锥均匀采样 ---
                % 因为第1阶是由光源几何决定的，不是由散射相函数决定的
                theta_i = (beta_T / 2) * rand;
                
                % 权重计算 (保持不变)
                f_Theta = sin(theta_i) / (1 - cos(beta_T/2)); 
                O_star_i = f_Theta * (beta_T/2);
            else
                % --- 第>1阶：使用 MCS 中的反变换法求散射角 (重要性采样) ---
                
                % 生成随机数
                xi_mu = rand;
                
                % 调用数值求解函数 (从 MCS.m 移植)
                % 注意：这里返回的是 cos(theta)
                mu_s = solve_phase_function(xi_mu, w_ray, w_mie, gamma, g, f);
                theta_i = acos(mu_s);
                
                % 权重计算 (核心修改)
                % 原理：因为我们现在的采样分布 q(theta) 完全等于物理相函数分布 p(theta)
                % 所以权重因子 p(theta)/q(theta) = 1。
                % 只需要保留能量衰减因子 (Albedo)
                O_star_i = (ks/ke) * 1.0; 
            end
            
            % 采样散射方位角（均匀分布 0-2pi）
            phi_i = 2 * pi * rand;
            
            % 累乘权重
            O_star = O_star * O_star_i;

            % 更新光子方向
            mu = update_direction(mu, theta_i, phi_i);
            
            % 更新光子位置
            pos = pos + d_i * mu;
            d_total = d_total + d_i;
            
            % 3. 接收检测 (Detection)
            % 这部分逻辑不变，因为我们需要计算光子“强行”折向接收器的概率
            
            % 计算到接收器的几何参数
            r_length = r_R / tan(beta_R);
            pos_R = r_length * (-mu_R); % 接收圆锥顶点
            d_to_rx = norm(pos - pos_R);
            d_to_rx_real = norm(pos); % 到原点(接收面中心)的距离
            
            cos_phi_r = dot(mu_R, pos - pos_R) / d_to_rx;
            
            % 检查是否在接收器视场内
            if cos_phi_r >= cos(beta_R/2)
                % 计算强行散射向接收器的角度 theta_n
                % 注意：这里的 theta_n 是光子当前方向 mu 与 "光子指向接收器方向" 的夹角
                vec_to_rx = -pos; % 指向原点
                mu_n_target = vec_to_rx / norm(vec_to_rx);
                cos_theta_n = dot(mu, mu_n_target);
                
                % 限制精度误差
                cos_theta_n = max(min(cos_theta_n, 1), -1);
                theta_n = acos(cos_theta_n);
                
                % 计算该角度下的相位函数值 f_Theta_n
                % 注意：这里必须显式计算，因为最后一步是“强行检测”，不是随机采样
                f_Theta_n = w_ray * f_Theta_Rayleigh(theta_n, gamma) + ...
                            w_mie * f_Theta_Mie(theta_n, g, f);
                
                % 计算接收立体角
                Omega_r = A_r / (d_to_rx_real^2) * cos_phi_r; % 近似公式
                
                % 计算接收概率 (Detection Probability)
                % 公式：p_d = Albedo * T(d) * P(theta) * Omega
                % 注意：这里需要除以 2pi 或 4pi 归一化因子，取决于 f_Theta 的定义
                % 原代码中的 f_Theta_Rayleigh/Mie 已经归一化积分 = 1 (对d_mu积分)
                % 对应的立体角积分概率是 f_Theta * Omega / (2*pi) ? 
                % 验证：P(theta) dOmega = P(cos_theta) dOmega / (2pi) 如果 P 对 mu 归一化为 1
                
                p_d = (ks/ke) * exp(-ke * d_to_rx_real) * min(1, f_Theta_n * Omega_r / (2*pi));

                % 计算最终贡献
                P_cur = p_d * O_star;
                
                % 累加结果
                P_hat(i) = P_hat(i) + P_cur;
                
                % 记录时间响应
                t_total = (d_total + d_to_rx_real) / c;
                bin_idx = floor(t_total / (T_bins(2)-T_bins(1))) + 1;
                % 注意：这里的bin_idx逻辑我微调了一下以匹配通用写法
                if bin_idx <= length(T_bins)-1 && bin_idx >= 1
                     h_hat_i(i, bin_idx) = h_hat_i(i, bin_idx) + P_cur;
                end
            end
        end
    end
    
    % 归一化
    P_hat = P_hat / N;
    P = sum(P_hat);
    h_hat_i = h_hat_i / N;
    h_hat = sum(h_hat_i, 1);
    
    % 归一化IRF
    delta_T = T_bins(2)-T_bins(1);
    h_hat = h_hat / delta_T;
    for i = 1:n_max
        h_hat_i(i, :) = h_hat_i(i, :) / delta_T;
    end
end

%% ================== 以下为辅助函数 ==================

% 1. 更新方向函数 (保持不变)
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

% 2. 数值求解相位函数采样 (从 MCS.m 移植并优化)
% 返回的是 scattering cosine (mu), 不是 angle
function mu_val = solve_phase_function(xi_mu, w_ray, w_mie, gamma, g, f)
    % 定义相位函数 P(mu)，这里 P(mu) 对 mu 在 [-1,1] 积分为 1
    % Rayleigh 部分
    P_ray = @(mu) (3*(1+3*gamma+(1-gamma)*mu.^2)) / (8*(1+2*gamma)); % 注意分母是 8 不是 16pi，因为去掉了方位角积分
    % Mie 部分 (Double Henyey-Greenstein)
    P_mie = @(mu) (1-g^2)/2 * (1./((1+g^2-2*g*mu).^(3/2)) + f*(0.5*(3*mu.^2-1))./((1+g^2).^(3/2)));
    
    % 混合相位函数
    P_total = @(mu) w_ray * P_ray(mu) + w_mie * P_mie(mu);
    
    % 使用 fzero 求解累积分布函数 (CDF) 的逆
    % equation: CDF(mu) - xi = 0
    % CDF(mu) = integral_{-1}^{mu} P_total(t) dt
    
    % 优化：为了提高速度，可以将 integral 替换为解析积分公式（如果可能），
    % 或者对于 fzero 保持现状。此处保持 MCS.m 的逻辑，但注意积分常数。
    
    % 定义目标函数
    target_fun = @(m) integral(P_total, -1, m) - xi_mu;
    
    % 求解
    % 注意：fzero 可能比较慢。在高性能MCI中，通常会预计算CDF表然后插值 lookup。
    % 这里为了保持与 MCS 一致性，直接解。
    options = optimset('Display','off', 'TolX', 1e-4); 
    try
        mu_val = fzero(target_fun, [-1 1], options);
    catch
        % 如果解算失败（极少数情况），返回前向散射
        mu_val = 1;
    end
    
    % 限制范围
    mu_val = max(min(mu_val, 1), -1);
end

% 3. 辅助函数：Rayleigh PDF 值计算 (用于检测概率)
function val = f_Theta_Rayleigh(theta, gamma)
    val = 3 * (1 + 3*gamma + (1-gamma)*cos(theta)^2)/ (8 * (1 + 2*gamma));
end

% 4. 辅助函数：Mie PDF 值计算 (用于检测概率)
function val = f_Theta_Mie(theta, g, f)
    val = (1 - g^2)/2 * (1/(1 + g^2 - 2*g*cos(theta))^(3/2) + f * 0.5 * (3*cos(theta)^2 - 1) / (1 + g^2)^(3/2));
end

%% 调用程序
% 参数设置（根据论文中的Table I）
param=struct();
param.r = 100; % 通信距离 (m)

param.theta_T = deg2rad(45); 
param.phi_T = 2*pi-deg2rad(90); % 发射器角度

param.theta_R = deg2rad(45); 
param.phi_R = deg2rad(90); % 接收器角度

param.beta_T = deg2rad(17); % 发射束发散角
param.beta_R = deg2rad(30); % 接收器视场角

param.A_r = 1.77*1e-4; % 接收器面积 (m^2)
param.ka = 0.802e-3; % 吸收系数 (1/m)
param.kr = 0.266e-3; % Rayleigh散射系数 (1/m)
param.km = 0.284e-3; % Mie散射系数 (1/m)
param.gamma = 0.017; % Rayleigh散射参数
param.g = 0.72; % Mie散射参数
param.f = 0.5; % Mie散射参数
param.h = 0.9; % HG函数参数（水下用，这里UV通信未使用）
param.N = 1e5; % 采样点数
param.n_max = 3; % 最大散射阶数

% 设置时间槽（根据接收时间范围调整）
param.c = 2.997046e8;
param.t_max = 2e-6; % 最大模拟时间 (s) (1 us)
param.dt = 5e-9; % 时间分辨率 (s) (1 ns)
param.T_bins = 0:param.dt:param.t_max;

% 运行MSI_PIS模型
[P, P_hat, h_hat, h_hat_i] = MSI_PIS_Function(param);

% 计算路径损耗 (dB)
PL = 10 * log10(1 / P);
PL1= 10*log10( 1 / trapz(param.T_bins, h_hat_i(1,:)) );
% 显示结果
fprintf('总路径损失: %.2f dB\n',PL);
fprintf('1阶散射路径损失: %.2f dB\n',PL1);
% 绘制IRF
figure;
plot(param.T_bins(1:end), h_hat_i(1:3,:),param.T_bins(1:end),h_hat);
xlabel('Time (s)');
ylabel('Impulse Response (W/m^2)');
title('Channel Impulse Response (IRF)');
legend('1st order', '2nd order', '3rd order','total');
grid on;

