function [P, P_hat, h_hat, h_hat_i] = MSI_PIS_Function(param)
    % 输入参数:
    % r: 收发器之间的距离 (m)
    % theta_T, phi_T: 发射器的天顶角和方位角 (rad)
    % beta_T: 发射光束发散角 (rad)
    % theta_R, phi_R: 接收器的天顶角和方位角 (rad)
    % beta_R: 接收器视场角 (rad)
    % A_r: 接收器面积 (m^2)
    % ka: 吸收系数 (1/m)
    % kr: Rayleigh散射系数 (1/m)
    % km: Mie散射系数 (1/m)
    % gamma, g, f: Rayleigh和Mie散射的模型参数
    % h: HG函数的模型参数（水下通信用）
    % N: 采样点数
    % n_max: 最大散射阶数
    % T_bins: 时间槽边界向量 (s)
    
    % 输出参数:
    % P: 总接收概率
    % P_hat: 各散射阶的接收概率 (1 x n_max)
    % h_hat: 总IRF (1 x length(T_bins)-1)
    % h_hat_i: 各散射阶的IRF (n_max x length(T_bins)-1)

    % 参数设置
    r = param.r; % 通信距离 (m)
    theta_T = param.theta_T; 
    phi_T = param.phi_T; % 发射器角度
    theta_R = param.theta_R; 
    phi_R = param.phi_R; % 接收器角度
    beta_T = param.beta_T; % 发射束发散角
    beta_R = param.beta_R; % 接收器视场角
    A_r = param.A_r; % 接收器面积 (m^2)
    ka = param.ka; % 吸收系数 (1/m)
    kr = param.kr; % Rayleigh散射系数 (1/m)
    km = param.km; % Mie散射系数 (1/m)
    gamma = param.gamma; % Rayleigh散射参数
    g = param.g; % Mie散射参数
    f = param.f; % Mie散射参数
    h = param.h; % HG函数参数（水下用，这里UV通信未使用）
    N = param.N; % 采样点数
    n_max = param.n_max; % 最大散射阶数
    c = param.c; % 光速 (m/s)
    t_max = param.t_max; % 最大接收时间（5倍传播时间）
    T_bins = param.T_bins; % 100个时间槽
    r_R = sqrt(A_r / pi); % 接收器半径 (m)
    
    % 计算散射系数和消光系数
    ks = kr + km;
    ke = ka + ks;

    
    % 计算方向余弦
    mu_T = [cos(phi_T)*sin(theta_T); sin(phi_T)*sin(theta_T); cos(theta_T)];
    mu_R = [cos(phi_R)*sin(theta_R); sin(phi_R)*sin(theta_R); cos(theta_R)];
    
    % 初始化
    P_hat = zeros(1, n_max);
    h_hat_i = zeros(n_max, length(T_bins));
    
    % 主循环：对每个采样点
    for k = 1:N
        % 初始化光子位置和方向
        pos = [0; r; 0]; % 初始位置（接收器在原点，发射器在(0,r,0)）
        mu = mu_T; % 初始方向
        d_total = 0; % 总传播距离
        O_star = 1; % 目标函数值
        
        % 对每个散射阶
        for i = 1:n_max
            % 采样传播距离（指数分布 - 重要性采样）
            d_i = -log(1 - rand) / ke;
            
            % 采样散射天顶角（均匀分布 - 部分重要性采样）
            if i == 1
                theta_i = (beta_T / 2) * rand;
            else
                theta_i = pi * rand;
            end
            
            % 采样散射方位角（均匀分布）
            phi_i = 2 * pi * rand;
            
            % 计算目标函数的乘积因子
            if i == 1
                O_star_i = pi * sin(theta_i) / (1 - cos(beta_T/2));
            else
                % 计算散射天顶角的PDF（Rayleigh+Mie）
                f_Theta = (kr/ks) * f_Theta_Rayleigh(theta_i, gamma) + (km/ks) * f_Theta_Mie(theta_i, g, f);
                O_star_i = (ks/ke) * pi * f_Theta;
            end
            O_star = O_star * O_star_i;
            
            % 更新光子方向
            mu = update_direction(mu, theta_i, phi_i);
            
            % 更新光子位置
            pos = pos + d_i * mu;
            d_total = d_total + d_i;
            
            % 计算到接收器的距离和角度
            d_to_rx = norm(pos);
            cos_phi_r = dot(mu_R, pos) / d_to_rx;
            
            % 检查是否在接收器视场内
            if cos_phi_r >= cos(beta_R/2)
                % 计算第n次散射的天顶角（相对于接收器方向）
                theta_n = acos(dot(-mu, pos) / d_to_rx);
                
                % 计算散射相函数
                f_Theta_n = (kr/ks) * f_Theta_Rayleigh(theta_n, gamma) + (km/ks) * f_Theta_Mie(theta_n, g, f);
                
                % 修改：使用正确的接收器立体角计算公式
                Omega_r = A_r / d_to_rx^2 * cos_phi_r;  % 使用准确的接收器立体角计算
                
                % 修改：使用正确的检测概率计算公式
                p_d = ks/ke * exp(-ke * d_to_rx) * f_Theta_n * Omega_r * cos_phi_r;
                
                % 计算当前目标函数值（包括检测概率）
                P_cur = p_d * O_star;
                
                % 更新接收概率
                P_hat(i) = P_hat(i) + P_cur;
                
                % 计算接收时间
                t_total = (d_total + d_to_rx) / c;
                
                % 找到对应的时间槽
                bin_idx = find(t_total >= T_bins(1:end-1) & t_total < T_bins(2:end), 1);
                if ~isempty(bin_idx)
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
    
    % 归一化IRF（除以时间间隔和接收面积）
    delta_T = T_bins(2)-T_bins(1);
    h_hat = h_hat / (delta_T);
    for i = 1:n_max
        h_hat_i(i, :) = h_hat_i(i, :) / (delta_T);
    end
end