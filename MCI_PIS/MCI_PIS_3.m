function [P, P_hat, h_hat, h_hat_i] = MCI_PIS_Function_Modified(param)
    % 输入参数:
    % param.model_type: 'UV' (混合Rayleigh/Mie) 或 'HG' (Henyey-Greenstein)
    % 其他参数同原代码...
    
    % 参数设置
    r = param.r; 
    theta_T = param.theta_T; 
    phi_T = param.phi_T; 
    theta_R = param.theta_R; 
    phi_R = param.phi_R; 
    beta_T = param.beta_T; 
    beta_R = param.beta_R; 
    A_r = param.A_r; 
    ka = param.ka; 
    kr = param.kr; 
    km = param.km; 
    gamma = param.gamma; 
    g = param.g; 
    f = param.f; 
    h_HG = param.h; % HG函数的非对称因子 g
    N = param.N; 
    n_max = param.n_max; 
    c = param.c; 
    T_bins = param.T_bins; 
    
    % 默认模型类型为 UV (兼容旧代码)
    if ~isfield(param, 'model_type')
        model_type = 'UV';
    else
        model_type = param.model_type;
    end

    % 计算散射系数和消光系数
    ks = kr + km;
    ke = ka + ks;
    
    % 计算方向余弦
    mu_T = [cos(phi_T)*sin(theta_T); sin(phi_T)*sin(theta_T); cos(theta_T)];
    mu_R = [cos(phi_R)*sin(theta_R); sin(phi_R)*sin(theta_R); cos(theta_R)];
    
    % 初始化
    P_hat = zeros(1, n_max);
    h_hat_i = zeros(n_max, length(T_bins));
    
    % 主循环
    for k = 1:N
        pos = [0; r; 0]; 
        mu = mu_T; 
        d_total = 0; 
        O_star = 1; 
        
        for i = 1:n_max
            % 1. 采样传播距离
            d_i = -log(1 - rand) / ke;
            
            % 2. 采样散射天顶角 (保持原代码的采样策略：第1次均匀分布于波束内，之后均匀分布于0-pi)
            if i == 1
                theta_i = (beta_T / 2) * rand;
            else
                theta_i = pi * rand;
            end
            
            % 3. 采样散射方位角
            phi_i = 2 * pi * rand;
            
            % 4. 计算目标函数的乘积因子
            if i == 1
                f_Theta = sin(theta_i) / (1 - cos(beta_T/2)); 
                O_star_i = f_Theta * (beta_T/2);
            else
                % === 修改点：根据模型选择相函数 PDF ===
                if strcmp(model_type, 'HG')
                    % HG 相函数
                    f_Theta = f_Theta_HG(theta_i, h_HG);
                else
                    % 原 UV 混合相函数
                    f_Theta = (kr/ks) * f_Theta_Rayleigh(theta_i, gamma) + (km/ks) * f_Theta_Mie(theta_i, g, f);
                end
                
                % 由于 theta_i 是在 [0, pi] 上均匀采样，其 proposal PDF 为 1/pi
                % 权重 weight = f_target / f_proposal = f_Theta / (1/pi) = f_Theta * pi
                O_star_i = (ks/ke) * f_Theta * pi; 
            end
            O_star = O_star * O_star_i;
            
            % 更新状态
            mu = update_direction(mu, theta_i, phi_i);
            pos = pos + d_i * mu;
            d_total = d_total + d_i;
            
            % 计算接收参数
            d_to_rx = norm(pos);
            cos_phi_r = dot(mu_R, pos) / d_to_rx;
            
            % 检查视场角
            if cos_phi_r >= cos(beta_R/2)
                theta_n = acos(dot(-mu, pos) / d_to_rx);
                
                % === 修改点：计算第 n 次散射的概率密度 ===
                if strcmp(model_type, 'HG')
                    f_Theta_n = f_Theta_HG(theta_n, h_HG);
                else
                    f_Theta_n = (kr/ks) * f_Theta_Rayleigh(theta_n, gamma) + (km/ks) * f_Theta_Mie(theta_n, g, f);
                end
                
                Omega_r = A_r / d_to_rx^2 * cos_phi_r; 
                p_d = ks/ke * exp(-ke * d_to_rx)^2 * cos_phi_r * min(1, f_Theta_n * Omega_r / (2*pi*sin(theta_n)));
                
                P_cur = p_d * O_star;
                
                P_hat(i) = P_hat(i) + P_cur;
                
                t_total = (d_total + d_to_rx) / c;
                bin_idx = find(t_total >= T_bins(1:end-1) & t_total < T_bins(2:end), 1);
                if ~isempty(bin_idx)
                    h_hat_i(i, bin_idx) = h_hat_i(i, bin_idx) + P_cur;
                end
            end
        end
    end
    
    % 归一化处理
    P_hat = P_hat / N;
    P = sum(P_hat);
    h_hat_i = h_hat_i / N;
    h_hat = sum(h_hat_i, 1);
    
    delta_T = T_bins(2)-T_bins(1);
    h_hat = h_hat / delta_T;
    for i = 1:n_max
        h_hat_i(i, :) = h_hat_i(i, :) / delta_T;
    end
end

% === 辅助函数区域 ===

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

function f = f_Theta_Rayleigh(theta, gamma)
    f = 3 * (1 + 3*gamma + (1-gamma)*cos(theta)^2) * sin(theta) / (8 * (1 + 2*gamma));
end

function f = f_Theta_Mie(theta, g, f_param)
    f = (1 - g^2)/2 * (1/(1 + g^2 - 2*g*cos(theta))^(3/2) + f_param * 0.5 * (3*cos(theta)^2 - 1) / (1 + g^2)^(3/2)) * sin(theta);
end

% === 新增：HG相函数的PDF ===
function f = f_Theta_HG(theta, g)
    % Henyey-Greenstein function PDF f_Theta(theta)
    % P(theta) = (1-g^2) / (4*pi*(1+g^2-2g*cos(theta))^1.5)
    % f_Theta(theta) = \int_0^2pi P(theta) sin(theta) dphi = 2*pi * P(theta) * sin(theta)
    % 结果: f = (1-g^2) / (2*(1+g^2-2g*cos(theta))^1.5) * sin(theta)
    
    numerator = 1 - g^2;
    denominator = 2 * (1 + g^2 - 2*g*cos(theta)).^(1.5);
    f = (numerator ./ denominator) .* sin(theta);
end

% ==============================
% 主脚本部分：参数设置与调用
% ==============================

% 初始化参数结构体
param = struct();

% --- 通用几何参数 ---
param.r = 30;              % 通信距离 (m) (水下距离通常较短)
param.theta_T = deg2rad(20); 
param.phi_T = 0;           % 发射器角度
param.theta_R = deg2rad(20); 
param.phi_R = pi;          % 接收器角度 (对准发射端)
param.beta_T = deg2rad(17); 
param.beta_R = deg2rad(45); 
param.A_r = 1e-4;          % 接收器面积

% --- 核心修改：模型选择 ---
% 可选: 'UV' (大气,使用Mie/Rayleigh) 或 'HG' (水下,使用Henyey-Greenstein)
param.model_type = 'HG';   

if strcmp(param.model_type, 'HG')
    % --- HG模型参数 (水下环境) ---
    % 典型水下参数 (例如: 清澈海水)
    param.ka = 0.05;       % 吸收系数 (1/m)
    param.ks = 0.20;       % 散射系数 (1/m) (HG模式下直接指定ks)
    param.kr = 0; 
    param.km = param.ks; % 保持兼容性，kr设为0
    
    % Henyey-Greenstein 非对称因子 g
    % g = 0: 各向同性散射
    % g -> 1: 强前向散射 (水下通常 0.8 ~ 0.95)
    param.h = 0.924;         
    
    % 下面参数在 HG 模式下不使用，但为了防止报错保留默认值
    param.gamma = 0; param.g = 0; param.f = 0; 
else
    % --- UV模型参数 (大气环境) ---
    param.ka = 0.802e-3; 
    param.kr = 0.266e-3; 
    param.km = 0.284e-3; 
    param.gamma = 0.017; 
    param.g = 0.72; 
    param.f = 0.5; 
    param.h = 0; 
end

% --- 仿真控制参数 ---
param.N = 1e5;            % 采样点数 (调试时可减少)
param.n_max = 4;          % 最大散射阶数
param.c = 2.25e8;         % 水中光速 (approx c/1.33)

% 时间相关
param.t_max = 300e-9;     % 300 ns
param.dt = 1e-9;          % 1 ns
param.T_bins = 0:param.dt:param.t_max;

% 运行模型
fprintf('当前模型模式: %s\n', param.model_type);
[P, P_hat, h_hat, h_hat_i] = MCI_PIS_Function_Modified(param);

% 结果计算与绘图
PL = 10 * log10(1 / P);
if ~isempty(h_hat_i)
    PL1 = 10*log10(1 / trapz(param.T_bins, h_hat_i(1,:)));
else
    PL1 = Inf;
end

fprintf('总路径损失: %.2f dB\n', PL);
fprintf('1阶散射路径损失: %.2f dB\n', PL1);

figure;
plot(param.T_bins, h_hat, 'k-', 'LineWidth', 0.5); hold on;
colors = {'r--', 'g--', 'b--', 'm--'};
for i = 1:min(param.n_max, 4)
    plot(param.T_bins, h_hat_i(i,:), colors{i}, 'DisplayName', sprintf('%d-order', i));
end
xlabel('Time (s)');
ylabel('Impulse Response (W/m^2)');
title(['Channel Impulse Response (' param.model_type ' Mode)']);
legend('Total', '1st', '2nd', '3rd', '4th');
grid on;