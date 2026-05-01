%% 水下无线光通信信道MC-MPS模型仿真
% 参考文献: Wen et al., "Modeling and performance analysis of underwater 
% wireless optical absorption, scattering, and turbulence channels employing 
% Monte Carlo-multiple phase screens", Applied Optics, Vol. 62, No. 26, 2023.

clc;
clear;
close all;

%% 1. 参数设置
% 系统参数
lambda = 532e-9;          % 光波长 (m)
r_T = 0.01;               % 发射机束腰半径 (m)
alpha = deg2rad(20);      % 发散角 (rad)
r_R = 5*0.0254;           % 接收机孔径半径 (m), 5英寸转换为米
beta_R = deg2rad(30);        % 接收机视场角 (rad)
L = 10;                   % 传输距离 (m)
N_p = 1e5;                % 光子包数量
N = 256;                  % 相位屏采样点数
D = 1.0;                  % 相位屏边长 (m)
delta = D/N;              % 相位屏网格间距 (m)
u_threshold = 1e-6;        % 光强度阈值

% 海水参数
water_type = 'coastal';   % 海水类型: 'clear', 'coastal', 'harbor'
switch water_type
    case 'clear'
        a = 0.1125;       % 吸收系数 (m^-1)
        b = 0.0375;       % 散射系数 (m^-1)
    case 'coastal'
        a = 0.18;         % 吸收系数 (m^-1)
        b = 0.22;         % 散射系数 (m^-1)
    case 'harbor'
        a = 0.3723;       % 吸收系数 (m^-1)
        b = 1.8177;       % 散射系数 (m^-1)
end
c = a + b;                % 衰减系数 (m^-1)
W0 = b/c;                 % 散射反照率
g = 0.924;                % 各向异性因子

% 湍流参数
turbulence_type = 'weak'; % 湍流强度: 'none', 'weak', 'strong'
T_avg = 20;               % 平均温度 (°C)
S_avg = 35;               % 平均盐度 (ppt)
switch turbulence_type
    case 'none'
        epsilon = 0;      % 湍流动能耗散率 (m^2/s^3)
        chi_T = 0;        % 均方温度耗散率 (K^2/s^3)
    case 'weak'
        epsilon = 1e-3;   % 湍流动能耗散率 (m^2/s^3)
        chi_T = 1e-7;     % 均方温度耗散率 (K^2/s^3)
    case 'strong'
        epsilon = 1e-10;  % 湍流动能耗散率 (m^2/s^3)
        chi_T = 1e-5;     % 均方温度耗散率 (K^2/s^3)
end
eta = 1e-3;               % 湍流内尺度 (m)
H = -20;                   % 温度和盐度随深度变化比率 (°C·ppt^-1)

% 计算其他湍流相关参数
k = 2*pi/lambda;          % 波数
d = min(delta^2*N/lambda, 2); % 相位屏间距 (m)

% 计算TEOS-10相关参数
[c_p, sigma_T, mu, rho] = calculate_teos10_params(T_avg, S_avg);%这个函数是简化版
Pr = mu*c_p/sigma_T;      % 普朗特数
Sc = mu^2/(5.954e-15*(T_avg+273.15)*rho); % 施密特数

% 计算湍流引起的衰减,对应公式10
if epsilon > 0 && chi_T > 0
    gamma = calculate_turbulence_attenuation(lambda, epsilon, chi_T, T_avg, S_avg, H, eta);
else
    gamma = 0;
end

%% 2. 初始化
% 计算焦距
f = -2*r_T/alpha;

% 初始化接收统计
received_power = 0;
received_photons = 0;
received_intensity = [];
received_time = [];
path_loss = 0;

% 初始化光子包数组
photons = struct('position', [], 'direction', [], 'weight', [], 'field', [], 'active', true(N_p,1));

% 初始化所有光子包
for i = 1:N_p
    % 初始位置 (高斯分布)
    phi_0 = 2*pi*rand();
    r_0 = r_T*sqrt(-log(1-rand()));
    x_0 = r_0*cos(phi_0);
    y_0 = r_0*sin(phi_0);
    z_0 = 0;
    photons(i).position = [x_0, y_0, z_0];
    
    % 初始方向
    theta_0 = -r_0/f;
    photons(i).direction = [sin(theta_0)*cos(phi_0), sin(theta_0)*sin(phi_0), -cos(theta_0)];
    
    % 初始权重
    photons(i).weight = 1;
    
    % 初始光场
    w = r_T; % 模场半径,注意这里用到的公式和正规的不一样，正规的是exp(-2r^2/w^2),这里用的是exp(-r^2/r_T^2)
    A0 = 1;  % 初始幅度
    photons(i).field = A0*exp(-(x_0^2+y_0^2)/w^2);
end

%% 3. 生成多重相位屏
if epsilon > 0 && chi_T > 0
    % 计算OTOPS功率谱
    [kappa_x, kappa_y, Phi_n] = calculate_OTOPS(T_avg, S_avg, epsilon, chi_T, H, eta, lambda, N, D);
    
    % 生成相位屏
    phase_screens = generate_phase_screens(Phi_n, kappa_x, kappa_y, d, k, N, D, L);
    num_screens = size(phase_screens, 3);
else
    phase_screens = [];
    num_screens = 0;
end

%% 4. 蒙特卡洛模拟
fprintf('开始MC-MPS仿真...\n');
tic;

for i = 1:N_p
    if mod(i, 10000) == 0
        fprintf('已处理 %d/%d 光子包\n', i, N_p);
    end
    
    % 当前光子包
    photon = photons(i);
    
    % 如果光子包不活跃，跳过
    if ~photon.active
        continue;
    end
    
    % 光子包传播
    while photon.active
        % 计算步长
        l_i = -log(rand())/c;
        
        % 更新位置
        new_position = photon.position + photon.direction * l_i;
        
        % 检查是否到达接收平面
        if new_position(3) >= L
            % 计算到达接收平面的确切位置
            t = (L - photon.position(3)) / photon.direction(3);
            final_position = photon.position + photon.direction * t;
            
            % 检查是否在接收孔径内
            distance_from_center = sqrt(final_position(1)^2 + final_position(2)^2);
            if distance_from_center <= r_R
                % 检查是否在接收视场内
                if photon.direction(3) >= cos(beta_R/2)
                    % 光子包被接收
                    received_photons = received_photons + 1;
                    received_power = received_power + photon.weight;
                    received_intensity(end+1) = abs(photon.field)^2 * photon.weight;
                    received_time(end+1) = t;
                end
            end
            
            % 停止传播
            photon.active = false;
            break;
        end
        
        % 检查是否超出边界
        if abs(new_position(1)) > 10*D || abs(new_position(2)) > 10*D
            photon.active = false;
            break;
        end
        
        % 更新位置
        photon.position = new_position;
        
        % 处理湍流效应（通过相位屏）
        if num_screens > 0
            % 计算通过相位屏的数量
            screens_passed = floor(l_i / d);
            
            % 通过每个相位屏
            for s = 1:min(screens_passed, num_screens)
                % 获取当前相位屏
                phase_screen = phase_screens(:,:,s);
                
                % 计算光子在相位屏上的位置
                x_idx = round((photon.position(1) + D/2) / delta) + 1;
                y_idx = round((photon.position(2) + D/2) / delta) + 1;
                
                % 确保索引在有效范围内
                x_idx = max(1, min(N, x_idx));
                y_idx = max(1, min(N, y_idx));
                
                % 应用相位屏
                phase_shift = phase_screen(y_idx, x_idx);
                photon.field = photon.field * exp(1j * phase_shift);
            end
        end
        
        % 处理散射
        if rand() < W0
            % 计算散射角
            phi_i = 2*pi*rand();
            zeta_theta = rand();
            
            % 使用Henyey-Greenstein函数计算散射极角
            if g == 0
                cos_theta_i = 1 - 2*zeta_theta;
            else
                cos_theta_i = (1 + g^2 - ((1 - g^2)/(1 - g + 2*g*zeta_theta))^2) / (2*g);
            end
            theta_i = acos(cos_theta_i);
            
            % 更新方向
            if abs(photon.direction(3)) < 0.99999
                sin_theta = sqrt(1 - photon.direction(3)^2);
                temp = [photon.direction(1)/sin_theta, photon.direction(2)/sin_theta, 0;
                        photon.direction(1)*photon.direction(3)/sin_theta, ...
                        photon.direction(2)*photon.direction(3)/sin_theta, -sin_theta;
                        photon.direction(2), -photon.direction(1), 0];
                
                new_dir_local = [sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos_theta_i]';
                photon.direction = temp * new_dir_local;
            else
                photon.direction = [sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos_theta_i * sign(photon.direction(3))]';
            end
            
            % 归一化方向向量
            photon.direction = photon.direction / norm(photon.direction);
        end
        
        % 处理吸收
        photon.weight = photon.weight * exp(-a * l_i);
        
        % 更新光场（考虑吸收和湍流衰减）
        photon.field = photon.field * sqrt(exp(-a * l_i)) * exp(-gamma * l_i / 2);
        
        % 检查光子包权重是否低于阈值
        if photon.weight < u_threshold || abs(photon.field)^2 < u_threshold
            photon.active = false;
            break;
        end
    end
    
    % 更新光子包状态
    photons(i) = photon;
end

%% 5. 计算性能指标
% 计算路径损耗
if received_photons > 0
    path_loss = -10*log10(received_power / N_p); % dB
else
    path_loss = Inf; % 无限大损耗
end

fprintf('仿真完成！耗时 %.2f 秒\n', toc);
fprintf('接收光子数: %d/%d\n', received_photons, N_p);
fprintf('路径损耗: %.2f dB\n', path_loss);

%% 6. 结果分析与可视化
% 接收光强度PDF
if ~isempty(received_intensity)
    figure;
    histogram(received_intensity, 50, 'Normalization', 'pdf');
    title('接收光强度概率密度函数');
    xlabel('归一化光强度');
    ylabel('概率密度');
    grid on;
    
    % 计算闪烁指数
    I_mean = mean(received_intensity);
    I_mean_sq = mean(received_intensity.^2);
    SI = (I_mean_sq - I_mean^2) / I_mean^2;
    fprintf('闪烁指数: %.4f\n', SI);
    
    % 尝试对数正态分布拟合（仅适用于弱湍流）
    if strcmp(turbulence_type, 'weak')
        x = linspace(min(received_intensity), max(received_intensity), 100);
        params = lognfit(received_intensity);
        y = lognpdf(x, params(1), params(2));
        hold on;
        plot(x, y, 'r-', 'LineWidth', 2);
        legend('仿真数据', '对数正态拟合');
    end
end

% 信道冲激响应
if ~isempty(received_time)
    figure;
    histogram(received_time*1e9, 50, 'Normalization', 'probability');
    title('信道冲激响应');
    xlabel('时间 (ns)');
    ylabel('归一化幅度');
    grid on;
end

%% 辅助函数
% 计算TEOS-10相关参数
function [c_p, sigma_T, mu, rho] = calculate_teos10_params(T, S)
    % 简化版本，实际应用中应使用TEOS-10工具箱
    c_p = 4000;      % 比热容 (J/kg/K)
    sigma_T = 0.6;   % 热导率 (W/m/K)
    mu = 1e-3;       % 动态粘度 (Pa·s)
    rho = 1025;      % 密度 (kg/m^3)
end

% 计算湍流引起的衰减
function gamma = calculate_turbulence_attenuation(lambda, epsilon, chi_T, T_avg, S_avg, H, eta)
    k = 2*pi/lambda;
    
    % 计算OTOPS功率谱参数
    a1 = 1.779e-4; a2 = -1.05e-6; a3 = 1.6e-8;
    a4 = -2.02e-6; a5 = 1.155e-2; a6 = -4.23e-3;
    
    A = a2*S_avg + 2*a3*T_avg*S_avg + 2*a4*T_avg + a6/lambda;
    B = a1 + a2*T_avg + a4*T_avg^2 + a5/lambda;
    
    % 计算温度、盐度和协方差谱
    kappa = linspace(0, 2*k, 1000);
    Phi_T = calculate_spectrum(kappa, epsilon, chi_T, eta, T_avg, S_avg, 'T');
    Phi_S = calculate_spectrum(kappa, epsilon, chi_T, eta, T_avg, S_avg, 'S');
    Phi_TS = calculate_spectrum(kappa, epsilon, chi_T, eta, T_avg, S_avg, 'TS');
    
    % 计算折射率起伏功率谱
    Phi_n = A^2*Phi_T + B^2*Phi_S + 2*A*B*Phi_TS;
    
    % 计算衰减
    integrand = Phi_n .* kappa;
    gamma = 4*pi^2*k^2 * trapz(kappa, integrand);
end

% 计算温度、盐度或协方差谱
function Phi_M = calculate_spectrum(kappa, epsilon, chi_M, eta, T_avg, S_avg, M)
    % 计算Prandtl数和Schmidt数
    [c_p, sigma_T, mu, rho] = calculate_teos10_params(T_avg, S_avg);
    Pr = mu*c_p/sigma_T;
    Sc = mu^2/(5.954e-15*(T_avg+273.15)*rho);
    
    % 计算c_M
    switch M
        case 'T'
            c_M = 0.072*(4/3)*0.72*Pr^(-1);
        case 'S'
            c_M = 0.072*(4/3)*0.72*Sc^(-1);
        case 'TS'
            c_M = 0.072*(4/3)*0.72*(Pr+Sc)/(2*Pr*Sc);
    end
    
    % 计算功率谱
    Phi_M = 0.72*chi_M*epsilon^(-1/3)/(4*pi) .* kappa.^(-11/3) .* ...
            exp(-176.90*(kappa*eta).^2/c_M^(0.96)) .* ...
            (1 + 21.61*(kappa*eta).^0.61.*c_M^0.02 - 18.18*(kappa*eta).^0.55.*c_M^0.04);
end

% 计算OTOPS功率谱
function [kappa_x, kappa_y, Phi_n] = calculate_OTOPS(T_avg, S_avg, epsilon, chi_T, H, eta, lambda, N, D)
    % 计算空间频率
    kappa_x = 2*pi*(-N/2:N/2-1)/D;
    kappa_y = 2*pi*(-N/2:N/2-1)/D;
    [KX, KY] = meshgrid(kappa_x, kappa_y);
    kappa = sqrt(KX.^2 + KY.^2);
    
    % 计算OTOPS参数
    a1 = 1.779e-4; a2 = -1.05e-6; a3 = 1.6e-8;
    a4 = -2.02e-6; a5 = 1.155e-2; a6 = -4.23e-3;
    
    A = a2*S_avg + 2*a3*T_avg*S_avg + 2*a4*T_avg + a6/lambda;
    B = a1 + a2*T_avg + a4*T_avg^2 + a5/lambda;
    
    % 计算温度、盐度和协方差谱
    Phi_T = calculate_spectrum(kappa, epsilon, chi_T, eta, T_avg, S_avg, 'T');
    Phi_S = calculate_spectrum(kappa, epsilon, chi_T, eta, T_avg, S_avg, 'S');
    
    % 计算温度-盐度相关参数
    %两个参数都应该从TOES——10里获得
    alpha_c = 2e-4;  % 热膨胀系数    
    beta_c = 7.6e-4; % 盐收缩系数
    omega = alpha_c*H/beta_c;
    
    if abs(omega) >= 1
        d_r = abs(omega) + sqrt(abs(omega)*(abs(omega)-1));
    elseif abs(omega) >= 0.5
        d_r = 1.85*abs(omega) - 0.85;
    else
        d_r = 0.15*abs(omega);
    end
    
    chi_S = chi_T*d_r/H^2;
    chi_TS = chi_T*(1+d_r)/(2*H);
    
    Phi_TS = calculate_spectrum(kappa, epsilon, chi_TS, eta, T_avg, S_avg, 'TS');
    
    % 计算折射率起伏功率谱
    Phi_n = A^2*Phi_T + B^2*Phi_S + 2*A*B*Phi_TS;
end

% 生成相位屏
function phase_screens = generate_phase_screens(Phi_n, kappa_x, kappa_y, d, k, N, D, L)
    % 计算相位屏数量
    num_screens = floor(L/d);
    
    % 初始化相位屏数组
    phase_screens = zeros(N, N, num_screens);
    
    % 生成每个相位屏
    for s = 1:num_screens
        % 生成随机复数高斯场
        cn = (randn(N) + 1j*randn(N)) / sqrt(2);
        
        % 计算相位扰动功率谱
        F_phi = 2*pi^2*k^2*d*Phi_n;
        
        % 应用功率谱
        field = fftshift(sqrt(F_phi)) .* cn;
        
        % 逆傅里叶变换得到相位屏
        phase_screen = real(ifft2(ifftshift(field)));
        
        % 存储相位屏
        phase_screens(:,:,s) = phase_screen;
    end
end
