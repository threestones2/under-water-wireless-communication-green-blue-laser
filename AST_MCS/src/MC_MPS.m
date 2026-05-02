%% 水下无线光通信混合仿真: MC-PIS (散射) + MPS (湍流) + Eq.10 (衍射衰减)
% 
% 功能:
% 1. 生成 OTOPS 海洋湍流相位屏 (次谐波补偿版)
% 2. 蒙特卡洛光子追踪 (Monte Carlo Ray Tracing)
% 3. 相位屏交互: 使用相位梯度法模拟波前畸变 (几何光学近似)
% 4. 传输衰减: 包含吸收(a)、散射(b) 以及 Eq.10 定义的湍流衍射衰减(alpha_turb)
% 5. 重要性采样: 使用下一事件估计 (Next Event Estimation) 加速收敛
%
% 参考: Wen et al. (2023), Yuan et al. (2020), Schmidt (2010)

clc; clear; close all;

%% ================= 1. 参数设置 =================
% --- 系统参数 ---
lambda_nm = 532;            % 波长 (nm)
lambda = lambda_nm * 1e-9;  % 波长 (m)
k_wave = 2*pi/lambda;       % 波数
L = 100;                     % 总传输距离 (m)
N_screens = 5;              % 相位屏数量 (等间距分布)
D_screen = 1.0;             % 相位屏边长 (m)
N_grid = 256;               % 相位屏网格数 (256x256)

% --- 发射与接收 ---
w0 = 0.01;                  % 发射束腰半径 (m)
div_angle = 1e-3;           % 初始发散角 (rad)
rx_aperture = 0.5;         % 接收孔径直径 (m)
rx_fov = 30 * pi/180;       % 接收视场角 (rad)

% --- 海水信道参数 (Coastal) ---
coef_a = 0.179;             % 吸收系数 (1/m)
coef_b = 0.219;             % 散射系数 (1/m)
coef_c = coef_a + coef_b;   % 衰减系数
g_HG = 0.924;               % Henyey-Greenstein 散射各向异性因子

% --- 湍流参数 (OTOPS) ---
T_avg = 20; S_avg = 35;     % 温盐度
epsilon = 1e-5;             % 动能耗散率 (m^2/s^3)
chi_T = 1e-7;               % 温度方差耗散率 (K^2/s)
eta = 1e-3;                 % 内尺度 (m)
H_ratio = -20;              % 温盐比

% --- MC 仿真参数 ---
N_packets = 1e5;            % 光子包数量 (建议 >1e5 以获得平滑结果)
max_scatter = 3;           % 最大散射次数

%% ================= 2. 预计算: 湍流衰减因子 (Eq.10) & 相位屏 =================
fprintf('正在计算湍流参数与相位屏...\n');

% 2.1 计算 Eq.10 中的湍流衰减系数 alpha_turb
% 公式: alpha = 4 * pi^2 * k^2 * integral(Phi_n * k, 0, 2k)
% 我们使用数值积分计算
dk_int = 0.1; 
k_int = 0.1:dk_int:(2*k_wave);
[Phi_val, ~] = calc_otops_spectrum(k_int, epsilon, chi_T, eta, T_avg, S_avg, lambda_nm, H_ratio);
integrand = Phi_val .* k_int;
integral_val = sum(integrand) * dk_int;
alpha_turb = 4 * pi^2 * k_wave^2 * integral_val;

fprintf('  > 海水衰减系数 c = %.4f /m\n', coef_c);
fprintf('  > 湍流衍射衰减 alpha (Eq.10) = %.4e /m\n', alpha_turb);

% 2.2 生成多张相位屏
% 每张屏幕代表 L/N_screens 距离的累积相位
dz_screen = L / N_screens;
screens = struct('data', [], 'x', [], 'y', [], 'grad_x', [], 'grad_y', []);
dx = D_screen / N_grid;
x_axis = (-N_grid/2 : N_grid/2-1) * dx;

for i = 1:N_screens
    % 调用相位屏生成函数 (次谐波补偿版)
    phi_map = generate_phase_screen_SH_Corrected(D_screen, N_grid, lambda, dz_screen, ...
                                                 T_avg, S_avg, epsilon, chi_T, eta, H_ratio);
    
    % 预计算梯度 (用于光子偏折)
    [gx, gy] = gradient(phi_map, dx);
    % 物理偏折角 theta = grad_phi / k
    screens(i).data = phi_map;
    screens(i).grad_x = gx / k_wave; 
    screens(i).grad_y = gy / k_wave;
    screens(i).z_pos = i * dz_screen;
    screens(i).x_axis = x_axis;
end

%% ================= 3. 主循环: Monte Carlo 光子追踪 =================
fprintf('开始 MC 仿真 (N=%d)...\n', N_packets);
Rx_Power_Total = 0;
Rx_Power_Ballistic = 0;
Rx_Power_Scattered = 0;

tic;
for p = 1:N_packets
    % --- 初始化光子 ---
    % 空间: 高斯分布
    r = w0 * sqrt(-0.5 * log(rand())); 
    phi_pos = 2*pi*rand();
    pos = [r*cos(phi_pos), r*sin(phi_pos), 0]; % x, y, z
    
    % 角度: 高斯发散角
    theta = div_angle * sqrt(-log(rand()));
    phi_dir = 2*pi*rand();
    % 方向余弦 (dir_x, dir_y, dir_z)
    sin_t = sin(theta);
    dir = [sin_t*cos(phi_dir), sin_t*sin(phi_dir), cos(theta)];
    
    weight = 1.0;
    current_z = 0;
    screen_idx = 1;
    
    % --- 路径追踪 ---
    for s = 1:max_scatter
        % 1. 确定步长
        % 采样散射步长 (仅由 scattering 决定，吸收和湍流处理在权重中)
        dist_to_scatter = -log(rand()) / coef_b;
        
        % 检查是否会先撞到相位屏
        if screen_idx <= N_screens
            dist_to_screen = screens(screen_idx).z_pos - pos(3);
        else
            dist_to_screen = Inf;
        end
        
        % 检查是否会先到达接收面
        dist_to_rx = L - pos(3);
        
        % 选取最小距离
        [step_dist, event_type] = min([dist_to_scatter, dist_to_screen, dist_to_rx]);
        % event_type: 1=散射, 2=相位屏, 3=接收面
        
        % 2. 移动光子
        pos = pos + dir * step_dist;
        
        % 3. 应用连续衰减 (Beer-Lambert + Eq.10 湍流衰减)
        % 注意: 散射衰减由事件采样处理，这里只处理吸收和湍流衍射损耗
        % 但在加权 MC 中，通常步长由 b 采样，权重乘 exp(-a*step)。
        % 此处我们加上 Eq.10: weight *= exp(-(a + alpha_turb) * step)
        path_loss = coef_a + alpha_turb;
        weight = weight * exp(-path_loss * step_dist);
        
        if weight < 1e-9, break; end
        
        % 4. 处理事件
        if event_type == 3 % 到达接收面
            % 检查是否在孔径和视场内
            r_spot = sqrt(pos(1)^2 + pos(2)^2);
            angle_inc = acos(dir(3)); % 入射角 (相对于 z 轴)
            
            if r_spot <= (rx_aperture/2) && angle_inc <= rx_fov
                Rx_Power_Total = Rx_Power_Total + weight;
                if s == 1 % 未发生散射 (只经历了相位屏偏折)
                    Rx_Power_Ballistic = Rx_Power_Ballistic + weight;
                else
                    Rx_Power_Scattered = Rx_Power_Scattered + weight;
                end
            end
            break; % 结束该光子
            
        elseif event_type == 2 % 撞到相位屏
            % 湍流相互作用: 更新方向 (参考 ang_spec_prop 的几何等效)
            % 读取梯度
            % 将物理坐标映射到网格索引
            c_idx = round((pos(1) - x_axis(1))/dx) + 1;
            r_idx = round((pos(2) - x_axis(1))/dx) + 1;
            
            if c_idx > 0 && c_idx <= N_grid && r_idx > 0 && r_idx <= N_grid
                d_theta_x = screens(screen_idx).grad_x(r_idx, c_idx);
                d_theta_y = screens(screen_idx).grad_y(r_idx, c_idx);
                
                % 更新方向向量 (小角度近似: dir_new = dir + [dx, dy, 0])
                % 更精确的做法是旋转坐标系
                dir(1) = dir(1) + d_theta_x;
                dir(2) = dir(2) + d_theta_y;
                % 归一化
                dir = dir / norm(dir);
            end
            screen_idx = screen_idx + 1;
            % 光子继续传输，不改变 s (不计为散射事件)
            s = s - 1; 
            
        elseif event_type == 1 % 发生散射
            % --- 重要性采样 (Next Event Estimation) ---
            % 估计直接散射到接收机的概率
            vec_to_rx = [0, 0, L] - pos;
            dist_rx = norm(vec_to_rx);
            dir_to_rx = vec_to_rx / dist_rx;
            
            % 计算散射角余弦
            cos_theta_s = dot(dir, dir_to_rx);
            
            % 计算 HG 相函数概率密度 p(cos_theta)
            p_hg = (1 - g_HG^2) / (4*pi * (1 + g_HG^2 - 2*g_HG*cos_theta_s)^1.5);
            
            % 接收机立体角 approx
            solid_angle = (pi * (rx_aperture/2)^2 * abs(dir_to_rx(3))) / (dist_rx^2);
            
            % 路径上的衰减 (到接收机)
            % 注意: 这一段是虚拟路径，需要乘上所有的衰减 (c + alpha)
            loss_to_rx = exp(-(coef_c + alpha_turb) * dist_rx);
            
            % 视场角检查
            angle_at_rx = acos(-dir_to_rx(3)); % 接收机看来的入射角
            if angle_at_rx <= rx_fov
                % 累加估计功率
                P_est = weight * p_hg * solid_angle * loss_to_rx;
                Rx_Power_Scattered = Rx_Power_Scattered + P_est;
                Rx_Power_Total = Rx_Power_Total + P_est;
            end
            
            % --- 真实散射 (继续游走) ---
            % 使用 HG 函数抽样新的方向
            dir = scatter_Henyey_Greenstein(dir, g_HG);
            
            % 能量权重不再减小 (因为 step 已经处理了 absorbtion)
            % 但如果是 implicit capture, 权重应该保持。
            % 注意: Next Event Estimation 是额外的估计，不消耗光子本体。
        end
    end
    
    if mod(p, N_packets/10) == 0
        fprintf('  已完成 %.0f%%\n', p/N_packets*100);
    end
end
toc;

%% ================= 4. 结果分析 =================
Rx_Power_Total = Rx_Power_Total / N_packets;
Rx_Power_Ballistic = Rx_Power_Ballistic / N_packets;
Rx_Power_Scattered = Rx_Power_Scattered / N_packets;

fprintf('\n=== 仿真结果 (L = %.1f m) ===\n', L);
fprintf('总归一化接收功率: %.6e\n', Rx_Power_Total);
fprintf('  - 准弹道光 (Turbulence Only): %.6e\n', Rx_Power_Ballistic);
fprintf('  - 散射光 (Scattered): %.6e\n', Rx_Power_Scattered);
loss_dB = 10*log10(Rx_Power_Total);
fprintf('链路损耗: %.2f dB\n', loss_dB);


%% ================= 5. 辅助函数 =================
function dir_new = scatter_Henyey_Greenstein(dir_old, g)
    % HG 相函数散射方向更新
    % 采样天顶角 theta
    rnd = rand();
    if g == 0
        cos_theta = 2*rnd - 1;
    else
        temp = (1 - g^2) / (1 - g + 2*g*rnd);
        cos_theta = (1 + g^2 - temp^2) / (2*g);
    end
    sin_theta = sqrt(1 - cos_theta^2);
    phi = 2*pi*rand();
    
    % 旋转坐标系 (Rodrigues / Local frame)
    u = dir_old(1); v = dir_old(2); w = dir_old(3);
    if abs(w) > 0.99999
        dir_new = [sin_theta*cos(phi), sin_theta*sin(phi), sign(w)*cos_theta];
    else
        factor = sqrt(1 - w^2);
        dir_new = zeros(1,3);
        dir_new(1) = sin_theta * (u * w * cos(phi) - v * sin(phi)) / factor + u * cos_theta;
        dir_new(2) = sin_theta * (v * w * cos(phi) + u * sin(phi)) / factor + v * cos_theta;
        dir_new(3) = -sin_theta * cos(phi) * factor + w * cos_theta;
    end
end
%%
function [Phi_n, Phi_Hill_Handle] = calc_otops_spectrum(k_vals, eps, chi, eta, T, S, lam_nm, Hr)
    % 计算 OTOPS 功率谱密度 
    % 系数
    a1=1.779e-4; a2=-1.05e-6; a3=1.6e-8; a4=-2.02e-6; a5=1.155e-2; a6=-4.23e-3;
    A = a2*S + 2*a3*T*S + 2*a4*T + a6/lam_nm; 
    B = a1 + a2*T + a3*T^2 + a5/lam_nm;
    
    rho=1025; mu=1e-3; cp=4182; 
    Pr=mu*cp/0.6; Sc=mu^2/(5.954e-15*(T+273.15)*rho);
    
    alpha_c=2.6e-4; beta_c=7.6e-4;
    R_rho=alpha_c*abs(Hr)/beta_c;
    if R_rho>=1, dr=R_rho+sqrt(R_rho*(R_rho-1)); else, dr=1.85*R_rho-0.85; end
    
    chi_S = chi*dr/Hr^2; chi_TS = chi*(1+dr)/(2*Hr);
    
    % Hill Spectrum Function
    f_Hill = @(K, chi_x, c_x) 0.033 * chi_x * eps^(-1/3) .* (K.^2 + 1e-20).^(-11/6) .* ...
        exp(-176.9*(K*eta).^2 .* c_x^(0.96)) .* ...
        (1 + 21.61*(K*eta).^(0.61).*c_x^(0.02) - 18.18*(K*eta).^(0.55).*c_x^(0.04));
    
    c_T=0.072^(4/3)/Pr; c_S=0.072^(4/3)/Sc; c_TS=0.072^(4/3)*(Pr+Sc)/(2*Pr*Sc);
    
    Phi_n = A^2*f_Hill(k_vals, chi, c_T) + B^2*f_Hill(k_vals, chi_S, c_S) + ...
            2*A*B*f_Hill(k_vals, chi_TS, c_TS);
    
    Phi_Hill_Handle = f_Hill; % 返回句柄以备用
end
