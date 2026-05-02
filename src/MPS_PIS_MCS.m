%% UV Communication Model Comparison: MCI_PIS vs. MC_MPS vs. MCS
% 目的: 对比三种蒙特卡洛策略在相同物理环境下的计算结果与效率
% 1. MCI_PIS: 经典概率积分 (每个散射点均估计)
% 2. MC_MPS:  修正框架 (同 PIS 物理核)
% 3. MCS:     视场截断法 (仅在进入 FOV 后估计并终止，复现 Ding 2009 逻辑)

clc; clear; close all;

%% ================= 1. 全局参数定义 (Standardized) =================
param = struct();
% --- 几何 ---
param.r = 100;              % 距离 (m)
param.theta_T = deg2rad(45);% Tx 仰角
param.phi_T = 2*pi-deg2rad(90); % Tx 方位 (指向 -y)
param.theta_R = deg2rad(45);% Rx 仰角
param.phi_R = deg2rad(90);  % Rx 方位 (指向 +y)
param.beta_T = deg2rad(17); % 发射发散角
param.beta_R = deg2rad(30); % 接收 FOV
param.A_r = 1.77e-4;        % 接收面积

% --- 介质 (UV 266nm) ---
param.ka = 0.802e-3;        
param.kr = 0.266e-3;        
param.km = 0.284e-3;        
param.ks = param.kr + param.km; 
param.ke = param.ka + param.ks; 
param.albedo = param.ks / param.ke;
param.gamma = 0.017;        
param.g = 0.72; 
param.f = 0.5;

% --- 仿真控制 ---
param.c = 2.997046e8;       
param.dt = 1e-8;            
param.t_max = 2e-6; 
param.T_bins = 0:param.dt:param.t_max;
param.N = 1e5;              % 光子数
param.n_max = 10;           % 最大散射阶数 (MCS 需要多跳才能进 FOV)
param.bound = 5 * param.r;  % 空间边界 (防止 MCS 跑飞)

% --- 预计算矢量 ---
mu_T = [cos(param.phi_T)*sin(param.theta_T); sin(param.phi_T)*sin(param.theta_T); cos(param.theta_T)];
mu_R = [cos(param.phi_R)*sin(param.theta_R); sin(param.phi_R)*sin(param.theta_R); cos(param.theta_R)];
Tx_Pos = [0; param.r; 0];   
Rx_Pos = [0; 0; 0];         

fprintf('=== 开始仿真对比 (N=%d) ===\n', param.N);

%% ================= 2. 算法 1: MCI_PIS (基准) =================
h_pis = zeros(1, length(param.T_bins));
total_E_pis = 0;
tic;
for k = 1:param.N
    pos = Tx_Pos; mu = mu_T; d_total = 0; 
    % 初始权重 (修正版)
    theta_i = (param.beta_T / 2) * rand;
    weight = (param.beta_T/2) * sin(theta_i) / (1 - cos(param.beta_T/2)); 
    mu = update_dir(mu, theta_i, 2*pi*rand);

    for i = 1:param.n_max
        d = -log(1-rand)/param.ke; pos = pos + d*mu; d_total = d_total + d;
        
        % PIS 估计
        v2r = Rx_Pos - pos; d2r = norm(v2r); dir2r = v2r/d2r;
        cos_phi = dot(mu_R, -dir2r); % 注意方向: 接收机主轴 vs 入射光(反向)
        % 修正: MCI_PIS 代码中常把 pos 当作源, 所以检测 pos 是否在 FOV
        % 此时 v2r 是 pos 指向 Rx. 接收机主轴 mu_R 和 (Rx-pos) 的夹角?
        % 不, 应该是 mu_R 和 (pos-Rx) 的夹角. 
        % 让我们保持之前的逻辑: dot(mu_R, pos-Rx)/dist 
        cos_fov = dot(mu_R, -dir2r); % pos在Rx看来得在FOV里
        
        if cos_fov >= cos(param.beta_R/2)
            theta_s = acos(dot(mu, dir2r)); % 散射角
            p_val = phase_fun(theta_s, param);
            % 修正版概率公式
            p_d = (param.ks/param.ke) * exp(-param.ke*d2r) * p_val * (param.A_r/d2r^2) * cos_fov;
            E = weight * p_d;
            
            total_E_pis = total_E_pis + E;
            bin = floor((d_total + d2r)/param.c/param.dt)+1;
            if bin<=length(h_pis), h_pis(bin) = h_pis(bin) + E; end
        end
        
        % 散射
        weight = weight * param.albedo;
        if weight < 1e-9, break; end
        mu = scatter_dir(mu, param);
    end
end
t_pis = toc;
PL_pis = 10*log10(param.N / total_E_pis);

%% ================= 3. 算法 2: MCS (视场截断 + 能量截断) =================
% 特征: 只有当光子物理移动进入 FOV 区域时, 才计算概率并停止
h_mcs = zeros(1, length(param.T_bins));
total_E_mcs = 0;
tic;
for k = 1:param.N
    pos = Tx_Pos; mu = mu_T; d_total = 0;
    theta_i = (param.beta_T / 2) * rand;
    weight = (param.beta_T/2) * sin(theta_i) / (1 - cos(param.beta_T/2)); 
    mu = update_dir(mu, theta_i, 2*pi*rand);
    
    for i = 1:param.n_max
        d = -log(1-rand)/param.ke; 
        pos_new = pos + d*mu; 
        
        % [MCS 特有] 边界检查
        if max(abs(pos_new)) > param.bound, break; end
        
        d_total = d_total + d;
        pos = pos_new;
        
        % [MCS 特有] 检查是否在接收机 FOV 锥体内部
        v2r = pos - Rx_Pos; d2r = norm(v2r);
        % 判断 pos 是否在 Rx 的视场锥内: dot(mu_R, pos-Rx) / dist >= cos(FOV/2)
        if d2r > 1e-6
            cos_fov = dot(mu_R, v2r) / d2r;
            
            if cos_fov >= cos(param.beta_R/2)
                % === 光子进入 FOV, 触发强制探测 (Forced Detection) ===
                % 计算从当前位置 pos 直连 Rx 的概率 (类似 PIS, 但只做一次)
                % 注意: 这里假设光子在 pos 处再次散射指向 Rx 中心
                % Ding 2009: "The arrival probability... is calculated... The photon is then terminated"
                
                theta_s = acos(dot(mu, -v2r/d2r)); % 散射向 Rx
                p_val = phase_fun(theta_s, param);
                
                % 此时不需要再乘 exp(-ke*d2r) ? 
                % 论文逻辑: 光子实际上已经到了 pos. "概率"是指它下一跳刚好打中 detector 的概率.
                % 所以需要: P(scatter to Rx) * P(transmission) * SolidAngle
                
                prob_hit = p_val * (param.A_r / d2r^2) * cos_fov * exp(-param.ke * d2r);
                
                % MCS 能量 = 当前光子权重 * 吸收衰减(已在weight里?) * 探测概率
                % 注意: 上面的循环里我们还没乘 albedo (生存概率).
                % 在 MCS 中, 移动 d 后, 权重应乘 exp(-ka * d) ?
                % 通常 Monte Carlo: d 由 -log(rand)/ke 采样, 则不需显式乘 exp(-ke*d).
                % 但需乘 albedo (ks/ke) 在散射时.
                % 这里是"入场即算", 尚未发生第 i 次散射. 
                % 所以 E = weight * prob_hit.
                
                % 此外需乘 param.ks/param.ke ? 
                % 如果看作是"在 pos 处发生了散射指向 Rx", 则需要乘 albedo.
                E = weight * (param.ks/param.ke) * prob_hit;
                
                total_E_mcs = total_E_mcs + E;
                t_arr = (d_total + d2r)/param.c;
                bin = floor(t_arr/param.dt)+1;
                if bin<=length(h_mcs), h_mcs(bin) = h_mcs(bin) + E; end
                
                break; % MCS: 贡献一次后立即终止
            end
        end
        
        % 若未被接收, 继续散射
        weight = weight * param.albedo;
        
        % [MCS 特有] 能量截断
        if weight < 1e-9, break; end
        
        mu = scatter_dir(mu, param);
    end
end
t_mcs = toc;
PL_mcs = 10*log10(param.N / total_E_mcs);

%% ================= 4. 结果对比 =================
fprintf('\n--- 算法性能对比 ---\n');
fprintf('MCI_PIS: Path Loss = %.2f dB | Time = %.4f s\n', PL_pis, t_pis);
fprintf('MCS    : Path Loss = %.2f dB | Time = %.4f s\n', PL_mcs, t_mcs);
fprintf('差异    : %.2f dB\n', abs(PL_pis - PL_mcs));

% 绘图
figure('Position',[200,200,900,400]);
subplot(1,2,1);
plot(param.T_bins*1e6, h_pis/(param.N*param.dt), 'r-', 'LineWidth', 1.5); hold on;
plot(param.T_bins*1e6, h_mcs/(param.N*param.dt), 'b--', 'LineWidth', 1.0);
xlabel('Time (\mus)'); ylabel('Impulse Response (W/m^2/s)');
legend('MCI PIS', 'MCS (Ding)');
title('CIR Comparison'); grid on;

subplot(1,2,2);
bar([t_pis, t_mcs]);
set(gca, 'XTickLabel', {'PIS', 'MCS'});
ylabel('Simulation Time (s)');
title('Computational Cost');

%% ================= 5. 辅助函数 =================
function val = phase_fun(theta, p)
    % 标准单位立体角相函数 (integral = 1)
    p_ray = 3*(1+3*p.gamma+(1-p.gamma)*cos(theta)^2)/(16*pi*(1+2*p.gamma));
    p_mie = (1-p.g^2)/(4*pi) * (1./(1+p.g^2-2*p.g*cos(theta)).^1.5 + ...
             p.f*0.5*(3*cos(theta)^2-1)./(1+p.g^2).^1.5);
    val = (p.kr/p.ks)*p_ray + (p.km/p.ks)*p_mie;
end

function mu_new = update_dir(mu, theta, phi)
    % 罗德里格斯旋转或坐标变换
    ux = mu(1); uy = mu(2); uz = mu(3);
    if abs(uz) < 0.99999
        inv = 1/sqrt(1-uz^2);
        mu_new = [sin(theta)*(ux*uz*cos(phi)-uy*sin(phi))*inv + ux*cos(theta);
                  sin(theta)*(uy*uz*cos(phi)+ux*sin(phi))*inv + uy*cos(theta);
                  -sin(theta)*cos(phi)/inv + uz*cos(theta)];
    else
        % 接近 Z 轴 (Singularity)
        mu_new = [sin(theta)*cos(phi); sin(theta)*sin(phi); sign(uz)*cos(theta)];
    end
    mu_new = mu_new / norm(mu_new);
end

function mu_new = scatter_dir(mu, p)
    % 采样新的散射方向
    % 1. 决定 Rayleigh 或 Mie
    if rand < (p.kr/p.ks)
        % Rayleigh 近似 (简化: 均匀分布或Rejection Sampling, 这里用简化偶极子)
        th = acos(2*rand-1); % 粗略近似
    else
        % HG 近似
        g = p.g;
        if abs(g) < 1e-3
            th = acos(2*rand-1);
        else
            sq = (1-g^2)/(1-g+2*g*rand);
            th = acos((1+g^2 - sq^2)/(2*g));
        end
    end
    ph = 2*pi*rand;
    mu_new = update_dir(mu, th, ph);
end