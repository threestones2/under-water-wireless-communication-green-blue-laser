%%
%本代码尝试复现《Modeling of Non-Line-of-Sight Ultraviolet  Scattering Channels for Communication》论文的蒙特卡洛仿真
%Published in: IEEE Journal on Selected Areas in Communications ( Volume: 27, Issue: 9, December 2009)
%DOI:10.1109/JSAC.2009.091203
%MCS和MCI其实都比较像单次散射模型，进入FOV后就开始,MCS从能量角度限制，MCI从散射次数限制


function [h, path_loss] = monte_carlo_uv_channel(params)
    % 参数初始化（保持不变）
    lambda = params.lambda;
    k_s_ray = params.k_s_ray;
    k_s_mie = params.k_s_mie;
    k_a = params.k_a;
    k_s = k_s_ray + k_s_mie;
    k_e = k_s + k_a;
    
    phi1 = params.phi1;
    phi2 = params.phi2;
    theta1 = params.theta1;
    theta2 = params.theta2;
    psi1 = params.psi1;
    psi2 = params.psi2;
    r = params.r;
    rx_area = params.rx_area;
    rx = sqrt(rx_area/pi);
    
    gamma = params.gamma;
    g = params.g;
    f = params.f;
    c = params.c;% 光速
    
    num_photons = params.num_photons;
    max_scattering = params.max_scattering; % 新增：最大散射次数
    
    max_time = params.max_time;
    dt = params.dt;
    time_bins = 0:dt:max_time;
    h = zeros(size(time_bins));
    
    rx_pos = [0, r, 0];
    miu0 = [cos(psi1)*sin(theta1), sin(psi1)*sin(theta1), cos(theta1)];
    miur = [cos(psi2)*sin(theta2), sin(psi2)*sin(theta2), cos(theta2)];
    
    % 添加最大传播距离（基线距离的100倍）
    max_distance = 100 * r;

    for i = 1:num_photons
        photon = initialize_photon(miu0, phi1);
        scattering_count = 0; % 初始化散射计数器
        
        while scattering_count < max_scattering % 修改终止条件：基于散射次数
            % 计算步长
            delta_s = -log(rand) / k_s;
            
            % 更新位置
            new_pos = photon.pos + photon.direction * delta_s;
            photon.total_distance = photon.total_distance + delta_s;
            
            % 检查最大传播距离
            if photon.total_distance > max_distance
                break;
            end
            
            % 更新生存概率（只考虑吸收）
            photon.survival_prob = photon.survival_prob * exp(-k_a * delta_s);
            
            % 更新当前位置
            photon.pos = new_pos;
            
            % 检查是否在接收视场内
            if is_in_rx_fov(photon.pos, rx_pos, miur, phi2)
                % 计算到接收器的距离
                d = norm(photon.pos - rx_pos);
                
                % 更新总距离
                photon.total_distance = photon.total_distance + d;
                arrival_time = photon.total_distance / c;
                
                % 计算接收概率（修正计算）
                p1n = calculate_arrival_probability(photon.direction, photon.pos, rx_pos, ...
                    k_s_ray, k_s_mie, gamma, g, f, rx_area, miur);
                
                % 计算传播损失
                p2n = exp(-k_e * d);
                
                % 确保概率值在合理范围内
                arrival_prob = min([photon.survival_prob * p1n * p2n, 1.0]);
                
                % 记录到冲激响应
                bin_index = floor(arrival_time / dt) + 1;
                if bin_index <= length(h)
                    h(bin_index) = h(bin_index) + arrival_prob;
                end
                
                % 光子被接收，结束传播
                break;
            end
            
            % 如果未被接收，计算新的散射方向
            photon.direction = calculate_scattering_direction(photon.direction, ...
                k_s_ray, k_s_mie, gamma, g, f);
            
            % 增加散射计数器
            scattering_count = scattering_count + 1;
        end
    end
    
    % 归一化处理
    h = h / (num_photons * dt);
    
    % 计算路径损失
    path_loss = 1 / trapz(time_bins, h);
end

% ========== 参数设置和调用示例 ==========
% 设置参数（基于论文Table I和实验条件）
params = struct();
params.lambda = 260; % 波长 (nm)
params.k_s_ray = 0.266e-3; % Rayleigh散射系数 (m^-1) 从表I转换
params.k_s_mie = 0.284e-3; % Mie散射系数 (m^-1)
params.k_a = 0.802e-3; % 吸收系数 (m^-1)

params.phi1 = deg2rad(17); % 发射光束全宽发散角 (rad)
params.phi2 = deg2rad(30); % 接收视场角 (rad)

params.theta1 = deg2rad(45); % 发射天顶角 (rad)
params.theta2 = deg2rad(45); % 接收天顶角 (rad)

params.psi1=deg2rad(90); % 光源方位角(rad)
params.psi2=2*pi-deg2rad(90); % 接收机方位角(rad)

params.r = 100; % 基线距离 (m)
params.rx_area = 1.77e-4; % 接收面积 (m^2) (1.77 cm^2)
params.gamma = 0.017; % Rayleigh散射参数
params.g = 0.72; % Mie散射不对称参数
params.f = 0.5; % Mie散射参数
params.num_photons = 1e5; % 光子数量
params.max_scattering = 3; % 新增：最大散射次数
params.max_time = 2e-6; % 最大模拟时间 (s) (10 us)
params.dt = 2e-8; % 时间分辨率 (s) (1 ns)
params.c = 2.997046e8; % 光速 (m/s)

% 运行蒙特卡洛模拟
[h, path_loss] = monte_carlo_uv_channel(params);

% 显示结果
fprintf('路径损失: %.2f dB\n', 10*log10(path_loss));

% 绘制冲激响应
time_axis = (0:length(h)-1) * params.dt;
figure;
plot(time_axis * 1e6, h); % 时间单位转换为us
xlabel('时间 (\mus)');
ylabel('强度 (W/m^2)');
title('NLOS紫外散射信道冲激响应');
grid on;



%%
% 初始化光子(根据光源束散角phi1和光源顶角theta1,方位角2π内均匀分布)
function photon = initialize_photon(miu0,phi1)
    photon.pos = [0, 0, 0]; % 发射器位置
    photon.survival_prob = 1.0;%初始生存概率为1
    photon.total_distance = 0;%初始走过的路程为0
    
    % 生成初始方向（在光束内均匀分布）(公式1)
    xi_theta = rand;
    xi_psi = rand;
    U = acos(1 - xi_theta * (1 - cos(phi1/2)));
    psi_ini = 2 * pi * xi_psi;
    % theta0=U+pi/2-deg2rad(70);
    % psi0=psi_ini;
    % photon.direction=[sin(theta0)*cos(psi0), sin(theta0)*sin(psi0), cos(theta0)];
    photon.direction=rotate_direction(miu0,U,psi_ini);
end

% 检查点是否在接收视场内
function in_fov = is_in_rx_fov(pos, rx_pos, miur, phi2)
    % 计算从接收器到点的向量
    vec = pos - rx_pos;
    distance = norm(vec);

    %计算向量和接收器主轴方向向量之间的cos值
    cos_theta_n=dot(vec,miur)/distance;
    
    % 检查是否在视场角内,余弦值大代表角度小
    in_fov = cos_theta_n>=cos(phi2/2);
end

% 修正接收概率计算函数，这已经算是最后一次散射了，但MCS仿真思路是接收概率还不够小就可以继续散射
function p1n = calculate_arrival_probability(photon_direction, pos, rx_pos, ...
    k_s_ray, k_s_mie, gamma, g, f, rx_area, miur)
    
    rn = rx_pos - pos;
    distance_rn = norm(rn);
    
    % 防止除以零,此时可以认为被接收，所以p1n=1
    if distance_rn < 1e-6
        p1n = 1;
        return;
    end
    
    % 归一化向量
    rn_unit = rn / distance_rn;
    
    %cos_theta_s：朝接收器中心需要偏转的俯仰角吧（我用惯导里的术语了）
    cos_theta_s = dot(rn_unit, photon_direction);
    cos_theta_omiga = dot(rn_unit, miur);
    
    % 确保余弦值在有效范围内
    cos_theta_s = max(min(cos_theta_s, 1), -1);
    cos_theta_omiga = max(min(cos_theta_omiga, 1), -1);
    
    % 计算散射相函数（修正计算）
    w_ray = k_s_ray / (k_s_ray + k_s_mie);
    w_mie = 1 - w_ray;
    
    % Rayleigh相位函数
    P_ray = (3*(1+3*gamma+(1-gamma)*cos_theta_s^2)) / (16*pi*(1+2*gamma));
    
    % Mie相位函数（修正计算）
    hg_term = (1-g^2) / (4*pi) / (1+g^2-2*g*cos_theta_s)^(3/2);
    correction_term = (1-g^2)/(4*pi) *f * (0.5*(3*cos_theta_s^2-1)) / ((1+g^2)^(3/2));
    P_mie = hg_term + correction_term;
    
    % 组合相位函数（确保非负）
    P = w_ray * P_ray + w_mie * P_mie;
    P = max(P, 0);  % 确保非负
    
    % 计算立体角（修正计算）
    solid_angle = rx_area * abs(cos_theta_omiga) / (distance_rn^2);
    solid_angle = min(solid_angle, pi);  % 限制最大立体角

    % 好像缺了个sin_theta_s
    %sin_theta_s=sqrt(1-cos_theta_s^2);
    sin_theta_s=1;
    
    % 计算接收概率
    p1n = P * solid_angle/sin_theta_s;%前面sin_theta_s=1，所以相当于没除
    
    % 确保概率值在合理范围内
    p1n = min(max(p1n, 0), 1.0);
end

% 计算散射方向 (公式4-8)
function new_direction = calculate_scattering_direction(direction, k_s_ray, k_s_mie, gamma, g, f)
    % 计算相位函数权重
    w_ray = k_s_ray / (k_s_ray + k_s_mie);
    w_mie = 1 - w_ray;
    
    % 生成散射角（根据相位函数）(公式7)
    xi_mu = rand;
    mu = solve_phase_function(xi_mu, w_ray, w_mie, gamma, g, f);
    theta_s = acos(mu);
    
    % 方位角均匀分布
    psi_s = 2 * pi * rand;
    
    % 将散射方向转换到全局坐标系 (公式8)
    new_direction = rotate_direction(direction, theta_s, psi_s);
end

% 解相位函数方程（数值解）(公式5-7)
function mu = solve_phase_function(xi_mu, w_ray, w_mie, gamma, g, f)
    % 定义相位函数 (公式4-6)
    P = @(mu) (w_ray * (3*(1+3*gamma+(1-gamma)*mu.^2))/(16*pi*(1+2*gamma)) + ...
               w_mie * (1-g^2)/(4*pi) * (1./((1+g^2-2*g*mu).^(3/2)) + ...
               f*(0.5*(3*mu.^2-1))/((1+g^2).^(3/2))));
    
    % 数值求解 (公式7)
    options = optimset('Display','off');
    mu = fzero(@(m) 2*pi*integral(@(u) P(u), -1, m) - xi_mu, [-1 1], options);
end

% 旋转方向向量（根据散射角和方位角）(公式8)
function new_dir = rotate_direction(dir, theta_s, psi_s)
    % 当前方向余弦
    mu_x = dir(1);
    mu_y = dir(2);
    mu_z = dir(3);
    
    % 处理数值稳定性问题
    denom = sqrt(1 - mu_z^2);
    if denom < 1e-10
        % 当原方向接近z轴时的特殊处理
        if mu_z > 0
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), cos(theta_s)];
        else
            new_dir = [sin(theta_s)*cos(psi_s), sin(theta_s)*sin(psi_s), -cos(theta_s)];
        end
    else
        % 常规计算
        sin_theta = sin(theta_s);
        cos_theta = cos(theta_s);
        cos_psi = cos(psi_s);
        sin_psi = sin(psi_s);
        
        new_dir_x = sin_theta/denom * (mu_x*mu_z*cos_psi - mu_y*sin_psi) + mu_x*cos_theta;
        new_dir_y = sin_theta/denom * (mu_y*mu_z*cos_psi + mu_x*sin_psi) + mu_y*cos_theta;
        new_dir_z = -sin_theta*cos_psi*denom + mu_z*cos_theta;
        
        new_dir = [new_dir_x, new_dir_y, new_dir_z];
    end
    
    % 归一化
    new_dir = new_dir / norm(new_dir);
end


