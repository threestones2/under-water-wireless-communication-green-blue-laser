    % ========== 1. 参数设置 ==========
    params = struct();
    params.lambda = 514; 
    params.c = 2.25e8;   
    
    % 清澈海水参数
    params.c_att = 0.1514; 
    params.k_a = 0.114;    
    params.k_s = 0.0374;   
    
    % Ding (Mix) 相函数
    params.w_mie = 0.8440; 
    params.gamma = 0.017;  
    params.g_mie = 0.9814; 
    params.f_mie = 0.49;   
    
    % 几何配置
    params.r = 100;        
    params.phi1 = 1e-3;    
    params.phi2 = deg2rad(30); 
    params.rx_radius = 0.1;    
    
    params.tx_pos = [0, 0, 0];
    params.tx_dir = [0, 1, 0]; 
    params.rx_pos = [0, params.r, 0];
    params.rx_dir = [0, -1, 0]; 
    
    % 仿真控制
    params.num_photons = 1e5;   
    params.max_scattering = 200; 
    params.max_time = 1e-6;     
    params.dt = 1e-9;           
    
    % 预计算 CDF
    fprintf('正在预计算相函数 CDF 表...\n');
    [cdf_table, theta_table] = precompute_phase_function_cdf(params);
    
    % ========== 2. 主循环 ==========
    fprintf('开始仿真 (距离: %dm)... \n', params.r);
    
    time_bins = 0:params.dt:params.max_time;
    
    % [修改A]：初始化三个数组，分别记录总能量、直射、散射
    h_total = zeros(size(time_bins));
    h_direct = zeros(size(time_bins));
    h_scatter = zeros(size(time_bins));
    
    detected_count = 0;
    rx_area = pi * params.rx_radius^2;
    cos_fov_half = cos(params.phi2/2);
    max_dist_limit = 500; 

    tic;
    for i = 1:params.num_photons
        photon = initialize_photon(params.tx_pos, params.tx_dir, params.phi1);
        scattering_order = 0;
        
        while scattering_order <= params.max_scattering
            delta_s = -log(rand) / params.k_s; 
            
            [is_hit, dist_to_rx] = check_intersection(photon.pos, photon.dir, delta_s, ...
                params.rx_pos, params.rx_dir, params.rx_radius, cos_fov_half);
            
            if is_hit
                total_dist = photon.total_dist + dist_to_rx;
                arrival_time = total_dist / params.c;
                final_energy = photon.weight * exp(-params.k_a * dist_to_rx);
                
                bin_idx = floor(arrival_time / params.dt) + 1;
                if bin_idx <= length(h_total)
                    % [修改B]：根据散射阶数分类记录
                    h_total(bin_idx) = h_total(bin_idx) + final_energy;
                    
                    if scattering_order == 0
                        % 0阶散射 = 直射 (Ballistic)
                        h_direct(bin_idx) = h_direct(bin_idx) + final_energy;
                    else
                        % >0阶散射 = 散射 (Scattered)
                        h_scatter(bin_idx) = h_scatter(bin_idx) + final_energy;
                    end
                end
                
                detected_count = detected_count + 1;
                break; 
            end
            
            photon.pos = photon.pos + photon.dir * delta_s;
            photon.total_dist = photon.total_dist + delta_s;
            photon.weight = photon.weight * exp(-params.k_a * delta_s);
            
            if photon.weight < 1e-20 || photon.total_dist > max_dist_limit
                break;
            end
            
            xi = rand;
            theta_s = interp1(cdf_table, theta_table, xi, 'linear');
            phi_s = 2 * pi * rand;
            photon.dir = rotate_direction(photon.dir, theta_s, phi_s);
            
            scattering_order = scattering_order + 1;
        end
        
        if mod(i, params.num_photons/10) == 0
            fprintf('进度: %.0f%% \n', i/params.num_photons*100);
        end
    end
    toc;
    
    % ========== 3. 结果处理与绘图 (重点修改) ==========
    % 归一化
    norm_factor = params.num_photons * params.dt * rx_area;
    h_total = h_total / norm_factor;
    h_direct = h_direct / norm_factor;
    h_scatter = h_scatter / norm_factor;
    
    total_rx_energy = sum(h_total) * params.dt * rx_area * params.num_photons; % 还原回能量比例
    if total_rx_energy > 0
        path_loss_db = 10*log10(1 / sum(h_total * params.dt * rx_area)); % 简化的PL计算
    else
        path_loss_db = Inf;
    end
    
    fprintf('Path Loss: %.2f dB\n', path_loss_db);
    
    % [修改C]：使用 Semilogy 对数图来展示巨大差异
    figure;
    t_ns = time_bins * 1e9;
    
    % 为了防止 log(0) 出错，将 0 值替换为极小值 (仅用于绘图)
    min_val = 1e-25; 
    h_total(h_total==0) = min_val;
    h_direct(h_direct==0) = min_val;
    h_scatter(h_scatter==0) = min_val;
    
    % 绘图
    semilogy(t_ns, h_total, 'k-', 'LineWidth', 0.5, 'DisplayName', 'Total Response');
    hold on;
    semilogy(t_ns, h_direct, 'r--', 'LineWidth', 0.5, 'DisplayName', 'Ballistic (LOS)');
    semilogy(t_ns, h_scatter, 'b:', 'LineWidth', 1.0, 'DisplayName', 'Scattered (NLOS)');
    
    xlabel('Time (ns)');
    ylabel('Impulse Response (W/m^2) [Log Scale]');
    title(['UV/Optical Channel Components (r=' num2str(params.r) 'm)']);
    legend('show', 'Location', 'northeast');
    grid on;
    
    % 设置坐标轴范围
    t_los = params.r / params.c * 1e9; % 直射到达时间 (~444ns)
    xlim([t_los - 20, t_los + 100]); % 聚焦直射附近
    
    % 自动调整Y轴范围，确保能看到散射分量
    max_val = max(h_total);
    ylim([max_val * 1e-8, max_val * 10]); % 显示 8 个数量级的动态范围

% ========== 辅助函数 (保持不变，请确保包含) ==========
function photon = initialize_photon(pos, main_dir, divergence)
    photon.pos = pos; photon.total_dist = 0; photon.weight = 1.0;
    xi_theta = rand; xi_psi = rand;
    cos_val = 1 - xi_theta * (1 - cos(divergence/2));
    theta_gen = acos(cos_val); psi_gen = 2 * pi * xi_psi;
    dir_local = [sin(theta_gen)*cos(psi_gen), sin(theta_gen)*sin(psi_gen), cos(theta_gen)];
    up = [0, 0, 1]; if abs(dot(main_dir, up)) > 0.99, up = [1, 0, 0]; end
    u = cross(up, main_dir); u = u / norm(u); v = cross(main_dir, u);
    photon.dir = dir_local(1)*u + dir_local(2)*v + dir_local(3)*main_dir;
    photon.dir = photon.dir / norm(photon.dir);
end

function [is_hit, dist] = check_intersection(pos, dir, step_len, rx_c, rx_n, rx_r, cos_fov)
    is_hit = false; dist = 0;
    denom = dot(rx_n, dir);
    if denom > -1e-6, return; end 
    vec = rx_c - pos; t = dot(vec, rx_n) / denom;
    if t <= 1e-9 || t > step_len, return; end
    p_int = pos + t * dir;
    if norm(p_int - rx_c) > rx_r, return; end
    if -denom < cos_fov, return; end
    is_hit = true; dist = t;
end

function new_dir = rotate_direction(dir, theta, psi)
    u = dir(1); v = dir(2); w = dir(3); denom = sqrt(1 - w^2);
    if denom < 1e-10
        if w > 0, new_dir = [sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta)];
        else, new_dir = [sin(theta)*cos(psi), sin(theta)*sin(psi), -cos(theta)]; end
    else
        st = sin(theta); ct = cos(theta); sp = sin(psi); cp = cos(psi);
        nx = st/denom * (u*w*cp - v*sp) + u*ct;
        ny = st/denom * (v*w*cp + u*sp) + v*ct;
        nz = -st*cp*denom + w*ct;
        new_dir = [nx, ny, nz];
    end
    new_dir = new_dir / norm(new_dir);
end

function [cdf, theta] = precompute_phase_function_cdf(p)
    N = 2000; theta = linspace(0, pi, N); cos_theta = cos(theta);
    k_b = p.k_s; k_r = k_b * (1 - p.w_mie); k_m = k_b * p.w_mie;
    pdf_vals = zeros(1, N);
    for i = 1:N
        ct = cos_theta(i);
        p_ray = 3 * (1 + 3*p.gamma + (1-p.gamma)*ct^2) / (16 * pi * (1 + 2*p.gamma));
        term1 = 1 / (1 + p.g_mie^2 - 2*p.g_mie*ct)^1.5;
        term2 = p.f_mie * 0.5 * (3*ct^2 - 1) / (1 + p.g_mie^2)^1.5;
        p_mie = (1 - p.g_mie^2)/2 * (term1 + term2) / (2*pi); 
        pdf_vals(i) = (k_r * p_ray + k_m * p_mie) / k_b;
    end
    integrand = pdf_vals .* 2 * pi .* sin(theta);
    cdf = cumtrapz(theta, integrand);
    cdf = cdf / cdf(end); cdf(1) = 0; cdf(end) = 1;
    [cdf, unique_idx] = unique(cdf); theta = theta(unique_idx);
end