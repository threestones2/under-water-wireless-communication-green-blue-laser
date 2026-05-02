% mc_func_Berrocal.m

% 功能说明:
% 蒙特卡洛模拟函数，用于追踪光子在散射介质中的传播。
% 它计算被接收器探测到的光子的统计信息，如总功率、数量、位置和传播距离。

% 输入参数:
% num_photons: 模拟中使用的光子总数
% scattering_events: (在此版本中未使用) 每次散射的事件数
% c: 介质的总衰减系数 (c = a + b)
% a: 介质的吸收系数
% receiver_z: 接收器平面所在的Z轴坐标
% cdf_scatter: 散射相函数的累积分布函数 (CDF) 查找表
% angle: 与CDF对应的散射角数组
% init_angle, init_angle2: (在此版本中未使用) 发射器的初始指向角
% diverg: (在此版本中未使用) 光束的发散角

% 输出参数:
% total_time: 模拟运行的总时间
% total_rec_power: 接收到的总功率 (所有被接收光子的权重之和)
% total_rec_packets: 接收到的光子总数
% rec_loc_final: 被接收光子的最终位置和入射方向 [x, y, ux, uy, uz]
% total_rec_dist: 每个被接收光子的总传播距离
% rec_weights: 每个被接收光子的最终权重

function [total_time,total_rec_power,total_rec_packets,rec_loc_final,total_rec_dist,rec_weights] = ...
    mc_func_Berrocal1(num_photons,scattering_events,c,a,receiver_z,...
    cdf_scatter,angle,init_angle,init_angle2,diverg)

% --- 初始化阶段 ---

useLimits = 'true'; % 是否启用模拟区域边界限制

if strcmp(useLimits,'true')
    xLimMax = 0.005000;  % X轴正向边界
    xLimMin = -0.005;    % X轴负向边界
    yLimMax = 0.005;     % Y轴正向边界
    yLimMin = -0.005;    % Y轴负向边界
    zLimMax = 0.01;      % Z轴正向边界 (代码逻辑中未使用)
    zLimMin = -0.005;    % Z轴负向边界
end

% 设置光子状态矩阵
% 列定义: 1:x, 2:y, 3:z, 4:ux, 5:uy, 6:uz, 7:weight, 8:status
% status: 1 = 存活, 0 = 已被接收(探测到), -1 = 已被终止(飞出边界或轮盘赌淘汰)
photon = zeros(num_photons,8);
photon(:,7) = ones(num_photons,1);  % 初始化所有光子权重为 1
photon(:,8) = ones(num_photons,1);  % 初始化所有光子为存活状态

totaldist = zeros(num_photons,1);   % 记录每个光子走过的总距离
rec_dist = zeros(num_photons,1);    % 记录被接收光子在接收点的总传播距离
rec_loc = zeros(num_photons,5);     % 临时存储被接收光子的位置 [x, y] 和方向 [ux, uy, uz]
total_rec_packets = 0;              % 计数器：穿过探测器平面的光子总数
total_rec_power = 0;                % 累加器：穿过探测器平面的总功率

% 预计算常量以提高循环效率
prob_of_survival = ((c-a)/c);       % 单次散射反照率 (b/c), 即每次散射后光子存活的概率
rouletteConst = 10;                 % 轮盘赌常数. 低权重光子有 1/rouletteConst 的概率存活
rouletteConst_inv = 1/rouletteConst; % 轮盘赌常数的倒数
inv_c = 1/c;                        % 衰减系数的倒数, 用于计算光程
% inv_b = 1/(c-a);                  % (在此版本中未使用)

min_power = 1e-4;                   % 光子权重的阈值, 低于此值将进行轮盘赌

max_uz = 0.99999;                   % 用于处理方向计算中 uz 接近 1 的特殊情况, 避免除零

tic; % 开始计时

% --- 设置光子初始位置和方向 ---
[x,y] = startingDistribution(num_photons); % 调用一个外部函数来生成光束横截面上的初始位置分布

% 将光束设置为沿Z轴正方向传播
photon(:,1) = x.*xLimMax;           % 初始 x 坐标
photon(:,2) = y.*yLimMax;           % 初始 y 坐标
photon(:,3) = 0;                    % 初始 z 坐标 (假设发射器在 z=0 平面)
photon(:,6) = ones(num_photons,1);  % 初始方向 uz = 1

photonsRemaining = num_photons;     % 剩余待处理的光子数量

clear x y % 释放不再需要的变量内存

% --- 主循环 ---
% 持续循环直到所有光子都被处理完毕 (被接收或被终止)
while photonsRemaining > 0
    % 为每个光子批量生成该迭代步中所需的随机数, 提高效率
    rand_array = rand(num_photons,3);         % 生成 N行3列 的随机数矩阵
    rand_array(:,3) = rand_array(:,3).*2.*pi; % 将第3列随机数映射到 [0, 2*pi] 用于方位角 phi

    % 遍历每一个光子, 计算其新的位置和状态
    for i = 1:num_photons
        
        % 只处理仍然存活的光子
        if (photon(i,8) == 1)
        
            % 根据Beer-Lambert定律的逆变换采样, 随机生成一个光程 r
            r = -inv_c*log(rand_array(i,1));

            % --- 生成散射角 theta ---
            % 使用二分查找在累积分布函数(CDF)上找到对应的散射角
            minIndex = 2;
            maxIndex = length(cdf_scatter);
            
            % 核心的二分查找循环
            while maxIndex >= minIndex
                midIndex = minIndex + ceil((maxIndex - minIndex)/2);
                if rand_array(i,2) > cdf_scatter(midIndex)
                    minIndex = midIndex + 1;
                elseif rand_array(i,2) < cdf_scatter(midIndex)
                    maxIndex = midIndex - 1;
                else % 正好命中
                    break;
                end
            end
            % 确定最终的索引区间
            midIndex = minIndex + ceil((maxIndex - minIndex)/2);
            k = midIndex;

            % 通过线性插值得到更精确的散射角 theta
            theta = angle(k-1) + (rand_array(i,2) - cdf_scatter(k-1))*(angle(k) - angle(k-1))/(cdf_scatter(k) - cdf_scatter(k-1));
            % 从预生成的随机数中获取方位角 phi
            phi = rand_array(i,3);

            % 根据上一步的方向向量计算本次移动的位移
            x_step = r*photon(i,4);  % x 方向位移
            y_step = r*photon(i,5);  % y 方向位移
            z_step = r*photon(i,6);  % z 方向位移
                
            % --- 判断光子是否穿过接收器平面 ---
            if ((photon(i,3) + z_step) >= receiver_z)
          
                if photon(i,6) ~= 0 % 确保光子在Z轴上有速度分量
                    % 计算到接收器平面的Z轴距离
                    z_dist_rec_intersection = receiver_z - photon(i,3);
                    % 根据相似三角形原理计算在接收平面上的x, y交点
                    y_dist_rec_intersection = z_dist_rec_intersection*photon(i,5)/photon(i,6);
                    x_dist_rec_intersection = z_dist_rec_intersection*photon(i,4)/photon(i,6);
                else
                    % 理论上不应发生, Z方向速度为0的光子无法穿越Z平面
                    disp('光子Z轴速度为0, 如何穿过接收平面?');
                    disp(z_step); disp(photon(i,3));
                end

                % 到达接收器交点的欧几里得距离
                dist_to_rec = z_dist_rec_intersection / photon(i,6);

                % 记录在接收平面上的最终信息
                rec_loc(i,1) = photon(i,1) + x_dist_rec_intersection;   % 接收点 x 坐标
                rec_loc(i,2) = photon(i,2) + y_dist_rec_intersection;   % 接收点 y 坐标
                rec_loc(i,3) = photon(i,4);                             % 入射方向 ux
                rec_loc(i,4) = photon(i,5);                             % 入射方向 uy
                rec_loc(i,5) = photon(i,6);                             % 入射方向 uz

                total_rec_packets = total_rec_packets + 1;              % 接收到的光子数+1
                total_rec_power = total_rec_power + photon(i,7);        % 累加接收到的功率
                
                rec_dist(i) = totaldist(i) + dist_to_rec;               % 记录该光子的总传播距离
                photon(i,8) = 0;                                        % 标记光子为“已接收”
                photonsRemaining = photonsRemaining - 1;                % 剩余光子数-1
                
                % 更新该光子的总传播距离
                totaldist(i) = totaldist(i) + dist_to_rec;

            else % 如果光子没有到达探测器, 则移动它并进行散射
                
                % 更新光子位置
                photon(i,1) = photon(i,1) + x_step;
                photon(i,2) = photon(i,2) + y_step;
                photon(i,3) = photon(i,3) + z_step;
                
                % 检查光子是否飞出模拟区域边界
                if ((photon(i,1) > xLimMax) || (photon(i,1) < xLimMin) || (photon(i,2) > yLimMax) || (photon(i,2) < yLimMin) || (photon(i,3) < zLimMin))
                    photon(i,8) = -1;                           % 标记为“已终止”
                    photonsRemaining = photonsRemaining - 1;    % 剩余光子数-1
                else 	% 如果光子仍在边界内
               
                    % 因吸收而衰减权重
                    photon(i,7) = photon(i,7)*prob_of_survival;
                    
                    % 轮盘赌: 处理低权重光子
                    if  photon(i,7) < min_power
                        if rand() > rouletteConst_inv % 以 (1 - 1/m) 的概率终止光子
                            photon(i,8) = -1;                               
                            photonsRemaining = photonsRemaining - 1;
                        else % 以 1/m 的概率让光子存活, 并将其权重放大 m 倍
                            photon(i,7) = photon(i,7)*rouletteConst;
                        end
                    end
                    
                    % --- 更新散射后的方向 ---
                    if abs(photon(i,6)) > max_uz % 特殊情况：如果光子几乎沿Z轴传播
                        photon(i,4) = sin(theta)*cos(phi);
                        photon(i,5) = sin(theta)*sin(phi);
                        photon(i,6) = sign(photon(i,6))*cos(theta);
                    else % 一般情况：使用标准的旋转公式
                        sqrt_uz = sqrt(1 - photon(i,6)^2);
                        old_ux = photon(i,4);
                        old_uy = photon(i,5);
                        old_uz = photon(i,6);
                        % Henyey-Greenstein散射相函数的方向更新公式
                        photon(i,4) = (sin(theta)/sqrt_uz)*(old_ux*old_uz*cos(phi) - old_uy*sin(phi)) + old_ux*cos(theta);   % new ux
                        photon(i,5) = (sin(theta)/sqrt_uz)*(old_uy*old_uz*cos(phi) + old_ux*sin(phi)) + old_uy*cos(theta);   % new uy
                        photon(i,6) = (-sin(theta)*cos(phi))*sqrt_uz + old_uz*cos(theta);                                    % new uz
                    end
                end
                % 更新光子已传播的总距离
                totaldist(i) = totaldist(i) + r;
            end
        end
    end
end

total_time = toc; % 停止计时

% --- 后处理阶段 ---
% 将临时存储的接收数据整理到最终的输出数组中
% 预分配最终数组大小可以避免动态增长, 速度更快
rec_loc_final = ones(total_rec_packets,5);
j = 1;
for i = 1:num_photons
    if (photon(i,8) == 0) % 如果光子是被接收的
       rec_loc_final(j,:) = rec_loc(i,:); % 记录接收位置和角度
       j = j + 1;
    end
end

% 错误检查
if ((j-1) ~= total_rec_packets)
    disp('错误! 最终接收到的光子总数与计数器不符。');
    disp(sprintf('j = %d and total_rec_packets = %d',j, total_rec_packets));
end

j = 1;
total_rec_dist = zeros(total_rec_packets,1);
rec_weights = zeros(total_rec_packets,1);
total_rec_power = 0; % 重新精确计算总接收功率
for i = 1:num_photons    
   if (photon(i,8) == 0)
       total_rec_dist(j) = rec_dist(i);   % 存储每个被接收光子的传播距离
       rec_weights(j) = photon(i,7);      % 存储每个被接收光子的最终权重
       j = j + 1;       
   end
end

% 再次进行错误检查
if ((j-1) ~= total_rec_packets)
    disp('错误! 最终记录的距离条目数与计数器不符。');
    disp(sprintf('j = %d and total_rec_packets = %d',j, total_rec_packets));
end

end
