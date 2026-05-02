% 经验 Petzold 相函数，源自论文:
% THEORETICAL AND EMPIRICAL PHASE FUNCTIONS FOR MONTE CARLO CALCULATIONS OF
% LIGHT SCATTERING IN SEAWATER by Vladimir I. Haltrin

% 定义散射角度 theta (单位: 度)
% 在 0-10 度范围内使用高分辨率，在 10-180 度范围内使用较低分辨率
theta = [[0:0.01:10] [10.1:0.1:180]]; 
% 备选方案1: 使用均匀间隔定义 theta
% theta = [0:180/5000:180];
% 备选方案2: 使用弧度制定义 theta
% theta = [0:pi/5000:pi];
 

% % 港口水体参数
% c = 2.190;
% a = 0.366;
% % 近岸水体参数 (当前启用)
% c = 0.22 + 0.179;
% a = 0.179;
% % 清澈水体参数
% c = 0.0374 + 0.114;
% a = 0.114;
% % Petzold 清澈水体模型 1
c = 0.199;
a = 0.082;
% % Petzold 近岸水体模型 1
% c = 0.470;
% a = 0.195;


% 计算散射系数 b, b = c - a (衰减 = 吸收 + 散射)
b = c-a;
% 计算单次散射反照率 albedo (散射系数 / 衰减系数, b/c, 无量纲)
albedo = (c-a)/c;


% 根据 Haltrin 经验公式计算中间系数
q = 2.598 + 17.748*sqrt(b) - 16.722*b + 5.932*b*sqrt(b);
k1 = 1.188 - 0.688*albedo;
k2 = 0.1*(3.07 - 1.90*albedo);
k3 = 0.01*(4.58 - 3.02*albedo);
k4 = 0.001*(3.24 - 2.25*albedo);
k5 = 0.0001*(0.84 - 0.61*albedo);


% 计算经验相函数
phase_func = (4*pi/b)*exp(q*(1 + ((-1)^1*k1*theta.^(1/2)) ...
                                + ((-1)^2*k2*theta.^(2/2)) ...
                                + ((-1)^3*k3*theta.^(3/2)) ...
                                + ((-1)^4*k4*theta.^(4/2)) ...
                                + ((-1)^5*k5*theta.^(5/2))));
                            
% 计算体积散射函数 beta
beta = (b/(4*pi)).*phase_func;
% 对相函数进行加权（通常用于绘图或某些计算），此处未使用
%phase_func_adj = phase_func.*sind(theta);

%%
% 归一化方程校验：积分结果应接近 1
% 0.5 * integral(p(theta)*sin(theta)*d(theta)) from 0 to pi 应该等于 1
normalization_eq = 0.5.*trapz(theta,(phase_func.*sin(theta*pi/180)));
% 备选方案: 使用 sind 函数，功能相同
% normalization_eq = 0.5.*trapz(theta,(phase_func.*sind(theta)))
  
%%
% 计算累积分布函数 (CDF)，用于蒙特卡洛抽样
cdf_phase_func = cumtrapz(theta,phase_func.*sind(theta));
% 将 CDF 归一化到 [0, 1] 区间
cdf_phase_func = cdf_phase_func ./ max(cdf_phase_func);
% 通过插值计算 CDF 的反函数，生成一个查找表
% 这允许我们通过一个 [0,1] 上的均匀随机数来查找对应的散射角 theta
cdf_phase_func_interp = interp1(cdf_phase_func,theta,[0:0.001:1]);

% 绘图命令 (取消注释以启用)
hold on;    
% 在对数-对数坐标系下绘制相函数，这是可视化相函数的常用方法
plot(theta,phase_func);
