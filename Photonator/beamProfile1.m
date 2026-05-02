function [x,y,ux,uy,uz] = beamProfile1(n,beamWaist,diverg,~)

% 原始光束函数有误，已于2011年10月31日更新
% 更新方法：从x/y平面的高斯半径光线分布开始（无初始发散）
% 通过薄透镜光线矩阵变换进行发散：(x1 = x2, theta2 = (-1/f)x1 + theta1)
% 选择1/f值使得当光束腰斑处的光线通过发散透镜时产生所需发散角

%beamWaist：光束腰斑半径（光束最窄处的半径）
%diverg：光束发散角（弧度）

% 生成[0,1]均匀分布随机数
randVals = rand(n,1);

% 计算等效薄透镜焦距倒数（负值表示发散透镜）
invF = -diverg/beamWaist;

% 生成符合高斯光束强度分布的半径（概率密度与exp(-r^2)成正比）
radius = beamWaist * sqrt(-log(1-randVals)); 

% 通过薄透镜变换计算极化偏折角
divAng = -invF .* radius;         

% 生成均匀分布的方位角 [0, 2π]
phiAng = (2*pi) * rand(n,1);       	

% 预计算方位角的三角函数值（提高效率）
cosPhi = cos(phiAng);
sinPhi = sin(phiAng);

% 计算方向向量的z分量（与光轴夹角为divAng）
uz = cos(divAng);
% 计算垂直于光轴的平面内分量幅值
sin_uz = sqrt(1-uz.^2);

% 计算方向向量的x和y分量
ux = sin_uz .* cosPhi;
uy = sin_uz .* sinPhi;

% 计算光线在腰斑平面上的位置坐标
x = radius .* cosPhi;
y = radius .* sinPhi;

% ====== 以下为被注释的错误旧实现 ======
% sd = 1/e光束半径，即E(b) = (1/e)E(0)
% muD = cos(diverg);
% randVals = rand(n,1);
% ang = (2*pi) * rand(n,1);
% cosPhi = cos(ang);
% sinPhi = sin(ang);
% x = r .* cosPhi;
% y = r .* sinPhi;
% uz = 1 - rand(n,1).*(1 - muD);  % 在[1, muD]区间随机采样，[0,acos(muD)]上各向同性分布？ <- 错误方法！

% ====== 分析工具（函数外调用） ======
% 使用 hist3([x y],[100 100]) 可查看光线在xy平面的二维分布直方图
% 正确结果应表现为中心亮、边缘暗的高斯分布

end
