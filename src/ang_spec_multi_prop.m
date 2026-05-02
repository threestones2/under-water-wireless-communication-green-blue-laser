function [xn yn Uout] = ang_spec_multi_prop(Uin, wvl, delta1, deltan, z, t)

% Uin：源光场
% wvl：波数
% delta1：源平面网格大小
% deltan：

% --- 1. 初始化 ---
N = size(Uin, 1); 
[nx ny] = meshgrid((-N/2 : 1 : N/2 - 1)); 
k = 2*pi/wvl; 

% --- 2. 超高斯吸收边界 (Super-Gaussian Absorbing Boundary) ---
% 作用：因为使用了 FFT，数值计算隐含周期性边界。为了防止光束扩算到网格
% 边缘后“卷绕”回另一侧（Aliasing），在边缘强制让光强衰减为 0。
% 这里的 sg 是一个边缘为 0，中间为 1 的方形/圆形窗口。
nsq = nx.^2 + ny.^2; 
w = 0.47*N;  
sg = exp(-nsq.^8/w^16); clear('nsq', 'w');  

% --- 3. 纵向坐标与网格缩放设置 ---
z = [0 z]; % 在距离向量前加一个 0，方便计算 delta_z
n = length(z);  
Delta_z = z(2:n) - z(1:n-1); % 每一步的传播距离

% 核心：计算每一步的网格大小 (Variable Grid Spacing)
% 这种算法允许网格从 delta1 (源) 线性变化到 deltan (终)
% 这样可以动态适应发散或聚焦的光束，保持采样率。
alpha = z / z(n); 
delta = (1-alpha) * delta1 + alpha * deltan;  
m = delta(2:n) ./ delta(1:n-1); % 计算每一步的放大倍率 m

% 源平面物理坐标
x1 = nx * delta(1);  
y1 = ny * delta(1);  
r1sq = x1.^2 + y1.^2; 

% --- 4. 预补偿相位 Q1 ---
% 由于坐标系在变大/变小，相当于引入了额外的球面波因子。
% Q1 用于在传播开始前预先抵消一部分这种几何效应。
Q1 = exp(i*k/2*(1-m(1))/Delta_z(1)*r1sq);  

% 【关键点】应用初始相位屏 t(:,:,1)
Uin = Uin .* Q1 .* t(:,:,1);  

% --- 5. 分步传播循环 (Split-Step Loop) ---
for idx = 1 : n-1  
    % 计算当前切片的频率坐标
    deltaf = 1 / (N*delta(idx));  
    fX = nx * deltaf;  
    fY = ny * deltaf;  
    fsq = fX.^2 + fY.^2;  
    
    Z = Delta_z(idx); 
    
    % 传递函数 Q2 (含缩放因子 m)
    % 这是菲涅尔衍射传递函数的缩放版本
    Q2 = exp(-i*pi^2*2*Z/m(idx)/k*fsq);  
    
    % 【核心传播步骤】
    % 1. Uin / m(idx): 幅度归一化，保持能量守恒
    % 2. ft2 -> Q2 -> ift2: 真空传播 (Vacuum Propagation)
    % 3. t(:,:,idx+1): 叠加到达新平面后的湍流相位 (Phase Screen)
    % 4. sg: 吸收掉碰到边缘的光 (Boundary)
    Uin = sg .* t(:,:,idx+1) ...  
        .* ift2(Q2 ...  
        .* ft2(Uin / m(idx), delta(idx)), deltaf);  
end  

% --- 6. 后处理 ---
% 观察平面坐标
xn = nx * delta(n);  
yn = ny * delta(n);  
rnsq = xn.^2 + yn.^2; 

% 后补偿相位 Q3
% 去除由于网格缩放残留的球面相位因子，得到正确的物理波前。
Q3 = exp(i*k/2*(m(n-1)-1)/(m(n-1)*Z)*rnsq);  
Uout = Q3 .* Uin;