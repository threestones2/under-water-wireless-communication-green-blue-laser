function [x2, y2, Uout] = ang_spec_prop_circular_gaussian(Uin, wvl, d1, d2, Dz, w0, D_ap)
% function [x2, y2, Uout] = ang_spec_prop_circular_gaussian(Uin, wvl, d1, d2, Dz, w0, D_ap)
% 
% 修改说明：
% 1. 增加了 w0 (高斯束腰半径) 和 D_ap (圆形光阑直径) 作为输入参数。
% 2. 在传播前对 Uin 施加高斯强度分布和圆形截断。
%
% Uin: 输入光场 (通常输入 ones(N) 或初始平面波即可)
% wvl: 波长
% d1: 源平面网格间距
% d2: 观察平面网格间距
% Dz: 传播距离
% w0: 高斯光束束腰半径 (设为 [] 或 Inf 则不应用高斯)
% D_ap: 圆形光阑直径 (设为 [] 或 Inf 则不应用光阑)

    N = size(Uin, 1);   % 假设是正方形网格
    k = 2*pi/wvl;       % 光波数

    % 源平面坐标 (Source-plane coordinates)
    [x1, y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);
    r1sq = x1.^2 + y1.^2;

    % ==========================================================
    % [新增部分] 光束整形：高斯波束 + 圆形光阑
    % ==========================================================
    
    % 1. 应用高斯光束轮廓 (Gaussian Beam Profile)
    % 公式: exp( -(x^2+y^2) / w0^2 )
    if ~isempty(w0) && ~isinf(w0)
        gaussian_profile = exp(-r1sq / w0^2);
        Uin = Uin .* gaussian_profile;
    end

    % 2. 应用圆形光阑 (Circular Aperture)
    % 如果点在半径 D_ap/2 之外，置为 0
    if ~isempty(D_ap) && ~isinf(D_ap)
        mask = (r1sq <= (D_ap/2)^2);
        Uin = Uin .* mask; % 也就是将圆形外的部分“切掉”
    end
    % ==========================================================

    % 源平面空间频率 (Spatial frequencies of source plane)
    df1 = 1 / (N*d1);
    [fX, fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
    fsq = fX.^2 + fY.^2;

    % 缩放参数 (Scaling parameter)
    m = d2/d1;

    % 观察平面坐标 (Observation-plane coordinates)
    [x2, y2] = meshgrid((-N/2 : 1 : N/2 - 1) * d2);
    r2sq = x2.^2 + y2.^2;

    % 二次相位因子 (Quadratic phase factors)
    % 对应公式 6.65 中的 Q1, Q2, Q3
    Q1 = exp(1i * k/2 * (1-m)/Dz * r1sq);
    Q2 = exp(-1i * pi^2 * 2 * Dz/m/k * fsq);
    Q3 = exp(1i * k/2 * (m-1)/(m * Dz) * r2sq);

    % 计算传播后的场 (Compute the propagated field)
    % ft2 和 ift2 是本书自定义的傅里叶变换函数（含 fftshift）
    Uout = Q3 .* ift2(Q2 .* ft2(Q1 .* Uin / m, d1), df1);

end