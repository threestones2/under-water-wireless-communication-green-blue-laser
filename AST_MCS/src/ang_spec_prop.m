    % Listing 6.5: MATLAB function for the angular-spectrum method
    % 文件名: ang_spec_prop.m
    
    function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
    % function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
    
        N = size(Uin, 1);   % 假设是正方形网格
        k = 2*pi/wvl;       % 光波数
    
        % 源平面坐标 (Source-plane coordinates)
        [x1, y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);
        r1sq = x1.^2 + y1.^2;
    
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