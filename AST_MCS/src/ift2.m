function g = ift2(F_CFT_approx, delta_f)
    N=size(F_CFT_approx,1);
    delta_d=1/(N*delta_f);% 空间分辨率
    F_DFT = F_CFT_approx / delta_d^2; % 把 CFT 近似值转换为 DFT 傅里叶系数
    g = fftshift(ifft2(ifftshift(F_DFT))); % MATLAB 的 ifft2 已经包含了 1/N^2 因子
end
