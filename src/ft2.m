% ft2.m
function F = ft2(g, delta)
% F = ft2(g, delta)
%   计算2D傅里叶变换
%   g: 输入的二维空间域函数
%   delta: 空间采样间隔 (例如，delta_x = delta_y = delta)
%   F: 输出的二维频域函数 (中心化)

% 确保输入是双精度
g = double(g);

% 将零频率分量移到矩阵的左上角，以便fft2正确处理
% 然后进行2D FFT
% 再将零频率分量移回中心
F = fftshift(fft2(ifftshift(g)));

% 根据离散傅里叶变换的定义，通常需要一个尺度因子
% 对于连续傅里叶变换的近似，这个因子是 (delta * delta)
F = F * delta^2; 

end
