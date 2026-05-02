g = 0.924;
theta = 0:0.0001*2*pi:pi;
HG=(1 - g^2) / 4./pi./(1 + g^2 - 2*g*cos(theta)).^(3/2).*sin(theta);
plot(theta/pi, HG);
% 使用梯形法计算积分
integral_value = 2*pi*trapz(theta, HG);
disp(['积分结果: ', num2str(integral_value)]);

