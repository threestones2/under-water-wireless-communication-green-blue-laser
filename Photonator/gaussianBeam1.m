% beam_histogram.m
%
% 一个 MATLAB 程序
%   由 Steven L. Jacques 于 1998年12月编写
%   由 AI 助手于 2025年10月使用 histogram 函数进行现代化修改
%
% 本程序生成与激光束相关的二维高斯光束的
% 径向概率密度函数 p(r) 和相对辐照度 E(r) [单位: mm^-2]。
%
% 程序通过模拟10,000个随机光子进行统计，
% 并将统计直方图与理论解析表达式进行比较。
%
% 程序将生成两个图像：p(r) 和 E(r) 的对比图。

% --- 初始化环境 ---
close all;  % 关闭所有已打开的图形窗口
clear;      % 清除工作空间中的所有变量

% --- 参数定义 ---
b = 1;                       % 定义高斯光束的 1/e 辐照度半径为 1 [mm] (可任意选择)
N = 10000;                   % 用于模拟的随机光子总数

% --- 蒙特卡洛模拟核心 ---
% 1. 生成 N 个在 [0, 1] 区间内均匀分布的随机数
rnd1 = rand(1, N);
% 2. 使用“逆变换采样”方法，将均匀分布的随机数转换为符合高斯径向概率分布的 N 个光子落点半径 r1
%    推导公式: r = b * sqrt(-log(1 - u)), 其中 u 是[0,1]的随机数, log()是自然对数ln()
r1 = b * sqrt(-log(1 - rnd1));

% --- 数据统计与处理 (使用 histogram 函数) ---
nb = 30;                     % 设置直方图的“箱子”数量

% 创建一个直方图对象，并直接进行“概率密度函数(pdf)”归一化。
% 'Normalization','pdf' 会自动完成概率密度计算，使得整个直方图面积为1。
% 这完全等价于原代码中繁琐的手动归一化步骤 (n/N/dx)。
h = histogram(r1, nb, 'Normalization', 'pdf');

% 从直方图对象 h 中提取模拟结果
p_r_sim = h.Values;          % 模拟得到的径向概率密度 p(r)
edges = h.BinEdges;          % 获取每个箱子的边界位置
% 通过箱子边界计算每个箱子的中心位置，用于绘图
x_centers = edges(1:end-1) + h.BinWidth/2;

% 从模拟得到的 p(r) 反算出模拟的辐照度 E(r)
% 理论关系: p(r) = E(r) * 2*pi*r  =>  E(r) = p(r) / (2*pi*r)
E_r_sim = p_r_sim ./ (2 * pi * x_centers);

% --- 理论曲线计算 ---
% 创建一个平滑的、用于绘制理论曲线的半径数组 r
r_theory = (0:0.05*b:3*b);
% 计算理论上的辐照度 E(r)，这是一个标准的二维高斯分布
% 除以 (pi*b^2) 是为了归一化，使得 E(r) 在整个二维平面积分后总能量为1
E_theory = exp(-r_theory.^2/b^2) / (pi*b^2);
% 计算理论上的径向概率密度 p(r)
p_theory = E_theory .* 2 .* pi .* r_theory;

% --- 绘图：比较模拟结果与理论曲线 ---

% 图1：径向概率密度 p(r)
figure;
plot(x_centers, p_r_sim, 'bo', 'DisplayName', '模拟数据'); % 'bo' 表示蓝色圆圈
hold on;
plot(r_theory, p_theory, 'r-', 'LineWidth', 1.5, 'DisplayName', '理论曲线'); % 'r-' 表示红色实线
hold off;
title(['径向概率密度 p(r) (b = ' num2str(b) ' mm)']);
xlabel('半径 r [mm]');
ylabel('概率密度 p(r)');
legend; % 显示图例
grid on; % 显示网格

% 图2：相对辐照度 E(r)
figure;
plot(x_centers, E_r_sim, 'bo', 'DisplayName', '模拟数据'); % 'bo' 表示蓝色圆圈
hold on;
plot(r_theory, E_theory, 'r-', 'LineWidth', 1.5, 'DisplayName', '理论曲线'); % 'r-' 表示红色实线
hold off;
title(['相对辐照度 E(r) (b = ' num2str(b) ' mm)']);
xlabel('半径 r [mm]');
ylabel('概率密度 p(r)');
legend; % 显示图例
grid on; % 显示网格
