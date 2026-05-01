% 两种蒙特卡洛衰减计算方法的计算复杂度与耗时对比验证
clear; clc;

% ================= 1. 参数设置 =================
a = 0.366;          % 吸收系数 (m^-1)
b = 1.829;          % 散射系数 (m^-1)
c = a + b;          % 总衰减系数 (m^-1)
albedo = b / c;     % 反照率

L_target = 10;      % 设定目标传输总光程 (m)
N_photons = 1e6;    % 仿真光子总数

fprintf('--- 物理参数 ---\n');
fprintf('吸收系数 a: %.3f, 散射系数 b: %.3f, 衰减系数 c: %.3f\n', a, b, c);
fprintf('目标光程 L: %.2f m, 光子数: 1e6\n\n', L_target);

% ================= 2. 方法一：解耦法 (基于散射系数 b 游走) =================
tic;
weight_m1 = zeros(1, N_photons);
steps_m1_total = 0; % 记录总步数

for i = 1:N_photons
    path_len = 0;
    steps = 0;
    while true
        % 按照【散射系数 b】抽取步长
        step = -log(1 - rand()) / b; 
        if path_len + step >= L_target
            break;
        end
        path_len = path_len + step;
        steps = steps + 1;
    end
    steps_m1_total = steps_m1_total + steps;
    % 解析计算吸收衰减
    weight_m1(i) = exp(-a * L_target); 
end
time_m1 = toc;

% ================= 3. 方法二：耦合统计法 (基于总衰减系数 c 游走) =================
tic;
weight_m2 = zeros(1, N_photons);
steps_m2_total = 0; % 记录总步数

for i = 1:N_photons
    path_len = 0;
    steps = 0;
    while true
        % 按照【总衰减系数 c】抽取步长
        step = -log(1 - rand()) / c; 
        if path_len + step >= L_target
            break;
        end
        path_len = path_len + step;
        steps = steps + 1;
    end
    steps_m2_total = steps_m2_total + steps;
    % 统计计算概率衰减
    weight_m2(i) = albedo^steps;
end
time_m2 = toc;

% ================= 4. 耗时与步数统计结果 =================
fprintf('--- 耗时与运算量分析 ---\n');
fprintf('【方法一】(b-步长) 总耗时: %.4f 秒, 模拟总循环步数: %d\n', time_m1, steps_m1_total);
fprintf('【方法二】(c-步长) 总耗时: %.4f 秒, 模拟总循环步数: %d\n', time_m2, steps_m2_total);
fprintf('耗时比 (Method2/Method1): %.4f\n', time_m2 / time_m1);
fprintf('步数比 (Method2/Method1): %.4f\n', steps_m2_total / steps_m1_total);
fprintf('理论步数比 (c/b): %.4f\n', c / b);