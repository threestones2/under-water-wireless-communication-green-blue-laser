%% exp10_4: 三算法对比分析与绘图
%  比较 WCI-MC (无湍流) vs SA-MCS (半解析) vs MCS Physical (参考真值)
%  评估收敛速度、精度、计算效率

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= 加载数据 =================
data_wci  = load('data_exp10_WCIMC_None.mat');
data_sa   = load('data_exp10_SA_MCS.mat');
data_phys = load('data_exp10_MCS_Physical.mat');

dist_cell = data_wci.dist_cell;
water_types = data_wci.water_types;
num_W = length(water_types);
N_arr = data_wci.N_packets_arr;
num_N = length(N_arr);

fprintf('=== exp10_4: 三算法对比分析 ===\n');
fprintf('参考真值: MCS Physical (N = %d)\n', data_phys.N_packets);
fprintf('WCI-MC / SA-MCS 光子数扫描: [%s]\n\n', num2str(N_arr));

% ================= 创建图形 =================
N_labels = {'10^3', '10^4', '10^5', '10^6', '10^7'};
algo_names = {'WCI-MC', 'SA-MCS', 'Physical (Ref)'};
bar_colors = {[0.2 0.4 0.8], [0.8 0.3 0.3], [0.3 0.3 0.3]};

for w_idx = 1:num_W
    dists = dist_cell{w_idx};
    num_D = length(dists);

    % --- 图1: 路径损耗柱状图 ---
    figure('Position', [50, 100, 320 * num_D, 500]);
    for d_idx = 1:num_D
        subplot(1, num_D, d_idx);
        L = dists(d_idx);
        PL_ref = data_phys.PL_Cell{w_idx}(d_idx);

        % 构建柱状图数据: [N x 3] — WCI / SA / Ref
        PL_data = zeros(num_N, 3);
        PL_data(:, 1) = -data_wci.PL_Cell{w_idx}(d_idx, :);
        PL_data(:, 2) = -data_sa.PL_Cell{w_idx}(d_idx, :);
        PL_data(:, 3) = -PL_ref * ones(1, num_N);  % 参考值扩展

        b = bar(PL_data, 'grouped');
        for a = 1:3
            b(a).FaceColor = bar_colors{a};
            b(a).EdgeColor = 'none';
        end

        set(gca, 'XTickLabel', N_labels);
        xlabel('光子数 N'); ylabel('Path Loss (dB)');
        title(sprintf('%s  L=%dm', water_types{w_idx}, L));
        if d_idx == num_D, legend(algo_names, 'Location', 'best'); end
        grid on;
    end
    sgtitle(sprintf('%s — 路径损耗对比 (柱状图)', water_types{w_idx}));
    saveas(gcf, sprintf('fig_PL_bar_%s.png', water_types{w_idx}));

    % --- 图2: 绝对误差柱状图 ---
    figure('Position', [50, 100, 320 * num_D, 500]);
    for d_idx = 1:num_D
        subplot(1, num_D, d_idx);
        L = dists(d_idx);
        PL_ref = data_phys.PL_Cell{w_idx}(d_idx);

        err_data = zeros(num_N, 2);
        err_data(:, 1) = abs(data_wci.PL_Cell{w_idx}(d_idx, :) - PL_ref);
        err_data(:, 2) = abs(data_sa.PL_Cell{w_idx}(d_idx, :) - PL_ref);

        b = bar(err_data, 'grouped');
        for a = 1:2
            b(a).FaceColor = bar_colors{a};
            b(a).EdgeColor = 'none';
        end

        set(gca, 'XTickLabel', N_labels);
        xlabel('光子数 N'); ylabel('|PL - PL_{ref}| (dB)');
        title(sprintf('%s  L=%dm', water_types{w_idx}, L));
        if d_idx == num_D, legend(algo_names(1:2), 'Location', 'best'); end
        grid on;
    end
    sgtitle(sprintf('%s — 绝对误差 vs 参考真值', water_types{w_idx}));
    saveas(gcf, sprintf('fig_err_bar_%s.png', water_types{w_idx}));

    % --- 图3: 计算时间柱状图 ---
    figure('Position', [50, 100, 320 * num_D, 500]);
    for d_idx = 1:num_D
        subplot(1, num_D, d_idx);
        L = dists(d_idx);

        time_data = zeros(num_N, 2);
        time_data(:, 1) = data_wci.Time_Cell{w_idx}(d_idx, :);
        time_data(:, 2) = data_sa.Time_Cell{w_idx}(d_idx, :);

        b = bar(time_data, 'grouped');
        for a = 1:2
            b(a).FaceColor = bar_colors{a};
            b(a).EdgeColor = 'none';
        end

        set(gca, 'XTickLabel', N_labels);
        xlabel('光子数 N'); ylabel('计算时间 (s)');
        title(sprintf('%s  L=%dm', water_types{w_idx}, L));
        if d_idx == num_D, legend(algo_names(1:2), 'Location', 'best'); end
        grid on;
    end
    sgtitle(sprintf('%s — 计算时间对比', water_types{w_idx}));
    saveas(gcf, sprintf('fig_time_bar_%s.png', water_types{w_idx}));
end

% ================= 汇总对比表 =================
fprintf('\n===== 汇总对比 (N = 1e5) =====\n');
n5_idx = find(N_arr == 1e5, 1);
if isempty(n5_idx), n5_idx = 3; end

fprintf('%-14s | %5s | %8s | %8s | %7s | %7s\n', ...
    '水质', '距离', 'Ref(dB)', 'WCI(dB)', 'SA(dB)', 'WCI(s)', 'SA(s)');
fprintf('%s\n', repmat('-', 1, 75));

for w_idx = 1:num_W
    dists = dist_cell{w_idx};
    for d_idx = 1:length(dists)
        L = dists(d_idx);
        PL_ref = data_phys.PL_Cell{w_idx}(d_idx);
        PL_wci = data_wci.PL_Cell{w_idx}(d_idx, n5_idx);
        PL_sa  = data_sa.PL_Cell{w_idx}(d_idx, n5_idx);
        T_wci  = data_wci.Time_Cell{w_idx}(d_idx, n5_idx);
        T_sa   = data_sa.Time_Cell{w_idx}(d_idx, n5_idx);
        fprintf('%-14s | %5d | %8.2f | %8.2f | %7.2f | %7.2f | %7.2f\n', ...
            water_types{w_idx}, L, -PL_ref, -PL_wci, -PL_sa, T_wci, T_sa);
    end
end

% ================= 效率加速比 =================
fprintf('\n===== SA-MCS vs WCI-MC 加速比 (N = 1e5) =====\n');
fprintf('%-14s | %5s | %10s | %10s\n', '水质', '距离', 'Time_WCI', 'Speedup');
fprintf('%s\n', repmat('-', 1, 50));

for w_idx = 1:num_W
    dists = dist_cell{w_idx};
    for d_idx = 1:length(dists)
        T_wci = data_wci.Time_Cell{w_idx}(d_idx, n5_idx);
        T_sa  = data_sa.Time_Cell{w_idx}(d_idx, n5_idx);
        speedup = T_wci / T_sa;
        fprintf('%-14s | %5d | %9.2fs | %9.1fx\n', ...
            water_types{w_idx}, dists(d_idx), T_wci, speedup);
    end
end

fprintf('\n所有图片已保存至当前文件夹。\n');
