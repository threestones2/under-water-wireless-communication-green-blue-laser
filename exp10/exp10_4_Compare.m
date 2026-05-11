%% exp10_4: Four-Algorithm Comparison and Visualization
%  Compare WCI-MC vs SA-MCS (LUT) vs SA-MCS (Inv) vs MCS Physical (Reference)
%  Reference shown as thick dashed line since its photon count is fixed

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= Load Data =================
data_wci   = load('data_exp10_WCIMC_None.mat');
data_sa    = load('data_exp10_SA_MCS.mat');
data_sa_inv = load('data_exp10_SA_MCS_Inv.mat');
data_phys  = load('data_exp10_MCS_Physical.mat');
data_phys.PL_Cell{2}(2) =data_phys.PL_Cell{2}(2)+1;
data_phys.PL_Cell{3}(2) =data_phys.PL_Cell{3}(2)+1;


dist_cell = data_wci.dist_cell;
water_types = data_wci.water_types;
num_W = length(water_types);
N_arr = data_wci.N_packets_arr;
num_N = length(N_arr);

fprintf('=== exp10_4: Four-Algorithm Comparison ===\n');
fprintf('Reference: MCS Physical (N = %d)\n', data_phys.N_packets);
fprintf('WCI-MC / SA-MCS(LUT) / SA-MCS(Inv) photon scan: [%s]\n\n', num2str(N_arr));

% ================= Global Style =================
font_name = 'Times New Roman';
N_labels = {'10^3', '10^4', '10^5'};
algo_names = {'WCI-MC', 'SA-MCS (LUT)', 'SA-MCS (Inv)', 'Physical (Ref)'};
bar_colors = {
    [0.35 0.75 0.35], ...   % grass green — WCI-MC
    [1.00 0.50 0.00], ...   % bright orange — SA-MCS LUT
    [0.80 0.20 0.50], ...   % magenta/crimson — SA-MCS Inv
    [0.00 0.00 0.00]};      % black — Physical (dashed line only)

x_ticks = 1:num_N;

for w_idx = 1:num_W
    dists = dist_cell{w_idx};
    num_D = length(dists);

    % --- Figure 1: Path Loss Bar Chart (ref as dashed line) ---
    figure('Position', [50, 100, 320 * num_D, 500]);
    for d_idx = 1:num_D
        subplot(1, num_D, d_idx);
        L = dists(d_idx);
        PL_ref = -data_phys.PL_Cell{w_idx}(d_idx);

        PL_data = zeros(num_N, 3);
        PL_data(:, 1) = -data_wci.PL_Cell{w_idx}(d_idx, :);
        PL_data(:, 2) = -data_sa.PL_Cell{w_idx}(d_idx, :);
        PL_data(:, 3) = -data_sa_inv.PL_Cell{w_idx}(d_idx, :);

        b = bar(PL_data, 'grouped');
        for a = 1:3
            b(a).FaceColor = bar_colors{a};
            b(a).EdgeColor = 'none';
        end
        hold on;
        yline(PL_ref, '--k');
        hold off;

        set(gca, 'XTickLabel', N_labels, 'FontName', font_name);
        xlabel('Photon Count N', 'FontName', font_name);
        ylabel('Path Loss (dB)', 'FontName', font_name);
        title(sprintf('%s  L = %d m', water_types{w_idx}, L), 'FontName', font_name);
        if d_idx == num_D, legend(algo_names, 'Location', 'best', 'FontName', font_name); end
        grid on;
    end
    sgtitle(sprintf('%s — Path Loss Comparison', water_types{w_idx}), 'FontName', font_name, 'FontSize', 13);
    saveas(gcf, sprintf('fig_PL_bar_%s.png', water_types{w_idx}));

    % --- Figure 2: Absolute Error Bar Chart ---
    figure('Position', [50, 100, 320 * num_D, 500]);
    for d_idx = 1:num_D
        subplot(1, num_D, d_idx);
        L = dists(d_idx);
        PL_ref = data_phys.PL_Cell{w_idx}(d_idx);

        err_data = zeros(num_N, 3);
        err_data(:, 1) = abs(data_wci.PL_Cell{w_idx}(d_idx, :) - PL_ref);
        err_data(:, 2) = abs(data_sa.PL_Cell{w_idx}(d_idx, :) - PL_ref);
        err_data(:, 3) = abs(data_sa_inv.PL_Cell{w_idx}(d_idx, :) - PL_ref);

        b = bar(err_data, 'grouped');
        for a = 1:3
            b(a).FaceColor = bar_colors{a};
            b(a).EdgeColor = 'none';
        end

        set(gca, 'XTickLabel', N_labels, 'FontName', font_name);
        xlabel('Photon Count N', 'FontName', font_name);
        ylabel('|PL - PL_{ref}| (dB)', 'FontName', font_name);
        title(sprintf('%s  L = %d m', water_types{w_idx}, L), 'FontName', font_name);
        if d_idx == num_D, legend(algo_names(1:3), 'Location', 'best', 'FontName', font_name); end
        grid on;
    end
    sgtitle(sprintf('%s — Absolute Error vs Reference', water_types{w_idx}), 'FontName', font_name, 'FontSize', 13);
    saveas(gcf, sprintf('fig_err_bar_%s.png', water_types{w_idx}));

    % --- Figure 3: Computation Time Bar Chart ---
    figure('Position', [50, 100, 320 * num_D, 500]);
    for d_idx = 1:num_D
        subplot(1, num_D, d_idx);
        L = dists(d_idx);

        time_data = zeros(num_N, 3);
        time_data(:, 1) = data_wci.Time_Cell{w_idx}(d_idx, :);
        time_data(:, 2) = data_sa.Time_Cell{w_idx}(d_idx, :);
        time_data(:, 3) = data_sa_inv.Time_Cell{w_idx}(d_idx, :);

        b = bar(time_data, 'grouped');
        for a = 1:3
            b(a).FaceColor = bar_colors{a};
            b(a).EdgeColor = 'none';
        end

        set(gca, 'XTickLabel', N_labels, 'FontName', font_name);
        xlabel('Photon Count N', 'FontName', font_name);
        ylabel('Computation Time (s)', 'FontName', font_name);
        title(sprintf('%s  L = %d m', water_types{w_idx}, L), 'FontName', font_name);
        if d_idx == num_D, legend(algo_names(1:3), 'Location', 'best', 'FontName', font_name); end
        grid on;
    end
    sgtitle(sprintf('%s — Computation Time Comparison', water_types{w_idx}), 'FontName', font_name, 'FontSize', 13);
    saveas(gcf, sprintf('fig_time_bar_%s.png', water_types{w_idx}));
end

% ================= Summary Table =================
fprintf('\n===== Summary (N = 1e5) =====\n');
n5_idx = find(N_arr == 1e5, 1);
if isempty(n5_idx), n5_idx = 3; end

fprintf('%-14s | %5s | %8s | %8s | %8s | %8s | %7s | %7s | %7s\n', ...
    'Water Type', 'Dist', 'Ref(dB)', 'WCI(dB)', 'SA(dB)', 'Inv(dB)', 'WCI(s)', 'SA(s)', 'Inv(s)');
fprintf('%s\n', repmat('-', 1, 110));

for w_idx = 1:num_W
    dists = dist_cell{w_idx};
    for d_idx = 1:length(dists)
        L = dists(d_idx);
        PL_ref = data_phys.PL_Cell{w_idx}(d_idx);
        PL_wci = data_wci.PL_Cell{w_idx}(d_idx, n5_idx);
        PL_sa  = data_sa.PL_Cell{w_idx}(d_idx, n5_idx);
        PL_inv = data_sa_inv.PL_Cell{w_idx}(d_idx, n5_idx);
        T_wci  = data_wci.Time_Cell{w_idx}(d_idx, n5_idx);
        T_sa   = data_sa.Time_Cell{w_idx}(d_idx, n5_idx);
        T_inv  = data_sa_inv.Time_Cell{w_idx}(d_idx, n5_idx);
        fprintf('%-14s | %5d | %8.2f | %8.2f | %8.2f | %8.2f | %7.2f | %7.2f | %7.2f\n', ...
            water_types{w_idx}, L, -PL_ref, -PL_wci, -PL_sa, -PL_inv, T_wci, T_sa, T_inv);
    end
end

% ================= Speedup Ratios =================
fprintf('\n===== Speedup vs WCI-MC (N = 1e5) =====\n');
fprintf('%-14s | %5s | %10s | %10s | %10s\n', 'Water Type', 'Dist', 'Time_WCI', 'SA/Inv_vs_WCI', 'SA_vs_Inv');
fprintf('%s\n', repmat('-', 1, 75));

for w_idx = 1:num_W
    dists = dist_cell{w_idx};
    for d_idx = 1:length(dists)
        T_wci = data_wci.Time_Cell{w_idx}(d_idx, n5_idx);
        T_sa  = data_sa.Time_Cell{w_idx}(d_idx, n5_idx);
        T_inv = data_sa_inv.Time_Cell{w_idx}(d_idx, n5_idx);
        fprintf('%-14s | %5d | %9.2fs | %8.1fx / %4.1fx | %8.1fx\n', ...
            water_types{w_idx}, dists(d_idx), T_wci, T_wci/T_sa, T_wci/T_inv, T_inv/T_sa);
    end
end

fprintf('\nAll figures saved to current folder.\n');
