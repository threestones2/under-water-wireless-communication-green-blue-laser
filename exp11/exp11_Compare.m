%% exp11_Compare: 四算法对比 — Jerlov IB 路径衰减 / 计算时间 vs 传输距离
%  VRT-MC (No Turb.) / VRT-MC (With Turb.) / SA-MCS (LUT) / MCS Physical (No Turb.)
%  加载 exp11_1_VRTMC_None / exp11_2_VRTMC / exp11_3_SA_MCS_LUT / exp11_6_MCS_Physical_None 的仿真结果

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= Load Data =================
data_vrt_none  = load('data_exp11_VRTMC_None.mat');
data_vrt_turb  = load('data_exp11_VRTMC.mat');
data_sa        = load('data_exp11_SA_MCS_LUT.mat');
data_phys_none = load('data_exp11_MCS_Physical_None.mat');

dist_arr = data_vrt_none.dist_arr;
PL_vrt_none   = -data_vrt_none.PL_arr;
PL_vrt_turb   = -data_vrt_turb.PL_arr;
PL_sa         = -data_sa.PL_arr;
PL_phys_none  = -data_phys_none.PL_arr;
T_vrt_none    = data_vrt_none.Time_arr;
T_vrt_turb    = data_vrt_turb.Time_arr;
T_sa          = data_sa.Time_arr;
T_phys_none   = data_phys_none.Time_arr;
N_packets = data_vrt_none.N_packets;

fprintf('=== exp11_Compare: Four-Algorithm Comparison ===\n');
fprintf('Water: %s | N = %d | Distances: [%s] m\n\n', ...
    data_vrt_none.water_type, N_packets, num2str(dist_arr));

% ================= Style =================
font_name = 'Times New Roman';
blue_light    = [0.40 0.65 0.90];   % 浅蓝 — SA-MCS (LUT)
teal_light    = [0.30 0.70 0.65];   % 浅青 — VRT-MC (No Turbulence)
teal_dark     = [0.00 0.45 0.45];   % 深青 — VRT-MC (With Turbulence)
orange_warm   = [0.90 0.50 0.10];   % 暖橙 — MCS Physical (No Turbulence)

algo_labels = {'SA-MCS (LUT)', 'VRT-MC (No Turb.)', ...
               'VRT-MC (With Turb.)', 'MCS Physical (No Turb.)'};
markers = {'d-', 's-', 's--', 'p-'};
colors = {blue_light, teal_light, teal_dark, orange_warm};

% ================= Figure 1: Path Loss vs Distance =================
figure('Position', [100, 100, 960, 600]);
hold on;

plot(dist_arr, PL_sa, markers{1}, ...
    'Color', colors{1}, 'LineWidth', 1.2, ...
    'MarkerFaceColor', colors{1}, 'MarkerSize', 6, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_vrt_none, markers{2}, ...
    'Color', colors{2}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{2}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_vrt_turb, markers{3}, ...
    'Color', colors{3}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{3}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_phys_none, markers{4}, ...
    'Color', colors{4}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{4}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');

hold off;

set(gca, 'FontName', font_name, 'FontSize', 11);
xlabel('Transmission Distance (m)', 'FontName', font_name, 'FontSize', 12);
ylabel('Path Loss (dB)', 'FontName', font_name, 'FontSize', 12);
title(sprintf('Jerlov IB — Path Loss vs Distance  (N = 10^{%d})', log10(N_packets)), ...
    'FontName', font_name, 'FontSize', 13);
legend(algo_labels, 'Location', 'northwest', 'FontName', font_name, 'FontSize', 9);
grid on; box on;

saveas(gcf, 'fig_exp11_PL_vs_Distance.png');
fprintf('Figure saved: fig_exp11_PL_vs_Distance.png\n');

% ================= Figure 2: Computation Time vs Distance =================
figure('Position', [100, 100, 960, 600]);
hold on;

plot(dist_arr, T_sa, markers{1}, ...
    'Color', colors{1}, 'LineWidth', 1.2, ...
    'MarkerFaceColor', colors{1}, 'MarkerSize', 6, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, T_vrt_none, markers{2}, ...
    'Color', colors{2}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{2}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, T_vrt_turb, markers{3}, ...
    'Color', colors{3}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{3}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, T_phys_none, markers{4}, ...
    'Color', colors{4}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{4}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');

hold off;

set(gca, 'FontName', font_name, 'FontSize', 11);
xlabel('Transmission Distance (m)', 'FontName', font_name, 'FontSize', 12);
ylabel('Computation Time (s)', 'FontName', font_name, 'FontSize', 12);
title(sprintf('Jerlov IB — Computation Time vs Distance  (N = 10^{%d})', log10(N_packets)), ...
    'FontName', font_name, 'FontSize', 13);
legend(algo_labels, 'Location', 'northwest', 'FontName', font_name, 'FontSize', 9);
grid on; box on;

saveas(gcf, 'fig_exp11_Time_vs_Distance.png');
fprintf('Figure saved: fig_exp11_Time_vs_Distance.png\n');

% ================= Figure 3: Excess Loss — VRT-MC Turbulence Contribution =================
PL_excess_vrt = PL_vrt_turb - PL_vrt_none;

figure('Position', [100, 100, 700, 500]);

b = bar(dist_arr, PL_excess_vrt(:), 'FaceColor', teal_dark, 'EdgeColor', 'none');

set(gca, 'FontName', font_name, 'FontSize', 11);
xlabel('Transmission Distance (m)', 'FontName', font_name, 'FontSize', 12);
ylabel('Excess Path Loss (dB)', 'FontName', font_name, 'FontSize', 12);
title(sprintf('Jerlov IB — VRT-MC Turbulence-Induced Excess Loss  (N = 10^{%d})', log10(N_packets)), ...
    'FontName', font_name, 'FontSize', 13);
grid on; box on;

for i = 1:length(dist_arr)
    if abs(PL_excess_vrt(i)) > 0.01
        text(dist_arr(i), PL_excess_vrt(i) + 0.04, sprintf('%.1f', PL_excess_vrt(i)), ...
            'FontName', font_name, 'FontSize', 9, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

saveas(gcf, 'fig_exp11_ExcessLoss.png');
fprintf('Figure saved: fig_exp11_ExcessLoss.png\n');

% ================= Console Summary =================
fprintf('\n===== Full Summary (N = %d) =====\n', N_packets);
fprintf('%-6s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'Dist', 'SA(dB)', 'V0(dB)', 'VTurb(dB)', 'P0(dB)', ...
    'SA(s)', 'V0(s)', 'VTurb(s)', 'P0(s)');
fprintf('%s\n', repmat('-', 1, 100));
for i = 1:length(dist_arr)
    fprintf('%-6d | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f\n', ...
        dist_arr(i), PL_sa(i), PL_vrt_none(i), PL_vrt_turb(i), PL_phys_none(i), ...
        T_sa(i), T_vrt_none(i), T_vrt_turb(i), T_phys_none(i));
end

fprintf('\nDone.\n');
