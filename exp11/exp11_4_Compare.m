%% exp11_4: 六算法对比 — Jerlov IB 路径衰减 / 计算时间 vs 传输距离
%  WCI-MC (无湍流) vs WCI-MC (Snell-Only) vs WCI-MC (Full Turb.) vs
%  MCS Physical (No Turb.) vs MCS Physical (Snell-Only) vs SA-MCS (LUT)
%  加载 exp11_1 / exp11_2 / exp11_2b / exp11_3 / exp11_5 / exp11_6 的仿真结果

clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% ================= Load Data =================
data_none      = load('data_exp11_WCIMC_None.mat');
data_snell     = load('data_exp11_WCIMC_Snell.mat');
data_full      = load('data_exp11_WCIMC_Strong.mat');
data_phys_none = load('data_exp11_MCS_Physical_None.mat');
data_phys_turb = load('data_exp11_MCS_Physical_Turb.mat');
data_sa        = load('data_exp11_SA_MCS_LUT.mat');

dist_arr = data_none.dist_arr;
PL_none       = -data_none.PL_arr;
PL_snell      = -data_snell.PL_arr;
PL_full       = -data_full.PL_arr;
PL_phys_none  = -data_phys_none.PL_arr;
PL_phys_turb  = -data_phys_turb.PL_arr;
PL_phys_std   = data_phys_turb.PL_std_arr;     % 物理接收 Snell 帧间标准差
N_turb_phys   = data_phys_turb.N_turb;         % 湍流实现次数
PL_sa         = -data_sa.PL_arr;
T_none        = data_none.Time_arr;
T_snell       = data_snell.Time_arr;
T_full        = data_full.Time_arr;
T_phys_none   = data_phys_none.Time_arr;
T_phys_turb   = data_phys_turb.Time_arr;
T_sa          = data_sa.Time_arr;
N_packets = data_none.N_packets;

fprintf('=== exp11_4: Six-Algorithm Comparison ===\n');
fprintf('Water: %s | N = %d | Distances: [%s] m\n', ...
    data_none.water_type, N_packets, num2str(dist_arr));
fprintf('MCS Physical (Snell): N_turb = %d 帧相位屏平均\n\n', N_turb_phys);

% ================= Style =================
font_name = 'Times New Roman';
blue_light    = [0.40 0.65 0.90];   % 浅蓝 — SA-MCS (LUT)
blue_medium   = [0.15 0.40 0.75];   % 中蓝 — WCI-MC (No Turbulence)
blue_dark     = [0.05 0.20 0.55];   % 深蓝 — WCI-MC (Snell-Only)
teal_dark     = [0.00 0.45 0.45];   % 深青 — WCI-MC (Full Turbulence)
orange_warm   = [0.90 0.50 0.10];   % 暖橙 — MCS Physical (No Turbulence)
red_brick     = [0.70 0.15 0.15];   % 砖红 — MCS Physical (Snell-Only)

algo_labels = {'SA-MCS (LUT)', 'WCI-MC (No Turb.)', 'WCI-MC (Snell-Only)', ...
               'WCI-MC (Full Turb.)', 'MCS Physical (No Turb.)', 'MCS Physical (Snell)'};
markers = {'d-', 'o-', '^--', 's-', 'p-', 'v-.'};
colors = {blue_light, blue_medium, blue_dark, teal_dark, orange_warm, red_brick};

% ================= Figure 1: Path Loss vs Distance =================
figure('Position', [100, 100, 960, 600]);
hold on;

plot(dist_arr, PL_sa, markers{1}, ...
    'Color', colors{1}, 'LineWidth', 1.2, ...
    'MarkerFaceColor', colors{1}, 'MarkerSize', 6, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_none, markers{2}, ...
    'Color', colors{2}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{2}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_snell, markers{3}, ...
    'Color', colors{3}, 'LineWidth', 1.5, ...
    'MarkerFaceColor', colors{3}, 'MarkerSize', 8, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_full, markers{4}, ...
    'Color', colors{4}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{4}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_phys_none, markers{5}, ...
    'Color', colors{5}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{5}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, PL_phys_turb, markers{6}, ...
    'Color', colors{6}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{6}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
errorbar(dist_arr, PL_phys_turb, PL_phys_std, 'LineStyle', 'none', ...
    'Color', colors{6}, 'LineWidth', 1.2, 'CapSize', 5);

hold off;

set(gca, 'FontName', font_name, 'FontSize', 11);
xlabel('Transmission Distance (m)', 'FontName', font_name, 'FontSize', 12);
ylabel('Path Loss (dB)', 'FontName', font_name, 'FontSize', 12);
title(sprintf('Jerlov IB — Path Loss vs Distance  (N = 10^{%d})', log10(N_packets)), ...
    'FontName', font_name, 'FontSize', 13);
legend(algo_labels, 'Location', 'northwest', 'FontName', font_name, 'FontSize', 7);
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
plot(dist_arr, T_none, markers{2}, ...
    'Color', colors{2}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{2}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, T_snell, markers{3}, ...
    'Color', colors{3}, 'LineWidth', 1.5, ...
    'MarkerFaceColor', colors{3}, 'MarkerSize', 8, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, T_full, markers{4}, ...
    'Color', colors{4}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{4}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, T_phys_none, markers{5}, ...
    'Color', colors{5}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{5}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');
plot(dist_arr, T_phys_turb, markers{6}, ...
    'Color', colors{6}, 'LineWidth', 1.8, ...
    'MarkerFaceColor', colors{6}, 'MarkerSize', 9, ...
    'MarkerEdgeColor', 'none');

hold off;

set(gca, 'FontName', font_name, 'FontSize', 11);
xlabel('Transmission Distance (m)', 'FontName', font_name, 'FontSize', 12);
ylabel('Computation Time (s)', 'FontName', font_name, 'FontSize', 12);
title(sprintf('Jerlov IB — Computation Time vs Distance  (N = 10^{%d})', log10(N_packets)), ...
    'FontName', font_name, 'FontSize', 13);
legend(algo_labels, 'Location', 'northwest', 'FontName', font_name, 'FontSize', 7);
grid on; box on;

saveas(gcf, 'fig_exp11_Time_vs_Distance.png');
fprintf('Figure saved: fig_exp11_Time_vs_Distance.png\n');

% ================= Figure 3: Excess Loss — Turbulence Contribution =================
PL_excess_snell     = PL_snell - PL_none;         % WCI-MC Snell-Only vs WCI-MC None
PL_excess_full      = PL_full - PL_none;          % Full Turb. vs WCI-MC None
PL_excess_phys_turb = PL_phys_turb - PL_phys_none; % Physical Snell vs Physical None

figure('Position', [100, 100, 960, 500]);

b = bar(dist_arr, [PL_excess_snell(:), PL_excess_full(:), PL_excess_phys_turb(:)], 'grouped');
b(1).FaceColor = blue_dark;   b(1).EdgeColor = 'none';
b(2).FaceColor = teal_dark;   b(2).EdgeColor = 'none';
b(3).FaceColor = red_brick;   b(3).EdgeColor = 'none';

set(gca, 'FontName', font_name, 'FontSize', 11);
xlabel('Transmission Distance (m)', 'FontName', font_name, 'FontSize', 12);
ylabel('Excess Path Loss (dB)', 'FontName', font_name, 'FontSize', 12);
title(sprintf('Jerlov IB — Turbulence-Induced Excess Loss  (N = 10^{%d})', log10(N_packets)), ...
    'FontName', font_name, 'FontSize', 13);
legend({'WCI-MC Snell (refraction)', 'Full Turb. (refr. + spreading)', ...
    'Physical Snell (refraction)'}, ...
    'Location', 'northwest', 'FontName', font_name, 'FontSize', 8);
grid on; box on;

for i = 1:length(dist_arr)
    if abs(PL_excess_snell(i)) > 0.01
        text(dist_arr(i) - 0.38, PL_excess_snell(i) + 0.04, sprintf('%.1f', PL_excess_snell(i)), ...
            'FontName', font_name, 'FontSize', 7, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    if abs(PL_excess_full(i)) > 0.01
        text(dist_arr(i), PL_excess_full(i) + 0.04, sprintf('%.1f', PL_excess_full(i)), ...
            'FontName', font_name, 'FontSize', 7, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    if abs(PL_excess_phys_turb(i)) > 0.01
        text(dist_arr(i) + 0.38, PL_excess_phys_turb(i) + 0.04, sprintf('%.1f', PL_excess_phys_turb(i)), ...
            'FontName', font_name, 'FontSize', 7, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

saveas(gcf, 'fig_exp11_ExcessLoss.png');
fprintf('Figure saved: fig_exp11_ExcessLoss.png\n');

% ================= Console Summary =================
fprintf('\n===== Full Summary (N = %d) =====\n', N_packets);
fprintf('%-6s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'Dist', 'SA(dB)', 'W0(dB)', 'WSnl(dB)', 'WFul(dB)', 'P0(dB)', 'PSnl(dB)', ...
    'SA(s)', 'W0(s)', 'WSnl(s)', 'WFul(s)', 'P0(s)', 'PSnl(s)');
fprintf('%s\n', repmat('-', 1, 158));
for i = 1:length(dist_arr)
    fprintf('%-6d | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f\n', ...
        dist_arr(i), PL_sa(i), PL_none(i), PL_snell(i), PL_full(i), PL_phys_none(i), PL_phys_turb(i), ...
        T_sa(i), T_none(i), T_snell(i), T_full(i), T_phys_none(i), T_phys_turb(i));
end

fprintf('\nDone.\n');
