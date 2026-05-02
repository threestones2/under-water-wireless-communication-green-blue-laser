%% Exp2_PDF: Joint Plotting Script for 20m and 50m Data
% 功能: 读取硬盘上的 20m 与 50m 仿真数据文件，同窗完成 PDF 统计并结合对数正态拟合绘图
clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

% 检查数据文件是否存在
if ~isfile('data_exp2_PDF_20m.mat') || ~isfile('data_exp2_PDF_50m.mat')
    error('未找到数据文件，请先运行 20m 和 50m 的仿真脚本。');
end

% ================= 数据加载 =================
data_20m = load('data_exp2_PDF_20m.mat', 'Intensity_Samples', 'water_types');
samples_20m = squeeze(data_20m.Intensity_Samples(1, 1, :));

data_50m = load('data_exp2_PDF_50m.mat', 'Intensity_Samples');
samples_50m = squeeze(data_50m.Intensity_Samples(1, 1, :));

% 获取标题用水质环境字符串
water_type_str = data_20m.water_types{1};

% ================= 联合 PDF 统计分析与对数正态分布拟合 =================
% 采用指定图窗配置
figure('Name', ['Joint Intensity PDF - ', water_type_str], 'Color', 'w', ...
       'Position', [100, 100, 640, 480]);
hold on; grid on; box on;

% 绘图样式与基础参数定义
distances = [20, 50];
samples_cell = {samples_20m, samples_50m};

% 更新为指定的绿色和黄色
plot_colors = {'#77AC30', '#EDB120'}; % 绿 (20m), 黄 (50m)
plot_markers = {'o', 'd'};            % 圆圈 (20m), 菱形 (50m)

for d_idx = 1:2
    samples = samples_cell{d_idx};
    dist_val = distances(d_idx);
    c = plot_colors{d_idx};
    m = plot_markers{d_idx};
    
    mean_I = mean(samples);
    I_norm = samples / mean_I;
    
    % 计算闪烁指数 (Scintillation Index) 用于计算对数正态方差
    var_I = var(samples);
    SI = var_I / (mean_I^2);
    
    % 提取概率密度数据并绘制仿真散点
    num_bins = 20; 
    [pdf_counts, edges] = histcounts(I_norm, num_bins, 'Normalization', 'pdf');
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    
    % 绘制散点并在图例中加入闪烁指数 (SI) 数值
    plot(bin_centers, pdf_counts, m, 'Color', c, ...
         'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'none', ...
         'DisplayName', sprintf('simulation (L=%dm, SI=%.4f)', dist_val, SI));
         
    % 计算并绘制理论对数正态分布 (Lognormal PDF) 连续曲线
    I_axis = linspace(0.01, max(I_norm)*1.2, 500);
    if SI > 0
        sigma2_ln = log(SI + 1);
        mu_ln = -0.5 * sigma2_ln; 
        PDF_theo = (1 ./ (I_axis .* sqrt(sigma2_ln * 2 * pi))) .* exp( - (log(I_axis) - mu_ln).^2 ./ (2 * sigma2_ln) );
        
        plot(I_axis, PDF_theo, '-', 'Color', c, 'LineWidth', 1.5, ...
             'DisplayName', sprintf('lognormal (L=%dm)', dist_val));
    end
end

% ================= 图表格式与坐标轴规范设置 =================
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'LineWidth', 1);
xlabel('Normalized Intensity $I / \langle I \rangle$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Probability Density Function (PDF)', 'FontSize', 13, 'FontName', 'Times New Roman');
title(['Received Intensity PDF: ', water_type_str], 'FontSize', 14, 'FontName', 'Times New Roman');
legend('show', 'Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');