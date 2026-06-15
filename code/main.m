%% 多种群IA vs 标准IA 对比实验 (50维Rastrigin, 30次运行)
clear; clc; close all;

%% 参数设置
D = 50;  N_total = 100;  T = 500;  BOUNDS = [-5.12, 5.12];  N_RUNS = 30;

fprintf('MPIA vs IA 对比实验\n');
fprintf('问题: %d维Rastrigin | 种群: %d | 迭代: %d | 运行: %d次\n', D, N_total, T, N_RUNS);

params_ia = struct('clone_rate',10, 'mutate_scale',0.5, 'replace_rate',0.05);
params_mpia = struct('K',5, 'clone_rate',10, 'mutate_scale',0.5, ...
    'cauchy_scale',0.3, 'replace_rate',0.05, 'migrate_T',10, 'migrate_M',2);

%% 30次独立运行：标准IA 
fprintf('\n>>> 运行标准免疫算法 (IA) ...\n');

ia_results    = zeros(N_RUNS, 1);      
ia_conv_best  = zeros(T, N_RUNS);       
ia_conv_avg   = zeros(T, N_RUNS);       

tic;
for run = 1:N_RUNS
    [~, best_f, cb, ca] = IA(D, N_total, T, BOUNDS, params_ia);
    ia_results(run)      = best_f;
    ia_conv_best(:, run) = cb;
    ia_conv_avg(:, run)  = ca;
    fprintf('  IA Run %2d/%d  最优适应度: %.4f\n', run, N_RUNS, best_f);
end
ia_time = toc;
fprintf('  IA 总耗时: %.2f 秒 | 平均每次: %.2f 秒\n', ia_time, ia_time/N_RUNS);

%% 30次独立运行：多种群IA 
fprintf('\n>>> 运行多种群免疫算法 (MPIA) ...\n');

mpia_results    = zeros(N_RUNS, 1);
mpia_conv_best  = zeros(T, N_RUNS);
mpia_conv_avg   = zeros(T, N_RUNS);

tic;
for run = 1:N_RUNS
    [~, best_f, cb, ca, ~] = MPIA(D, N_total, T, BOUNDS, params_mpia);
    mpia_results(run)      = best_f;
    mpia_conv_best(:, run) = cb;
    mpia_conv_avg(:, run)  = ca;
    fprintf('  MPIA Run %2d/%d  最优适应度: %.4f\n', run, N_RUNS, best_f);
end
mpia_time = toc;
fprintf('  MPIA 总耗时: %.2f 秒 | 平均每次: %.2f 秒\n', mpia_time, mpia_time/N_RUNS);

%% 统计结果汇总
fprintf(' 统计结果汇总\n');
fprintf('%-12s %12s %12s %12s %12s %12s\n', ...
    '算法', '最优值', '最差值', '中位数', '均值', '标准差');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-12s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '标准IA', ...
    min(ia_results), max(ia_results), median(ia_results), ...
    mean(ia_results), std(ia_results));
fprintf('%-12s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '多种群IA', ...
    min(mpia_results), max(mpia_results), median(mpia_results), ...
    mean(mpia_results), std(mpia_results));

% 改进幅度
ia_mean  = mean(ia_results);
mpia_mean = mean(mpia_results);
improvement_best = (min(ia_results) - min(mpia_results)) / min(ia_results) * 100;
improvement_mean = (ia_mean - mpia_mean) / ia_mean * 100;
fprintf('\n改进幅度:\n');
fprintf('  最优值提升: %.2f%%\n', improvement_best);
fprintf('  均值提升:   %.2f%%\n', improvement_mean);

%% Wilcoxon秩和检验
try
    [p_ranksum, h_ranksum] = ranksum(ia_results, mpia_results);
catch
    [p_ranksum, h_ranksum] = manual_ranksum(ia_results, mpia_results, 0.05);
end
fprintf('\nWilcoxon检验 (α=0.05): p = %.6f\n', p_ranksum);
if h_ranksum == 1
    fprintf('  结论: 拒绝H0，两种算法存在显著性差异 ★\n');
else
    fprintf('  结论: 不能拒绝H0，两种算法无显著性差异\n');
end

% 计算均值和标准差
ia_best_mean   = mean(ia_conv_best, 2);
ia_best_std    = std(ia_conv_best, 0, 2);
mpia_best_mean = mean(mpia_conv_best, 2);
mpia_best_std  = std(mpia_conv_best, 0, 2);

ia_avg_mean   = mean(ia_conv_avg, 2);
ia_avg_std    = std(ia_conv_avg, 0, 2);
mpia_avg_mean = mean(mpia_conv_avg, 2);
mpia_avg_std  = std(mpia_conv_avg, 0, 2);

%% 收敛曲线图
x_vec = (1:T)';  x_fill = [x_vec; flipud(x_vec)];

figure('Position', [50, 50, 1400, 550], 'Color', 'w');
subplot(1,2,1); hold on;
fill(x_fill, [ia_best_mean+ia_best_std; flipud(ia_best_mean-ia_best_std)], 'b', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
fill(x_fill, [mpia_best_mean+mpia_best_std; flipud(mpia_best_mean-mpia_best_std)], 'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
h1 = plot(x_vec, ia_best_mean, 'b-', 'LineWidth', 1.8);
h2 = plot(x_vec, mpia_best_mean, 'r-', 'LineWidth', 1.8);
xlabel('迭代次数', 'FontSize', 11); ylabel('最优适应度', 'FontSize', 11);
title(sprintf('最优适应度收敛曲线 (D=%d)', D), 'FontSize', 12, 'FontWeight', 'bold');
legend([h1,h2], {'标准IA','多种群IA'}, 'Location', 'northeast', 'FontSize', 10);
grid on; box on; hold off;

subplot(1,2,2); hold on;
fill(x_fill, [ia_avg_mean+ia_avg_std; flipud(ia_avg_mean-ia_avg_std)], 'b', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
fill(x_fill, [mpia_avg_mean+mpia_avg_std; flipud(mpia_avg_mean-mpia_avg_std)], 'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
h3 = plot(x_vec, ia_avg_mean, 'b-', 'LineWidth', 1.8);
h4 = plot(x_vec, mpia_avg_mean, 'r-', 'LineWidth', 1.8);
xlabel('迭代次数', 'FontSize', 11); ylabel('平均适应度', 'FontSize', 11);
title(sprintf('平均适应度收敛曲线 (D=%d)', D), 'FontSize', 12, 'FontWeight', 'bold');
legend([h3,h4], {'标准IA','多种群IA'}, 'Location', 'northeast', 'FontSize', 10);
grid on; box on; hold off;
sgtitle('MPIA vs IA — 50维Rastrigin优化对比', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, '收敛曲线对比.png');
fprintf('图片已保存: 收敛曲线对比.png\n');

%% 箱线图
figure('Position', [50, 50, 650, 500], 'Color', 'w'); hold on;
data_all = [ia_results, mpia_results];
for idx = 1:2
    d = data_all(:, idx);
    q1 = prctile_basic(d, 25);  q2 = prctile_basic(d, 50);  q3 = prctile_basic(d, 75);
    iqr_v = q3 - q1;
    wl = max(min(d), q1-1.5*iqr_v);  wh = min(max(d), q3+1.5*iqr_v);
    c = [0.3 0.5 1.0; 1.0 0.4 0.4];  bw = 0.25;  xc = idx;
    plot([xc xc], [wl q1], 'k-', 'LineWidth', 1.2);
    plot([xc xc], [q3 wh], 'k-', 'LineWidth', 1.2);
    plot([xc-bw/2 xc+bw/2], [wl wl], 'k-', 'LineWidth', 1);
    plot([xc-bw/2 xc+bw/2], [wh wh], 'k-', 'LineWidth', 1);
    rectangle('Position', [xc-bw/2, q1, bw, iqr_v], 'FaceColor', c(idx,:), 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot([xc-bw/2 xc+bw/2], [q2 q2], 'k-', 'LineWidth', 1.8);
    out = d(d < wl | d > wh);
    if ~isempty(out), scatter(xc*ones(size(out)), out, 20, 'k', 'x', 'LineWidth', 1.2); end
    plot(xc, mean(d), 's', 'MarkerSize', 9, 'MarkerFaceColor', c(idx,:)*0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'标准IA', '多种群IA'});
ylabel('最优适应度', 'FontSize', 12);
title(sprintf('30次运行最优值分布 (D=%d)', D), 'FontSize', 13, 'FontWeight', 'bold');
grid on; hold off;
saveas(gcf, '箱线图对比.png');
fprintf('图片已保存: 箱线图对比.png\n');

%% 收敛前期放大图（前100代）
figure('Position', [50, 50, 650, 500], 'Color', 'w'); hold on;
plot(1:100, ia_best_mean(1:100), 'b-', 'LineWidth', 1.8);
plot(1:100, mpia_best_mean(1:100), 'r-', 'LineWidth', 1.8);
xlabel('迭代次数', 'FontSize', 11); ylabel('最优适应度', 'FontSize', 11);
title('收敛前期对比（前100代）', 'FontSize', 12, 'FontWeight', 'bold');
legend({'标准IA','多种群IA'}, 'Location', 'northeast', 'FontSize', 10);
grid on; box on; hold off;
saveas(gcf, '收敛前期对比.png');
fprintf('图片已保存: 收敛前期对比.png\n');

%% 输出最终结果
fprintf('\n实验完成\n');
fprintf('标准IA:\n');
fprintf('  最优: %.4f | 均值: %.4f | 标准差: %.4f\n', ...
    min(ia_results), mean(ia_results), std(ia_results));
fprintf('多种群IA:\n');
fprintf('  最优: %.4f | 均值: %.4f | 标准差: %.4f\n', ...
    min(mpia_results), mean(mpia_results), std(mpia_results));
fprintf('  均值改进: %.2f%% | p-value: %.6f\n', improvement_mean, p_ranksum);

%%  Wilcoxon 秩和检验
function [p, h] = manual_ranksum(x, y, alpha)
    if nargin < 3, alpha = 0.05; end
    x = x(:); y = y(:);
    n1 = length(x); n2 = length(y);
    combined = [x; y];
    group = [ones(n1,1); 2*ones(n2,1)]; 
    [sorted_vals, sort_idx] = sort(combined);
    ranks = zeros(size(combined));
    tie_ranks = zeros(size(sorted_vals));
    i = 1;
    while i <= length(sorted_vals)
        j = i;
        while j < length(sorted_vals) && sorted_vals(j+1) == sorted_vals(i)
            j = j + 1;
        end
        tie_ranks(i:j) = mean(i:j);
        i = j + 1;
    end
    ranks(sort_idx) = tie_ranks;
    R1 = sum(ranks(group == 1));
    U1 = R1 - n1*(n1+1)/2;
    U2 = n1*n2 - U1;
    U = min(U1, U2);
    EU = n1*n2/2;
    [~, ~, freq] = unique(combined);
    tj = histcounts(freq, 1:max(freq)+1);
    tj = tj(tj > 1);
    tie_correction = sum(tj.^3 - tj) / ((n1+n2)*((n1+n2)-1));
    sigma_U = sqrt(n1*n2/12 * ((n1+n2+1) - tie_correction));
    z = (abs(U - EU) - 0.5) / sigma_U;
    p = 2 * (1 - (1 + erf(abs(z)/sqrt(2))) / 2);
    h = double(p < alpha);
end

%% 计算百分位数
function pct = prctile_basic(data, percentile)
    sorted = sort(data(:));
    n = length(sorted);
    pos = (percentile/100) * (n - 1) + 1;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        pct = sorted(lo);
    else
        pct = sorted(lo) + (pos - lo) * (sorted(hi) - sorted(lo));
    end
end

