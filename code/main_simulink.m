%% MPIA vs IA vs ZN — DC电机PID参数优化对比实验
%  被控对象: G(s) = 10/(s²+10s)  (直流电机简化模型)
%  优化目标: 最小化ITAE (积分时间绝对误差)
clear; clc; close all;

fprintf('══════════════════════════════════════════\n');
fprintf('  DC电机PID优化 — MPIA vs IA vs ZN\n');
fprintf('  被控对象: G(s)=10/(s²+10s)\n');
fprintf('══════════════════════════════════════════\n\n');

%% ── 环境检测 ──────────────────────────────
fprintf('>>> 环境检测 ...\n');
simulink_available = false;
cst_available = false;
try
    v = ver('simulink');
    if ~isempty(v), simulink_available = true; end
end
try
    v = ver('control');
    if ~isempty(v), cst_available = true; end
end

if simulink_available
    fprintf('  ✓ Simulink 可用\n');
    % 尝试构建Simulink模型
    try
        build_simulink_model();
        fprintf('  ✓ Simulink模型已构建\n');
    catch ME
        fprintf('  ⚠ Simulink模型构建失败: %s\n', ME.message);
        simulink_available = false;
    end
end
if cst_available
    fprintf('  ✓ Control System Toolbox 可用\n');
end
if ~simulink_available && ~cst_available
    fprintf('  ○ 使用手动仿真模式（RK4数值积分）\n');
end
fprintf('\n');

%% ── 参数配置 ──────────────────────────────
D = 3;                  % 优化维度: Kp, Ki, Kd
N_total = 50;           % 总种群大小
T = 100;                % 最大迭代次数
BOUNDS = [0, 100];      % PID参数搜索范围（统一上下界）
N_RUNS = 30;            % 独立运行次数

params_ia = struct('clone_rate', 10, 'mutate_scale', 0.5, 'replace_rate', 0.05);
params_mpia = struct('K', 5, 'clone_rate', 10, 'mutate_scale', 0.5, ...
    'cauchy_scale', 0.3, 'replace_rate', 0.05, 'migrate_T', 10, 'migrate_M', 2);

%% ── 基准: Ziegler-Nichols整定 ─────────────────
fprintf('>>> Ziegler-Nichols 整定 ...\n');
[K_zn, fit_zn] = zn_pid_tune();
fprintf('  ZN结果: Kp=%.4f, Ki=%.4f, Kd=%.4f, ITAE=%.4f\n\n', ...
    K_zn(1), K_zn(2), K_zn(3), fit_zn);

%% ── 标准IA优化PID (30次独立运行) ────────────
fprintf('>>> 运行标准免疫算法 (IA) ...\n');
ia_results    = zeros(N_RUNS, 1);
ia_K_best     = zeros(N_RUNS, 3);
ia_conv_best  = zeros(T, N_RUNS);
ia_conv_avg   = zeros(T, N_RUNS);

tic;
for run = 1:N_RUNS
    [best_K, best_f, cb, ca] = PID_IA(N_total, T, BOUNDS, params_ia);
    ia_results(run)      = best_f;
    ia_K_best(run, :)    = best_K;
    ia_conv_best(:, run) = cb;
    ia_conv_avg(:, run)  = ca;
    fprintf('  IA Run %2d/%d  K=[%6.2f,%6.2f,%6.2f]  ITAE=%.4f\n', ...
        run, N_RUNS, best_K(1), best_K(2), best_K(3), best_f);
end
ia_time = toc;
fprintf('  IA 总耗时: %.2f秒 | 平均: %.2f秒\n\n', ia_time, ia_time/N_RUNS);

%% ── MPIA优化PID (30次独立运行) ──────────────
fprintf('>>> 运行多种群免疫算法 (MPIA) ...\n');
mpia_results    = zeros(N_RUNS, 1);
mpia_K_best     = zeros(N_RUNS, 3);
mpia_conv_best  = zeros(T, N_RUNS);
mpia_conv_avg   = zeros(T, N_RUNS);

tic;
for run = 1:N_RUNS
    [best_K, best_f, cb, ca, ~] = PID_MPIA(N_total, T, BOUNDS, params_mpia);
    mpia_results(run)      = best_f;
    mpia_K_best(run, :)    = best_K;
    mpia_conv_best(:, run) = cb;
    mpia_conv_avg(:, run)  = ca;
    fprintf('  MPIA Run %2d/%d  K=[%6.2f,%6.2f,%6.2f]  ITAE=%.4f\n', ...
        run, N_RUNS, best_K(1), best_K(2), best_K(3), best_f);
end
mpia_time = toc;
fprintf('  MPIA 总耗时: %.2f秒 | 平均: %.2f秒\n\n', mpia_time, mpia_time/N_RUNS);

%% ── 统计汇总 ────────────────────────────────
fprintf('══════════════ 统计结果汇总 ══════════════\n');
fprintf('%-14s %12s %12s %12s %12s %12s\n', ...
    '算法', '最优值', '最差值', '中位数', '均值', '标准差');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-14s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    'ZN整定', fit_zn, fit_zn, fit_zn, fit_zn, 0);
fprintf('%-14s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '标准IA', ...
    min(ia_results), max(ia_results), median(ia_results), ...
    mean(ia_results), std(ia_results));
fprintf('%-14s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '多种群IA', ...
    min(mpia_results), max(mpia_results), median(mpia_results), ...
    mean(mpia_results), std(mpia_results));

% 改进幅度
ia_mean_val   = mean(ia_results);
mpia_mean_val = mean(mpia_results);
impr_best = (min(ia_results) - min(mpia_results)) / min(ia_results) * 100;
impr_mean = (ia_mean_val - mpia_mean_val) / ia_mean_val * 100;
impr_zn   = (fit_zn - mpia_mean_val) / fit_zn * 100;
fprintf('\n改进幅度:\n');
fprintf('  MPIA vs IA 最优提升:   %.2f%%\n', impr_best);
fprintf('  MPIA vs IA 均值提升:   %.2f%%\n', impr_mean);
fprintf('  MPIA vs ZN 均值提升:   %.2f%%\n', impr_zn);

% 最优PID参数
[~, ia_best_idx] = min(ia_results);
[~, mpia_best_idx] = min(mpia_results);
fprintf('\n最优PID参数:\n');
fprintf('  ZN:      Kp=%.4f, Ki=%.4f, Kd=%.4f\n', K_zn(1), K_zn(2), K_zn(3));
fprintf('  IA最优:  Kp=%.4f, Ki=%.4f, Kd=%.4f\n', ...
    ia_K_best(ia_best_idx, 1), ia_K_best(ia_best_idx, 2), ia_K_best(ia_best_idx, 3));
fprintf('  MPIA最优: Kp=%.4f, Ki=%.4f, Kd=%.4f\n', ...
    mpia_K_best(mpia_best_idx, 1), mpia_K_best(mpia_best_idx, 2), mpia_K_best(mpia_best_idx, 3));

%% ── Wilcoxon秩和检验 ────────────────────────
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

%% ── 绘图1: 阶跃响应对比 ──────────────────────
fprintf('\n>>> 生成对比图 ...\n');

figure('Position', [50, 50, 1400, 500], 'Color', 'w');

% 选择最优参数分别仿真阶跃响应
[t_zn, y_zn] = simulate_step_response(K_zn(1), K_zn(2), K_zn(3));
K_ia_best  = ia_K_best(ia_best_idx, :);
K_mpia_best = mpia_K_best(mpia_best_idx, :);
[t_ia, y_ia] = simulate_step_response(K_ia_best(1), K_ia_best(2), K_ia_best(3));
[t_mpia, y_mpia] = simulate_step_response(K_mpia_best(1), K_mpia_best(2), K_mpia_best(3));

subplot(1,2,1); hold on;
plot(t_zn, y_zn, 'g-', 'LineWidth', 1.8);
plot(t_ia, y_ia, 'b--', 'LineWidth', 1.8);
plot(t_mpia, y_mpia, 'r-', 'LineWidth', 1.8);
plot([t_zn(1) t_zn(end)], [1 1], 'k:', 'LineWidth', 0.8);
xlabel('时间 (s)', 'FontSize', 11);
ylabel('输出 y(t)', 'FontSize', 11);
title('阶跃响应对比', 'FontSize', 12, 'FontWeight', 'bold');
legend({sprintf('ZN (ITAE=%.2f)', fit_zn), ...
    sprintf('IA (ITAE=%.2f)', min(ia_results)), ...
    sprintf('MPIA (ITAE=%.2f)', min(mpia_results))}, ...
    'Location', 'southeast', 'FontSize', 10);
grid on; box on; hold off;

% 阶跃响应局部放大（超调区域）
subplot(1,2,2); hold on;
xlim_ms = [0.05, 0.5];
idx_zn  = t_zn >= xlim_ms(1) & t_zn <= xlim_ms(2);
idx_ia  = t_ia >= xlim_ms(1) & t_ia <= xlim_ms(2);
idx_mpia = t_mpia >= xlim_ms(1) & t_mpia <= xlim_ms(2);
plot(t_zn(idx_zn), y_zn(idx_zn), 'g-', 'LineWidth', 1.8);
plot(t_ia(idx_ia), y_ia(idx_ia), 'b--', 'LineWidth', 1.8);
plot(t_mpia(idx_mpia), y_mpia(idx_mpia), 'r-', 'LineWidth', 1.8);
plot(xlim_ms, [1 1], 'k:', 'LineWidth', 0.8);
xlabel('时间 (s)', 'FontSize', 11);
ylabel('输出 y(t)', 'FontSize', 11);
title('阶跃响应对比（超调区域放大）', 'FontSize', 12, 'FontWeight', 'bold');
legend({'ZN', 'IA', 'MPIA'}, 'Location', 'northeast', 'FontSize', 10);
grid on; box on; hold off;
sgtitle('DC电机PID控制 — 阶跃响应对比 (ZN vs IA vs MPIA)', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, '阶跃响应对比.png');
fprintf('  图片已保存: 阶跃响应对比.png\n');

%% ── 绘图2: 收敛曲线 ──────────────────────────
ia_best_mean   = mean(ia_conv_best, 2);
ia_best_std    = std(ia_conv_best, 0, 2);
mpia_best_mean = mean(mpia_conv_best, 2);
mpia_best_std  = std(mpia_conv_best, 0, 2);

figure('Position', [50, 50, 1400, 550], 'Color', 'w');
x_vec = (1:T)';
x_fill = [x_vec; flipud(x_vec)];

subplot(1,2,1); hold on;
fill(x_fill, [ia_best_mean+ia_best_std; flipud(ia_best_mean-ia_best_std)], ...
    'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill(x_fill, [mpia_best_mean+mpia_best_std; flipud(mpia_best_mean-mpia_best_std)], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
h1 = plot(x_vec, ia_best_mean, 'b-', 'LineWidth', 1.8);
h2 = plot(x_vec, mpia_best_mean, 'r-', 'LineWidth', 1.8);
xlabel('迭代次数', 'FontSize', 11);
ylabel('最优ITAE', 'FontSize', 11);
title('PID优化 — 最优ITAE收敛曲线', 'FontSize', 12, 'FontWeight', 'bold');
legend([h1, h2], {'标准IA', '多种群IA'}, 'Location', 'northeast', 'FontSize', 10);
grid on; box on; hold off;

ia_avg_mean   = mean(ia_conv_avg, 2);
ia_avg_std    = std(ia_conv_avg, 0, 2);
mpia_avg_mean = mean(mpia_conv_avg, 2);
mpia_avg_std  = std(mpia_conv_avg, 0, 2);

subplot(1,2,2); hold on;
fill(x_fill, [ia_avg_mean+ia_avg_std; flipud(ia_avg_mean-ia_avg_std)], ...
    'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill(x_fill, [mpia_avg_mean+mpia_avg_std; flipud(mpia_avg_mean-mpia_avg_std)], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
h3 = plot(x_vec, ia_avg_mean, 'b-', 'LineWidth', 1.8);
h4 = plot(x_vec, mpia_avg_mean, 'r-', 'LineWidth', 1.8);
xlabel('迭代次数', 'FontSize', 11);
ylabel('平均ITAE', 'FontSize', 11);
title('PID优化 — 平均ITAE收敛曲线', 'FontSize', 12, 'FontWeight', 'bold');
legend([h3, h4], {'标准IA', '多种群IA'}, 'Location', 'northeast', 'FontSize', 10);
grid on; box on; hold off;
sgtitle('PID参数优化收敛曲线 — IA vs MPIA', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'PID收敛曲线对比.png');
fprintf('  图片已保存: PID收敛曲线对比.png\n');

%% ── 绘图3: ITAE箱线图 ────────────────────────
figure('Position', [50, 50, 650, 500], 'Color', 'w'); hold on;
data_all = [ia_results, mpia_results];
for idx = 1:2
    d = data_all(:, idx);
    q1 = prctile_basic(d, 25);  q2 = prctile_basic(d, 50);  q3 = prctile_basic(d, 75);
    iqr_v = q3 - q1;
    wl = max(min(d), q1 - 1.5*iqr_v);
    wh = min(max(d), q3 + 1.5*iqr_v);
    c = [0.3 0.5 1.0; 1.0 0.4 0.4];
    bw = 0.25;  xc = idx;
    plot([xc xc], [wl q1], 'k-', 'LineWidth', 1.2);
    plot([xc xc], [q3 wh], 'k-', 'LineWidth', 1.2);
    plot([xc-bw/2 xc+bw/2], [wl wl], 'k-', 'LineWidth', 1);
    plot([xc-bw/2 xc+bw/2], [wh wh], 'k-', 'LineWidth', 1);
    rectangle('Position', [xc-bw/2, q1, bw, iqr_v], ...
        'FaceColor', c(idx,:), 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot([xc-bw/2 xc+bw/2], [q2 q2], 'k-', 'LineWidth', 1.8);
    out = d(d < wl | d > wh);
    if ~isempty(out)
        scatter(xc*ones(size(out)), out, 20, 'k', 'x', 'LineWidth', 1.2);
    end
    plot(xc, mean(d), 's', 'MarkerSize', 9, ...
        'MarkerFaceColor', c(idx,:)*0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'标准IA', '多种群IA'});
ylabel('ITAE值', 'FontSize', 12);
title(sprintf('PID优化 — %d次运行ITAE分布', N_RUNS), 'FontSize', 13, 'FontWeight', 'bold');
grid on; hold off;
saveas(gcf, 'PID箱线图对比.png');
fprintf('  图片已保存: PID箱线图对比.png\n');

%% ── 绘图4: PID参数散点图 ─────────────────────
figure('Position', [50, 50, 900, 700], 'Color', 'w');

% Kp-Ki平面
subplot(2,2,1); hold on;
scatter(ia_K_best(:,1), ia_K_best(:,2), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(mpia_K_best(:,1), mpia_K_best(:,2), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
plot(K_zn(1), K_zn(2), 'gs', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Kp', 'FontSize', 11); ylabel('Ki', 'FontSize', 11);
title('Kp-Ki 参数分布', 'FontSize', 12, 'FontWeight', 'bold');
legend({'IA', 'MPIA', 'ZN'}, 'Location', 'best', 'FontSize', 9);
grid on; box on; hold off;

% Kp-Kd平面
subplot(2,2,2); hold on;
scatter(ia_K_best(:,1), ia_K_best(:,3), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(mpia_K_best(:,1), mpia_K_best(:,3), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
plot(K_zn(1), K_zn(3), 'gs', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Kp', 'FontSize', 11); ylabel('Kd', 'FontSize', 11);
title('Kp-Kd 参数分布', 'FontSize', 12, 'FontWeight', 'bold');
legend({'IA', 'MPIA', 'ZN'}, 'Location', 'best', 'FontSize', 9);
grid on; box on; hold off;

% Ki-Kd平面
subplot(2,2,3); hold on;
scatter(ia_K_best(:,2), ia_K_best(:,3), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(mpia_K_best(:,2), mpia_K_best(:,3), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
plot(K_zn(2), K_zn(3), 'gs', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Ki', 'FontSize', 11); ylabel('Kd', 'FontSize', 11);
title('Ki-Kd 参数分布', 'FontSize', 12, 'FontWeight', 'bold');
legend({'IA', 'MPIA', 'ZN'}, 'Location', 'best', 'FontSize', 9);
grid on; box on; hold off;

% ITAE vs 参数关系
subplot(2,2,4); hold on;
scatter(ia_results, ia_K_best(:,1), 20, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
scatter(mpia_results, mpia_K_best(:,1), 20, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('ITAE', 'FontSize', 11); ylabel('Kp', 'FontSize', 11);
title('ITAE vs Kp', 'FontSize', 12, 'FontWeight', 'bold');
legend({'IA', 'MPIA'}, 'Location', 'best', 'FontSize', 9);
grid on; box on; hold off;

sgtitle('PID参数优化 — 参数分布与相关性分析', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'PID参数分布.png');
fprintf('  图片已保存: PID参数分布.png\n');

%% ── 最终报告 ─────────────────────────────────
fprintf('\n══════════════════ 实验完成 ══════════════════\n');
fprintf('DC电机模型: G(s) = 10/(s^2 + 10s)\n');
fprintf('优化维度: D=%d | 种群: %d | 迭代: %d | 运行: %d次\n\n', ...
    D, N_total, T, N_RUNS);

fprintf('ZN整定:\n');
fprintf('  Kp=%.4f, Ki=%.4f, Kd=%.4f → ITAE=%.4f\n\n', ...
    K_zn(1), K_zn(2), K_zn(3), fit_zn);

fprintf('标准IA (均值):\n');
fprintf('  Kp=%.2f±%.2f, Ki=%.2f±%.2f, Kd=%.2f±%.2f\n', ...
    mean(ia_K_best(:,1)), std(ia_K_best(:,1)), ...
    mean(ia_K_best(:,2)), std(ia_K_best(:,2)), ...
    mean(ia_K_best(:,3)), std(ia_K_best(:,3)));
fprintf('  ITAE: 最优=%.4f, 均值=%.4f±%.4f\n\n', ...
    min(ia_results), mean(ia_results), std(ia_results));

fprintf('多种群IA (均值):\n');
fprintf('  Kp=%.2f±%.2f, Ki=%.2f±%.2f, Kd=%.2f±%.2f\n', ...
    mean(mpia_K_best(:,1)), std(mpia_K_best(:,1)), ...
    mean(mpia_K_best(:,2)), std(mpia_K_best(:,2)), ...
    mean(mpia_K_best(:,3)), std(mpia_K_best(:,3)));
fprintf('  ITAE: 最优=%.4f, 均值=%.4f±%.4f\n\n', ...
    min(mpia_results), mean(mpia_results), std(mpia_results));

fprintf('性能提升:\n');
fprintf('  MPIA vs IA:  %.2f%% (均值) | p=%.6f\n', impr_mean, p_ranksum);
fprintf('  MPIA vs ZN:  %.2f%% (均值)\n', impr_zn);

%% ── Simulink 展示最优结果 ────────────────────
if simulink_available && exist('DC_Motor_Model.slx', 'file')
    fprintf('\n>>> 在 Simulink 中展示最优 PID 结果 ...\n');

    % 用最优 MPIA 参数
    K_best = mpia_K_best(mpia_best_idx, :);
    assignin('base', 'Kp_sim', K_best(1));
    assignin('base', 'Ki_sim', K_best(2));
    assignin('base', 'Kd_sim', K_best(3));
    assignin('base', 'N_sim', 100);

    fprintf('  最优PID参数已写入workspace:\n');
    fprintf('    Kp_sim=%.4f, Ki_sim=%.4f, Kd_sim=%.4f\n', K_best(1), K_best(2), K_best(3));
    fprintf('  打开 Simulink 模型，点击 Run 查看 Scope\n');

    % 打开模型（如果没打开）
    if ~bdIsLoaded('DC_Motor_Model')
        open_system('DC_Motor_Model');
    end

    % 自动跑一次仿真
    try
        simOut = sim('DC_Motor_Model', 'StopTime', '2');
        y_data = simOut.get('yout');
        if ~isempty(y_data)
            t_sim = y_data.time;
            y_sim = y_data.signals.values;
            fprintf('  Simulink仿真完成: y(∞)=%.4f, max(y)=%.4f\n', y_sim(end), max(y_sim));

            % 添加 Simulink 响应到已存在的阶跃响应图上
            figure(1); subplot(1,2,1); hold on;
            plot(t_sim, y_sim, 'm-.', 'LineWidth', 2);
            legend off;
            lh = findobj(gca, 'Type', 'Line');
            lg = legend(lh([end, end-1, end-2, end-3]), ...
                {'ZN','IA','MPIA','Simulink(MPIA)'}, ...
                'Location', 'southeast', 'FontSize', 10);
        end
    catch ME
        fprintf('  Simulink仿真出错: %s\n', ME.message);
    end
end

fprintf('\n══════════════════ 全部完成 ══════════════════\n');


%% ═══════════════ 辅助函数 ═══════════════════

%% Ziegler-Nichols PID整定
function [K_zn, fit_zn] = zn_pid_tune()
% ZN整定: 数值搜索临界增益和周期，然后应用ZN-PID公式
% 对于G(s)=10/(s²+10s)，用P-only控制搜索"临界"振荡点

    % 搜索Ku：找到产生约40%超调的Kp作为"临界"增益近似
    Kp_test = logspace(-2, 1, 50);  % 0.01 ~ 10
    best_overshoot_diff = inf;
    best_Kp = 1;
    best_Tu = 0.5;

    for i = 1:length(Kp_test)
        [t, y] = simulate_step_response(Kp_test(i), 0, 0);
        y_ss = y(end);
        y_max = max(y);
        os = (y_max - y_ss) / max(y_ss, 1e-6) * 100;

        % 找第一个峰值时间作为Tu/2的估计
        if os > 10
            [~, peak_idx] = findpeaks(y);
            if ~isempty(peak_idx) && length(peak_idx) >= 2
                Tu_est = t(peak_idx(2)) - t(peak_idx(1));
            elseif ~isempty(peak_idx)
                Tu_est = t(peak_idx(1)) * 2;
            else
                Tu_est = 0.5;
            end

            diff = abs(os - 40);  % 目标: 40%超调
            if diff < best_overshoot_diff
                best_overshoot_diff = diff;
                best_Kp = Kp_test(i);
                best_Tu = Tu_est;
            end
        end
    end

    Ku = best_Kp;
    Tu = best_Tu;

    % ZN-PID公式
    Kp_zn = 0.6 * Ku;
    Ki_zn = 1.2 * Ku / max(Tu, 0.01);
    Kd_zn = 0.075 * Ku * Tu;

    K_zn = [Kp_zn, Ki_zn, Kd_zn];
    fit_zn = pid_fitness(K_zn);
end

%% 阶跃响应仿真（用于绘图）
function [t, y] = simulate_step_response(Kp, Ki, Kd)
% 使用与pid_fitness一致的仿真方式
    dt = 0.002;
    t_end = 2.0;
    t = (0:dt:t_end)';
    n_steps = length(t);

    x1 = 0;  x2 = 0;
    int_e = 0;  e_prev = 1;
    y = zeros(n_steps, 1);
    y(1) = x1;

    for k = 2:n_steps
        e = 1 - x1;
        de = (e - e_prev) / dt;
        u = Kp*e + Ki*int_e + Kd*de;
        u = max(min(u, 1000), -1000);

        [x1, x2] = rk4_step_pid(@motor_dyn, [x1; x2], u, dt);
        int_e = int_e + e * dt;
        e_prev = e;
        y(k) = x1;
    end
end

function dxdt = motor_dyn(state, u)
    dxdt = [state(2); -10*state(2) + 10*u];
end

function [x1, x2] = rk4_step_pid(dyn_func, state, u, dt)
    k1 = dyn_func(state, u);
    k2 = dyn_func(state + 0.5*dt*k1, u);
    k3 = dyn_func(state + 0.5*dt*k2, u);
    k4 = dyn_func(state + dt*k3, u);
    ns = state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x1 = ns(1);  x2 = ns(2);
end

%% Wilcoxon秩和检验（手动实现，兼容无Stats Toolbox）
function [p, h] = manual_ranksum(x, y, alpha)
    if nargin < 3, alpha = 0.05; end
    x = x(:); y = y(:);
    n1 = length(x); n2 = length(y);
    combined = [x; y];
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
    R1 = sum(ranks(1:n1));
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

%% 百分位数计算
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
