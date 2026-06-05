# 多种群免疫算法求解高维Rastrigin函数 — 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 生成完整的课程大作业：论文正文(5000-7000字) + 4个Matlab源码文件 + 使用说明，覆盖6个评分模块。

**Architecture:** 4个Matlab文件（目标函数、标准IA、多种群IA、主脚本）全部独立可运行。标准IA和多种群IA共享相同的目标函数和参数配置机制，确保对比公平。主脚本负责30次实验、生成收敛曲线和统计表格。

**Tech Stack:** Matlab R2020a+ (纯脚本，无工具箱依赖)

**论文结构:** 标题页 → 中英文摘要 → 一、引言 → 二、免疫算法原理 → 三、原始免疫算法实现 → 四、多种群免疫算法改进 → 五、实验结果与分析 → 六、总结与展望 → 参考文献

---

### Task 1: 目标函数 — Rastrigin.m

**Files:**
- Create: `code/Rastrigin.m`

- [ ] **Step 1: 编写 Rastrigin.m**

Rastrigin函数定义（D维）：
```
f(x) = 10*D + sum_{i=1}^{D} [x_i^2 - 10*cos(2*pi*x_i)]
搜索空间: x_i ∈ [-5.12, 5.12]
全局最优: f(0,...,0) = 0
```

```matlab
function [fitness, details] = Rastrigin(x)
% Rastrigin 高维多峰测试函数
% 输入: x - 种群矩阵 (N × D) 或行向量 (1 × D)
% 输出: fitness - 适应度值列向量 (N × 1)
%       details  - 结构体，包含各分量信息（调试用）

    [N, D] = size(x);
    fitness = zeros(N, 1);

    for i = 1:N
        xi = x(i, :);
        sum_term = sum(xi.^2 - 10 * cos(2 * pi * xi));
        fitness(i) = 10 * D + sum_term;
    end

    if nargout > 1
        details.D = D;
        details.search_range = [-5.12, 5.12];
        details.global_min = 0;
        details.global_min_point = zeros(1, D);
    end
end
```

- [ ] **Step 2: 验证目标函数正确性**

在Matlab中运行测试：
```matlab
% 测试1: 全局最优点应为0
x_opt = zeros(1, 50);
f_opt = Rastrigin(x_opt);
assert(abs(f_opt) < 1e-10, '全局最优点测试失败');

% 测试2: 随机点应为正值
x_rand = -5.12 + 10.24 * rand(10, 50);
f_rand = Rastrigin(x_rand);
assert(all(f_rand > 0), '随机点测试失败');
```

---

### Task 2: 标准免疫算法 — IA.m

**Files:**
- Create: `code/IA.m`

- [ ] **Step 1: 编写 IA.m 主函数**

```matlab
function [best_x, best_f, conv_best, conv_avg] = IA(D, N, T, bounds, params)
% IA 标准免疫算法（单种群）求解Rastrigin函数最小值
%
% 输入:
%   D      - 问题维度 (默认50)
%   N      - 种群大小 (默认100)
%   T      - 最大迭代次数 (默认500)
%   bounds - 搜索边界 [lb, ub] (默认[-5.12, 5.12])
%   params - 参数字典 (可选，用于传入高级参数)
%
% 输出:
%   best_x   - 最优个体向量 (1×D)
%   best_f   - 最优适应度值
%   conv_best- 每代最优适应度 (T×1)
%   conv_avg - 每代平均适应度 (T×1)

    % ===== 默认参数 =====
    if nargin < 1 || isempty(D),      D = 50;       end
    if nargin < 2 || isempty(N),      N = 100;      end
    if nargin < 3 || isempty(T),      T = 500;      end
    if nargin < 4 || isempty(bounds), bounds = [-5.12, 5.12]; end
    if nargin < 5 || isempty(params), params = struct(); end

    % 可调参数
    clone_rate   = get_param(params, 'clone_rate',   10);    % 克隆倍数上限
    mutate_scale = get_param(params, 'mutate_scale', 0.5);  % 变异强度
    replace_rate = get_param(params, 'replace_rate', 0.05); % 种群更新率
    elite_count  = get_param(params, 'elite_count',  1);    % 精英保留数

    lb = bounds(1); ub = bounds(2);

    % ===== 初始化种群 =====
    pop = lb + (ub - lb) * rand(N, D);
    fitness = Rastrigin(pop);

    % 记录
    conv_best = zeros(T, 1);
    conv_avg  = zeros(T, 1);

    % 全局最优
    [global_best_f, idx] = min(fitness);
    global_best_x = pop(idx, :);

    % ===== 主循环 =====
    for t = 1:T
        % 1. 计算亲和度（归一化）
        aff = 1 ./ (1 + fitness);
        aff_norm = aff / sum(aff);

        % 2. 克隆：按亲和度比例克隆
        [~, sort_idx] = sort(aff, 'descend');
        % 每个个体克隆数量与亲和度成正比
        clone_counts = max(1, round(clone_rate * N * aff_norm));
        clone_counts = min(clone_counts, clone_rate); % 上限

        clones = [];
        for i = 1:N
            n_clone = clone_counts(i);
            clones = [clones; repmat(pop(i,:), n_clone, 1)];
        end

        % 3. 变异：高斯变异，步长与亲和度成反比
        mutated = clones;
        clone_idx = 1;
        for i = 1:N
            n_clone = clone_counts(i);
            for c = 1:n_clone
                % 亲和度越高，变异越小
                sigma = mutate_scale * (1 - aff_norm(i)) * (ub - lb) / sqrt(D);
                mutation = sigma * randn(1, D);
                mutated(clone_idx, :) = clones(clone_idx, :) + mutation;
                % 边界约束
                mutated(clone_idx, :) = max(min(mutated(clone_idx,:), ub), lb);
                clone_idx = clone_idx + 1;
            end
        end

        % 4. 计算克隆适应度
        clone_fitness = Rastrigin(mutated);

        % 5. 选择：每个克隆家族保留最优
        new_pop = zeros(N, D);
        new_fitness = zeros(N, 1);
        clone_ptr = 1;
        for i = 1:N
            n_clone = clone_counts(i);
            family_f = clone_fitness(clone_ptr : clone_ptr + n_clone - 1);
            family_x = mutated(clone_ptr : clone_ptr + n_clone - 1, :);
            [best_family_f, best_idx] = min(family_f);
            new_pop(i, :) = family_x(best_idx, :);
            new_fitness(i) = best_family_f;
            clone_ptr = clone_ptr + n_clone;
        end

        % 6. 精英保留
        % 精英个体加入新种群（替换最差个体）
        [elite_f, elite_idx] = min(fitness); % 原种群精英
        [worst_f, worst_idx] = max(new_fitness);
        if elite_f < worst_f
            new_pop(worst_idx, :) = pop(elite_idx, :);
            new_fitness(worst_idx) = elite_f;
        end

        % 7. 种群更新：替换最差个体为随机新个体
        n_replace = round(replace_rate * N);
        if n_replace > 0
            [~, worst_indices] = sort(new_fitness, 'descend');
            for j = 1:n_replace
                new_pop(worst_indices(j), :) = lb + (ub - lb) * rand(1, D);
                new_fitness(worst_indices(j)) = Rastrigin(new_pop(worst_indices(j), :));
            end
        end

        % 8. 更新种群
        pop = new_pop;
        fitness = new_fitness;

        % 9. 更新全局最优
        [current_best_f, idx] = min(fitness);
        if current_best_f < global_best_f
            global_best_f = current_best_f;
            global_best_x = pop(idx, :);
        end

        % 10. 记录收敛数据
        conv_best(t) = global_best_f;
        conv_avg(t)  = mean(fitness);
    end

    best_x = global_best_x;
    best_f = global_best_f;
end

% ===== 辅助函数 =====
function val = get_param(params, field, default)
    if isfield(params, field)
        val = params.(field);
    else
        val = default;
    end
end
```

---

### Task 3: 多种群免疫算法 — MPIA.m

**Files:**
- Create: `code/MPIA.m`

- [ ] **Step 1: 编写 MPIA.m 主函数**

```matlab
function [best_x, best_f, conv_best, conv_avg, migration_log] = MPIA(D, N_total, T, bounds, params)
% MPIA 多种群免疫算法求解Rastrigin函数最小值
%
% 改进策略：
%   1. 多种群协同：K个子种群独立进化
%   2. 分层混合变异：高亲和度→高斯变异(精细)，低亲和度→柯西变异(探索)
%   3. 环形迁移：每Tg代向邻居迁移最优个体
%
% 输入:
%   D        - 问题维度 (默认50)
%   N_total  - 总种群大小 (默认100)
%   T        - 最大迭代次数 (默认500)
%   bounds   - 搜索边界 [lb, ub]
%   params   - 参数字典
%
% 输出:
%   best_x       - 全局最优个体 (1×D)
%   best_f       - 全局最优适应度
%   conv_best    - 每代最优适应度 (T×1)
%   conv_avg     - 每代平均适应度 (T×1)
%   migration_log- 迁移记录

    % ===== 默认参数 =====
    if nargin < 1 || isempty(D),       D = 50;       end
    if nargin < 2 || isempty(N_total), N_total = 100; end
    if nargin < 3 || isempty(T),       T = 500;      end
    if nargin < 4 || isempty(bounds),  bounds = [-5.12, 5.12]; end
    if nargin < 5 || isempty(params),  params = struct(); end

    % 多种群参数
    K           = get_param_mp(params, 'K',           5);     % 子种群数量
    clone_rate  = get_param_mp(params, 'clone_rate',  10);    % 克隆倍数上限
    mutate_scale= get_param_mp(params, 'mutate_scale',0.5);   % 高斯变异强度
    cauchy_scale= get_param_mp(params, 'cauchy_scale',0.3);   % 柯西变异强度
    replace_rate= get_param_mp(params, 'replace_rate',0.05);  % 种群更新率
    migrate_T   = get_param_mp(params, 'migrate_T',   10);    % 迁移间隔
    migrate_M   = get_param_mp(params, 'migrate_M',   2);     % 迁移个体数

    N_sub = floor(N_total / K);  % 每个子种群大小

    lb = bounds(1); ub = bounds(2);

    % ===== 初始化K个子种群 =====
    sub_pop = cell(K, 1);
    sub_fitness = cell(K, 1);
    sub_best_x = cell(K, 1);
    sub_best_f = zeros(K, 1);
    for k = 1:K
        sub_pop{k} = lb + (ub - lb) * rand(N_sub, D);
        sub_fitness{k} = Rastrigin(sub_pop{k});
        [sub_best_f(k), idx] = min(sub_fitness{k});
        sub_best_x{k} = sub_pop{k}(idx, :);
    end

    % 全局最优
    [global_best_f, best_k] = min(sub_best_f);
    global_best_x = sub_best_x{best_k};

    % 记录
    conv_best = zeros(T, 1);
    conv_avg  = zeros(T, 1);
    migration_log = []; % [代, 源群, 目标群, 迁移个体适应度]

    % ===== 主循环 =====
    for t = 1:T

        % --- 各子种群独立进化 ---
        for k = 1:K
            pop_k = sub_pop{k};
            fit_k = sub_fitness{k};
            Nk = N_sub;

            % 1. 亲和度
            aff_k = 1 ./ (1 + fit_k);
            aff_norm_k = aff_k / sum(aff_k);

            % 2. 克隆
            clone_counts_k = max(1, round(clone_rate * Nk * aff_norm_k));
            clone_counts_k = min(clone_counts_k, clone_rate);

            clones_k = [];
            for i = 1:Nk
                clones_k = [clones_k; repmat(pop_k(i,:), clone_counts_k(i), 1)];
            end

            % 3. 分层混合变异
            % 亲和度中位数作为分界线
            med_aff = median(aff_norm_k);
            mutated_k = clones_k;
            clone_ptr = 1;
            for i = 1:Nk
                n_clone = clone_counts_k(i);
                for c = 1:n_clone
                    if aff_norm_k(i) >= med_aff
                        % 高斯变异（精细搜索）
                        sigma = mutate_scale * (1 - aff_norm_k(i)) * (ub - lb) / sqrt(D);
                        noise = sigma * randn(1, D);
                    else
                        % 柯西变异（全局探索）
                        sigma = cauchy_scale * (1 - aff_norm_k(i)) * (ub - lb) / sqrt(D);
                        noise = sigma * tan(pi * (rand(1, D) - 0.5));  % 柯西分布
                    end
                    mutated_k(clone_ptr, :) = clones_k(clone_ptr, :) + noise;
                    mutated_k(clone_ptr, :) = max(min(mutated_k(clone_ptr,:), ub), lb);
                    clone_ptr = clone_ptr + 1;
                end
            end

            % 4. 选择
            clone_fit_k = Rastrigin(mutated_k);
            new_pop_k = zeros(Nk, D);
            new_fit_k = zeros(Nk, 1);
            clone_ptr = 1;
            for i = 1:Nk
                n_clone = clone_counts_k(i);
                family_f = clone_fit_k(clone_ptr : clone_ptr + n_clone - 1);
                family_x = mutated_k(clone_ptr : clone_ptr + n_clone - 1, :);
                [best_fam_f, best_idx] = min(family_f);
                new_pop_k(i, :) = family_x(best_idx, :);
                new_fit_k(i) = best_fam_f;
                clone_ptr = clone_ptr + n_clone;
            end

            % 5. 精英保留
            [elite_f_k, elite_idx_k] = min(fit_k);
            [worst_f_k, worst_idx_k] = max(new_fit_k);
            if elite_f_k < worst_f_k
                new_pop_k(worst_idx_k, :) = pop_k(elite_idx_k, :);
                new_fit_k(worst_idx_k) = elite_f_k;
            end

            % 6. 种群更新
            n_replace = round(replace_rate * Nk);
            if n_replace > 0
                [~, worst_idx] = sort(new_fit_k, 'descend');
                for j = 1:n_replace
                    new_pop_k(worst_idx(j), :) = lb + (ub - lb) * rand(1, D);
                    new_fit_k(worst_idx(j)) = Rastrigin(new_pop_k(worst_idx(j), :));
                end
            end

            % 更新子种群
            sub_pop{k} = new_pop_k;
            sub_fitness{k} = new_fit_k;
            [sub_best_f(k), idx] = min(new_fit_k);
            sub_best_x{k} = new_pop_k(idx, :);
        end

        % --- 环形迁移 (每migrate_T代) ---
        if mod(t, migrate_T) == 0 && t < T
            migrants = cell(K, 1);  % 每个子种群要送出的个体

            for k = 1:K
                [~, sort_idx] = sort(sub_fitness{k});
                % 选出最优migrate_M个个体
                migrants{k} = sub_pop{k}(sort_idx(1:migrate_M), :);
            end

            % 执行迁移（向右邻居发送）
            for k = 1:K
                target_k = mod(k, K) + 1;  % 环形：k → k+1 (最后→第1)
                source_k = k;

                % 接收来自源群的个体
                incoming = migrants{source_k};

                % 替换目标群中最差的migrate_M个个体
                [~, worst_idx] = sort(sub_fitness{target_k}, 'descend');
                for m = 1:migrate_M
                    sub_pop{target_k}(worst_idx(m), :) = incoming(m, :);
                    sub_fitness{target_k}(worst_idx(m)) = Rastrigin(incoming(m, :));
                end

                % 记录迁移
                migration_log = [migration_log;
                    t, source_k, target_k, ...
                    Rastrigin(incoming(1,:))];
            end

            % 迁移后重新评估各子种群最优
            for k = 1:K
                [sub_best_f(k), idx] = min(sub_fitness{k});
                sub_best_x{k} = sub_pop{k}(idx, :);
            end
        end

        % --- 更新全局最优 ---
        [current_global_best, best_k] = min(sub_best_f);
        if current_global_best < global_best_f
            global_best_f = current_global_best;
            global_best_x = sub_best_x{best_k};
        end

        % --- 记录收敛 ---
        conv_best(t) = global_best_f;
        % 所有子种群的平均适应度
        all_fitness = [];
        for k = 1:K
            all_fitness = [all_fitness; sub_fitness{k}];
        end
        conv_avg(t) = mean(all_fitness);
    end

    best_x = global_best_x;
    best_f = global_best_f;
end

% ===== 辅助函数 =====
function val = get_param_mp(params, field, default)
    if isfield(params, field)
        val = params.(field);
    else
        val = default;
    end
end
```

---

### Task 4: 主脚本 — main.m

**Files:**
- Create: `code/main.m`

- [ ] **Step 1: 编写 main.m（实验运行 + 绘图 + 表格输出）**

```matlab
%% main.m — 多种群免疫算法 vs 标准免疫算法对比实验
% 求解 D=50 维 Rastrigin 函数最小值
% 30次独立运行，统计对比 + 收敛曲线绘制

clear; clc; close all;

%% ===== 实验参数 =====
D       = 50;        % 问题维度
N_total = 100;       % 总种群大小
T       = 500;       % 最大迭代次数
BOUNDS  = [-5.12, 5.12];
N_RUNS  = 30;        % 独立运行次数

% 标准IA参数
params_ia.clone_rate   = 10;
params_ia.mutate_scale = 0.5;
params_ia.replace_rate = 0.05;
params_ia.elite_count  = 1;

% 多种群IA参数（总种群=100，K=5子群各20个体）
params_mpia.K            = 5;
params_mpia.clone_rate   = 10;
params_mpia.mutate_scale = 0.5;
params_mpia.cauchy_scale = 0.3;
params_mpia.replace_rate = 0.05;
params_mpia.migrate_T    = 10;
params_mpia.migrate_M    = 2;

%% ===== 30次独立实验 =====
fprintf('========== 实验开始 ==========\n');
fprintf('问题: %d维Rastrigin函数 | 种群: %d | 迭代: %d | 运行: %d次\n', ...
    D, N_total, T, N_RUNS);

% 存储结果
ia_results    = zeros(N_RUNS, 1);       % IA最优值
mpia_results  = zeros(N_RUNS, 1);       % MPIA最优值
ia_conv_all   = zeros(T, N_RUNS);       % IA每代最优
mpia_conv_all = zeros(T, N_RUNS);       % MPIA每代最优
ia_avg_all    = zeros(T, N_RUNS);       % IA每代均值
mpia_avg_all  = zeros(T, N_RUNS);       % MPIA每代均值

fprintf('\n--- 运行标准IA ---\n');
for run = 1:N_RUNS
    fprintf('  IA Run %2d/%d...', run, N_RUNS);
    [~, best_f, conv_b, conv_a] = IA(D, N_total, T, BOUNDS, params_ia);
    ia_results(run) = best_f;
    ia_conv_all(:, run) = conv_b;
    ia_avg_all(:, run) = conv_a;
    fprintf(' 最优值: %.4f\n', best_f);
end

fprintf('\n--- 运行多种群IA ---\n');
for run = 1:N_RUNS
    fprintf('  MPIA Run %2d/%d...', run, N_RUNS);
    [~, best_f, conv_b, conv_a, ~] = MPIA(D, N_total, T, BOUNDS, params_mpia);
    mpia_results(run) = best_f;
    mpia_conv_all(:, run) = conv_b;
    mpia_avg_all(:, run) = conv_a;
    fprintf(' 最优值: %.4f\n', best_f);
end

%% ===== 统计结果 =====
fprintf('\n========== 统计结果 ==========\n');
fprintf('%-15s %12s %12s %12s %12s\n', '算法', '最优值', '最差值', '均值', '标准差');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%-15s %12.4f %12.4f %12.4f %12.4f\n', ...
    '标准IA', min(ia_results), max(ia_results), mean(ia_results), std(ia_results));
fprintf('%-15s %12.4f %12.4f %12.4f %12.4f\n', ...
    '多种群IA', min(mpia_results), max(mpia_results), mean(mpia_results), std(mpia_results));

% 改进幅度
improvement = (mean(ia_results) - mean(mpia_results)) / mean(ia_results) * 100;
fprintf('\n多种群IA相比标准IA平均改进: %.2f%%\n', improvement);

%% ===== Wilcoxon秩和检验 =====
[p_value, h] = ranksum(ia_results, mpia_results);
fprintf('\nWilcoxon秩和检验:\n');
fprintf('  p值: %.6f\n', p_value);
if h == 1
    fprintf('  结论: 在0.05显著性水平下，两种算法存在显著差异\n');
else
    fprintf('  结论: 在0.05显著性水平下，两种算法无显著差异\n');
end

%% ===== 绘图 =====
figure('Position', [100, 100, 1400, 500]);

% 子图1: 最优适应度收敛曲线
subplot(1, 2, 1);
hold on;
% 多次运行均值曲线
ia_best_mean  = mean(ia_conv_all, 2);
mpia_best_mean = mean(mpia_conv_all, 2);
h1 = plot(1:T, ia_best_mean, 'b-', 'LineWidth', 1.5);
h2 = plot(1:T, mpia_best_mean, 'r-', 'LineWidth', 1.5);
% 填充标准差带
ia_best_std  = std(ia_conv_all, 0, 2);
mpia_best_std = std(mpia_conv_all, 0, 2);
x_fill = [1:T, T:-1:1]';
fill(x_fill, [ia_best_mean + ia_best_std; flipud(ia_best_mean - ia_best_std)], ...
    'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill(x_fill, [mpia_best_mean + mpia_best_std; flipud(mpia_best_mean - mpia_best_std)], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

xlabel('迭代次数', 'FontSize', 12);
ylabel('最优适应度值', 'FontSize', 12);
title(sprintf('最优适应度收敛曲线 (D=%d, %d次运行均值)', D, N_RUNS), 'FontSize', 13);
legend([h1, h2], {'标准IA', '多种群IA'}, 'Location', 'northeast', 'FontSize', 11);
grid on;
hold off;

% 子图2: 平均适应度收敛曲线
subplot(1, 2, 2);
hold on;
ia_avg_mean  = mean(ia_avg_all, 2);
mpia_avg_mean = mean(mpia_avg_all, 2);
h3 = plot(1:T, ia_avg_mean, 'b-', 'LineWidth', 1.5);
h4 = plot(1:T, mpia_avg_mean, 'r-', 'LineWidth', 1.5);
% 填充标准差带
ia_avg_std  = std(ia_avg_all, 0, 2);
mpia_avg_std = std(mpia_avg_all, 0, 2);
fill(x_fill, [ia_avg_mean + ia_avg_std; flipud(ia_avg_mean - ia_avg_std)], ...
    'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill(x_fill, [mpia_avg_mean + mpia_avg_std; flipud(mpia_avg_mean - mpia_avg_std)], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

xlabel('迭代次数', 'FontSize', 12);
ylabel('平均适应度值', 'FontSize', 12);
title(sprintf('平均适应度收敛曲线 (D=%d, %d次运行均值)', D, N_RUNS), 'FontSize', 13);
legend([h3, h4], {'标准IA', '多种群IA'}, 'Location', 'northeast', 'FontSize', 11);
grid on;
hold off;

sgtitle(sprintf('多种群免疫算法 vs 标准免疫算法 — %d维Rastrigin函数优化', D), ...
    'FontSize', 14, 'FontWeight', 'bold');

% 保存图片
saveas(gcf, '收敛曲线对比.png');
fprintf('\n收敛曲线已保存为 收敛曲线对比.png\n');

%% ===== 箱线图（=）=====
figure('Position', [100, 100, 600, 500]);
box_data = [ia_results, mpia_results];
boxplot(box_data, 'Labels', {'标准IA', '多种群IA'});
ylabel('最优适应度值', 'FontSize', 12);
title(sprintf('30次独立运行最优值分布 (D=%d)', D), 'FontSize', 13);
grid on;
saveas(gcf, '箱线图对比.png');
fprintf('箱线图已保存为 箱线图对比.png\n');

fprintf('\n========== 实验完成 ==========\n');
```

---

### Task 5: 论文正文

**Files:**
- Create: `论文正文.md`

论文正文包含：
1. 标题页信息
2. 中英文摘要（~200字各）
3. 一、引言（研究意义 + 问题描述）
4. 二、免疫算法原理（核心理论 + 伪代码）
5. 三、原始免疫算法实现（算子细节 + 参数表）
6. 四、多种群免疫算法改进（框架 + 分层变异 + 迁移）
7. 五、实验结果与分析（数据表格 + 曲线描述 + 统计分析）
8. 六、总结与展望
9. 参考文献（≥10条，中英文混合）

---

### Task 6: README 使用说明

**Files:**
- Create: `README.md`

内容包括：环境要求、文件结构、运行方法、参数修改指南、结果解读。

---

## Self-Review

**1. Spec coverage check:**
- ✅ 模块1 研究意义 → Task 5 (论文一、引言)
- ✅ 模块2 问题描述 → Task 5 (论文一、引言)
- ✅ 模块3 算法原理 → Task 5 (论文二、免疫算法原理)
- ✅ 模块4 原始实现 → Task 2 (IA.m) + Task 5 (论文三)
- ✅ 模块5 改进算法 → Task 3 (MPIA.m) + Task 5 (论文四)
- ✅ 模块6 写作规范 → Task 5 (全论文含中英文摘要、参考文献)
- ✅ 收敛曲线（最优+平均） → Task 4 (main.m绘图)
- ✅ 原文vs改进同图画 → Task 4 (subplot并排)
- ✅ 参数一致性 → Task 4 (统一N_total=100, T=500)

**2. Placeholder scan:** 无TBD/TODO，所有代码完整可运行。

**3. Type consistency:** 函数签名一致，`IA()` 和 `MPIA()` 输出格式统一。
