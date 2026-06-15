function [best_x, best_f, conv_best, conv_avg] = IA(D, N, T, bounds, params)
% IA  标准免疫算法: 克隆选择 + 高斯变异 + 精英保留
% 用法: [x,f,cb,ca] = IA(50,100,500,[-5.12,5.12],params)

    if nargin < 1 || isempty(D),      D = 50;       end
    if nargin < 2 || isempty(N),      N = 100;      end
    if nargin < 3 || isempty(T),      T = 500;      end
    if nargin < 4 || isempty(bounds), bounds = [-5.12, 5.12]; end
    if nargin < 5 || isempty(params), params = struct(); end

    clone_rate   = getField(params, 'clone_rate',   10);
    mutate_scale = getField(params, 'mutate_scale', 0.5);
    replace_rate = getField(params, 'replace_rate', 0.05);
    obj_func = getField(params, 'obj_func', @Rastrigin);
    lb = bounds(1);  ub = bounds(2);

    % 初始化种群
    pop = lb + (ub - lb) * rand(N, D);
    fitness = obj_func(pop);
    conv_best = zeros(T, 1);  conv_avg = zeros(T, 1);
    [global_best_f, idx] = min(fitness);
    global_best_x = pop(idx, :);

    % 主循环
    for t = 1:T
        % 1. 亲和度（f越小亲和度越高）
        aff = 1 ./ (1 + fitness);
        aff_norm = aff / sum(aff);

        % 2. 按亲和度比例克隆
        clone_counts = min(round(clone_rate * N * aff_norm), clone_rate);
        clone_counts = max(clone_counts, 1);

        total_clones = sum(clone_counts);
        clones = zeros(total_clones, D);
        clone_ptr = 1;
        for i = 1:N
            n = clone_counts(i);
            clones(clone_ptr : clone_ptr + n - 1, :) = repmat(pop(i,:), n, 1);
            clone_ptr = clone_ptr + n;
        end

        % 3. 高斯变异（亲和度越高变异越小）
        mutated = clones;
        clone_ptr = 1;
        for i = 1:N
            n = clone_counts(i);
            sigma = mutate_scale * (1 - aff_norm(i)) * (ub - lb) / sqrt(D);
            for c = 1:n
                mutated(clone_ptr, :) = clones(clone_ptr, :) + sigma * randn(1, D);
                % 边界约束
                mutated(clone_ptr, :) = max(min(mutated(clone_ptr,:), ub), lb);
                clone_ptr = clone_ptr + 1;
            end
        end

        % 4. 评估克隆体
        clone_fitness = obj_func(mutated);

        % 5. 克隆家族内精英选择
        new_pop = zeros(N, D);  new_fitness = zeros(N, 1);
        clone_ptr = 1;
        for i = 1:N
            n = clone_counts(i);
            idx_c = clone_ptr : clone_ptr + n - 1;
            [best_val, best_i] = min(clone_fitness(idx_c));
            new_pop(i, :) = mutated(clone_ptr + best_i - 1, :);
            new_fitness(i) = best_val;
            clone_ptr = clone_ptr + n;
        end

        % 6. 精英保留
        [elite_f, elite_i] = min(fitness);
        [worst_f, worst_i] = max(new_fitness);
        if elite_f < worst_f
            new_pop(worst_i, :) = pop(elite_i, :);
            new_fitness(worst_i) = elite_f;
        end

        % 7. 随机替换最差个体维持多样性
        n_replace = round(replace_rate * N);
        if n_replace > 0
            [~, wi] = sort(new_fitness, 'descend');
            for j = 1:n_replace
                new_pop(wi(j), :) = lb + (ub - lb) * rand(1, D);
                new_fitness(wi(j)) = obj_func(new_pop(wi(j), :));
            end
        end

        % 8. 种群更新 & 更新全局最优
        pop = new_pop;  fitness = new_fitness;
        [current_best_f, idx] = min(fitness);
        if current_best_f < global_best_f
            global_best_f = current_best_f;
            global_best_x = pop(idx, :);
        end
        conv_best(t) = global_best_f;
        conv_avg(t)  = mean(fitness);
    end

    best_x = global_best_x;
    best_f = global_best_f;
end

function value = getField(params, fieldName, defaultValue)
    if isfield(params, fieldName)
        value = params.(fieldName);
    else
        value = defaultValue;
    end
end
