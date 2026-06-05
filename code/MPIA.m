function [best_x, best_f, conv_best, conv_avg, migration_log] = MPIA(D, N_total, T, bounds, params)
% MPIA  多种群免疫算法: K子种群 + 分层混合变异 + 环形迁移
% 用法: [x,f,cb,ca,log] = MPIA(50,100,500,[-5.12,5.12],params)

    if nargin < 1 || isempty(D),       D = 50;       end
    if nargin < 2 || isempty(N_total), N_total = 100; end
    if nargin < 3 || isempty(T),       T = 500;      end
    if nargin < 4 || isempty(bounds),  bounds = [-5.12, 5.12]; end
    if nargin < 5 || isempty(params),  params = struct(); end

    K            = getField(params, 'K',            5);
    clone_rate   = getField(params, 'clone_rate',   10);
    mutate_scale = getField(params, 'mutate_scale', 0.5);
    cauchy_scale = getField(params, 'cauchy_scale', 0.3);
    replace_rate = getField(params, 'replace_rate', 0.05);
    migrate_T    = getField(params, 'migrate_T',    10);
    migrate_M    = getField(params, 'migrate_M',    2);

    N_sub = floor(N_total / K);
    lb = bounds(1);  ub = bounds(2);

    % 初始化K个子种群
    sub_pop = cell(K,1);  sub_fitness = cell(K,1);
    sub_best_x = cell(K,1);  sub_best_f = zeros(K,1);
    for k = 1:K
        sub_pop{k} = lb + (ub - lb) * rand(N_sub, D);
        sub_fitness{k} = Rastrigin(sub_pop{k});
        [sub_best_f(k), idx] = min(sub_fitness{k});
        sub_best_x{k} = sub_pop{k}(idx, :);
    end
    [global_best_f, best_k] = min(sub_best_f);
    global_best_x = sub_best_x{best_k};

    conv_best = zeros(T, 1);  conv_avg = zeros(T, 1);
    migration_log = zeros(0, 4);

    % 主循环
    for t = 1:T
        % --- 各子种群独立免疫进化 ---
        for k = 1:K
            pop_k = sub_pop{k};  fit_k = sub_fitness{k};  Nk = N_sub;

            % 亲和度 + 克隆
            aff_norm_k = (1 ./ (1 + fit_k));
            aff_norm_k = aff_norm_k / sum(aff_norm_k);
            clone_cnt = max(1, min(round(clone_rate * Nk * aff_norm_k), clone_rate));

            total = sum(clone_cnt);
            clones_k = zeros(total, D);  ptr = 1;
            for i = 1:Nk
                n = clone_cnt(i);
                clones_k(ptr:ptr+n-1, :) = repmat(pop_k(i,:), n, 1);
                ptr = ptr + n;
            end

            % 分层混合变异: 亲和度>=中位数→高斯, <中位数→柯西
            med_aff = median(aff_norm_k);
            mutated_k = clones_k;  ptr = 1;
            for i = 1:Nk
                n = clone_cnt(i);
                for c = 1:n
                    if aff_norm_k(i) >= med_aff
                        sigma = mutate_scale * (1-aff_norm_k(i)) * (ub-lb)/sqrt(D);
                        noise = sigma * randn(1, D);
                    else
                        sigma = cauchy_scale * (1-aff_norm_k(i)) * (ub-lb)/sqrt(D);
                        noise = sigma * tan(pi*(rand(1,D)-0.5));  % 柯西采样
                    end
                    mutated_k(ptr, :) = max(min(mutated_k(ptr,:)+noise, ub), lb);
                    ptr = ptr + 1;
                end
            end

            % 评估 + 家族选择 + 精英保留 + 随机更新
            clone_fit = Rastrigin(mutated_k);
            new_pop_k = zeros(Nk, D);  new_fit_k = zeros(Nk, 1);  ptr = 1;
            for i = 1:Nk
                n = clone_cnt(i);
                [best_val, best_i] = min(clone_fit(ptr:ptr+n-1));
                new_pop_k(i, :) = mutated_k(ptr+best_i-1, :);
                new_fit_k(i) = best_val;  ptr = ptr + n;
            end
            [elite_f, elite_i] = min(fit_k);
            [worst_f, worst_i] = max(new_fit_k);
            if elite_f < worst_f
                new_pop_k(worst_i,:) = pop_k(elite_i,:);  new_fit_k(worst_i) = elite_f;
            end
            n_rep = round(replace_rate * Nk);
            if n_rep > 0
                [~, wi] = sort(new_fit_k, 'descend');
                for j = 1:n_rep
                    new_pop_k(wi(j),:) = lb+(ub-lb)*rand(1,D);
                    new_fit_k(wi(j)) = Rastrigin(new_pop_k(wi(j),:));
                end
            end
            sub_pop{k} = new_pop_k;  sub_fitness{k} = new_fit_k;
            [sub_best_f(k), idx] = min(new_fit_k);
            sub_best_x{k} = new_pop_k(idx, :);
        end

        % --- 环形迁移 ---
        if mod(t, migrate_T) == 0 && t < T
            migrants = cell(K, 1);
            for k = 1:K
                [~, si] = sort(sub_fitness{k});
                migrants{k} = sub_pop{k}(si(1:migrate_M), :);
            end
            for k = 1:K
                dst = mod(k, K) + 1;
                incoming = migrants{k};
                [~, wi] = sort(sub_fitness{dst}, 'descend');
                for m = 1:migrate_M
                    sub_pop{dst}(wi(m), :) = incoming(m, :);
                    sub_fitness{dst}(wi(m)) = Rastrigin(incoming(m, :));
                end
                migration_log(end+1,:) = [t, k, dst, Rastrigin(incoming(1,:))]; %#ok<AGROW>
            end
            for k = 1:K
                [sub_best_f(k), idx] = min(sub_fitness{k});
                sub_best_x{k} = sub_pop{k}(idx, :);
            end
        end

        % 更新全局最优
        [cur_best, best_k] = min(sub_best_f);
        if cur_best < global_best_f
            global_best_f = cur_best;  global_best_x = sub_best_x{best_k};
        end
        conv_best(t) = global_best_f;
        all_fit = cell2mat(sub_fitness);
        conv_avg(t) = mean(all_fit);
    end

    best_x = global_best_x;  best_f = global_best_f;
end

function value = getField(params, fieldName, defaultValue)
    if isfield(params, fieldName), value = params.(fieldName); else, value = defaultValue; end
end
