function [best_K, best_fit, conv_best, conv_avg] = PID_IA(N_total, T, bounds, params)
% PID_IA  使用标准免疫算法优化PID参数
%   用法: [K, fit, cb, ca] = PID_IA(50, 100, [0,100], params)
%
%   输入:
%     N_total — 种群大小（默认50）
%     T       — 最大迭代次数（默认100）
%     bounds  — 参数搜索范围 [lb, ub]（默认[0, 100]）
%     params  — IA参数结构体（可选，clone_rate等）
%   输出:
%     best_K   — 最优PID参数 [Kp, Ki, Kd]
%     best_fit — 最优ITAE值
%     conv_best — 收敛曲线（最优）
%     conv_avg  — 收敛曲线（平均）

    if nargin < 1 || isempty(N_total), N_total = 50; end
    if nargin < 2 || isempty(T),       T = 100;      end
    if nargin < 3 || isempty(bounds),  bounds = [0, 100]; end
    if nargin < 4 || isempty(params),  params = struct(); end

    D = 3;  % PID三参数

    % 设置目标函数为pid_fitness
    params.obj_func = @pid_fitness;

    % 默认IA参数
    if ~isfield(params, 'clone_rate'),   params.clone_rate = 10;   end
    if ~isfield(params, 'mutate_scale'), params.mutate_scale = 0.5; end
    if ~isfield(params, 'replace_rate'), params.replace_rate = 0.05; end

    [best_K, best_fit, conv_best, conv_avg] = ...
        IA(D, N_total, T, bounds, params);
end
