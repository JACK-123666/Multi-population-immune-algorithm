function fitness = pid_fitness(x)
% pid_fitness  PID控制器的适应度函数（ITAE指标）
%   输入: x — N×3矩阵，每行 [Kp, Ki, Kd]
%   输出: fitness — N×1向量，越小越好
%
%   被控对象: 直流电机 G(s) = 10/(s²+10s)
%   控制目标: 单位阶跃响应 — ITAE最小化 + 超调惩罚
%
%   自动检测环境: Simulink > Control System Toolbox > 手动仿真

    N = size(x, 1);
    fitness = zeros(N, 1);

    % 检测可用仿真后端
    sim_backend = detect_backend();

    for i = 1:N
        Kp = x(i, 1);
        Ki = x(i, 2);
        Kd = x(i, 3);

        % 负参数惩罚（PID参数不应为负）
        if Kp < 0 || Ki < 0 || Kd < 0
            fitness(i) = 1e10 + abs(Kp) + abs(Ki) + abs(Kd);
            continue;
        end

        try
            switch sim_backend
                case 'simulink'
                    fitness(i) = simulate_simulink(Kp, Ki, Kd);
                case 'cst'
                    fitness(i) = simulate_cst(Kp, Ki, Kd);
                case 'manual'
                    fitness(i) = simulate_manual(Kp, Ki, Kd);
            end
        catch
            fitness(i) = 1e10;  % 仿真失败，给大惩罚
        end
    end
end

%% 环境检测
function backend = detect_backend()
% 检测可用的仿真后端 (优先级: Simulink > manual > CST)
    backend = 'manual';  % 默认手动仿真(RK4), 最可靠

    % 检测 Simulink (需要模型文件存在)
    try
        v = ver('simulink');
        if ~isempty(v) && exist('DC_Motor_Model.slx', 'file')
            backend = 'simulink';
            return;
        end
    catch
    end

    % 手动仿真已验证可靠，CST仅作备选
    % Control System Toolbox 对含积分器的系统有 step() 兼容性问题
end

%% Control System Toolbox 仿真
function itae = simulate_cst(Kp, Ki, Kd)
% 使用 tf + feedback + step 进行仿真
    s = tf('s');

    % 被控对象: G(s) = 10/(s²+10s)
    G = 10 / (s^2 + 10*s);

    % PID控制器: C(s) = Kp + Ki/s + Kd*s
    C = Kp + Ki/s + Kd*s;

    % 闭环传递函数
    T = feedback(C * G, 1);

    % 仿真时间
    t_end = 2.0;
    t = linspace(0, t_end, 500)';

    % 阶跃响应（用 lsim 避免 step() 的积分器初始化问题）
    u = ones(size(t));
    [y, t] = lsim(T, u, t);

    % 计算ITAE = ∫ t * |error| dt
    error = 1 - y;
    itae_raw = trapz(t, t .* abs(error));

    % 超调量
    y_ss = y(end);
    y_max = max(y);
    overshoot = max(0, (y_max - y_ss) / max(y_ss, 1e-6) * 100);

    % 超调惩罚
    penalty = 1.0;
    if overshoot > 30
        penalty = 1.0 + 0.1 * (overshoot - 30);
    end
    if overshoot > 50
        penalty = penalty + 0.5 * (overshoot - 50);
    end

    % 稳态误差惩罚
    ss_error = abs(1 - y_ss);
    if ss_error > 0.02
        penalty = penalty + ss_error * 100;
    end

    % 上升时间惩罚（太慢不好）
    try
        rise_idx = find(y >= 0.9 * y_ss, 1);
        if ~isempty(rise_idx)
            rise_time = t(rise_idx);
        else
            rise_time = t_end;
        end
    catch
        rise_time = t_end;
    end
    if rise_time > 1.0
        penalty = penalty + (rise_time - 1.0) * 2;
    end

    itae = itae_raw * penalty;
end

%% 手动仿真（无工具箱备选）
function itae = simulate_manual(Kp, Ki, Kd)
% 使用四阶龙格-库塔法手动仿真闭环系统
% 状态空间: G(s)=10/(s²+10s) → x1'=x2, x2'=-10*x2+10*u
% PID: u = Kp*e + Ki*∫e*dt + Kd*de/dt

    dt = 0.002;
    t_end = 2.0;
    t = 0:dt:t_end;
    n_steps = length(t);

    % 状态变量 [x1; x2; integral_error]
    x1 = 0;      % 输出 y
    x2 = 0;      % 输出导数
    int_e = 0;   % 误差积分
    e_prev = 1;  % 上一步误差

    y = zeros(n_steps, 1);
    y(1) = x1;

    for k = 2:n_steps
        % 当前误差
        e = 1 - x1;

        % PID 输出（用后向差分近似微分）
        de = (e - e_prev) / dt;
        u = Kp * e + Ki * int_e + Kd * de;

        % 限幅（防止过大控制量）
        u = max(min(u, 1000), -1000);

        % RK4 积分一步
        [x1, x2] = rk4_step(@motor_dynamics, [x1; x2], u, dt);

        % 更新积分和误差记忆
        int_e = int_e + e * dt;
        e_prev = e;

        y(k) = x1;
    end

    % 计算ITAE
    error = 1 - y;
    itae_raw = trapz(t, t .* abs(error));

    % 超调量
    y_ss = y(end);
    y_max = max(y);
    overshoot = max(0, (y_max - y_ss) / max(y_ss, 1e-6) * 100);

    % 惩罚项
    penalty = 1.0;
    if overshoot > 30
        penalty = 1.0 + 0.1 * (overshoot - 30);
    end
    if overshoot > 50
        penalty = penalty + 0.5 * (overshoot - 50);
    end

    ss_error = abs(1 - y_ss);
    if ss_error > 0.02
        penalty = penalty + ss_error * 100;
    end

    try
        rise_idx = find(y >= 0.9 * y_ss, 1);
        if ~isempty(rise_idx)
            rise_time = t(rise_idx);
        else
            rise_time = t_end;
        end
    catch
        rise_time = t_end;
    end
    if rise_time > 1.0
        penalty = penalty + (rise_time - 1.0) * 2;
    end

    itae = itae_raw * penalty;
end

%% 直流电机动力学
function dxdt = motor_dynamics(state, u)
    x1 = state(1);
    x2 = state(2);

    dx1 = x2;
    dx2 = -10 * x2 + 10 * u;

    dxdt = [dx1; dx2];
end

%% 四阶龙格-库塔积分
function [x1_new, x2_new] = rk4_step(dyn_func, state, u, dt)
    k1 = dyn_func(state, u);
    k2 = dyn_func(state + 0.5*dt*k1, u);
    k3 = dyn_func(state + 0.5*dt*k2, u);
    k4 = dyn_func(state + dt*k3, u);

    new_state = state + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    x1_new = new_state(1);
    x2_new = new_state(2);
end

%% Simulink 仿真（预留接口）
function itae = simulate_simulink(Kp, Ki, Kd)
% 使用 Simulink 模型 DC_Motor_Model.slx 仿真
% 如果模型不存在，降级到 Control System Toolbox
    model_name = 'DC_Motor_Model';

    try
        % 检查模型是否已加载
        if ~bdIsLoaded(model_name)
            load_system(model_name);
        end

        % 将参数写入工作区
        assignin('base', 'Kp_sim', Kp);
        assignin('base', 'Ki_sim', Ki);
        assignin('base', 'Kd_sim', Kd);

        % 运行仿真
        simOut = sim(model_name, 'SrcWorkspace', 'base', ...
                     'StopTime', '2', 'SaveOutput', 'on');

        % 提取输出
        y = simOut.get('yout');
        if isempty(y)
            error('Simulink输出为空');
        end
        y_data = y.signals.values;
        t_data = y.time;

        % 计算ITAE
        error = 1 - y_data;
        itae_raw = trapz(t_data, t_data .* abs(error));

        y_ss = y_data(end);
        y_max = max(y_data);
        overshoot = max(0, (y_max - y_ss) / max(y_ss, 1e-6) * 100);

        penalty = 1.0;
        if overshoot > 30
            penalty = 1.0 + 0.1 * (overshoot - 30);
        end
        if overshoot > 50
            penalty = penalty + 0.5 * (overshoot - 50);
        end

        ss_error = abs(1 - y_ss);
        if ss_error > 0.02
            penalty = penalty + ss_error * 100;
        end

        itae = itae_raw * penalty;

    catch ME
        % Simulink 仿真失败，静默降级到 CST
        persistent sim_warned;
        if isempty(sim_warned)
            fprintf('  [提示] Simulink仿真不可用，使用CST纯脚本模式\n');
            sim_warned = true;
        end
        itae = simulate_cst(Kp, Ki, Kd);
    end
end
