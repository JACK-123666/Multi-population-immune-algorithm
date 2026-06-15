function build_simulink_model()
% build_simulink_model  程序化构建 DC Motor PID 控制 Simulink 模型
%
%   需要: Simulink + Simulink Control Design Toolbox
%   生成: DC_Motor_Model.slx
%
%   模型结构:
%     Step → Sum → PID(s) → G(s)=10/(s²+10s) → To Workspace
%              ↑__________________________________|

    model_name = 'DC_Motor_Model';

    % 检查 Simulink 是否可用
    try
        v = ver('simulink');
        if isempty(v)
            error('Simulink 未安装');
        end
    catch
        error('Simulink 不可用，无法构建模型');
    end

    % 如果模型已打开，先关闭
    if bdIsLoaded(model_name)
        close_system(model_name, 0);
    end

    % 删除已有的模型文件
    model_file = [model_name, '.slx'];
    if exist(model_file, 'file')
        delete(model_file);
    end

    % 创建新模型
    new_system(model_name);
    open_system(model_name);

    % ── 添加模块 ──────────────────────────────
    % Step 输入
    add_block('simulink/Sources/Step', [model_name, '/Step']);
    set_param([model_name, '/Step'], ...
        'Time', '0', ...
        'Before', '0', ...
        'After', '1');

    % Sum (误差计算)
    add_block('simulink/Math Operations/Sum', [model_name, '/Sum']);
    set_param([model_name, '/Sum'], ...
        'Inputs', '|+-', ...
        'IconShape', 'round');

    % PID Controller
    add_block('simulink/Continuous/PID Controller', [model_name, '/PID']);
    set_param([model_name, '/PID'], ...
        'Controller', 'PID', ...
        'P', 'Kp_sim', ...
        'I', 'Ki_sim', ...
        'D', 'Kd_sim', ...
        'FilterCoefficient', 'N_sim');

    % 如果 PID Controller 不可用，使用 Transfer Fcn 构建自定义 PID
    % （这里用 PID Controller 模块，MATLAB R2010b+ 都支持）

    % Transfer Fcn (被控对象)
    add_block('simulink/Continuous/Transfer Fcn', [model_name, '/DC_Motor']);
    set_param([model_name, '/DC_Motor'], ...
        'Numerator', '[10]', ...
        'Denominator', '[1 10 0]');

    % Scope
    add_block('simulink/Sinks/Scope', [model_name, '/Scope']);
    set_param([model_name, '/Scope'], ...
        'NumInputPorts', '2');

    % To Workspace (输出数据)
    add_block('simulink/Sinks/To Workspace', [model_name, '/yout']);
    set_param([model_name, '/yout'], ...
        'VariableName', 'yout', ...
        'SaveFormat', 'Structure With Time');

    add_block('simulink/Sinks/To Workspace', [model_name, '/uout']);
    set_param([model_name, '/uout'], ...
        'VariableName', 'uout', ...
        'SaveFormat', 'Structure With Time');

    % ── 连接模块 ──────────────────────────────
    % Step → Sum(+)
    add_line(model_name, 'Step/1', 'Sum/1', 'autorouting', 'smart');

    % Sum → PID → DC_Motor
    add_line(model_name, 'Sum/1', 'PID/1', 'autorouting', 'smart');
    add_line(model_name, 'PID/1', 'DC_Motor/1', 'autorouting', 'smart');

    % DC_Motor → yout, Scope
    add_line(model_name, 'DC_Motor/1', 'yout/1', 'autorouting', 'smart');
    add_line(model_name, 'DC_Motor/1', 'Scope/1', 'autorouting', 'smart');

    % PID → uout
    add_line(model_name, 'PID/1', 'uout/1', 'autorouting', 'smart');

    % 反馈: DC_Motor → Sum(-)
    add_line(model_name, 'DC_Motor/1', 'Sum/2', 'autorouting', 'smart');

    % ── 连接 Scope 第二通道（控制信号）───────
    add_line(model_name, 'PID/1', 'Scope/2', 'autorouting', 'smart');

    % ── 设置模型参数 ──────────────────────────
    set_param(model_name, ...
        'Solver', 'ode45', ...
        'StopTime', '2', ...
        'MaxStep', '0.01', ...
        'SaveOutput', 'on', ...
        'SaveFormat', 'StructureWithTime');

    % ── 初始化工作区变量 ─────────────────────
    assignin('base', 'Kp_sim', 1);
    assignin('base', 'Ki_sim', 10);
    assignin('base', 'Kd_sim', 0.1);
    assignin('base', 'N_sim', 100);  % PID滤波器系数

    % ── 保存并打开模型 ──────────────────────────
    save_system(model_name, model_file);
    % 保持模型打开，让用户看到框图
    fprintf('  Simulink模型已构建并打开: %s\n', model_file);
    fprintf('  模型包含: Step→Sum→PID→DC_Motor→Scope\n');
    fprintf('  双击 Scope 查看阶跃响应波形\n');
end
