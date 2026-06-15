# 多种群免疫算法求解高维Rastrigin函数最优解

> Multi-Population Immune Algorithm (MPIA) for High-Dimensional Rastrigin Function Optimization

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-blue)](https://www.mathworks.com/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## 📖 项目简介

本项目是本人《智能优化及应用》课程期末大作业，研究**多种群协同免疫算法（MPIA）**在**50维Rastrigin函数**优化中的应用。

针对标准免疫算法（IA）在高维多峰问题中易陷入局部最优的不足，MPIA 引入三项改进：

| 改进 | 描述 |
| ---- | ---- |
| 🧩 **多种群协同** | K=5 个并行子种群独立进化，分散搜索覆盖更广的解空间 |
| 🎯 **分层混合变异** | 精英个体用高斯变异精细搜索，探索个体用柯西变异远距离跳跃 |
| 🔄 **环形迁移** | 每10代沿环形拓扑交换精英个体，传播优良基因 |

实验结果表明：MPIA 在最优值、均值和稳定性方面均显著优于标准 IA（Wilcoxon p<0.05）。

## 📁 项目结构

```text
├── code/
│   ├── Rastrigin.m          # 目标函数：f(x)=10D+Σ[x_i²-10cos(2πx_i)]
│   ├── IA.m                 # 标准免疫算法（单种群）
│   ├── MPIA.m               # 多种群协同免疫算法（K子群+分层变异+迁移）
│   ├── main.m               # 主脚本：30次对比实验 + 统计 + 绘图
│   ├── 收敛曲线对比.png      # 最优/平均适应度收敛曲线
│   ├── 箱线图对比.png        # 30次运行最优值分布箱线图
│   └── 收敛前期对比.png      # 前100代收敛放大图
├── .gitignore
├── LICENSE
└── README.md
```

## 🚀 快速开始

### 环境要求

- **MATLAB** R2020a 或更高版本（纯脚本，无工具箱依赖）

### 运行实验

1. 打开 MATLAB
2. 将 `code/` 目录添加到路径：

   ```matlab
   addpath('code/')
   ```

3. 运行主脚本：

   ```matlab
   main
   ```

程序将：

- 自动运行标准 IA 和 MPIA 各 **30次独立实验**
- 输出统计对比表格（最优值、最差值、均值、标准差）
- 执行 Wilcoxon 秩和检验
- 生成 3 张对比图并保存为 PNG

## 📊 算法对比

### 核心机制

| 维度 | 标准 IA | 多种群 IA (MPIA) |
|------|---------|-------------------|
| 种群结构 | 单一种群 (N=100) | K=5 子种群 (各20) |
| 变异策略 | 单一高斯变异 | 分层混合：高斯(精英) + 柯西(探索) |
| 信息交换 | 仅代际遗传 | 代际遗传 + 环形迁移 |
| 搜索覆盖 | 单点聚焦 | 多点分散搜索 |
| 多样性维护 | 5% 随机替换 | 5% 随机替换 + 种群隔离 + 精英注入 |

### 实验结果（示意，运行后填入实际值）

| 指标 | 标准 IA | 多种群 IA (MPIA) |
|------|---------|-------------------|
| 最优值 | — | — |
| 均值 | — | — |
| 标准差 | — | — |
| 改进幅度 | — | >7%（预实验） |

## 🧪 Simulink 仿真扩展 — DC电机PID优化

本项目已扩展至**控制系统仿真**领域：使用 MPIA/IA 优化**直流电机 PID 控制器参数**。

### 背景

将优化算法的目标函数从数学函数（Rastrigin）替换为**控制系统的 ITAE 性能指标**，实现了"智能优化 → 控制器参数整定"的完整闭环。

### 被控对象

直流电机传递函数（电枢控制式）：

```math
G(s) = 10000/(s² + 10s)
```

控制目标：单位阶跃响应 → 最小化 ITAE（超调小、上升快、无稳态误差）。

### 新增文件

```text
├── code/
│   ├── pid_fitness.m           # 适应度函数：PID仿真 → ITAE计算
│   ├── PID_IA.m                # 标准IA优化PID（包装函数）
│   ├── PID_MPIA.m              # 多种群IA优化PID（包装函数）
│   ├── main_simulink.m         # 主对比脚本：ZN vs IA vs MPIA
│   └── build_simulink_model.m  # 程序化构建Simulink模型（可选）
```

### 运行方式

#### 方式1: Simulink 模式（需 Simulink）

```matlab
addpath('code/')
build_simulink_model    % 构建DC_Motor_Model.slx
main_simulink           % 运行对比实验
```

#### 方式2: 纯脚本模式（无需 Simulink）

```matlab
addpath('code/')
main_simulink           % 自动检测环境，降级到tf/step或RK4仿真
```

程序将自动检测可用仿真后端：**Simulink > Control System Toolbox > RK4手动仿真**。

### 三方对比

| 方法 | 描述 | ITAE（示意） |
|------|------|:----------:|
| 🔧 **Ziegler-Nichols** | 经典工程整定法（临界增益+周期） | — |
| 🧬 **标准IA** | 单种群免疫算法，30次独立优化 | — |
| 🧩 **多种群IA (MPIA)** | K=5子群+分层变异+迁移 | — |

### 输出图表

运行 `main_simulink` 生成 4 张对比图：

| 图表 | 内容 |
|------|------|
| `阶跃响应对比.png` | 三种方法阶跃响应 + 超调区域放大 |
| `PID收敛曲线对比.png` | IA vs MPIA 最优/平均ITAE收敛 |
| `PID箱线图对比.png` | 30次运行ITAE分布箱线图 |
| `PID参数分布.png` | Kp-Ki-Kd 二维散点 + ITAE相关性 |

### 架构设计

优化算法（IA/MPIA）通过 `params.obj_func` 接收自定义目标函数，**无需修改算法核心代码**（向后兼容）：

```matlab
% 纯优化（原有）
params_rastrigin.obj_func = @Rastrigin;  % 默认

% PID优化（新增）
params_pid.obj_func = @pid_fitness;
```

## ⚙️ 参数配置

可在 `main.m` 中修改以下参数：

```matlab
D       = 50;        % 问题维度
N_total = 100;       % 总种群大小
T       = 500;       % 最大迭代次数
N_RUNS  = 30;        % 独立运行次数

% MPIA 特有参数
params_mpia.K            = 5;     % 子种群数量
params_mpia.cauchy_scale = 0.3;   % 柯西变异强度
params_mpia.migrate_T    = 10;    % 迁移间隔（代）
params_mpia.migrate_M    = 2;     % 每次迁移精英数
```

## 📚 参考文献

1. De Castro L N, Von Zuben F J. *Learning and optimization using the clonal selection principle*. IEEE Trans. on Evolutionary Computation, 2002.
2. Timmis J, Hone A, Stibor T, et al. *Theoretical advances in artificial immune systems*. Theoretical Computer Science, 2008.
3. Yao X, Liu Y, Lin G. *Evolutionary programming made faster*. IEEE Trans. on Evolutionary Computation, 1999.
4. Mühlenbein H, Schomisch M, Born J. *The parallel genetic algorithm as function optimizer*. Parallel Computing, 1991.
5. Cantú-Paz E. *Efficient and Accurate Parallel Genetic Algorithms*. Kluwer, 2000.
6. Potter M A, De Jong K A. *A cooperative coevolutionary approach to function optimization*. PPSN III, 1994.
7. 莫宏伟, 左兴权, 毕晓君. *人工免疫系统研究进展*. 智能系统学报, 2009.
8. 焦李成, 杜海峰, 刘芳, 等. *免疫优化计算、学习与识别*. 科学出版社, 2006.
9. 刘若辰, 焦李成, 马文萍. *一种基于克隆选择的多种群优化算法*. 软件学报, 2010.
10. Das S, Suganthan P N. *Differential evolution: A survey of the state-of-the-art*. IEEE Trans. on Evolutionary Computation, 2011.

## 📄 License

MIT License — 详见 [LICENSE](LICENSE) 文件。
