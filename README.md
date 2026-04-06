# quantum-mc

将原 MATLAB 量子跳跃（MCWF）逻辑迁移到 Julia + QuantumOptics.jl 的实现位于：

- `src/MCWFOptimized.jl`: 高性能 MCWF 核心（含三通道跳转、线程并行）。
- `scripts/verify_against_matlab.jl`: 与 MATLAB `.mat` 参考结果做一致性校验（目标 `1e-10`）。
- `scripts/benchmark.jl`: Julia 单机性能测试与 MATLAB `parfor` 耗时对比。

## 快速开始

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project -e 'using QuantumMC; using QuantumMC.MCWFOptimized; println(run_ensemble(MCWFParams(), RunConfig(trajectories=10_000)).count)'
```

## MATLAB 对标验证

准备 `matlab_reference.mat`（至少包含 `probability_ref` 向量）：

```bash
julia --project scripts/verify_against_matlab.jl matlab_reference.mat
```

## 性能测试

```bash
# 仅测试 Julia
julia --project -t auto scripts/benchmark.jl 200000

# 同时输入 MATLAB parfor 耗时（秒）计算 speedup
julia --project -t auto scripts/benchmark.jl 200000 18.5
```


## 详细文档

- `docs/代码结构与原理说明.md`：代码结构与原理详解。
- `docs/服务器运行与参数扫描指南.md`：服务器运行、参数扫描与结果导出说明。

## Jupyter Notebook 入口

- `notebooks/QuantumMC_入口.ipynb`：用于在 Jupyter 中直接运行单点仿真、参数扫描与结果导出（CSV/MAT）。
