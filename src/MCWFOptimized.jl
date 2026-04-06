module MCWFOptimized

using LinearAlgebra
using QuantumOptics
using Random
using Base.Threads

export MCWFParams, RunConfig, run_trajectory!, run_ensemble, build_quantumoptics_system

# =========================
# 参数结构体
# =========================

"""
物理参数与几何参数。
- 尽量与 MATLAB 原脚本变量名语义保持一致，方便逐项对照。
"""
Base.@kwdef struct MCWFParams
    v_z::Float64 = 290.0                    # 纵向速度
    hbar::Float64 = 1.0546e-34              # 约化普朗克常数
    tau::Float64 = 97.88e-9                 # 激发态寿命
    c_vac::Float64 = 2.99792458e8           # 真空光速
    nu_0::Float64 = 276736.600e9            # 跃迁频率
    m::Float64 = 4.0026 * 1.6606e-27        # 原子质量
    p::Int = 10                             # 动量离散窗口（-p...p）
    w0::Float64 = 1e-3                      # 光束腰
    cg::Float64 = 1 / sqrt(3)               # Clebsch-Gordan 系数
    dt_scale::Float64 = 0.1                 # dt = dt_scale * tau
    alpha::Float64 = 0.0                    # 附加速度偏置系数
    t2::Float64 = 0.67 / 290.0              # slit2 -> probe 飞行时间
    t3::Float64 = 1.53 / 290.0              # probe -> slit3 飞行时间
    d3::Float64 = 0.0003                    # slit3 宽度
    slit_center::Float64 = 0.0              # slit3 中心位置
end

"""
一次运行的控制参数。
"""
Base.@kwdef struct RunConfig
    i_sat::Float64 = 0.6                    # 饱和参数
    detuning::Float64 = 0.0                 # 失谐（单位 gamma）
    trajectories::Int = 10_000              # 轨迹数
    no_jump_step::Int = 90                  # no-jump 矩阵重建步长
    seed::Int = 1234                        # 随机种子
end

# =========================
# 线程本地缓存
# =========================

"""
每个线程独占一份缓存，避免多线程下重复分配和写冲突。
"""
struct TrajCache
    ψ::Vector{ComplexF64}                   # 当前态向量（2N）
    ψtmp::Vector{ComplexF64}                # 矩阵乘法临时向量
    ψg::Vector{ComplexF64}                  # 基态分量（N）
    nojump::Matrix{ComplexF64}              # no-jump 一步传播矩阵
    pgrid::Vector{Float64}                  # 每个离散动量通道的 ħk*i
end

# =========================
# 基础物理量
# =========================

@inline _γ(p::MCWFParams) = 1 / p.tau
@inline _k(p::MCWFParams) = 2π / (p.c_vac / p.nu_0)
@inline _vrec(p::MCWFParams) = p.hbar * _k(p) / p.m
@inline _num1(p::MCWFParams) = 2 * p.p + 1

"""
构建 QuantumOptics.jl 基矢对象（用于后续可观测量扩展）。
说明：主热循环为了性能采用手写内核，并不依赖通用主方程求解器。
"""
function build_quantumoptics_system(p::MCWFParams)
    n = _num1(p)
    mb = MomentumBasis(-p.p, p.p)
    ib = NLevelBasis(2)
    b = tensor(mb, ib)
    return (momentum_basis = mb, internal_basis = ib, basis = b, dimension = 2n)
end

"""
创建线程本地缓存。
"""
function _make_cache(p::MCWFParams)
    n = _num1(p)
    i1 = -p.p
    hk = p.hbar * _k(p)

    # 预计算离散动量通道对应的 ħk*i，避免在热点循环重复构造。
    pgrid = zeros(Float64, n)
    for j in 1:n
        mom_index = i1 + (j - 1)
        pgrid[j] = hk * mom_index
    end

    return TrajCache(
        zeros(ComplexF64, 2n),
        zeros(ComplexF64, 2n),
        zeros(ComplexF64, n),
        Matrix{ComplexF64}(I, 2n, 2n),
        pgrid,
    )
end

# =========================
# no-jump 矩阵更新
# =========================

"""
原地构造一步 no-jump 传播矩阵。

关键修复：
- 能量项使用 (ħk*i + m*v)^2/(2m) 的完整表达式（包含交叉项），
  避免遗漏 2*ħk*i*m*v 导致动力学偏差。
"""
@inline function _build_nojump!(M::Matrix{ComplexF64}, cache::TrajCache, p::MCWFParams,
                                cfg::RunConfig, v1::Float64, j02::Int, γ::Float64)
    n = _num1(p)
    fill!(M, 0)

    # 离散时间步
    dt = p.dt_scale * p.tau
    dt_hbar = dt / p.hbar

    # 失谐（单位转回频率）
    delta = cfg.detuning * γ

    # 原子横向速度偏置
    v_eff = v1 + p.alpha * p.v_z
    mv = p.m * v_eff

    # 高斯脉冲包络
    time = p.w0 / p.v_z
    sigma2 = (time / 2)^2
    t = j02 * dt
    omega = (1 / sqrt(2)) * p.cg * γ * sqrt(cfg.i_sat / 2)
    off = p.hbar * omega / 2 * exp(-((t - time / 2)^2) / sigma2)

    # 对角项 + 耦合项
    @inbounds for j in 1:n
        # 完整动能表达式（与 MATLAB 原始公式一致）
        kin = (cache.pgrid[j] + mv)^2 / (2p.m)

        # 激发态：含 detuning 与非厄米衰减项 -iħγ/2
        He = kin - p.hbar * delta - im * p.hbar * γ / 2

        # 基态：仅动能
        Hg = kin

        # Euler 一阶传播：I - i*dt/hbar*H
        M[j, j] = 1 - im * dt_hbar * He
        M[n + j, n + j] = 1 - im * dt_hbar * Hg
    end

    @inbounds for j in 1:(n - 1)
        # 激发态 <-> 基态 近邻动量耦合
        val = -im * dt_hbar * off
        M[j, n + j + 1] = val
        M[n + j + 1, j] = val
    end

    return M
end

# =========================
# 工具函数
# =========================

@inline function _normalize!(ψ::AbstractVector{ComplexF64})
    nrm = norm(ψ)
    if nrm > 0
        invn = inv(nrm)
        @inbounds @simd for i in eachindex(ψ)
            ψ[i] *= invn
        end
    end
end

"""
把激发块映射到基态块（等价于 MATLAB 中 matrix_jump_0 * vec）。
"""
@inline function _jump_to_ground!(ψ::Vector{ComplexF64}, ψg::Vector{ComplexF64}, n::Int)
    @inbounds @simd for i in 1:n
        ψg[i] = ψ[n + i]
    end
    _normalize!(ψg)

    fill!(ψ, 0)
    @inbounds @simd for i in 1:n
        ψ[n + i] = ψg[i]
    end
end

@inline _real_cuberoot(x::Float64) = x ≥ 0 ? cbrt(x) : -cbrt(-x)

"""
Recoil 分支角度公式（对应 MATLAB ep31）。
"""
@inline function _ep31(rand31::Float64)
    a = 3rand31 / 2
    b = sqrt(9rand31 * rand31 / 4 + 1)
    return acos(_real_cuberoot(a + b) + _real_cuberoot(a - b))
end

"""
Reset 分支角度公式（对应 MATLAB ep32）。

关键修复：
- 对可能为负的根号项，先转 Complex 再 sqrt，避免 DomainError/NaN。
"""
@inline function _ep32(rand32::Float64)
    w = ComplexF64((-1 + sqrt(3) * im) / 2)
    s1 = sqrt(ComplexF64(9rand32 * rand32 / 4 - 1))
    t1 = ComplexF64(-3rand32 / 2) + s1
    t2 = ComplexF64(-3rand32 / 2) - s1
    term = w^2 * t1^(1 / 3) + w * t2^(1 / 3)
    return acos(rand32 < 0 ? -abs(term) : abs(term))
end

# =========================
# 单轨迹
# =========================

"""
单条轨迹演化。
返回：
- detected: 该轨迹在 slit3 被探测到的概率贡献

三通道逻辑：
1) Quench: 终止
2) Recoil: 更新速度后继续
3) Reset: 更新后终止（与 MATLAB break 分支对齐）
"""
function run_trajectory!(cache::TrajCache, p::MCWFParams, cfg::RunConfig, rng::AbstractRNG)
    n = _num1(p)
    γ = _γ(p)
    vrec = _vrec(p)
    dt = p.dt_scale * p.tau
    steps = floor(Int, (p.w0 / p.v_z) / dt)

    # 初速度采样（正态分布）
    v1 = randn(rng) * (0.8 * p.v_z / 2190)
    x_initial = v1 * p.t2

    # 初态：只有中间动量通道的基态占据（对齐 MATLAB 索引）
    fill!(cache.ψ, 0)
    cache.ψ[(3n + 1) ÷ 2] = 1 + 0im
    fill!(cache.ψg, 0)

    # 用于控制 jump 后立即重建 no-jump 矩阵
    jump_num = 0
    jump_num1 = 0

    for j02 in 1:steps
        if (j02 - 1) % cfg.no_jump_step == 0 || jump_num1 < jump_num
            _build_nojump!(cache.nojump, cache, p, cfg, v1, j02, γ)
            jump_num1 = jump_num
        end

        # 跳跃概率：dp = γ * dt * ||e态||^2
        dp = γ * dt * sum(abs2, @view cache.ψ[1:n])

        ep1 = rand(rng)
        ep2 = rand(rng)
        ep4 = 2π * rand(rng)

        if ep1 > dp
            # 无跳跃：no-jump 演化
            mul!(cache.ψtmp, cache.nojump, cache.ψ)
            copyto!(cache.ψ, cache.ψtmp)
            _normalize!(cache.ψ)
            continue
        end

        # 有跳跃：三通道分支
        if ep2 > 2 / 3
            # Quench
            break
        elseif ep2 < 1 / 3
            # Recoil
            r31 = (8 / 3) * rand(rng) - 4 / 3
            θ = _ep31(r31)
            _jump_to_ground!(cache.ψ, cache.ψg, n)
            v1 += vrec * cos(ep4) * sin(θ)
            jump_num += 1
        else
            # Reset（对齐 MATLAB: jump 后 break）
            r32 = ((8 / 3) * rand(rng) - 4 / 3) / 2
            θ = _ep32(r32)
            _jump_to_ground!(cache.ψ, cache.ψg, n)
            v1 += vrec * cos(ep4) * sin(θ)
            break
        end
    end

    # 末态投影到 slit3 空间窗口统计
    detected = 0.0
    @inbounds for j in 1:n
        pfin = abs2(cache.ψg[j])
        y = x_initial + (v1 + (j - (n + 1) / 2) * vrec) * p.t3
        if y < p.slit_center + p.d3 / 2 && y > p.slit_center - p.d3 / 2
            detected += pfin
        end
    end

    return detected
end

# =========================
# 多轨迹并行
# =========================

"""
并行统计：
- 返回总探测计数 count
- 返回末态平均概率分布 probability
"""
function run_ensemble(p::MCWFParams, cfg::RunConfig)
    n = _num1(p)
    nt = nthreads()

    # 线程本地累加器
    counts = zeros(Float64, nt)
    probs = [zeros(Float64, n) for _ in 1:nt]

    # 线程本地随机数与缓存
    rngs = [Xoshiro(cfg.seed + i) for i in 1:nt]
    caches = [_make_cache(p) for _ in 1:nt]

    @threads for _ in 1:cfg.trajectories
        tid = threadid()
        cache = caches[tid]
        detected = run_trajectory!(cache, p, cfg, rngs[tid])
        counts[tid] += detected

        # 累加该轨迹的末态分布（来自线程本地 cache.ψg）
        @inbounds @simd for j in 1:n
            probs[tid][j] += abs2(cache.ψg[j])
        end
    end

    total_count = sum(counts)
    pmean = zeros(Float64, n)
    for tid in 1:nt
        @. pmean += probs[tid]
    end
    pmean ./= cfg.trajectories

    return (count = total_count, probability = pmean)
end

end
