#!/usr/bin/env julia
using Random
using MAT
using QuantumMC.MCWFOptimized

"""
用法:
  julia --project scripts/verify_against_matlab.jl matlab_reference.mat

mat 文件需包含变量 `probability_ref` (Vector{Float64})。
"""
function main()
    if length(ARGS) < 1
        error("请提供 MATLAB 参考结果 .mat 文件路径")
    end
    ref_path = ARGS[1]
    ref = matread(ref_path)
    @assert haskey(ref, "probability_ref") "mat 文件缺少 probability_ref"
    p_ref = Vector{Float64}(ref["probability_ref"])

    params = MCWFParams()
    cfg = RunConfig(trajectories=50_000, seed=20260406, detuning=0.0, i_sat=0.6)
    result = run_ensemble(params, cfg)
    p_julia = result.probability

    @assert length(p_ref) == length(p_julia) "参考分布长度不一致"
    err = maximum(abs.(p_julia .- p_ref))
    println("max |p_julia - p_matlab| = ", err)
    if err ≤ 1e-10
        println("PASS: 与 MATLAB 在 1e-10 内一致")
    else
        println("FAIL: 偏差超出 1e-10")
    end
end

main()
