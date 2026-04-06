#!/usr/bin/env julia
using BenchmarkTools
using QuantumMC.MCWFOptimized
using Base.Threads

function run_once(ntraj)
    params = MCWFParams()
    cfg = RunConfig(trajectories=ntraj, seed=42, detuning=0.0, i_sat=0.6)
    return run_ensemble(params, cfg)
end

function main()
    ntraj = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 200_000
    matlab_parfor_seconds = length(ARGS) > 1 ? parse(Float64, ARGS[2]) : NaN

    println("Julia threads = ", nthreads())
    bench = @benchmark run_once($ntraj)
    t = minimum(bench).time / 1e9

    println("Trajectories: ", ntraj)
    println("Best Julia runtime (s): ", round(t, digits=4))
    println("Estimated throughput (traj/s): ", round(ntraj / t, digits=2))

    if !isnan(matlab_parfor_seconds)
        println("MATLAB parfor runtime (s): ", matlab_parfor_seconds)
        println("Speedup (MATLAB / Julia): ", round(matlab_parfor_seconds / t, digits=2), "x")
    else
        println("未提供 MATLAB parfor 耗时，无法自动计算 speedup。")
    end
end

main()
