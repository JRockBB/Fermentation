using FermentModels, TOML, Dates, Test
include("../scripts/audit_cocoa_benchmark_bayes.jl")
const A = CocoaBayesAudit
const C,B = A.C,A.B
const history = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits")
const output = joinpath(history,"bayes_brwb3_budget_pilot1")
const comparison_path = joinpath(history,"bayes_brwb3_solver_budget1/settings.toml")
ispath(output) && error("Use a new experiment directory")
@test TOML.parsefile(comparison_path)["all_trajectories_and_gradients_identical"]
mkpath(output)
cp(@__FILE__,joinpath(output,"experiment.jl"))
settings = Dict("recorded_before_sampling_utc"=>string(now(UTC)),
    "purpose"=>"Bounded four-chain brwb3 comparison of whole-solve and output-interval budgets, with completed-call warmup/runtime records.",
    "changes_between_cases"=>"Only solver_maxiters_scope; each case starts fresh from the same per-chain seeds.",
    "samples_per_chain"=>25,"warmup_per_chain"=>50,"chains"=>4,"max_depth"=>8,
    "target_accept"=>.995,"trial_seed"=>141100,"abstol"=>1e-8,"reltol"=>1e-8,
    "max_sampling_seconds_per_case"=>600.,"solver_maxiters"=>10_000,
    "solver_max_step_attempts"=>500,"stop_on_divergence"=>false,
    "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
    "limitations"=>"Deliberately short adaptation and depth 8 are diagnostic controls, not the full source sampling protocol or a convergence test. Compilation and scheduling differ by case order, so elapsed times are not a solver speed benchmark. Soft deadlines are checked only at completed calls. No earlier draws are pooled. Accepted MCMC diagnostics are not expected at this sample count.",
    "numerical_comparison_sha256"=>B.digest(comparison_path),"script_sha256"=>B.digest(@__FILE__),
    "environment_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in ("Project.toml","Manifest-v1.10.toml")))
C.F.write_toml(joinpath(output,"experiment.toml"),settings)
reports = Dict{String,Any}()
for scope in (:solve,:output_interval)
    directory = joinpath(output,string(scope))
    try
        C.run_benchmark_bayes(;output=directory,scope=:all,trials=[:brwb3],samples=25,
            warmup=50,chains=4,target_accept=.995,max_depth=8,seed=141085,
            solver=:tsit5,abstol=1e-8,reltol=1e-8,solver_maxiters=:source,
            solver_maxiters_scope=scope,predictive_draws=50,chain_execution=:threads,
            max_sampling_seconds=600.,stop_on_divergence=false)
    catch err
        report_path = joinpath(directory,"brwb3/report.toml")
        isfile(report_path) && TOML.parsefile(report_path)["status"]=="sampling_stopped" || rethrow()
    end
    report = TOML.parsefile(joinpath(directory,"brwb3/report.toml"))
    reports[string(scope)] = report
    if report["status"] in ("posterior_diagnostics_failed","posterior_diagnostics_passed")
        audit = A.audit_bayes(directory;output=joinpath(output,"$(scope).audit.toml"),replay_per_chain=25)
        @test audit["integrity_passed"]
        @test audit["trials"]["brwb3"]["likelihood_replay_complete"]
        @test audit["solver_budget_scope"] == string(scope)
    end
    C.F.write_toml(joinpath(output,"reports.toml"),reports)
end
checks = Tuple[]
@testset "Bounded budget pilot checkpoint comparison" begin
    for chain in 1:4
        filename = "chain_$(lpad(chain,3,'0')).csv"
        a,b = [readlines(joinpath(output,scope,"brwb3",filename))[2:end] for scope in ("solve","output_interval")]
        common = min(length(a),length(b))
        identical = a[1:common] == b[1:common]
        max_difference = common==0 ? NaN : maximum(maximum(abs.(
            parse.(Float64,split(a[i],',')[3:27])-parse.(Float64,split(b[i],',')[3:27]))) for i in 1:common)
        push!(checks,(chain,length(a),length(b),common,identical,max_difference))
        for scope in ("solve","output_interval")
            runtime = [split(line,',';keepempty=true) for line in readlines(joinpath(output,scope,"brwb3/runtime.csv"))[2:end]]
            rows = filter(r->parse(Int,r[1])==chain,runtime)
            @test parse.(Int,getindex.(rows,2)) == collect(1:length(rows))
            @test all(>=(0),diff(parse.(Float64,getindex.(rows,5))))
            @test count(r->parse(Int,r[4])>0,rows) == (scope=="solve" ? length(a) : length(b))
        end
    end
end
B.write_csv(joinpath(output,"comparison.csv"),("chain","solve_draws","interval_draws","common_draws",
    "common_checkpoint_rows_identical","max_parameter_difference"),checks)
settings["completed_utc"] = string(now(UTC))
settings["all_common_draws_identical"] = all(r[5] for r in checks)
settings["fit_accepted"] = all(get(r,"diagnostics_passed",false) for r in values(reports))
settings["artifact_sha256"] = Dict(replace(relpath(joinpath(dir,f),output),'\\'=>'/')=>B.digest(joinpath(dir,f))
    for (dir,_,files) in walkdir(output) for f in files if f!="experiment.toml")
C.F.write_toml(joinpath(output,"experiment.toml"),settings)
println("Budget pilot completed. Accepted: ",settings["fit_accepted"])
