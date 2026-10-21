using FermentModels, TOML, Dates, Test
include(joinpath(pwd(),"scripts/diagnose_cocoa_runtime.jl"))
include(joinpath(pwd(),"scripts/audit_cocoa_benchmark_bayes.jl"))
const R,A = CocoaRuntimeDiagnostics,CocoaBayesAudit
const C,B = R.C,R.B
const history = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits")
const reference = joinpath(history,"bayes_brwb3_budget_pilot2")
const output = joinpath(history,"bayes_brwb3_runtime_pilot1")
ispath(output) && error("Use a new experiment directory")
previous = TOML.parsefile(joinpath(reference,"experiment.toml"))
reference_hashes = previous["artifact_sha256"]
@test all(B.digest(joinpath(reference,f))==hash for (f,hash) in reference_hashes)
@test all(B.digest(joinpath(B.ROOT,f))==hash for (f,hash) in previous["environment_sha256"])
previous_settings = TOML.parsefile(joinpath(reference,"output_interval/settings.toml"))
snapshots = Dict("src/core/bayes.jl"=>"bayes.jl","scripts/fit_cocoa_benchmark_bayes.jl"=>"fit_cocoa_benchmark_bayes.jl")
@test all(B.digest(haskey(snapshots,f) ? joinpath(B.ROOT,".julia-depot/p25_telemetry_before",snapshots[f]) :
    joinpath(B.ROOT,f))==hash for (f,hash) in previous_settings["implementation_sha256"])
mkpath(output)
cp(@__FILE__,joinpath(output,"experiment.jl"))
mkpath(joinpath(output,"reference_implementation"))
for file in values(snapshots)
    cp(joinpath(B.ROOT,".julia-depot/p25_telemetry_before",file),joinpath(output,"reference_implementation",file))
end
experiment = Dict("recorded_before_sampling_utc"=>string(now(UTC)),
    "purpose"=>"Record leapfrog work, step sizes, acceptance and energy diagnostics for the same bounded brwb3 pilot; audit raw-chain correspondence and replay selected state costs.",
    "changes_from_reference"=>"Additional transition statistics and independent raw-chain export only; same equations, priors, data, solver, controls and per-chain seeds. Fresh adaptation; no pooled draws.",
    "samples_per_chain"=>25,"warmup_per_chain"=>50,"chains"=>4,"max_depth"=>8,
    "target_accept"=>.995,"trial_seed"=>141100,"abstol"=>1e-8,"reltol"=>1e-8,
    "solver_maxiters_scope"=>"output_interval","solver_maxiters"=>10_000,"solver_max_step_attempts"=>500,
    "max_sampling_seconds"=>600.,"chain_execution"=>"threads","julia_threads"=>Threads.nthreads(),
    "julia_version"=>string(VERSION),"script_sha256"=>B.digest(@__FILE__),
    "reference_directory"=>replace(relpath(reference,B.ROOT),'\\'=>'/'),
    "reference_experiment_sha256"=>B.digest(joinpath(reference,"experiment.toml")),
    "reference_artifact_sha256"=>reference_hashes,
    "environment_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in ("Project.toml","Manifest-v1.10.toml")),
    "limitations"=>"Short adaptation and depth 8 are diagnostic controls; no convergence is expected. Wall times include scheduling, compilation and GC. Retained-state cost replays cannot reconstruct a transition path or diagnose a divergence causally. Unchanged acceptance gates.")
C.F.write_toml(joinpath(output,"experiment.toml"),experiment)
fit = joinpath(output,"fit")
try
    C.run_benchmark_bayes(;output=fit,scope=:all,trials=[:brwb3],samples=25,warmup=50,
        chains=4,target_accept=.995,max_depth=8,seed=141085,solver=:tsit5,abstol=1e-8,reltol=1e-8,
        solver_maxiters=:source,solver_maxiters_scope=:output_interval,
        predictive_draws=50,chain_execution=:threads,max_sampling_seconds=600.)
catch err
    experiment["status"]="failed_or_stopped"
    experiment["error"]=sprint(showerror,err)
    C.F.write_toml(joinpath(output,"experiment.toml"),experiment)
    rethrow()
end
@testset "Telemetry preserves the reference pilot" begin
    for file in ("parameters.csv","predictive.csv","predictive_selection.csv","likelihood_observations.csv",
            "log_likelihood.csv",("chain_$(lpad(c,3,'0')).csv" for c in 1:4)...)
        @test B.digest(joinpath(fit,"brwb3",file))==B.digest(joinpath(reference,"output_interval/brwb3",file))
    end
end
audit = A.audit_bayes(fit;output=joinpath(output,"fit.audit.toml"),replay_per_chain=25)
@test audit["integrity_passed"] && audit["trials"]["brwb3"]["likelihood_replay_complete"]
runtime = R.diagnose_runtime(fit;output=joinpath(output,"runtime_diagnostics"),repeats=3)
@test runtime["runtime_events_checked"]==300 && runtime["retained_draws_checked"]==100
@test all(B.digest(joinpath(reference,f))==hash for (f,hash) in reference_hashes)
report = TOML.parsefile(joinpath(fit,"brwb3/report.toml"))
experiment["status"]="complete"
experiment["completed_utc"]=string(now(UTC))
experiment["reference_scientific_exports_identical"]=true
experiment["fit_accepted"]=report["diagnostics_passed"]
experiment["scientific_audit_checks"]=audit["checks"]
experiment["runtime_events_checked"]=runtime["runtime_events_checked"]
experiment["retained_draws_checked"]=runtime["retained_draws_checked"]
experiment["profiled_states"]=runtime["profiled_states"]
experiment["artifact_sha256"]=Dict(replace(relpath(joinpath(dir,f),output),'\\'=>'/')=>B.digest(joinpath(dir,f))
    for (dir,_,files) in walkdir(output) for f in files if joinpath(dir,f)!=joinpath(output,"experiment.toml"))
C.F.write_toml(joinpath(output,"experiment.toml"),experiment)
println("Runtime pilot and audits complete. Sampling accepted: ",experiment["fit_accepted"])
