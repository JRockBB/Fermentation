using FermentModels, ForwardDiff, Statistics, TOML, Dates
include("../scripts/diagnose_cocoa_bayes.jl")
const D = CocoaBayesDiagnostics
const C,B = D.C,D.B
const parent = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits/bayes_brwb3_ghhp4_v2")
const output = joinpath(dirname(parent),"bayes_brwb3_solver_timing1")
ispath(output) && error("Use a new timing directory")
trial = only(filter(t->t.trial==:brwb3,B.benchmark_trials(;scope=:all)))
saved = D.saved_draws(parent,trial)
selected = D.select_draws(saved.rows;per_chain=1,tail_parameters=[:yc11])
initial = NamedTuple{D.NAMES}(ntuple(_->.5,25))
lower,upper = NamedTuple{D.NAMES}(ntuple(_->1e-5,25)),NamedTuple{D.NAMES}(ntuple(_->5.,25))
cases = [("saved",:tsit5,1e-6),("tighter",:tsit5,1e-8),("alternate_tighter",:rodas5p,1e-8)]
models = [cocoa_fit_model(;initial,lower,upper,solver,abstol=tol,reltol=tol,maxiters=10_000) for (_,solver,tol) in cases]
gradient(model,data,priors,z) = ForwardDiff.gradient(x->D.log_coordinate_objective(model,data,priors,x),z)
mkpath(output)
cp(@__FILE__,joinpath(output,"timing_script.jl"))
metadata = Dict("started_utc"=>string(now(UTC)),"julia_version"=>string(VERSION),
    "julia_threads"=>Threads.nthreads(),"cpu_name"=>Sys.CPU_NAME,"cpu_threads"=>Sys.CPU_THREADS,
    "parent"=>saved.provenance,"repeats"=>5,"per_chain_grid"=>1,"tail_parameters"=>["yc11"],
    "solver_maxiters"=>10_000,"script_sha256"=>B.digest(@__FILE__),
    "limitations"=>"Warm-compiled model-gradient timings under shared machine load, not end-to-end NUTS timings. Cases alternate within each repetition. Selected retained points do not cover all trajectories.",
    "environment_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in ("Project.toml","Manifest-v1.10.toml")),
    "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in
        ("scripts/diagnose_cocoa_bayes.jl","scripts/fit_cocoa_benchmark_bayes.jl","src/core/bayes.jl",C.F.LIKELIHOOD_IMPLEMENTATION...)))
C.F.write_toml(joinpath(output,"settings.toml"),metadata)
rows = Tuple[]
for point in selected
    z = log.([point.parameters[n] for n in D.NAMES])
    values = [gradient(model,trial.data,saved.priors,z) for model in models]
    all(v->all(isfinite,v),values) || error("Nonfinite warmup gradient")
    C.F.write_toml(joinpath(output,"chain_$(point.chain)_draw_$(point.draw).toml"),
        Dict("chain"=>point.chain,"draw"=>point.draw,"reasons"=>point.reasons,
            "parameters"=>Dict(string(n)=>point.parameters[n] for n in D.NAMES),
            "parameter_order"=>String.(collect(D.NAMES)),
            "gradient_log"=>Dict(label=>value for ((label,_,_),value) in zip(cases,values))))
    for repeat in 1:5, (index,(label,solver,tol)) in enumerate(cases)
        measured = @timed gradient(models[index],trial.data,saved.priors,z)
        isapprox(measured.value,values[index];rtol=1e-12,atol=1e-12) || error("Gradient changed")
        push!(rows,(point.chain,point.draw,label,solver,tol,repeat,measured.time,measured.bytes))
    end
    B.write_csv(joinpath(output,"timings.csv"),
        ("chain","draw","case","solver","tolerance","repeat","seconds","allocated_bytes"),rows)
    println("Timed chain ",point.chain," draw ",point.draw); flush(stdout)
end
summaries = [(label,solver,tol,minimum(r[7] for r in rows if r[3]==label),
    median([r[7] for r in rows if r[3]==label]),maximum(r[8] for r in rows if r[3]==label))
    for (label,solver,tol) in cases]
B.write_csv(joinpath(output,"summary.csv"),
    ("case","solver","tolerance","minimum_seconds","median_seconds","maximum_allocated_bytes"),summaries)
metadata["completed_utc"] = string(now(UTC))
metadata["points"] = length(selected)
metadata["timed_gradients"] = length(rows)
C.F.write_toml(joinpath(output,"settings.toml"),metadata)
foreach(println,summaries)
