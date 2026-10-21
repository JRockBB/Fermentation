using FermentModels, ForwardDiff, OrdinaryDiffEq, TOML, Dates, Test
isdefined(Main,:CocoaBayesDiagnostics) || include("../scripts/diagnose_cocoa_bayes.jl")
const D = CocoaBayesDiagnostics
const C,B = D.C,D.B
const history = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits")
const parent = joinpath(history,"bayes_brwb3_solver_timing1")
const output = joinpath(history,"bayes_brwb3_solver_budget1")
const old_code = joinpath(B.ROOT,".julia-depot/p25_budget_before")
ispath(output) && error("Use a new comparison directory")
source = TOML.parsefile(joinpath(parent,"settings.toml"))
checkpoint = joinpath(B.ROOT,source["parent"]["parent_directory"])
trial = only(filter(t->t.trial==:brwb3,B.benchmark_trials(;scope=:all)))
@test B.digest(joinpath(checkpoint,"settings.toml")) == source["parent"]["parent_settings_sha256"]
@test trial.source_sha256 == source["parent"]["source_sha256"]
for (file,hash) in source["parent"]["parent_files_sha256"]
    @test B.digest(joinpath(checkpoint,"brwb3",file)) == hash
end
for (file,hash) in source["implementation_sha256"]
    archived = joinpath(old_code,file)
    @test B.digest(isfile(archived) ? archived : joinpath(B.ROOT,file)) == hash
end
for (file,hash) in source["environment_sha256"]
    @test B.digest(joinpath(B.ROOT,file)) == hash
end
priors = C.benchmark_priors()
prior_record = TOML.parsefile(joinpath(checkpoint,"brwb3/priors.toml"))
@test all(all(prior_record[string(n)][string(k)]==v for (k,v) in pairs(priors[n].provenance)) for n in D.NAMES)
points = filter(f->startswith(f,"chain_") && endswith(f,".toml"),readdir(parent))
parent_hashes = Dict(f=>B.digest(joinpath(parent,f)) for f in ["settings.toml";points])
mkpath(output)
cp(@__FILE__,joinpath(output,"experiment.jl"))
cp(old_code,joinpath(output,"previous_implementation"))
settings = Dict("started_utc"=>string(now(UTC)),"purpose"=>"Compare continuous output-interval budgeting with the historical whole-solve mode at ten previously selected brwb3 retained points.",
    "parent_directory"=>replace(relpath(parent,B.ROOT),'\\'=>'/'),"parent_sha256"=>parent_hashes,
    "checkpoint_provenance"=>source["parent"],"script_sha256"=>B.digest(@__FILE__),
    "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
    "solver_maxiters"=>10_000,"solver_max_step_attempts"=>500,
    "limitations"=>"Selected retained states, including divergence flags and their predecessors, are not failed leapfrog states or a global accuracy/convergence certificate. The original completed-fit hash guards remain strict; old implementation snapshots validate the timing parent. No fit is relabelled.",
    "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in
        ("src/cocoa/Cocoa.jl","src/core/bayes.jl",C.F.LIKELIHOOD_IMPLEMENTATION...)),
    "environment_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in ("Project.toml","Manifest-v1.10.toml")))
C.F.write_toml(joinpath(output,"settings.toml"),settings)
rows = Tuple[]
@testset "Saved-point solver budget comparison" begin
    for file in points
        point = TOML.parsefile(joinpath(parent,file))
        @test point["parameter_order"] == String.(collect(D.NAMES))
        p = NamedTuple{D.NAMES}(Tuple(point["parameters"][string(n)] for n in D.NAMES))
        fields = split(readlines(joinpath(checkpoint,"brwb3","chain_$(lpad(point["chain"],3,'0')).csv"))[point["draw"]+1],',')
        @test parse.(Int,fields[1:2]) == [point["chain"],point["draw"]]
        @test parse.(Float64,fields[3:27]) == collect(values(p))
        lower,upper = NamedTuple{D.NAMES}(ntuple(_->1e-5,25)),NamedTuple{D.NAMES}(ntuple(_->5.,25))
        z = log.(collect(values(p)))
        for (label,solver,tol) in (("saved",:tsit5,1e-6),("tighter",:tsit5,1e-8),("alternate_tighter",:rodas5p,1e-8))
            initial = NamedTuple{D.NAMES}(ntuple(_->.5,25))
            models = [cocoa_fit_model(;initial,lower,upper,solver,abstol=tol,reltol=tol,
                maxiters=10_000,maxiters_scope=scope) for scope in (:solve,:output_interval)]
            predictions = [m.predict(p,trial.data) for m in models]
            gradients = [ForwardDiff.gradient(x->D.log_coordinate_objective(m,trial.data,priors,x),z) for m in models]
            nlls = [gaussian_nll(m,trial.data,p) for m in models]
            @test all(!isnothing,predictions) && all(v->all(isfinite,v),gradients)
            @test predictions[1] == predictions[2]
            @test gradients[1] == gradients[2]
            @test nlls[1] == nlls[2]
            @test gradients[1] ≈ point["gradient_log"][label] rtol=1e-12 atol=1e-12
            algorithm = solver==:tsit5 ? Tsit5() : Rodas5P()
            problem = cocoa_problem(trial.data.initial,(trial.data.t0,last(trial.data.times)),p)
            events = []
            sol = FermentModels.Cocoa._solve_output_intervals(problem,algorithm,trial.data.times;
                maxiters=10_000,abstol=tol,reltol=tol,on_interval=e->push!(events,e))
            @test vec(permutedims(Array(sol))) == predictions[2]
            cap = maximum(e.accepted for e in events)
            @test all(e->e.failure==:none && e.attempted>=e.accepted,events)
            @test !isnothing(FermentModels.Cocoa._solve_output_intervals(problem,algorithm,trial.data.times;
                maxiters=cap,abstol=tol,reltol=tol))
            @test isnothing(FermentModels.Cocoa._solve_output_intervals(problem,algorithm,trial.data.times;
                maxiters=cap-1,abstol=tol,reltol=tol))
            result = Dict("chain"=>point["chain"],"draw"=>point["draw"],"reasons"=>point["reasons"],
                "parameters"=>point["parameters"],"parameter_order"=>point["parameter_order"],
                "case"=>label,"solver"=>string(solver),"abstol"=>tol,"reltol"=>tol,
                "predictions_identical"=>true,"gradients_identical"=>true,"nll_identical"=>true,
                "nll"=>nlls[1],"gradient_log"=>gradients[1],"prediction_normalized"=>predictions[1],
                "minimum_accepted_step_cap"=>cap,
                "intervals"=>[Dict(string(k)=>v isa Symbol ? string(v) : v for (k,v) in pairs(e)) for e in events])
            C.F.write_toml(joinpath(output,replace(file,".toml"=>"_$label.toml")),result)
            push!(rows,(point["chain"],point["draw"],point["reasons"],label,solver,tol,cap,
                sum(e.accepted for e in events),sum(e.attempted for e in events),0.,0.,0.))
            B.write_csv(joinpath(output,"comparison.csv"),("chain","draw","reasons","case","solver",
                "tolerance","max_interval_accepted_steps","total_accepted_steps","total_attempted_steps",
                "delta_nll","max_prediction_difference","max_gradient_difference"),rows)
        end
        println("Compared ",file); flush(stdout)
    end
    @test length(rows)==30
    @test all(B.digest(joinpath(parent,f))==h for (f,h) in parent_hashes)
end
settings["completed_utc"] = string(now(UTC))
settings["cases"] = length(rows)
settings["gradient_evaluations"] = 2length(rows)
settings["all_trajectories_and_gradients_identical"] = true
C.F.write_toml(joinpath(output,"settings.toml"),settings)
