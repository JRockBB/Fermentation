using FermentModels, TOML, Dates, Test
include(joinpath(pwd(),"scripts/diagnose_cocoa_runtime.jl"))
include(joinpath(pwd(),"scripts/audit_cocoa_benchmark_bayes.jl"))
const R,A = CocoaRuntimeDiagnostics,CocoaBayesAudit
const D,C,B = R.D,R.C,R.B

function adaptation_pilot()
    history = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits")
    reference = joinpath(history,"bayes_brwb3_runtime_pilot1")
    numerical = joinpath(history,"bayes_brwb3_runtime_solvers1")
    output = joinpath(history,"bayes_brwb3_adaptation_pilot1")
    ispath(output) && error("Use a new experiment directory")
    previous = TOML.parsefile(joinpath(reference,"experiment.toml"))
    scan = TOML.parsefile(joinpath(numerical,"experiment.toml"))
    @test previous["status"]=="complete" && scan["status"]=="complete"
    @test scan["numerical_screen_passed"]["tsit5_1e8"]
    @test !any(values(scan["pilot_candidate_screen_passed"]))
    @test all(B.digest(joinpath(reference,f))==h for (f,h) in previous["artifact_sha256"])
    @test all(B.digest(joinpath(numerical,f))==h for (f,h) in scan["artifact_sha256"])
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in scan["implementation_sha256"])
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in scan["environment_sha256"])
    @test Threads.nthreads()==2
    mkpath(output)
    cp(@__FILE__,joinpath(output,"experiment.jl"))
    watchdog = joinpath(B.ROOT,".julia-depot/p25_adaptation_watchdog.ps1")
    cp(watchdog,joinpath(output,"watchdog.ps1"))
    metadata = Dict{String,Any}("status"=>"sampling","recorded_before_sampling_utc"=>string(now(UTC)),
        "purpose"=>"Bounded longer-adaptation diagnostic for brwb3 using the existing Tsit5 settings and telemetry.",
        "changes_from_reference"=>"250 discarded calls instead of 50; 100 retained draws per chain instead of 25; a 1800-second soft time budget and stop on retained divergence. Same four chain seeds, two threads, source priors, equations, data, target acceptance, depth limit and solver controls. Fresh adaptation; no pooled or resumed draws.",
        "samples_per_chain"=>100,"warmup_per_chain"=>250,"chains"=>4,"max_depth"=>8,
        "target_accept"=>.995,"base_seed"=>141085,"trial_seed"=>141100,
        "solver"=>"tsit5","abstol"=>1e-8,"reltol"=>1e-8,
        "solver_maxiters"=>10_000,"solver_maxiters_scope"=>"output_interval","solver_max_step_attempts"=>500,
        "max_sampling_seconds"=>1800.,"stop_on_divergence"=>true,"process_watchdog_seconds"=>2700,
        "chain_execution"=>"threads","julia_threads"=>Threads.nthreads(),"julia_version"=>string(VERSION),
        "reference_directory"=>replace(relpath(reference,B.ROOT),'\\'=>'/'),
        "reference_experiment_sha256"=>B.digest(joinpath(reference,"experiment.toml")),
        "reference_artifact_sha256"=>previous["artifact_sha256"],
        "numerical_comparison_directory"=>replace(relpath(numerical,B.ROOT),'\\'=>'/'),
        "numerical_comparison_experiment_sha256"=>B.digest(joinpath(numerical,"experiment.toml")),
        "numerical_comparison_artifact_sha256"=>scan["artifact_sha256"],
        "script_sha256"=>B.digest(@__FILE__),"watchdog_script_sha256"=>B.digest(watchdog),
        "environment_sha256"=>scan["environment_sha256"],
        "implementation_sha256"=>scan["implementation_sha256"],
        "diagnostic_plan"=>"For a completed fit, replay every likelihood draw and selected latent predictions; validate every runtime event against the raw retained sampler statistics; profile slow retained states and controls; compare full-chain diagnostics and matched first/last 25-draw blocks with the reference. For a stopped run, audit only saved partial checkpoints and completed runtime records; do not calculate joint convergence diagnostics or claim raw-chain correspondence.",
        "followup_rule"=>"A completed audited run with no retained divergences, lower depth-hit fraction and improved mixing can justify a separately declared longer run. Any divergence, incomplete run, or persistent poor mixing requires diagnosis before extending. Sampling acceptance always uses the existing R-hat/ESS/divergence/depth criteria.",
        "limitations"=>"This is a short diagnostic, not an expected converged fit. Depth 8 remains a diagnostic cap. Full runs have different retained lengths, so ESS is also reported per draw and 25-draw blocks are descriptive controls. Retained states do not reconstruct Hamiltonian paths. Warmup parameters and adaptation state are not saved. Soft limits act only at callback completion; an external watchdog may interrupt an in-flight call. Timings include compilation, scheduling and GC.")
    C.F.write_toml(joinpath(output,"experiment.toml"),metadata)
    fit = joinpath(output,"fit")
    try
        try
            C.run_benchmark_bayes(;output=fit,scope=:all,trials=[:brwb3],samples=100,warmup=250,
                chains=4,target_accept=.995,max_depth=8,seed=141085,solver=:tsit5,abstol=1e-8,reltol=1e-8,
                solver_maxiters=:source,solver_maxiters_scope=:output_interval,solver_max_step_attempts=500,
                predictive_draws=50,chain_execution=:threads,max_sampling_seconds=1800.,stop_on_divergence=true)
        catch err
            report_file = joinpath(fit,"brwb3/report.toml")
            report = isfile(report_file) ? TOML.parsefile(report_file) : Dict()
            if get(report,"status",nothing)=="sampling_stopped"
                metadata["status"]="sampling_stopped"
                metadata["stop_reason"]=report["stop_reason"]
                metadata["completed_sampler_calls"]=report["completed_sampler_calls"]
                metadata["partial_retained_draws"]=report["partial_retained_draws"]
                metadata["fit_accepted"]=false
                metadata["joint_diagnostics_available"]=false
            else
                rethrow()
            end
        end
        if metadata["status"]=="sampling"
            metadata["status"]="auditing"
            C.F.write_toml(joinpath(output,"experiment.toml"),metadata)
            audit = A.audit_bayes(fit;output=joinpath(output,"fit.audit.toml"),replay_per_chain=100)
            @test audit["integrity_passed"] && audit["trials"]["brwb3"]["likelihood_replay_complete"]
            runtime = R.diagnose_runtime(fit;output=joinpath(output,"runtime_diagnostics"),repeats=3)
            @test runtime["runtime_events_checked"]==1400 && runtime["retained_draws_checked"]==400
            report = TOML.parsefile(joinpath(fit,"brwb3/report.toml"))
            metadata["status"]="complete"
            metadata["fit_accepted"]=report["diagnostics_passed"]
            metadata["joint_diagnostics_available"]=true
            metadata["scientific_audit_checks"]=audit["checks"]
            metadata["runtime_events_checked"]=runtime["runtime_events_checked"]
            metadata["retained_draws_checked"]=runtime["retained_draws_checked"]
            metadata["profiled_states"]=runtime["profiled_states"]
            trial = only(filter(t->t.trial==:brwb3,B.benchmark_trials(;scope=:all)))
            comparisons = Tuple[]
            for (label,directory,indices) in (("reference_25",joinpath(reference,"fit"),1:25),
                    ("longer_full_100",fit,1:100),("longer_first_25",fit,1:25),("longer_last_25",fit,76:100))
                saved = D.saved_draws(directory,trial)
                count_per_chain = length(indices)
                rows = [r for r in saved.rows if r.draw in indices]
                @test length(rows)==4count_per_chain
                draws = NamedTuple{D.NAMES}(Tuple(reshape([r.parameters[n] for r in rows],count_per_chain,4) for n in D.NAMES))
                summary = FermentModels.CoreModels._posterior_summary(draws)
                rhat = maximum(v.rhat for v in values(summary))
                bulk,tail = minimum(v.ess_bulk for v in values(summary)),minimum(v.ess_tail for v in values(summary))
                divergences,depth_hits = count(r->r.divergent,rows),count(r->r.tree_depth>=8,rows)
                push!(comparisons,(label,count_per_chain,length(rows),rhat,bulk,tail,bulk/length(rows),tail/length(rows),
                    divergences,depth_hits,depth_hits/length(rows)))
            end
            B.write_csv(joinpath(output,"mixing_comparison.csv"),("case","draws_per_chain","total_draws","max_rhat",
                "min_ess_bulk","min_ess_tail","min_bulk_ess_per_draw","min_tail_ess_per_draw","divergences",
                "depth_hits","depth_hit_fraction"),comparisons)
        end
        @test all(B.digest(joinpath(reference,f))==h for (f,h) in previous["artifact_sha256"])
        @test all(B.digest(joinpath(numerical,f))==h for (f,h) in scan["artifact_sha256"])
        @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in metadata["implementation_sha256"])
        @test B.digest(joinpath(reference,"experiment.toml"))==metadata["reference_experiment_sha256"]
        @test B.digest(joinpath(numerical,"experiment.toml"))==metadata["numerical_comparison_experiment_sha256"]
    catch err
        metadata["status"]="failed"
        metadata["error"]=sprint(showerror,err)
        rethrow()
    finally
        metadata["finished_utc"]=string(now(UTC))
        metadata["artifact_sha256"]=Dict(replace(relpath(joinpath(dir,f),output),'\\'=>'/')=>B.digest(joinpath(dir,f))
            for (dir,_,files) in walkdir(output) for f in files if joinpath(dir,f)!=joinpath(output,"experiment.toml"))
        C.F.write_toml(joinpath(output,"experiment.toml"),metadata)
    end
    println("Adaptation diagnostic: ",metadata["status"],"; accepted: ",metadata["fit_accepted"])
end

@testset "Bounded longer-adaptation diagnostic integrity" begin
    adaptation_pilot()
end
