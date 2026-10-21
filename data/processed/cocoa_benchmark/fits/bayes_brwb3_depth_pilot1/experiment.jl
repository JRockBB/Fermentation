using FermentModels, TOML, Dates, Test, Statistics
include(joinpath(pwd(),"scripts/diagnose_cocoa_runtime.jl"))
include(joinpath(pwd(),"scripts/audit_cocoa_benchmark_bayes.jl"))
const R,A = CocoaRuntimeDiagnostics,CocoaBayesAudit
const D,C,B = R.D,R.C,R.B

function comparison_row(label,directory,trial)
    saved = D.saved_draws(directory,trial)
    settings = saved.settings
    count_per_chain = settings["samples_per_chain"]
    rows = saved.rows
    @test length(rows)==settings["chains"]*count_per_chain==400
    draws = NamedTuple{D.NAMES}(Tuple(reshape([r.parameters[n] for r in rows],count_per_chain,4) for n in D.NAMES))
    summary = FermentModels.CoreModels._posterior_summary(draws)
    rhat = maximum(v.rhat for v in values(summary))
    bulk,tail = minimum(v.ess_bulk for v in values(summary)),minimum(v.ess_tail for v in values(summary))
    events = R.read_runtime(joinpath(directory,"brwb3/runtime.csv"),settings,rows)
    retained = filter(r->r.draw>0,events)
    leapfrogs = sum(r.n_steps for r in retained)
    gaps = sum(r.duration for r in retained)
    elapsed = maximum(r.elapsed for r in events)
    depth = settings["max_depth"]
    divergences,depth_hits = count(r->r.divergent,rows),count(r->r.tree_depth>=depth,rows)
    return (;case=label,max_depth=depth,warmup_per_chain=settings["warmup_per_chain"],
        draws_per_chain=count_per_chain,total_draws=length(rows),max_rhat=rhat,min_ess_bulk=bulk,
        min_ess_tail=tail,divergences,depth_hits,depth_hit_fraction=depth_hits/length(rows),
        retained_leapfrog_steps=leapfrogs,retained_completion_gap_seconds=gaps,
        median_retained_leapfrog_steps=median([r.n_steps for r in retained]),
        min_bulk_ess_per_1000_leapfrogs=1000bulk/leapfrogs,min_tail_ess_per_1000_leapfrogs=1000tail/leapfrogs,
        sampling_elapsed_seconds=elapsed,min_bulk_ess_per_sampling_second=bulk/elapsed,
        min_tail_ess_per_sampling_second=tail/elapsed)
end

function depth_pilot()
    history = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits")
    reference = joinpath(history,"bayes_brwb3_adaptation_pilot1")
    numerical = joinpath(history,"bayes_brwb3_runtime_solvers1")
    output = joinpath(history,"bayes_brwb3_depth_pilot1")
    ispath(output) && error("Use a new experiment directory")
    previous = TOML.parsefile(joinpath(reference,"experiment.toml"))
    scan = TOML.parsefile(joinpath(numerical,"experiment.toml"))
    @test previous["status"]=="complete" && !previous["fit_accepted"] && scan["status"]=="complete"
    @test (previous["warmup_per_chain"],previous["samples_per_chain"],previous["chains"],previous["max_depth"])==(250,100,4,8)
    @test scan["numerical_screen_passed"]["tsit5_1e8"] && !any(values(scan["pilot_candidate_screen_passed"]))
    @test all(B.digest(joinpath(reference,f))==h for (f,h) in previous["artifact_sha256"])
    @test all(B.digest(joinpath(numerical,f))==h for (f,h) in scan["artifact_sha256"])
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in previous["implementation_sha256"])
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in previous["environment_sha256"])
    @test Threads.nthreads()==2
    mkpath(output)
    cp(@__FILE__,joinpath(output,"experiment.jl"))
    watchdog = joinpath(B.ROOT,".julia-depot/p25_depth_watchdog.ps1")
    cp(watchdog,joinpath(output,"watchdog.ps1"))
    metadata = Dict{String,Any}("status"=>"sampling","recorded_before_sampling_utc"=>string(now(UTC)),
        "purpose"=>"Bounded depth-9 comparison against the completed depth-8 longer-adaptation brwb3 pilot.",
        "changes_from_reference"=>"Maximum tree depth 9 instead of 8. Sampling soft limit 2700 instead of 1800 seconds and external watchdog 3600 instead of 2700 seconds to accommodate additional work. Same 250 discarded calls, 100 retained draws, four seeds, two threads, target acceptance, priors, equations, data and solver. Fresh adaptation; no resumed or pooled draws.",
        "samples_per_chain"=>100,"warmup_per_chain"=>250,"chains"=>4,"max_depth"=>9,
        "target_accept"=>.995,"base_seed"=>141085,"trial_seed"=>141100,
        "solver"=>"tsit5","abstol"=>1e-8,"reltol"=>1e-8,
        "solver_maxiters"=>10_000,"solver_maxiters_scope"=>"output_interval","solver_max_step_attempts"=>500,
        "max_sampling_seconds"=>2700.,"stop_on_divergence"=>true,"process_watchdog_seconds"=>3600,
        "chain_execution"=>"threads","julia_threads"=>Threads.nthreads(),"julia_version"=>string(VERSION),
        "reference_directory"=>replace(relpath(reference,B.ROOT),'\\'=>'/'),
        "reference_experiment_sha256"=>B.digest(joinpath(reference,"experiment.toml")),
        "reference_artifact_sha256"=>previous["artifact_sha256"],
        "reference_receipt_sha256"=>B.digest(reference*".process.json"),
        "reference_audit_sha256"=>B.digest(reference*".audit.json"),
        "reference_native_exit_code_available"=>false,
        "numerical_comparison_directory"=>replace(relpath(numerical,B.ROOT),'\\'=>'/'),
        "numerical_comparison_experiment_sha256"=>B.digest(joinpath(numerical,"experiment.toml")),
        "numerical_comparison_artifact_sha256"=>scan["artifact_sha256"],
        "script_sha256"=>B.digest(@__FILE__),"watchdog_script_sha256"=>B.digest(watchdog),
        "environment_sha256"=>previous["environment_sha256"],
        "implementation_sha256"=>previous["implementation_sha256"],
        "diagnostic_plan"=>"On completion replay all likelihood draws and selected latent predictions, validate runtime against raw retained sampler statistics, profile slow states and controls, and compare equal-length four-chain diagnostics and ESS per retained leapfrog with depth 8. On a soft stop validate only partial checkpoints and completed runtime events; no joint diagnostics or raw-chain correspondence claim.",
        "followup_rule"=>"A complete audited run with zero retained divergences, a lower depth-hit fraction, no worse maximum R-hat, and higher minimum bulk and tail ESS per retained leapfrog can justify a separately declared longer run. A failed screen or incomplete run requires further diagnosis. This exploratory screen does not change the existing sampling-acceptance criteria.",
        "limitations"=>"One seed set at two depth limits is a diagnostic comparison, not an estimate of reproducible efficiency. ESS per retained leapfrog excludes warmup cost and does not account for varying gradient cost. Wall times include compilation, scheduling and GC; both elapsed sampling time and leapfrog work are reported. A larger depth cap can lower the cap-hit fraction mechanically, so that metric alone is insufficient. Warmup states and Hamiltonian paths are unavailable. The reference has an explicitly unknown native process exit code but independently audited complete scientific outputs.")
    C.F.write_toml(joinpath(output,"experiment.toml"),metadata)
    fit = joinpath(output,"fit")
    try
        try
            C.run_benchmark_bayes(;output=fit,scope=:all,trials=[:brwb3],samples=100,warmup=250,
                chains=4,target_accept=.995,max_depth=9,seed=141085,solver=:tsit5,abstol=1e-8,reltol=1e-8,
                solver_maxiters=:source,solver_maxiters_scope=:output_interval,solver_max_step_attempts=500,
                predictive_draws=50,chain_execution=:threads,max_sampling_seconds=2700.,stop_on_divergence=true)
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
                metadata["longer_run_screen_passed"]=false
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
            trial = only(filter(t->t.trial==:brwb3,B.benchmark_trials(;scope=:all)))
            control = comparison_row("depth8",joinpath(reference,"fit"),trial)
            candidate = comparison_row("depth9",fit,trial)
            B.write_csv(joinpath(output,"depth_comparison.csv"),string.(keys(control)),[Tuple(control),Tuple(candidate)])
            metadata["longer_run_screen_passed"] = candidate.divergences==0 &&
                candidate.depth_hit_fraction<control.depth_hit_fraction && candidate.max_rhat<=control.max_rhat &&
                candidate.min_bulk_ess_per_1000_leapfrogs>control.min_bulk_ess_per_1000_leapfrogs &&
                candidate.min_tail_ess_per_1000_leapfrogs>control.min_tail_ess_per_1000_leapfrogs
            metadata["status"]="complete"
            metadata["fit_accepted"]=report["diagnostics_passed"]
            metadata["joint_diagnostics_available"]=true
            metadata["scientific_audit_checks"]=audit["checks"]
            metadata["runtime_events_checked"]=runtime["runtime_events_checked"]
            metadata["retained_draws_checked"]=runtime["retained_draws_checked"]
            metadata["profiled_states"]=runtime["profiled_states"]
        end
        @test all(B.digest(joinpath(reference,f))==h for (f,h) in previous["artifact_sha256"])
        @test all(B.digest(joinpath(numerical,f))==h for (f,h) in scan["artifact_sha256"])
        @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in metadata["implementation_sha256"])
        @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in metadata["environment_sha256"])
        @test B.digest(joinpath(reference,"experiment.toml"))==metadata["reference_experiment_sha256"]
        @test B.digest(joinpath(numerical,"experiment.toml"))==metadata["numerical_comparison_experiment_sha256"]
        @test B.digest(reference*".process.json")==metadata["reference_receipt_sha256"] && B.digest(reference*".audit.json")==metadata["reference_audit_sha256"]
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
    println("Depth diagnostic: ",metadata["status"],"; accepted: ",metadata["fit_accepted"],"; longer-run screen: ",metadata["longer_run_screen_passed"])
end

@testset "Bounded tree-depth diagnostic integrity" begin
    depth_pilot()
end
