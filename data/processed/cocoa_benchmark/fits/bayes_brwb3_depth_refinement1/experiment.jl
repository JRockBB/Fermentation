using FermentModels, TOML, Dates, Test, Statistics
include(joinpath(pwd(),"scripts/diagnose_cocoa_runtime.jl"))
include(joinpath(pwd(),"scripts/audit_cocoa_benchmark_bayes.jl"))
const R,A = CocoaRuntimeDiagnostics,CocoaBayesAudit
const D,C,B = R.D,R.C,R.B
const SAMPLES,WARMUP,CHAINS,PREDICTIVE = 1000,1000,4,500

function summarize_run(label,directory,trial)
    saved = D.saved_draws(directory,trial)
    settings,rows = saved.settings,saved.rows
    count_per_chain = settings["samples_per_chain"]
    @test length(rows)==settings["chains"]*count_per_chain && settings["chains"]==CHAINS
    draws = NamedTuple{D.NAMES}(Tuple(reshape([r.parameters[n] for r in rows],count_per_chain,CHAINS) for n in D.NAMES))
    summary = FermentModels.CoreModels._posterior_summary(draws)
    rhat = maximum(v.rhat for v in values(summary))
    bulk,tail = minimum(v.ess_bulk for v in values(summary)),minimum(v.ess_tail for v in values(summary))
    events = R.read_runtime(joinpath(directory,"brwb3/runtime.csv"),settings,rows)
    retained = filter(r->r.draw>0,events)
    leapfrogs = sum(r.n_steps for r in retained)
    elapsed = maximum(r.elapsed for r in events)
    depth = settings["max_depth"]
    divergences,depth_hits = count(r->r.divergent,rows),count(r->r.tree_depth>=depth,rows)
    return (;case=label,max_depth=depth,warmup_per_chain=settings["warmup_per_chain"],
        draws_per_chain=count_per_chain,total_draws=length(rows),max_rhat=rhat,min_ess_bulk=bulk,
        min_ess_tail=tail,min_bulk_ess_per_draw=bulk/length(rows),min_tail_ess_per_draw=tail/length(rows),
        divergences,depth_hits,depth_hit_fraction=depth_hits/length(rows),retained_leapfrog_steps=leapfrogs,
        min_bulk_ess_per_1000_leapfrogs=1000bulk/leapfrogs,min_tail_ess_per_1000_leapfrogs=1000tail/leapfrogs,
        sampling_elapsed_seconds=elapsed,min_bulk_ess_per_sampling_second=bulk/elapsed,
        min_tail_ess_per_sampling_second=tail/elapsed)
end

function summarize_block(label,saved,indices)
    rows = [r for r in saved.rows if r.draw in indices]
    n = length(indices)
    @test length(rows)==CHAINS*n
    draws = NamedTuple{D.NAMES}(Tuple(reshape([r.parameters[k] for r in rows],n,CHAINS) for k in D.NAMES))
    summary = FermentModels.CoreModels._posterior_summary(draws)
    depth_hits = count(r->r.tree_depth>=saved.settings["max_depth"],rows)
    return (;case=label,first_draw=first(indices),last_draw=last(indices),draws_per_chain=n,total_draws=length(rows),
        max_rhat=maximum(v.rhat for v in values(summary)),min_ess_bulk=minimum(v.ess_bulk for v in values(summary)),
        min_ess_tail=minimum(v.ess_tail for v in values(summary)),divergences=count(r->r.divergent,rows),
        depth_hits,depth_hit_fraction=depth_hits/length(rows))
end

function depth_refinement()
    history = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits")
    reference = joinpath(history,"bayes_brwb3_depth_pilot1")
    numerical = joinpath(history,"bayes_brwb3_runtime_solvers1")
    output = joinpath(history,"bayes_brwb3_depth_refinement1")
    ispath(output) && error("Use a new experiment directory")
    previous = TOML.parsefile(joinpath(reference,"experiment.toml"))
    scan = TOML.parsefile(joinpath(numerical,"experiment.toml"))
    @test previous["status"]=="complete" && previous["longer_run_screen_passed"] && !previous["fit_accepted"]
    @test (previous["warmup_per_chain"],previous["samples_per_chain"],previous["chains"],previous["max_depth"])==(250,100,4,9)
    @test scan["status"]=="complete" && scan["numerical_screen_passed"]["tsit5_1e8"]
    @test all(B.digest(joinpath(reference,f))==h for (f,h) in previous["artifact_sha256"])
    @test all(B.digest(joinpath(numerical,f))==h for (f,h) in scan["artifact_sha256"])
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in previous["implementation_sha256"])
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in previous["environment_sha256"])
    @test Threads.nthreads()==2
    mkpath(output)
    cp(@__FILE__,joinpath(output,"experiment.jl"))
    watchdog = joinpath(B.ROOT,".julia-depot/p25_depth_refinement_watchdog.ps1")
    cp(watchdog,joinpath(output,"watchdog.ps1"))
    metadata = Dict{String,Any}("status"=>"sampling","recorded_before_sampling_utc"=>string(now(UTC)),
        "purpose"=>"Bounded longer four-chain brwb3 convergence attempt at the screened depth-9 settings.",
        "changes_from_reference"=>"1000 discarded calls instead of 250 and 1000 retained draws instead of 100 per chain. Predictive selection 500 instead of 50 draws. Sampling soft limit 10800 seconds and external watchdog 12600 seconds. Same depth 9, four seeds, two threads, target acceptance, priors, equations, observations, solver and ODE controls. Fresh prior initialization and adaptation; no resumed or pooled draws.",
        "samples_per_chain"=>SAMPLES,"warmup_per_chain"=>WARMUP,"chains"=>CHAINS,"predictive_draws"=>PREDICTIVE,
        "max_depth"=>9,"target_accept"=>.995,"base_seed"=>141085,"trial_seed"=>141100,
        "solver"=>"tsit5","abstol"=>1e-8,"reltol"=>1e-8,
        "solver_maxiters"=>10_000,"solver_maxiters_scope"=>"output_interval","solver_max_step_attempts"=>500,
        "max_sampling_seconds"=>10800.,"stop_on_divergence"=>true,"process_watchdog_seconds"=>12600,
        "chain_execution"=>"threads","julia_threads"=>Threads.nthreads(),"julia_version"=>string(VERSION),
        "reference_directory"=>replace(relpath(reference,B.ROOT),'\\'=>'/'),
        "reference_experiment_sha256"=>B.digest(joinpath(reference,"experiment.toml")),
        "reference_artifact_sha256"=>previous["artifact_sha256"],
        "reference_receipt_sha256"=>B.digest(reference*".process.json"),
        "reference_audit_sha256"=>B.digest(reference*".audit.json"),
        "reference_native_exit_code_available"=>true,
        "numerical_comparison_directory"=>replace(relpath(numerical,B.ROOT),'\\'=>'/'),
        "numerical_comparison_experiment_sha256"=>B.digest(joinpath(numerical,"experiment.toml")),
        "numerical_comparison_artifact_sha256"=>scan["artifact_sha256"],
        "script_sha256"=>B.digest(@__FILE__),"watchdog_script_sha256"=>B.digest(watchdog),
        "environment_sha256"=>previous["environment_sha256"],"implementation_sha256"=>previous["implementation_sha256"],
        "acceptance_rule"=>"Use the existing gate: four complete chains, finite R-hat below 1.01 for every parameter, bulk and tail ESS at least 400 for every parameter, zero retained divergences and zero depth-limit hits. Integrity audits are additionally required; they cannot change sampling acceptance.",
        "diagnostic_plan"=>"On completion replay all 4000 likelihood draws and all 500 selected latent predictions, validate 8000 runtime events against 4000 raw retained statistics, profile slow retained states and controls, and compare full diagnostics plus descriptive first/last 100-draw blocks. On a soft stop audit partial checkpoints and completed runtime records only; no joint diagnostics or raw-chain correspondence claim.",
        "followup_rule"=>"An accepted and audited completed fit may advance to benchmark predictive assessment. A stopped or diagnostically failed fit remains unaccepted and requires diagnosis; do not extend or relax criteria merely because it ran longer.",
        "limitations"=>"The pilot and refinement differ in adaptation and retained lengths. ESS per draw and matched blocks are descriptive; one seed set does not establish stationary sampling or reproducible efficiency. ESS per retained leapfrog excludes warmup and variable gradient cost. Wall times include compilation, scheduling and GC. Warmup states, adaptation state, and Hamiltonian paths are not saved. A completed callback can overrun the soft limit; the watchdog can interrupt an in-flight call. Historical artifacts remain immutable.")
    C.F.write_toml(joinpath(output,"experiment.toml"),metadata)
    fit = joinpath(output,"fit")
    try
        try
            C.run_benchmark_bayes(;output=fit,scope=:all,trials=[:brwb3],samples=SAMPLES,warmup=WARMUP,
                chains=CHAINS,target_accept=.995,max_depth=9,seed=141085,solver=:tsit5,abstol=1e-8,reltol=1e-8,
                solver_maxiters=:source,solver_maxiters_scope=:output_interval,solver_max_step_attempts=500,
                predictive_draws=PREDICTIVE,chain_execution=:threads,max_sampling_seconds=10800.,stop_on_divergence=true)
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
            audit = A.audit_bayes(fit;output=joinpath(output,"fit.audit.toml"),replay_per_chain=SAMPLES)
            @test audit["integrity_passed"] && audit["trials"]["brwb3"]["likelihood_replay_complete"]
            runtime = R.diagnose_runtime(fit;output=joinpath(output,"runtime_diagnostics"),repeats=3)
            @test runtime["runtime_events_checked"]==CHAINS*(WARMUP+SAMPLES) && runtime["retained_draws_checked"]==CHAINS*SAMPLES
            report = TOML.parsefile(joinpath(fit,"brwb3/report.toml"))
            @test report["rhat_threshold"]==1.01 && report["minimum_ess"]==400
            trial = only(filter(t->t.trial==:brwb3,B.benchmark_trials(;scope=:all)))
            control = summarize_run("pilot",joinpath(reference,"fit"),trial)
            candidate = summarize_run("refinement",fit,trial)
            B.write_csv(joinpath(output,"refinement_comparison.csv"),string.(keys(control)),[Tuple(control),Tuple(candidate)])
            reference_saved,refined_saved = D.saved_draws(joinpath(reference,"fit"),trial),D.saved_draws(fit,trial)
            blocks = [summarize_block("pilot_100",reference_saved,1:100),
                summarize_block("refinement_first_100",refined_saved,1:100),
                summarize_block("refinement_last_100",refined_saved,(SAMPLES-99):SAMPLES)]
            B.write_csv(joinpath(output,"matched_blocks.csv"),string.(keys(first(blocks))),Tuple.(blocks))
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
    println("Depth refinement: ",metadata["status"],"; accepted: ",metadata["fit_accepted"])
end

if abspath(PROGRAM_FILE)==@__FILE__
    @testset "Bounded depth-9 refinement integrity" begin
        depth_refinement()
    end
end
