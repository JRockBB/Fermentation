using FermentModels, ForwardDiff, Statistics, TOML, Dates, Test
include(joinpath(pwd(),"scripts/diagnose_cocoa_runtime.jl"))
const R = CocoaRuntimeDiagnostics
const D,C,B = R.D,R.C,R.B

scaled_difference(a,b) = maximum(abs.(a-b)./max.(1.,abs.(b)))

function run_comparison()
    history = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits")
    parent = joinpath(history,"bayes_brwb3_runtime_pilot1")
    output = joinpath(history,"bayes_brwb3_runtime_solvers1")
    ispath(output) && error("Use a new experiment directory")
    previous = TOML.parsefile(joinpath(parent,"experiment.toml"))
    @test previous["status"]=="complete" && !previous["fit_accepted"]
    @test all(B.digest(joinpath(parent,f))==h for (f,h) in previous["artifact_sha256"])
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in previous["environment_sha256"])
    previous_profiles = TOML.parsefile(joinpath(parent,"runtime_diagnostics/settings.toml"))
    @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in previous_profiles["implementation_sha256"])
    trial = only(filter(t->t.trial==:brwb3,B.benchmark_trials(;scope=:all)))
    saved = D.saved_draws(joinpath(parent,"fit"),trial)
    events = R.read_runtime(joinpath(parent,"fit/brwb3/runtime.csv"),saved.settings,saved.rows)
    selected = R.select_runtime_draws(events,saved.rows)
    @test length(selected)==12
    cases = [("tsit5_1e8",:tsit5,1e-8),("tsit5_1e10",:tsit5,1e-10),
        ("rodas5p_1e8",:rodas5p,1e-8),("rodas5p_1e10",:rodas5p,1e-10)]
    # These are numerical diagnostic criteria, not source-model parameters or fit acceptance rules.
    gates = Dict("max_normalized_prediction_difference"=>1e-6,"max_absolute_nll_difference"=>1e-5,
        "max_scaled_gradient_difference"=>1e-4,"max_fd_scaled_error"=>1e-4,
        "max_expensive_state_time_ratio"=>.8,"max_median_paired_time_ratio"=>1.25,
        "minimum_faster_rounds_at_each_expensive_state"=>4)
    metadata = Dict("status"=>"running","recorded_before_evaluation_utc"=>string(now(UTC)),
        "purpose"=>"Compare Tsit5 and Rodas5P accuracy and warmed work at the 12 previously selected runtime states before considering another sampler pilot.",
        "parent_directory"=>replace(relpath(parent,B.ROOT),'\\'=>'/'),
        "parent_experiment_sha256"=>B.digest(joinpath(parent,"experiment.toml")),
        "parent_artifact_sha256"=>previous["artifact_sha256"],"script_sha256"=>B.digest(@__FILE__),
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"cpu_name"=>Sys.CPU_NAME,
        "cases"=>[Dict("name"=>label,"solver"=>string(solver),"abstol"=>tol,"reltol"=>tol) for (label,solver,tol) in cases],
        "solver_maxiters"=>10_000,"solver_maxiters_scope"=>"output_interval","solver_max_step_attempts"=>500,
        "profile_repeats"=>3,"gradient_timing_rounds"=>5,"fd_step"=>1e-4,
        "comparison_gates"=>gates,"expensive_states"=>[Dict("chain"=>4,"draw"=>d) for d in (24,25)],
        "gate_definition"=>"Every candidate must complete all 12 cases, agree with both 1e-10 solver references within every numerical gate, and pass both centered finite-difference step checks. The fine references must also agree and pass their finite-difference checks. A Rodas5P candidate may advance to a bounded pilot only if each expensive state is at least 20% faster by median gradient time, is faster in at least four of five matched rounds, and the median of all 12 paired time ratios is at most 1.25 versus Tsit5 1e-8. These are exploratory screening criteria, not MCMC acceptance.",
        "measurement_order"=>"Rotate case order by point for profiles/FD; warm all gradient callables, then rotate the first case by point and timing round. Preserve each time, allocation and GC measurement.",
        "limitations"=>"A small selected subset of an unaccepted fit is not a global accuracy certificate. Fine numerical solutions are comparison references, not exact truth. Retained states do not reconstruct NUTS paths. Timings depend on shared machine load; no end-to-end sampling speed claim. No priors, equations, data or fit acceptance gates change.",
        "diagnostics_accept_fit"=>false,
        "environment_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in ("Project.toml","Manifest-v1.10.toml")),
        "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in
            ("scripts/diagnose_cocoa_runtime.jl","scripts/diagnose_cocoa_bayes.jl",
             "scripts/fit_cocoa_benchmark_bayes.jl","src/core/bayes.jl",C.F.LIKELIHOOD_IMPLEMENTATION...)))
    @test saved.settings["solver_maxiters_scope"]=="output_interval"
    @test saved.settings["trial_solver_maxiters"]["brwb3"]==10_000
    @test saved.settings["solver_max_step_attempts"]==500
    mkpath(output)
    cp(@__FILE__,joinpath(output,"experiment.jl"))
    C.F.write_toml(joinpath(output,"experiment.toml"),metadata)
    B.write_csv(joinpath(output,"selection.csv"),("chain","draw","reasons"),
        [(r.chain,r.draw,r.reasons) for r in selected])
    @test read(joinpath(output,"selection.csv"))==read(joinpath(parent,"runtime_diagnostics/brwb3/selection.csv"))
    initial = NamedTuple{D.NAMES}(ntuple(_->.5,length(D.NAMES)))
    lower = NamedTuple{D.NAMES}(ntuple(_->1e-5,length(D.NAMES)))
    upper = NamedTuple{D.NAMES}(ntuple(_->5.,length(D.NAMES)))
    models = [cocoa_fit_model(;initial,lower,upper,solver,abstol=tol,reltol=tol,
        maxiters=10_000,maxiters_scope=:output_interval,max_step_attempts=500) for (_,solver,tol) in cases]
    reports = Dict{Tuple{Int,Int,String},Dict{String,Any}}()
    measurements = Tuple[]
    try
        @testset "Solver profiles match numerical diagnostics and saved states" begin
            for (point_index,row) in enumerate(selected)
                order = [mod1(point_index+k,length(cases)) for k in 0:length(cases)-1]
                for index in order
                    label,solver,tol = cases[index]
                    altered = (;saved...,settings=merge(saved.settings,Dict("solver"=>string(solver),"abstol"=>tol,"reltol"=>tol)))
                    profile = R.profile_point(trial,altered,row.parameters;repeats=3)
                    check = D.point_check(models[index],trial.data,saved.priors,row.parameters;fd_step=1e-4)
                    complete = profile["status"]==check["status"]=="complete"
                    if complete
                        @test profile["prediction_normalized"]==check["prediction_normalized"]
                        @test profile["nll"]==check["nll"]
                        @test profile["gradient_log"]==check["gradient_log"]
                        if index==1
                            old = TOML.parsefile(joinpath(parent,"runtime_diagnostics/brwb3/chain_$(row.chain)_draw_$(row.draw).toml"))
                            @test all(profile[k]==old[k] for k in ("prediction_normalized","nll","gradient_log",
                                "accepted_ode_steps","rejected_ode_steps","rhs_evaluations","parameters"))
                        end
                    end
                    report = Dict{String,Any}("chain"=>row.chain,"draw"=>row.draw,"reasons"=>row.reasons,
                        "case"=>label,"solver"=>string(solver),"abstol"=>tol,"reltol"=>tol,
                        "status"=>complete ? "complete" : "failed","profile"=>profile,"numerics"=>check)
                    reports[(row.chain,row.draw,label)] = report
                    C.F.write_toml(joinpath(output,"chain_$(row.chain)_draw_$(row.draw)_$label.toml"),report)
                end
                z = log.(collect(values(row.parameters)))
                functions = [let model=model
                    ()->ForwardDiff.gradient(x->D.log_coordinate_objective(model,trial.data,saved.priors,x),z)
                end for model in models]
                for index in order
                    reports[(row.chain,row.draw,cases[index][1])]["status"]=="complete" || continue
                    @test functions[index]() == reports[(row.chain,row.draw,cases[index][1])]["profile"]["gradient_log"]
                end
                for round in 1:5, k in 0:length(cases)-1
                    index = mod1(point_index+round+k-1,length(cases))
                    label = cases[index][1]
                    reports[(row.chain,row.draw,label)]["status"]=="complete" || continue
                    measured = @timed functions[index]()
                    @test measured.value==reports[(row.chain,row.draw,label)]["profile"]["gradient_log"]
                    push!(measurements,(row.chain,row.draw,label,round,k+1,measured.time,measured.bytes,measured.gctime))
                end
                B.write_csv(joinpath(output,"gradient_timings.csv"),
                    ("chain","draw","case","round","position","seconds","allocated_bytes","gc_seconds"),measurements)
                println("Compared chain ",row.chain," draw ",row.draw); flush(stdout)
            end
        end
        comparisons = Tuple[]
        accuracy = Dict(label=>true for (label,_,_) in cases)
        ratios = Dict(label=>Float64[] for (label,_,_) in cases)
        expensive_fast = Dict(label=>true for (label,_,_) in cases)
        for row in selected
            baseline_times = [m[6] for m in measurements if m[1]==row.chain && m[2]==row.draw && m[3]=="tsit5_1e8"]
            refs = [reports[(row.chain,row.draw,label)] for label in ("tsit5_1e10","rodas5p_1e10")]
            for (label,solver,tol) in cases
                report = reports[(row.chain,row.draw,label)]
                complete = report["status"]=="complete" && all(r["status"]=="complete" for r in refs)
                if !complete
                    accuracy[label]=false
                    expensive_fast[label]=false
                    push!(comparisons,(row.chain,row.draw,label,solver,tol,"failed",ntuple(_->missing,10)...,false))
                    continue
                end
                p = report["profile"]
                pred = maximum(maximum(abs.(p["prediction_normalized"]-r["profile"]["prediction_normalized"])) for r in refs)
                nll = maximum(abs(p["nll"]-r["profile"]["nll"]) for r in refs)
                grad = maximum(scaled_difference(p["gradient_log"],r["profile"]["gradient_log"]) for r in refs)
                fd = maximum(f["max_scaled_error"] for f in report["numerics"]["finite_differences"])
                ok = pred<=gates["max_normalized_prediction_difference"] && nll<=gates["max_absolute_nll_difference"] &&
                    grad<=gates["max_scaled_gradient_difference"] && fd<=gates["max_fd_scaled_error"]
                accuracy[label] &= ok
                times = [m[6] for m in measurements if m[1]==row.chain && m[2]==row.draw && m[3]==label]
                @test length(times)==5
                ratio = length(baseline_times)==5 ? median(times)/median(baseline_times) : Inf
                push!(ratios[label],ratio)
                faster = length(baseline_times)==5 ? count(times .< baseline_times) : 0
                if row.chain==4 && row.draw in (24,25)
                    expensive_fast[label] &= ratio<=gates["max_expensive_state_time_ratio"] &&
                        faster>=gates["minimum_faster_rounds_at_each_expensive_state"]
                end
                push!(comparisons,(row.chain,row.draw,label,solver,tol,"complete",
                    p["accepted_ode_steps"],p["rejected_ode_steps"],p["rhs_evaluations"],median(times),ratio,
                    pred,nll,grad,fd,faster,ok))
            end
        end
        B.write_csv(joinpath(output,"comparison.csv"),("chain","draw","case","solver","tolerance","status",
            "accepted_ode_steps","rejected_ode_steps","rhs_evaluations","median_gradient_seconds",
            "time_ratio_to_tsit5_1e8","max_prediction_difference_to_fine","max_absolute_nll_difference_to_fine",
            "max_scaled_gradient_difference_to_fine","max_fd_scaled_error","faster_rounds","numerical_screen_passed"),comparisons)
        reference_ok = accuracy["tsit5_1e10"] && accuracy["rodas5p_1e10"]
        candidates = Dict(label=>(reference_ok && accuracy[label] && expensive_fast[label] &&
            length(ratios[label])==12 && median(ratios[label])<=gates["max_median_paired_time_ratio"])
            for label in ("rodas5p_1e8","rodas5p_1e10"))
        metadata["numerical_screen_passed"] = accuracy
        metadata["fine_reference_screen_passed"] = reference_ok
        metadata["pilot_candidate_screen_passed"] = candidates
        metadata["median_paired_gradient_time_ratio"] = Dict(k=>median(v) for (k,v) in ratios if !isempty(v))
        metadata["profile_cases"] = length(reports)
        metadata["completed_cases"] = count(r->r["status"]=="complete",values(reports))
        metadata["timed_gradients"] = length(measurements)
        @test all(B.digest(joinpath(parent,f))==h for (f,h) in previous["artifact_sha256"])
        @test B.digest(joinpath(parent,"experiment.toml"))==metadata["parent_experiment_sha256"]
        @test all(B.digest(joinpath(B.ROOT,f))==h for (f,h) in metadata["implementation_sha256"])
        metadata["status"]="complete"
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
    println("Candidate screens: ",metadata["pilot_candidate_screen_passed"])
end

@testset "Runtime solver comparison integrity" begin
    run_comparison()
end
