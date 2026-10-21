include("../scripts/diagnose_cocoa_runtime.jl")

# Synthetic software fixture; no sampler performance or biological claim.
function runtime_diagnostic_fixture(directory)
    R = CocoaRuntimeDiagnostics
    trial = only(filter(t->t.trial==:ghhp4,R.B.benchmark_trials(;scope=:all)))
    path = joinpath(directory,"ghhp4")
    mkpath(path)
    settings = Dict("method"=>"bayes","scope"=>"all","trials"=>["ghhp4"],
        "samples_per_chain"=>3,"chains"=>1,"max_depth"=>12,"warmup_per_chain"=>2,
        "runtime_schema_version"=>2,"trial_seeds"=>Dict("ghhp4"=>141115),
        "solver"=>"tsit5","abstol"=>1e-6,"reltol"=>1e-6,
        "trial_solver_maxiters"=>Dict("ghhp4"=>10_000),
        "implementation_sha256"=>Dict(f=>R.B.digest(joinpath(R.B.ROOT,f)) for f in
            (R.C.F.LIKELIHOOD_IMPLEMENTATION...,"src/core/bayes.jl","scripts/fit_cocoa_benchmark_bayes.jl")))
    R.C.F.write_toml(joinpath(directory,"settings.toml"),settings)
    R.C.F.write_toml(joinpath(path,"input.toml"),Dict("trial"=>"ghhp4",
        "source_file"=>trial.source_file,"source_sha256"=>trial.source_sha256,
        "source_branch"=>trial.branch,"hidden_lactate"=>trial.hidden_lactate,
        "initial_normalized"=>trial.data.initial,"times_h"=>trial.data.times,
        "state_order"=>String.(collect(COCOA_STATES)),"seed"=>141115,"solver_maxiters"=>10_000))
    priors = R.C.benchmark_priors(;hidden_lactate=trial.hidden_lactate)
    R.C.F.write_toml(joinpath(path,"priors.toml"),Dict(string(n)=>
        Dict(string(k)=>v for (k,v) in pairs(priors[n].provenance)) for n in R.D.NAMES))
    R.C.F.write_toml(joinpath(path,"report.toml"),Dict("trial"=>"ghhp4",
        "status"=>"posterior_diagnostics_failed","retained_draws"=>3,
        "observations"=>105,"divergences"=>0,"depth_hits"=>0))
    R.B.write_csv(joinpath(path,"chain_001.csv"),
        ("chain","draw",string.(R.D.NAMES)...,"divergent","tree_depth"),
        [(1,i,ntuple(_->.5,25)...,false,3) for i in 1:3])
    durations = [.5,.7,1.,4.,2.]
    elapsed = cumsum(durations)
    rows = [(1,i,i==1 ? "initialization" : i==2 ? "warmup" : "sampling",max(0,i-2),
        elapsed[i],durations[i],
        (i==1 ? ntuple(_->missing,8) : (false,3,7,.1,.9,10.0+i,-.2,.3))...) for i in 1:5]
    R.B.write_csv(joinpath(path,"runtime.csv"),R.C.RUNTIME_HEADER,rows)
    R.B.write_csv(joinpath(path,"sampler_statistics.csv"),R.C.SAMPLER_HEADER,
        [(1,i-2,rows[i][7:end]...) for i in 3:5])
    return trial
end

@testset "Runtime records, corruption rejection, and slow-state selection" begin
    R = CocoaRuntimeDiagnostics
    mktempdir() do temporary
        parent = joinpath(temporary,"parent")
        trial = runtime_diagnostic_fixture(parent)
        saved = R.D.saved_draws(parent,trial)
        file = joinpath(parent,"ghhp4/runtime.csv")
        original = read(file,String)
        rows = R.read_runtime(file,saved.settings,saved.rows)
        @test length(rows)==5 && ismissing(first(rows).n_steps)
        @test getproperty.(rows,:draw)==[0,0,1,2,3]
        @test getproperty.(rows[2:end],:n_steps)==fill(7,4)
        summary = R.summarize_runtime(rows)
        @test getindex.(summary,2)==["initialization","warmup","sampling"]
        @test last(summary)[3]==3 && last(summary)[7]==21
        @test collect(last(summary)[4:6])≈[2.,3.8,4.]
        @test last(summary)[8]≈1/3
        selection = R.select_runtime_draws(rows,saved.rows)
        @test getproperty.(selection,:draw)==[1,2,3]
        @test getproperty.(selection,:reasons)==["preceding_slow_call","slow_retained_call","median_duration_control"]
        tied = [(;r...,duration=1.) for r in rows]
        @test getproperty.(R.select_runtime_draws(tied,saved.rows),:draw)==[1]
        for per_chain in (0,true,1.5)
            @test_throws ArgumentError R.select_runtime_draws(rows,saved.rows;per_chain)
        end
        # A missing value is not synthesized; nonfinite errors require divergence.
        for (line,column,value) in ((3,1,"2"),(3,2,"3"),(4,3,"warmup"),(4,4,"2"),
                (3,5,"0.1"),(3,6,"-1"),(4,7,"true"),(4,8,"2"),(3,9,"0"),
                (3,10,"NaN"),(3,11,"2"),(3,12,"Inf"),(3,13,"NaN"),(3,14,"Inf"),
                (2,9,"1"),(3,9,""))
            lines = split(chomp(original),'\n')
            fields = split(lines[line],',';keepempty=true)
            fields[column] = value
            lines[line] = join(fields,',')
            write(file,join(lines,'\n')*"\n")
            @test_throws ArgumentError R.read_runtime(file,saved.settings,saved.rows)
        end
        write(file,original)
        lines = split(chomp(original),'\n')
        fields = split(lines[3],',';keepempty=true)
        fields[7],fields[14] = "true","Inf" # warmup divergence with its original infinite error
        lines[3] = join(fields,',')
        write(file,join(lines,'\n')*"\n")
        @test isinf(R.read_runtime(file,saved.settings,saved.rows)[2].max_hamiltonian_energy_error)
        write(file,join(split(chomp(original),'\n')[1:end-1],'\n')*"\n")
        @test_throws ArgumentError R.read_runtime(file,saved.settings,saved.rows)
        write(file,original)
        @test_throws ArgumentError R.read_runtime(file,merge(saved.settings,Dict("runtime_schema_version"=>1)),saved.rows)
        @test_throws ArgumentError R.read_runtime(file,saved.settings,[saved.rows;first(saved.rows)])
        @test_throws ArgumentError R.read_runtime(file,saved.settings,saved.rows[1:2])
        sampler_file = joinpath(parent,"ghhp4/sampler_statistics.csv")
        sampler_original = read(sampler_file,String)
        for (before,after) in (("false,3,7","false,3,8"),("0.1,0.9","0.2,0.9"),
                ("0.9,13.0","0.8,13.0"),("13.0,-0.2","14.0,-0.2"),
                ("-0.2,0.3","-0.1,0.3"),("-0.2,0.3","-0.2,0.4"))
            write(sampler_file,replace(sampler_original,before=>after;count=1))
            @test_throws ArgumentError R.read_runtime(file,saved.settings,saved.rows)
        end
        write(sampler_file,sampler_original)
        # Initialization can also be retained when warmup is zero.
        zero_warmup = merge(saved.settings,Dict("warmup_per_chain"=>0))
        zero_draws = [(;r...,tree_depth=r.draw==1 ? 0 : r.tree_depth) for r in saved.rows]
        R.B.write_csv(file,R.C.RUNTIME_HEADER,
            [(1,i,i==1 ? "initialization" : "sampling",i,Float64(i),1.,
                (i==1 ? ntuple(_->missing,8) : (false,3,7,.1,.9,10.,-.2,.3))...) for i in 1:3])
        R.B.write_csv(sampler_file,R.C.SAMPLER_HEADER,
            [(1,i,(i==1 ? ntuple(_->missing,8) : (false,3,7,.1,.9,10.,-.2,.3))...) for i in 1:3])
        @test getproperty.(R.read_runtime(file,zero_warmup,zero_draws),:draw)==[1,2,3]
        counter = Ref(0)
        measured = R.measure_repeated(()->(counter[]+=1),3)
        @test measured.value==counter[]==4 # one unmeasured warm call
        @test length(measured.seconds)==length(measured.bytes)==length(measured.gc_seconds)==3
        @test all(>=(0),measured.seconds) && all(>=(0),measured.bytes) && all(>=(0),measured.gc_seconds)
        @test_throws ArgumentError R.measure_repeated(()->1,0)
    end
end

@testset "Saved runtime profile integration and provenance" begin
    R = CocoaRuntimeDiagnostics
    mktempdir() do temporary
        parent = joinpath(temporary,"parent")
        trial = runtime_diagnostic_fixture(parent)
        output = joinpath(temporary,"diagnostics")
        for options in ((per_chain=0,),(repeats=0,),(trials=[:unknown],))
            @test_throws ArgumentError R.diagnose_runtime(parent;output,options...)
            @test !ispath(output)
        end
        settings_file = joinpath(parent,"settings.toml")
        settings = read(settings_file,String)
        broken = TOML.parse(settings)
        broken["implementation_sha256"]["scripts/fit_cocoa_benchmark_bayes.jl"]="wrong"
        R.C.F.write_toml(settings_file,broken)
        @test_throws ArgumentError R.diagnose_runtime(parent;output)
        @test !ispath(output)
        write(settings_file,settings)
        result = R.diagnose_runtime(parent;output,repeats=2)
        @test result["status"]=="complete" && result["runtime_integrity_passed"] && !result["diagnostics_accept_fit"]
        @test result["runtime_events_checked"]==5 && result["retained_draws_checked"]==3 && result["profiled_states"]==3
        @test all(R.B.digest(joinpath(parent,p))==hash for (p,hash) in result["parent_sha256"])
        @test length(readlines(joinpath(output,"ghhp4/costs.csv")))==4
        saved = R.D.saved_draws(parent,trial)
        # Replay both accounting modes, using the same likelihood trajectory.
        saved_interval = (;saved...,settings=merge(saved.settings,Dict("solver_maxiters_scope"=>"output_interval")))
        interval = R.profile_point(trial,saved_interval,first(saved.rows).parameters;repeats=1)
        report = TOML.parsefile(joinpath(output,"ghhp4/chain_1_draw_1.toml"))
        @test interval["prediction_normalized"]==report["prediction_normalized"]
        @test interval["gradient_log"]==report["gradient_log"]
        @test report["status"]=="complete" && !report["diagnostics_accept_fit"]
        @test report["accepted_ode_steps"]>0 && report["rejected_ode_steps"]>=0 && report["rhs_evaluations"]>0
        @test length(report["gradient_seconds"])==length(report["ode_seconds"])==2
        @test all(>=(0),report["gradient_seconds"]) && all(>=(0),report["ode_seconds"])
        @test report["parameters"]==Dict(string(k)=>v for (k,v) in pairs(first(saved.rows).parameters))
        limited = (;saved_interval...,settings=merge(saved_interval.settings,
            Dict("trial_solver_maxiters"=>Dict("ghhp4"=>1))))
        failed = R.profile_point(trial,limited,first(saved.rows).parameters;repeats=1)
        @test failed["status"]=="failed" && haskey(failed,"error")
        @test !haskey(failed,"gradient_log") && !failed["diagnostics_accept_fit"]
        @test_throws ArgumentError R.diagnose_runtime(parent;output)
    end
end
