include("../scripts/audit_cocoa_benchmark_bayes.jl")
using Random: MersenneTwister, randn

# Synthetic checkpoints with real source metadata, not estimates or MCMC output.
function bayes_audit_fixture(directory;budget_scope=nothing)
    A = CocoaBayesAudit
    C,B = A.C,A.B
    names = A.D.NAMES
    trials = filter(t->t.trial in (:brwb3,:ghhp4),B.benchmark_trials(;scope=:all))
    settings = Dict("method"=>"bayes","scope"=>"all","trials"=>string.([t.trial for t in trials]),
        "time_unit"=>"h","source_concentration_unit"=>"mg/g pulp",
        "samples_per_chain"=>20,"chains"=>2,"max_depth"=>12,"predictive_draws"=>6,
        "trial_seeds"=>Dict(string(t.trial)=>141085+i-1 for (i,t) in enumerate(B.benchmark_trials(;scope=:all))),
        "solver"=>"tsit5","abstol"=>1e-6,"reltol"=>1e-6,
        "trial_solver_maxiters"=>Dict(string(t.trial)=>(t.hidden_lactate ? 1_000_000 : 10_000) for t in trials),
        "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in
            ("scripts/fit_cocoa_benchmark_bayes.jl","scripts/fit_cocoa_benchmark.jl",
             "scripts/cocoa_fit_output.jl","src/core/bayes.jl",C.F.LIKELIHOOD_IMPLEMENTATION...)))
    if !isnothing(budget_scope)
        settings["solver_maxiters_scope"] = string(budget_scope)
        settings["solver_max_step_attempts"] = 500
    end
    mkpath(directory)
    C.F.write_toml(joinpath(directory,"settings.toml"),settings)
    comparisons = Tuple[]
    rng = MersenneTwister(2531)
    for trial in trials
        id = string(trial.trial)
        path = joinpath(directory,id)
        mkpath(path)
        C.F.write_toml(joinpath(path,"input.toml"),Dict("trial"=>id,"source_file"=>trial.source_file,
            "source_sha256"=>trial.source_sha256,"source_branch"=>trial.branch,
            "hidden_lactate"=>trial.hidden_lactate,"initial_normalized"=>trial.data.initial,
            "times_h"=>trial.data.times,"state_order"=>String.(collect(COCOA_STATES)),
            "seed"=>settings["trial_seeds"][id],"solver_maxiters"=>settings["trial_solver_maxiters"][id]))
        if !isnothing(budget_scope)
            input = TOML.parsefile(joinpath(path,"input.toml"))
            input["solver_maxiters_scope"] = string(budget_scope)
            input["solver_max_step_attempts"] = 500
            C.F.write_toml(joinpath(path,"input.toml"),input)
        end
        priors = C.benchmark_priors(;hidden_lactate=trial.hidden_lactate)
        C.F.write_toml(joinpath(path,"priors.toml"),Dict(string(n)=>
            Dict(string(k)=>v for (k,v) in pairs(priors[n].provenance)) for n in names))
        draws = NamedTuple{names}(Tuple(.5 .+ .005 .* randn(rng,20,2) for _ in names))
        draws.sigma[1,1] = 8. # Bayesian support extends above the FitModel MLE template bound.
        summary = FermentModels.CoreModels._posterior_summary(draws)
        for chain in 1:2
            B.write_csv(joinpath(path,"chain_$(lpad(chain,3,'0')).csv"),
                ("chain","draw",String.(collect(names))...,"divergent","tree_depth"),
                [(chain,i,(draws[n][i,chain] for n in names)...,chain==1 && i==1,3) for i in 1:20])
        end
        B.write_csv(joinpath(path,"parameters.csv"),
            ("parameter","mean","median","lower_95","upper_95","rhat","ess_bulk","ess_tail"),
            [(n,(getproperty(summary[n],k) for k in (:mean,:median,:lower,:upper,:rhat,:ess_bulk,:ess_tail))...) for n in names])
        initial = NamedTuple{names}(ntuple(_->.5,25))
        lower,upper = NamedTuple{names}(ntuple(_->1e-5,25)),NamedTuple{names}(ntuple(_->5.,25))
        model = cocoa_fit_model(;initial,lower,upper,maxiters=settings["trial_solver_maxiters"][id],
            maxiters_scope=isnothing(budget_scope) ? :solve : budget_scope)
        predictive = posterior_predictive(model,(;draws),trial.data;draws=6,rng)
        result = (;draws,predictive)
        C.write_predictive(path,trial,result)
        C.write_log_likelihood(path,model,trial,result)
        report = Dict("trial"=>id,"status"=>"posterior_diagnostics_failed","diagnostics_passed"=>false,
            "retained_draws"=>40,"likelihood_draws"=>40,"observations"=>count(!ismissing,trial.data.values),
            "divergences"=>1,"depth_hits"=>0,"rhat_threshold"=>1.01,"minimum_ess"=>200,
            "max_rhat"=>maximum(r.rhat for r in values(summary)),
            "min_ess_bulk"=>minimum(r.ess_bulk for r in values(summary)),
            "min_ess_tail"=>minimum(r.ess_tail for r in values(summary)),
            "elapsed_seconds"=>0.,"direct_metric_comparison_available"=>false)
        C.F.write_toml(joinpath(path,"report.toml"),report)
        reference = only(filter(r->r.trial==trial.trial,B.published_baseline()))
        push!(comparisons,(trial.trial,trial.in_paper_benchmark,report["observations"],report["status"],
            report["max_rhat"],report["min_ess_bulk"],report["min_ess_tail"],1,0,40,0.,reference.success,
            reference.max_rhat,reference.min_bulk_ess,reference.min_tail_ess,reference.bma_weight,missing,missing))
    end
    B.write_csv(joinpath(directory,"comparison.csv"),
        ("trial","in_paper_benchmark","observations","status","max_rhat","min_ess_bulk","min_ess_tail",
         "divergences","depth_hits","retained_draws","elapsed_seconds","published_baseline_success",
         "published_max_rhat","published_min_bulk_ess","published_min_tail_ess","published_bma_weight",
         "psis_loo","published_psis_loo"),comparisons)
    return trials
end

@testset "Saved interval-budget replay" begin
    A = CocoaBayesAudit
    mktempdir() do temporary
        parent = joinpath(temporary,"parent")
        bayes_audit_fixture(parent;budget_scope=:output_interval)
        result = A.audit_bayes(parent;output=joinpath(temporary,"audit.toml"),replay_per_chain=20)
        @test result["solver_budget_scope"] == "output_interval"
        @test result["solver_max_step_attempts"] == 500
        @test all(t["likelihood_replay_complete"] for t in values(result["trials"]))
    end
end

function audit_csv_edit(text,line,column,value)
    lines = split(chomp(text),'\n')
    fields = split(lines[line],',';keepempty=true)
    fields[column] = value
    lines[line] = join(fields,',')
    return join(lines,'\n')*"\n"
end

@testset "Saved Bayesian artifact audit" begin
    A = CocoaBayesAudit
    @test A.matches(NaN,NaN) && !A.matches(NaN,1.)
    @test !A.matches([1.,2.],[1.])
    mktempdir() do temporary
        parent = joinpath(temporary,"parent")
        trials = bayes_audit_fixture(parent)
        output = joinpath(temporary,"audit.toml")
        before = Dict(replace(relpath(joinpath(dir,file),parent),'\\'=>'/')=>A.B.digest(joinpath(dir,file))
            for (dir,_,files) in walkdir(parent) for file in files)
        result = A.audit_bayes(parent;output)
        @test result["integrity_passed"]
        @test result["solver_budget_scope"] == "solve"
        @test result["checks"]>500
        @test Set(keys(result["trials"]))==Set(["brwb3","ghhp4"])
        for id in ("brwb3","ghhp4")
            report = result["trials"][id]
            @test report["integrity_passed"] && !report["sampling_diagnostics_passed"]
            @test report["likelihood_rows_checked"]==40
            @test report["likelihood_draws_replayed"]==6
            @test !report["likelihood_replay_complete"]
            @test report["latent_predictive_draws_replayed"]==6
            @test report["likelihood_values_checked"]==40*(id=="brwb3" ? 96 : 105)
        end
        @test TOML.parsefile(output)==result
        @test all(A.B.digest(joinpath(parent,f))==h for (f,h) in before)
        @test result["parent_sha256"]==before
        @test all(A.B.digest(joinpath(A.B.ROOT,f))==h for (f,h) in result["implementation_sha256"])
        subset = A.audit_bayes(parent;output=joinpath(temporary,"subset.toml"),trials=[:ghhp4],replay_per_chain=1)
        @test collect(keys(subset["trials"]))==["ghhp4"]
        @test subset["trials"]["ghhp4"]["likelihood_draws_replayed"]==2
        @test !haskey(subset["parent_sha256"],"brwb3/chain_001.csv")
        complete = A.audit_bayes(parent;output=joinpath(temporary,"complete.toml"),trials=[:ghhp4],replay_per_chain=99)
        @test complete["trials"]["ghhp4"]["likelihood_draws_replayed"]==40
        @test complete["trials"]["ghhp4"]["likelihood_replay_complete"]
        @test_throws ArgumentError A.audit_bayes(parent;output)
        invalid = joinpath(temporary,"invalid.toml")
        for kwargs in ((replay_per_chain=0,),(replay_per_chain=true,),(replay_per_chain=1.5,),
                (trials=[:unknown],),(trials=[:ghhp4,:ghhp4],))
            @test_throws ArgumentError A.audit_bayes(parent;output=invalid,kwargs...)
            @test !ispath(invalid)
        end
        toml_mutations = [
            ("settings.toml",d->(d["implementation_sha256"]["scripts/fit_cocoa_benchmark_bayes.jl"]="wrong")),
            ("settings.toml",d->(d["time_unit"]="d")),
            ("settings.toml",d->(d["source_concentration_unit"]="mmol/L")),
            ("settings.toml",d->(d["solver_maxiters_scope"]="output_interval")),
            ("ghhp4/input.toml",d->(d["solver_max_step_attempts"]=1)),
            ("ghhp4/report.toml",d->(d["max_rhat"]+=1)),
            ("ghhp4/report.toml",d->(d["minimum_ess"]=1)),
            ("ghhp4/report.toml",d->(d["diagnostics_passed"]=true)),
            ("ghhp4/report.toml",d->(d["status"]="posterior_diagnostics_passed")),
            ("ghhp4/report.toml",d->(d["likelihood_draws"]=39)),
            ("ghhp4/report.toml",d->(d["direct_metric_comparison_available"]=true))]
        for (file,mutate!) in toml_mutations
            path = joinpath(parent,file)
            original = read(path,String)
            try
                d = TOML.parse(original)
                mutate!(d)
                A.C.F.write_toml(path,d)
                @test_throws ArgumentError A.audit_bayes(parent;output=invalid,trials=[:ghhp4])
                @test !ispath(invalid)
            finally
                write(path,original)
            end
        end
        csv_mutations = [
            ("parameters.csv",t->audit_csv_edit(t,2,2,"99.0")),
            ("parameters.csv",t->audit_csv_edit(t,2,6,"NaN")),
            ("log_likelihood.csv",t->join(split(chomp(t),'\n')[1:end-1],'\n')*"\n"),
            ("log_likelihood.csv",t->t*last(split(chomp(t),'\n'))*"\n"),
            ("log_likelihood.csv",t->audit_csv_edit(t,3,2,"1")),
            ("log_likelihood.csv",t->audit_csv_edit(t,3,3,"NaN")),
            ("log_likelihood.csv",t->audit_csv_edit(t,2,3,"99.0")),
            ("log_likelihood.csv",t->audit_csv_edit(t,1,3,"wrong")),
            ("likelihood_observations.csv",t->audit_csv_edit(t,2,2,"2")),
            ("likelihood_observations.csv",t->audit_csv_edit(t,2,3,"999.0")),
            ("likelihood_observations.csv",t->audit_csv_edit(t,2,4,"L")),
            ("likelihood_observations.csv",t->audit_csv_edit(t,2,5,"99.0")),
            ("predictive_selection.csv",t->audit_csv_edit(t,2,2,"3")),
            ("predictive_selection.csv",t->audit_csv_edit(t,2,3,"21")),
            ("predictive_selection.csv",t->begin
                first_fields = split(split(t,'\n')[2],',')
                audit_csv_edit(audit_csv_edit(t,3,2,first_fields[2]),3,3,first_fields[3])
            end),
            ("predictive.csv",t->audit_csv_edit(t,2,1,"999.0")),
            ("predictive.csv",t->audit_csv_edit(t,2,2,"L")),
            ("predictive.csv",t->audit_csv_edit(t,2,3,"99.0")),
            ("predictive.csv",t->audit_csv_edit(t,2,4,"99.0")),
            ("predictive.csv",t->audit_csv_edit(t,2,8,"999.0")),
            ("predictive.csv",t->audit_csv_edit(t,2,10,"99.0")),
            ("predictive.csv",t->audit_csv_edit(t,2,11,"99.0")),
            ("predictive.csv",t->audit_csv_edit(t,2,12,"99.0")),
            # Hidden lactate occupies state-major rows 46:60 (CSV lines 47:61).
            ("predictive.csv",t->audit_csv_edit(t,47,8,"0.1")),
            ("predictive.csv",t->audit_csv_edit(t,47,12,"0.1"))]
        for (file,mutate) in csv_mutations
            path = joinpath(parent,"ghhp4",file)
            original = read(path,String)
            try
                write(path,mutate(original))
                @test_throws ArgumentError A.audit_bayes(parent;output=invalid,trials=[:ghhp4])
                @test !ispath(invalid)
            finally
                write(path,original)
            end
        end
        for (column,value) in ((2,"false"),(5,"99.0"),(12,"false"),(13,"99.0"),
                (16,"99.0"),(17,"1.0"))
            path = joinpath(parent,"comparison.csv")
            original = read(path,String)
            try
                write(path,audit_csv_edit(original,3,column,value))
                @test_throws ArgumentError A.audit_bayes(parent;output=invalid,trials=[:ghhp4])
                @test !ispath(invalid)
            finally
                write(path,original)
            end
        end
        # Full replay also checks finite likelihood corruption outside the default grid.
        path = joinpath(parent,"ghhp4","log_likelihood.csv")
        original = read(path,String)
        try
            write(path,audit_csv_edit(original,3,3,"99.0"))
            @test_throws ArgumentError A.audit_bayes(parent;output=invalid,trials=[:ghhp4],replay_per_chain=99)
            @test !ispath(invalid)
        finally
            write(path,original)
        end
        @test all(A.B.digest(joinpath(parent,f))==h for (f,h) in before)
        # Published failures and screened-out trials have no published numeric metrics.
        available = A.B.benchmark_trials(;scope=:all)
        references = Dict(r.trial=>r for r in A.B.published_baseline())
        absent = first(filter(t->!t.in_paper_benchmark,available))
        failed = first(filter(t->haskey(references,t.trial) && !references[t.trial].success,available))
        for trial in (absent,failed)
            comparison_directory = joinpath(temporary,string(trial.trial))
            mkpath(comparison_directory)
            report = TOML.parsefile(joinpath(parent,"ghhp4/report.toml"))
            report["observations"] = count(!ismissing,trial.data.values)
            saved = (settings=Dict("trials"=>[string(trial.trial)]),report=report)
            header = split(readline(joinpath(parent,"comparison.csv")),',')
            row = (trial.trial,trial.in_paper_benchmark,report["observations"],report["status"],
                report["max_rhat"],report["min_ess_bulk"],report["min_ess_tail"],
                report["divergences"],report["depth_hits"],report["retained_draws"],
                report["elapsed_seconds"],false,missing,missing,missing,missing,missing,missing)
            path = joinpath(comparison_directory,"comparison.csv")
            A.B.write_csv(path,header,[row])
            @test isnothing(A.audit_comparison!(A.Checks(0),comparison_directory,saved,trial))
            original = read(path,String)
            write(path,audit_csv_edit(original,2,16,"0.5"))
            @test_throws ArgumentError A.audit_comparison!(A.Checks(0),comparison_directory,saved,trial)
            write(path,audit_csv_edit(original,2,12,"true"))
            @test_throws ArgumentError A.audit_comparison!(A.Checks(0),comparison_directory,saved,trial)
        end
        A.C.F.write_toml(joinpath(parent,"aborted_run.toml"),Dict("status"=>"interrupted","fit_accepted"=>false))
        @test_throws ArgumentError A.audit_bayes(parent;output=invalid)
        @test !ispath(invalid)
    end
end
