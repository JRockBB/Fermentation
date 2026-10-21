include("../scripts/fit_cocoa_benchmark_bayes.jl")
using Distributions: Normal, Cauchy, truncated, logpdf, quantile
using Random

@testset "2022 benchmark prior provenance" begin
    C = CocoaBenchmarkBayes
    for hidden in (false,true)
        priors = C.benchmark_priors(;hidden_lactate=hidden)
        @test keys(priors) == (COCOA_PARAMETERS...,:sigma)
        source = read(joinpath(C.B.ROOT,priors.mu1.provenance.file),String)
        normals = collect(eachmatch(r"(\w+)\s*~\s*normal\(0\.5,0\.3\)",source))
        @test Set(Symbol(m.captures[1]) for m in normals) == Set(COCOA_PARAMETERS)
        @test occursin(r"sigma\s*~\s*cauchy\(0,1\)",source)
        @test count(r"real\s*<lower=0>\s*\w+\s*;",source) == 25
        @test occursin(hidden ? "1.0E6" : "1.0E4",source)
        @test all(p.provenance.source == "doi:10.1098/rsos.210274" for p in values(priors))
        @test priors.mu1.provenance.sha256 == C.B.digest(joinpath(C.B.ROOT,priors.mu1.provenance.file))
        @test logpdf(priors.mu1.distribution,.4) == logpdf(truncated(Normal(.5,.3);lower=0.),.4)
        @test logpdf(priors.sigma.distribution,8.) == logpdf(truncated(Cauchy(0.,1.);lower=0.),8.)
        @test minimum(priors.sigma.distribution) == 0 && maximum(priors.sigma.distribution) == Inf
        @test quantile(priors.sigma.distribution,.99) > 5 # MLE upper bounds do not constrain this prior
    end
end

@testset "2022 pointwise likelihood and missing-state exports" begin
    C = CocoaBenchmarkBayes
    trials = C.B.benchmark_trials(;scope=:all)
    names = (COCOA_PARAMETERS...,:sigma)
    initial = NamedTuple{names}(ntuple(_->0.5,25))
    lower = NamedTuple{names}(ntuple(_->1e-5,25))
    upper = NamedTuple{names}(ntuple(_->5.,25))
    model = FitModel((p,d)->fill(p.mu1,length(d.values));initial,lower,upper)
    for id in (:brwb3,:ghhp4)
        trial = only(filter(t->t.trial==id,trials))
        likelihood = C.pointwise_loglikelihood(model,trial.data,initial)
        @test -sum(likelihood.values) ≈ gaussian_nll(model,trial.data,initial)
        @test likelihood.values ≈ [logpdf(Normal(initial.mu1,initial.sigma),trial.data.values[i])
            for i in likelihood.indices]
        @test length(C.observation_rows(trial)) == count(!ismissing,trial.data.values)
        @test all(row[3] > trial.data.t0 for row in C.observation_rows(trial))
    end
    trial = only(filter(t->t.trial==:ghhp4,trials))
    hidden_nan = FitModel((p,d)->[ismissing(v) ? NaN : p.mu1 for v in d.values];initial,lower,upper)
    @test C.pointwise_loglikelihood(hidden_nan,trial.data,initial).values ==
        C.pointwise_loglikelihood(model,trial.data,initial).values
    failed = FitModel((p,d)->nothing;initial,lower,upper)
    @test_throws ErrorException C.pointwise_loglikelihood(failed,trial.data,initial)
    nonfinite = FitModel((p,d)->fill(NaN,length(d.values));initial,lower,upper)
    @test_throws ErrorException C.pointwise_loglikelihood(nonfinite,trial.data,initial)
    wrong_length = FitModel((p,d)->[0.5];initial,lower,upper)
    @test_throws DimensionMismatch C.pointwise_loglikelihood(wrong_length,trial.data,initial)
    @test_throws ArgumentError C.pointwise_loglikelihood(model,trial.data,merge(initial,(sigma=0.,)))
    @test_throws ArgumentError C.pointwise_loglikelihood(model,
        merge(trial.data,(values=fill(missing,length(trial.data.values)),)),initial)
    draws = NamedTuple{names}(Tuple(n==:sigma ? reshape([.2,1.,6.,8.],2,2) :
        fill(initial[n],2,2) for n in names))
    posterior = (;draws)
    predictive = posterior_predictive(model,posterior,trial.data;draws=4,rng=MersenneTwister(224))
    result = (;draws,predictive)
    mktempdir() do output
        @test C.write_log_likelihood(output,model,trial,result) == 4
        likelihood_rows = [split(line,',') for line in readlines(joinpath(output,"log_likelihood.csv"))[2:end]]
        @test [r[1:2] for r in likelihood_rows] == [["1","1"],["1","2"],["2","1"],["2","2"]]
        @test all(length(r)==107 for r in likelihood_rows)
        @test sum(parse.(Float64,last(likelihood_rows)[3:end])) ≈
            -gaussian_nll(model,trial.data,merge(initial,(sigma=8.,)))
        observation_rows = readlines(joinpath(output,"likelihood_observations.csv"))[2:end]
        @test length(observation_rows) == 105 && !any(r->occursin(",L,",r),observation_rows)
        @test split(readlines(joinpath(output,"log_likelihood.csv"))[1],',';limit=3)[3] ==
            join((row[1] for row in C.observation_rows(trial)),',')
        C.write_predictive(output,trial,result)
        rows = [split(line,',';keepempty=true) for line in readlines(joinpath(output,"predictive.csv"))[2:end]]
        lactate = filter(r->r[2]=="L",rows)
        @test length(rows) == 120 && length(lactate) == 15
        @test all(isempty(r[3]) && all(!isempty,r[4:7]) && all(isempty,r[8:end]) for r in lactate)
        measured = filter(r->r[2]!="L",rows)
        @test all(parse(Float64,r[12]) ≈ parse(Float64,r[4])*parse(Float64,r[10]) for r in measured)
        @test all(parse(Float64,r[11]) ≈ parse(Float64,r[3])*parse(Float64,r[10]) for r in measured)
        selection = [parse.(Int,split(line,',')) for line in
            readlines(joinpath(output,"predictive_selection.csv"))[2:end]]
        @test [(r[2]-1)*2+r[3] for r in selection] == predictive.draw_indices
    end
end

@testset "2022 Bayesian benchmark preflight and checkpoint integration" begin
    C = CocoaBenchmarkBayes
    mktempdir() do temporary
        output = joinpath(temporary,"invalid")
        for kwargs in ((samples=9,),(warmup=-1,),(chains=0,),(max_depth=0,),
                (seed=-1,),(target_accept=1.,),(solver_maxiters=0,),
                (predictive_draws=0,),(chain_execution=:unknown,),(trials=[:unknown],),
                (solver_maxiters_scope=:unknown,),(solver_max_step_attempts=0,),
                (max_sampling_seconds=0.,),(max_sampling_seconds=NaN,),(stop_on_divergence=1,))
            @test_throws ArgumentError C.run_benchmark_bayes(;output,kwargs...)
            @test !ispath(output)
        end
        output = joinpath(temporary,"smoke")
        comparison = C.run_benchmark_bayes(;output,scope=:all,trials=[:ghhp4],
            samples=10,warmup=5,chains=1,max_depth=2,target_accept=.8,
            predictive_draws=10,solver_maxiters=10_000)
        directory = joinpath(output,"ghhp4")
        report = TOML.parsefile(joinpath(directory,"report.toml"))
        settings = TOML.parsefile(joinpath(output,"settings.toml"))
        input = TOML.parsefile(joinpath(directory,"input.toml"))
        @test length(comparison) == 1 && only(comparison)[1] == :ghhp4
        @test !report["diagnostics_passed"] && report["status"] == "posterior_diagnostics_failed"
        @test !report["direct_metric_comparison_available"]
        @test report["retained_draws"] == report["likelihood_draws"] == 10
        @test report["observations"] == 105
        @test settings["trial_seeds"]["ghhp4"] == input["seed"] == 141085+30
        @test input["solver_maxiters"] == 10_000 && input["hidden_lactate"]
        @test input["solver_maxiters_scope"] == settings["solver_maxiters_scope"] == "output_interval"
        @test input["solver_max_step_attempts"] == settings["solver_max_step_attempts"] == 500
        runtime = [split(line,',';keepempty=true) for line in readlines(joinpath(directory,"runtime.csv"))[2:end]]
        @test length(runtime) == 15
        @test settings["runtime_schema_version"] == 2
        @test split(first(readlines(joinpath(directory,"runtime.csv"))),',') == collect(C.RUNTIME_HEADER)
        @test parse.(Int,getindex.(runtime,2)) == collect(1:15)
        @test first(runtime)[3] == "initialization" && all(isempty,first(runtime)[7:end])
        @test all(r->parse(Int,r[9])>0 && parse(Float64,r[10])>0,runtime[2:end])
        @test all(r->0<=parse(Float64,r[11])<=1,runtime[2:end])
        sampler_rows = [split(line,',';keepempty=true) for line in
            readlines(joinpath(directory,"sampler_statistics.csv"))[2:end]]
        @test length(sampler_rows)==report["sampler_statistics_draws"]==10
        @test all(raw[3:end]==event[7:end] for (raw,event) in zip(sampler_rows,runtime[6:end]))
        @test all(r->r[3]=="warmup" && r[4]=="0",runtime[2:5])
        @test parse.(Int,getindex.(runtime[6:end],4)) == collect(1:10)
        @test all(>=(0),diff(parse.(Float64,getindex.(runtime,5))))
        @test length(TOML.parsefile(joinpath(directory,"priors.toml"))) == 25
        @test length(readlines(joinpath(directory,"parameters.csv"))) == 26
        checkpoints = [split(line,',') for line in readlines(joinpath(directory,"chain_001.csv"))[2:end]]
        @test length(checkpoints) == 10 && parse.(Int,getindex.(checkpoints,2)) == collect(1:10)
        @test count(r->parse(Bool,r[end-1]),checkpoints) == report["divergences"]
        @test count(r->parse(Int,r[end])>=2,checkpoints) == report["depth_hits"]
        likelihood_rows = readlines(joinpath(directory,"log_likelihood.csv"))[2:end]
        @test length(likelihood_rows) == 10 && all(length(split(r,','))==107 for r in likelihood_rows)
        trial = only(filter(t->t.trial==:ghhp4,C.B.benchmark_trials(;scope=:all)))
        parameters = NamedTuple{(COCOA_PARAMETERS...,:sigma)}(Tuple(parse.(Float64,first(checkpoints)[3:end-2])))
        lower = NamedTuple{keys(parameters)}(ntuple(_->1e-5,25))
        upper = NamedTuple{keys(parameters)}(ntuple(_->5.,25))
        initial = NamedTuple{keys(parameters)}(ntuple(_->.5,25))
        model = cocoa_fit_model(;initial,lower,upper,maxiters=10_000)
        @test -sum(parse.(Float64,split(first(likelihood_rows),',')[3:end])) ≈
            gaussian_nll(model,trial.data,parameters)
        @test_throws ArgumentError C.run_benchmark_bayes(;output)
        stopped = joinpath(temporary,"stopped")
        @test_throws C.SamplingStopped C.run_benchmark_bayes(;output=stopped,scope=:all,
            trials=[:ghhp4],samples=10,warmup=5,chains=1,max_depth=2,target_accept=.8,
            predictive_draws=10,max_sampling_seconds=eps())
        stopped_path = joinpath(stopped,"ghhp4")
        stopped_report = TOML.parsefile(joinpath(stopped_path,"report.toml"))
        @test stopped_report["status"] == "sampling_stopped" && !stopped_report["diagnostics_passed"]
        @test stopped_report["completed_sampler_calls"] == [1]
        @test stopped_report["partial_retained_draws"] == [0]
        @test length(readlines(joinpath(stopped_path,"chain_001.csv"))) == 1
        @test !isfile(joinpath(stopped_path,"parameters.csv")) && !isfile(joinpath(stopped,"comparison.csv"))
    end
end
