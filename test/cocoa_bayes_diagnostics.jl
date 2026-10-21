include("../scripts/diagnose_cocoa_bayes.jl")
using Distributions: Normal, Cauchy, truncated, logpdf

@testset "Log-coordinate posterior numerical checks" begin
    D = CocoaBayesDiagnostics
    initial,lower,upper = (alpha=.5,sigma=.5),(alpha=.01,sigma=.01),(alpha=5.,sigma=5.)
    model = FitModel((p,d)->fill(p.alpha,length(d.values));initial,lower,upper)
    data = (values=Union{Missing,Float64}[.2,missing,.8,.6],sigma=:sigma)
    priors = (alpha=(distribution=truncated(Normal(.5,.3);lower=0.),),
        sigma=(distribution=truncated(Cauchy(0.,1.);lower=0.),))
    for p in ((alpha=.7,sigma=.4),(alpha=.7,sigma=8.))
        check = D.point_check(model,data,priors,p)
        y = collect(skipmissing(data.values))
        residual = p.alpha .- y
        expected = [p.alpha*sum(residual)/p.sigma^2 + p.alpha*(p.alpha-.5)/.3^2 - 1,
            length(y)-sum(abs2,residual)/p.sigma^2 + 2p.sigma^2/(1+p.sigma^2) - 1]
        target = -sum(logpdf(Normal(p.alpha,p.sigma),v) for v in y) -
            logpdf(priors.alpha.distribution,p.alpha) -
            logpdf(priors.sigma.distribution,p.sigma) - log(p.alpha) - log(p.sigma)
        @test check["status"] == "complete" && !check["diagnostics_accept_fit"]
        @test check["negative_log_density_log_coordinates"] ≈ target
        @test check["gradient_log"] ≈ expected
        @test check["parameter_order"] == ["alpha","sigma"]
        @test all(c["max_scaled_error"]<1e-7 for c in check["finite_differences"])
        @test check["nll"] ≈ -sum(logpdf(Normal(p.alpha,p.sigma),v) for v in y)
    end
    for predict in ((p,d)->nothing,(p,d)->fill(NaN,length(d.values)),
            (p,d)->[ismissing(v) ? NaN : p.alpha for v in d.values])
        failed = D.point_check(FitModel(predict;initial,lower,upper),data,priors,initial)
        @test failed["status"] == "failed" && haskey(failed,"error")
        @test !failed["diagnostics_accept_fit"] && !haskey(failed,"gradient_log")
    end
    malformed = FitModel((p,d)->[p.alpha];initial,lower,upper)
    @test_throws DimensionMismatch D.point_check(malformed,data,priors,initial)
    bug = FitModel((p,d)->error("programming error");initial,lower,upper)
    @test_throws ErrorException D.point_check(bug,data,priors,initial)
    @test_throws ArgumentError D.point_check(model,data,priors,(alpha=0.,sigma=.5))
    @test_throws ArgumentError D.point_check(model,data,priors,(alpha=.5,sigma=Inf))
    @test_throws ArgumentError D.point_check(model,data,priors,(alpha=.5,))
    @test_throws ArgumentError D.point_check(model,data,priors,initial;fd_step=0.)
end

@testset "Deterministic posterior diagnostic subset" begin
    D = CocoaBayesDiagnostics
    rows = [(chain=c,draw=i,parameters=(yc11=[.6,.1,.9,.4,.5][i],sigma=.2),
        divergent=(c==1 && i==4),tree_depth=3) for c in 1:2 for i in 1:5]
    selected = D.select_draws(rows;per_chain=1)
    @test [(r.chain,r.draw) for r in selected] == [(1,2),(1,3),(1,4),(2,2),(2,3)]
    @test selected[2].reasons == "grid;preceding_divergent;yc11_maximum"
    @test selected[3].reasons == "divergent_retained"
    @test D.select_draws(rows;per_chain=1) == selected
    @test length(D.select_draws(rows;per_chain=99)) == 10
    @test [(r.chain,r.draw) for r in D.select_draws(rows;per_chain=2,tail_parameters=Symbol[])] ==
        [(1,1),(1,3),(1,4),(1,5),(2,1),(2,5)]
    for kwargs in ((per_chain=0,),(per_chain=true,),(per_chain=1.5,),
            (tail_parameters=[:unknown],),(tail_parameters=[:yc11,:yc11],),
            (tail_parameters=["yc11"],))
        @test_throws ArgumentError D.select_draws(rows;kwargs...)
    end
    @test_throws ArgumentError D.select_draws(NamedTuple[])
end

# Synthetic checkpoint values exercise the reader and solvers, not fitted estimates.
function bayes_diagnostic_fixture(directory,trial)
    D = CocoaBayesDiagnostics
    C,B = D.C,D.B
    names = D.NAMES
    path = joinpath(directory,string(trial.trial))
    mkpath(path)
    settings = Dict("method"=>"bayes","scope"=>"all","trials"=>[string(trial.trial)],
        "samples_per_chain"=>3,"chains"=>1,"max_depth"=>12,
        "trial_seeds"=>Dict(string(trial.trial)=>141115),"solver"=>"tsit5",
        "abstol"=>1e-6,"reltol"=>1e-6,
        "trial_solver_maxiters"=>Dict(string(trial.trial)=>10_000),
        "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f)) for f in
            (C.F.LIKELIHOOD_IMPLEMENTATION...,"src/core/bayes.jl")))
    C.F.write_toml(joinpath(directory,"settings.toml"),settings)
    input = Dict("trial"=>string(trial.trial),"source_file"=>trial.source_file,
        "source_sha256"=>trial.source_sha256,"source_branch"=>trial.branch,
        "hidden_lactate"=>trial.hidden_lactate,"initial_normalized"=>trial.data.initial,
        "times_h"=>trial.data.times,"state_order"=>String.(collect(COCOA_STATES)),
        "seed"=>141115,"solver_maxiters"=>10_000)
    C.F.write_toml(joinpath(path,"input.toml"),input)
    priors = C.benchmark_priors(;hidden_lactate=trial.hidden_lactate)
    C.F.write_toml(joinpath(path,"priors.toml"),Dict(string(n)=>
        Dict(string(k)=>v for (k,v) in pairs(priors[n].provenance)) for n in names))
    C.F.write_toml(joinpath(path,"report.toml"),Dict("trial"=>string(trial.trial),
        "status"=>"posterior_diagnostics_failed","retained_draws"=>3,
        "observations"=>count(!ismissing,trial.data.values),"divergences"=>0,"depth_hits"=>0))
    B.write_csv(joinpath(path,"chain_001.csv"),
        ("chain","draw",String.(collect(names))...,"divergent","tree_depth"),
        [(1,i,ntuple(_->.5,length(names))...,false,3) for i in 1:3])
    return path
end

@testset "Checkpoint integrity and numerical report integration" begin
    D = CocoaBayesDiagnostics
    trial = only(filter(t->t.trial==:ghhp4,D.B.benchmark_trials(;scope=:all)))
    mktempdir() do temporary
        parent = joinpath(temporary,"parent")
        path = bayes_diagnostic_fixture(parent,trial)
        saved = D.saved_draws(parent,trial)
        @test length(saved.rows) == 3
        @test saved.provenance["source_sha256"] == trial.source_sha256
        @test saved.rows[2].chain == 1 && saved.rows[2].draw == 2
        original = read(joinpath(path,"chain_001.csv"),String)
        lines = readlines(joinpath(path,"chain_001.csv"))
        for text in (join(lines[1:3],'\n')*"\n",
                replace(original,"chain,draw"=>"draw,chain"),
                replace(original,"1,2,"=>"1,1,"),
                replace(original,"1,2,"=>"2,2,"),
                replace(original,"0.5"=>"0.0";count=1),
                replace(original,"0.5"=>"NaN";count=1),
                replace(original,"false,3"=>"true,3";count=1),
                replace(original,"false,3"=>"false,12";count=1),
                replace(original,"false,3"=>"false,13";count=1))
            write(joinpath(path,"chain_001.csv"),text)
            @test_throws ArgumentError D.saved_draws(parent,trial)
        end
        write(joinpath(path,"chain_001.csv"),original)
        mutations = [
            ("settings.toml",d->(d["implementation_sha256"]["src/core/bayes.jl"]="wrong")),
            ("settings.toml",d->(d["samples_per_chain"]=0)),
            ("ghhp4/input.toml",d->(d["source_sha256"]="wrong")),
            ("ghhp4/input.toml",d->(d["hidden_lactate"]=false)),
            ("ghhp4/input.toml",d->(d["seed"]=0)),
            ("ghhp4/priors.toml",d->(d["sigma"]["source"]="wrong")),
            ("ghhp4/report.toml",d->(d["status"]="sampling")),
            ("ghhp4/report.toml",d->(d["observations"]=120)),
            ("ghhp4/report.toml",d->(d["retained_draws"]=2))]
        invalid_output = joinpath(temporary,"invalid")
        for (file,mutate!) in mutations
            full = joinpath(parent,file)
            before = read(full,String)
            d = TOML.parse(before)
            mutate!(d)
            D.C.F.write_toml(full,d)
            @test_throws ArgumentError D.diagnose_bayes(parent;output=invalid_output)
            @test !ispath(invalid_output)
            write(full,before)
        end
        for kwargs in ((per_chain=0,),(tail_parameters=[:unknown],),
                (tolerance_factor=0.,),(tolerance_factor=1.,),(fd_step=NaN,),
                (trials=[:brwb3],))
            @test_throws ArgumentError D.diagnose_bayes(parent;output=invalid_output,kwargs...)
            @test !ispath(invalid_output)
        end
        output = joinpath(temporary,"diagnostics")
        rows = D.diagnose_bayes(parent;output,per_chain=1,tail_parameters=Symbol[])
        @test length(rows)==3 && all(r[10]=="complete" for r in rows)
        @test length(readlines(joinpath(output,"comparison.csv")))==4
        settings = TOML.parsefile(joinpath(output,"settings.toml"))
        @test !settings["diagnostics_accept_fit"] && settings["trials"]==["ghhp4"]
        @test settings["parent_settings_sha256"]==D.B.digest(joinpath(parent,"settings.toml"))
        @test read(joinpath(path,"chain_001.csv"),String)==original
        checks = [TOML.parsefile(joinpath(output,"ghhp4","chain_001_draw_000002_$case.toml"))
            for case in ("saved","tighter","alternate_tighter")]
        @test all(c["status"]=="complete" && !c["diagnostics_accept_fit"] for c in checks)
        @test all(length(c["prediction_normalized"])==120 && length(c["gradient_log"])==25 for c in checks)
        @test all(c["parameter_order"]==String.(collect(D.NAMES)) for c in checks)
        @test checks[1]["nll"] ≈ gaussian_nll(cocoa_fit_model(
            ;initial=saved.rows[2].parameters,
            lower=NamedTuple{D.NAMES}(ntuple(_->1e-5,25)),
            upper=NamedTuple{D.NAMES}(ntuple(_->5.,25))),trial.data,saved.rows[2].parameters)
        @test checks[1]["max_scaled_gradient_difference_from_saved"]==0.
        @test checks[3]["solver"]=="rodas5p" && checks[3]["abstol"]≈1e-8
        @test all(c["finite_differences"][2]["max_scaled_error"]<.01 for c in checks)
        @test all(abs(c["delta_nll_from_saved"])<.01 for c in checks)
        @test all(c["max_prediction_difference_from_saved"]<.001 for c in checks)
        @test_throws ArgumentError D.diagnose_bayes(parent;output)
    end
end
