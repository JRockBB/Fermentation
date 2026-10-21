include("../scripts/fit_cocoa_benchmark.jl")
using .CocoaBenchmarkFits
using .CocoaBenchmarkFits.CocoaBenchmarkData
isdefined(@__MODULE__,:mle_restart_fixture) || include("cocoa_benchmark_fixture.jl")

@testset "2022 benchmark provenance and observation extraction" begin
    B = CocoaBenchmarkData
    @test B.verify_benchmark() == 78
    @test B.BENCHMARK_STATES == COCOA_STATES
    all_trials = B.benchmark_trials(;scope=:all)
    benchmark = B.benchmark_trials()
    @test length(all_trials) == 32
    @test [r.trial for r in benchmark] == collect(B.BENCHMARK_TRIALS)
    @test count(r->r.hidden_lactate, all_trials) == 4
    @test Set(r.trial for r in all_trials if !r.in_paper_benchmark) ==
        Set((:cipb1,:mywb1,:mywb2,:ghhp6,:ghhp7,:brwb5,:brwb6,:brwb8,:ecpb2))
    for trial in all_trials
        d = trial.data
        @test size(trial.observations) == (length(d.times),8)
        @test count(!ismissing,d.values) == length(d.times)*(trial.hidden_lactate ? 7 : 8)
        @test all(v -> 0 <= v <= 1, skipmissing(d.values))
        @test all(>(0), diff([d.t0;d.times]))
        indices = findall(!ismissing,d.scales)
        @test vec(maximum(vcat(trial.physical_initial[indices]',trial.observations[:,indices]);dims=1)) == d.scales[indices]
        @test d.initial[indices] == trial.physical_initial[indices] ./ d.scales[indices]
    end
    for trial in filter(r->r.hidden_lactate, all_trials)
        @test ismissing(trial.physical_initial[4]) && ismissing(trial.data.scales[4])
        @test all(ismissing,trial.observations[:,4])
        @test trial.data.initial[4] == 0 # source's fixed hidden-state initial
        @test count(ismissing,trial.data.values) == 15
    end
    ghhp2 = only(filter(r->r.trial==:ghhp2,all_trials))
    @test ghhp2.physical_initial[[1,2,3,5,6,7,8]] ==
        [47.956,51.386,0.286,0,0.048539049,0.004044921,0.000906062]
    @test ghhp2.observations[2,[3,5,6,7,8]] ==
        [5.389,1.536,0.394540199,0.03287835,0.00736475]
    # An independent pre-existing parser checks the three overlapping 2018 trials.
    for (id,old_file) in ((:ghhp1,"camu.R"),(:brwb1,"p_box1.R"),(:brwb2,"p_box2.R"))
        path = joinpath(B.ROOT,"data","raw","cocoa","CocoaModel-master","simulation",old_file)
        old = CocoaSourceData.inspect_dataset(path)
        current = only(filter(r->r.trial==id, all_trials))
        @test current.physical_initial == old.initial
        @test current.observations == old.observations
        @test current.data.times == old.times
    end
    @test B.literals("1e-3, .2, 3.") == [0.001,0.2,3.]
    @test_throws ErrorException B.literals("1, system('unexpected')")
    @test_throws ArgumentError B.literals("1e999")
    @test_throws ErrorException B.literals("NA, 1")
    driver = read(joinpath(B.ROOT,B.RAW,B.ARCHIVE_DIR,"posteriors","post_sims2.R"),String)
    @test_throws ErrorException B.parse_driver(replace(driver,".Dim=c(15,7)"=>".Dim=c(14,7)");hidden_lactate=true)
    @test_throws ErrorException B.parse_driver(replace(driver,"6,12,18,24"=>"6,6,18,24");hidden_lactate=true)
    @test_throws ArgumentError B.benchmark_trials(;scope=:unknown)
end

@testset "2022 benchmark restart provenance and selection" begin
    F = CocoaBenchmarkFits
    B = CocoaBenchmarkData
    available = B.benchmark_trials(;scope=:all)
    selected = F.select_trials(available,[:ghhp2,:brwb3,:ghhp6])
    @test [t.trial for t in selected] == [:ghhp6,:brwb3,:ghhp2]
    @test F.select_trials(available,nothing) === available
    @test_throws ArgumentError F.select_trials(available,Symbol[])
    @test_throws ArgumentError F.select_trials(available,[:ghhp2,:ghhp2])
    @test_throws ArgumentError F.select_trials(available,[:unknown])
    @test_throws ArgumentError F.select_trials(available,["ghhp2"])
    trial = only(filter(t->t.trial==:ghhp2,available))
    historical = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits/mle_pilot")
    @test_throws ArgumentError F.saved_candidate(historical,trial)
    pilot = mle_restart_fixture(mktempdir(),F,trial)
    saved = F.saved_candidate(pilot,trial)
    @test saved.parameters.sigma > 0
    @test saved.provenance["previous_accepted"] == false
    @test length(saved.provenance["files_sha256"]) == 4
    @test saved.provenance["files_sha256"]["ghhp2/report.toml"] ==
        B.digest(joinpath(pilot,"ghhp2/report.toml"))
    names = keys(saved.parameters)
    lower = NamedTuple{names}(ntuple(_->1e-5,25))
    upper = NamedTuple{names}(ntuple(_->5.0,25))
    center = NamedTuple{names}(ntuple(_->0.5,25))
    interior = F.restart_initial(center,lower,upper)
    @test interior.initial == center
    @test isempty(interior.adjusted)
    at_bound = merge(center,(mu1=lower.mu1,mu2=upper.mu2))
    adjusted = F.restart_initial(at_bound,lower,upper)
    @test Set(adjusted.adjusted) == Set((:mu1,:mu2))
    @test all(n->lower[n] < adjusted.initial[n] < upper[n],names)
    @test all(n->adjusted.initial[n] == center[n],setdiff(names,(:mu1,:mu2)))
    restarted = F.restart_initial(saved.parameters,lower,upper)
    @test restarted.adjusted == [:yc6]
    @test 0 < saved.parameters.yc6 - restarted.initial.yc6 < 1e-6
    @test all(n->restarted.initial[n] == saved.parameters[n],setdiff(names,(:yc6,)))
    @test_throws ArgumentError F.restart_initial(merge(saved.parameters,(mu1=10.,)),lower,upper)
    @test_throws ArgumentError F.restart_initial(at_bound,lower,upper;margin=0.)
    # Tamper only disposable copies, preserving historical fit artifacts.
    mktempdir() do temporary
        fixture = joinpath(temporary,"fixture")
        mkpath(joinpath(fixture,"ghhp2"))
        paths = ("settings.toml","ghhp2/input.toml","ghhp2/report.toml","ghhp2/start_001.toml")
        for path in paths
            cp(joinpath(pilot,path),joinpath(fixture,path))
        end
        @test F.saved_candidate(fixture,trial).parameters == saved.parameters
        for (path,mutate!) in (
            ("settings.toml",d->(d["method"]="bayes")),
            ("settings.toml",d->(d["implementation_sha256"]["src/cocoa/baseline.jl"]="changed")),
            ("ghhp2/input.toml",d->(d["source_sha256"]="changed")),
            ("ghhp2/input.toml",d->(d["source_branch"]=99)),
            ("ghhp2/input.toml",d->reverse!(d["state_order"])),
            ("ghhp2/input.toml",d->(d["initial_normalized"][4]=1.)),
            ("ghhp2/input.toml",d->(d["times_h"][1]=1.)),
            ("ghhp2/report.toml",d->(d["trial"]="ghhp3")),
            ("ghhp2/report.toml",d->(d["candidate_start"]=0)),
            ("ghhp2/report.toml",d->delete!(d["candidate_parameters"],"mu1")),
            ("ghhp2/start_001.toml",d->(d["parameters"]["mu1"]=0.1)))
            copy = TOML.parsefile(joinpath(pilot,path))
            mutate!(copy)
            F.write_toml(joinpath(fixture,path),copy)
            @test_throws ArgumentError F.saved_candidate(fixture,trial)
            cp(joinpath(pilot,path),joinpath(fixture,path);force=true)
        end
        output = joinpath(temporary,"invalid")
        @test_throws ArgumentError F.run_benchmark(;output,scope=:all,trials=[:ghhp2],maxiters=0)
        @test !ispath(output)
        @test_throws ArgumentError F.run_benchmark(;output,trials=[:unknown],restart_from=pilot)
        @test !ispath(output)
        # A small real restart checks wiring, provenance, missing observations,
        # and seeds without asking the optimizer to converge.
        output = joinpath(temporary,"restart")
        rows = F.run_benchmark(;output,scope=:all,trials=[:ghhp2],restart_from=pilot,
            maxiters=1,outer_iterations=1,maxtime=5.,seed=1234)
        config = TOML.parsefile(joinpath(output,"settings.toml"))
        input = TOML.parsefile(joinpath(output,"ghhp2/input.toml"))
        first_start = TOML.parsefile(joinpath(output,"ghhp2/start_001.toml"))
        report = TOML.parsefile(joinpath(output,"ghhp2/report.toml"))
        @test length(rows) == 1 && first(rows)[1] == :ghhp2
        @test config["trials"] == ["ghhp2"]
        @test input["seed"] == config["trial_seeds"]["ghhp2"] ==
            1234 + findfirst(t->t.trial==:ghhp2,available) - 1
        @test input["restart"]["files_sha256"] == saved.provenance["files_sha256"]
        @test input["restart"]["adjusted_parameters"] == ["yc6"]
        @test first_start["initial"] == input["optimizer_initial"] ==
            Dict(string(k)=>v for (k,v) in pairs(restarted.initial))
        @test report["accepted"] == first_start["accepted"]
        @test report["direct_metric_comparison_available"] == false
        @test count(line->occursin(",L,,",line),
            readlines(joinpath(output,"ghhp2/predictions.csv"))) == length(trial.data.times)
    end
end

@testset "2022 benchmark mechanisms and published references" begin
    B = CocoaBenchmarkData
    variants = B.mechanism_variants()
    @test length(variants) == length(unique(r.code for r in variants)) == 32
    @test first(variants).kinetic_parameters == 24
    @test last(variants).kinetic_parameters == 43
    @test only(filter(r->r.code=="23",variants)).fitted_parameters == 34
    reference = B.published_baseline()
    @test count(r->r.success, reference) == 18
    @test only(filter(r->r.trial==:ecpt2,reference)).min_tail_ess == 2000
    @test only(filter(r->r.trial==:brwb7,reference)).mean_ess_fraction == .514
    @test only(filter(r->r.trial==:mywb3,reference)).mean_ess_fraction == .501
    @test only(filter(r->r.trial==:brwb4,reference)).min_tail_ess == 288
    @test ismissing(only(filter(r->r.trial==:ghhp2,reference)).bma_weight)
    mktempdir() do output
        prepared = B.prepare_benchmark(;output)
        @test length(prepared.trials) == 32
        @test length(readlines(joinpath(output,"datasets.csv"))) == 33
        ghhp2 = readlines(joinpath(output,"ghhp2.csv"))
        @test count(line->occursin(",L,,unobserved,",line),ghhp2) == 15
        @test any(line->occursin(",L,,unobserved_initial,mg/g pulp,,0.0",line),ghhp2)
        original = B.digest(joinpath(output,"ghhp2.csv"))
        B.prepare_benchmark(;output)
        @test B.digest(joinpath(output,"ghhp2.csv")) == original
        @test_throws ArgumentError CocoaBenchmarkFits.run_benchmark(;output)
    end
end

@testset "2022 hidden-lactate likelihood" begin
    B = CocoaBenchmarkData
    trial = only(filter(r->r.trial==:ghhp2,B.benchmark_trials()))
    initial = load_params(joinpath(B.ROOT,"params","cocoa_camu_scaled.toml")).values
    lower = NamedTuple{keys(initial)}(ntuple(_->1e-5,25))
    upper = NamedTuple{keys(initial)}(ntuple(_->5.,25))
    model = cocoa_fit_model(;initial,lower,upper,abstol=1e-9,reltol=1e-8)
    score = CocoaBenchmarkFits.training_scores(model,trial.data,initial)
    @test score.n == 105
    @test all(isfinite,score.prediction) && all(>=(0),score.prediction)
    @test score.nll ≈ 105log(initial.sigma) + 105log(2pi)/2 +
        sum(abs2,score.prediction[score.indices]-trial.data.values[score.indices])/(2initial.sigma^2)
    @test score.normalized_rmse ≈ sqrt(sum(abs2,score.prediction[score.indices]-
        trial.data.values[score.indices])/105)
    hidden_indices = findall(ismissing,trial.data.values)
    @test length(hidden_indices) == 15
    hidden_nan = FitModel((p,d)->begin
        values = copy(score.prediction)
        values[hidden_indices] .= NaN
        values
    end;initial,lower,upper)
    @test gaussian_nll(hidden_nan,trial.data,initial) == score.nll
end
