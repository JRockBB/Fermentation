include("../scripts/diagnose_cocoa_benchmark.jl")
isdefined(@__MODULE__,:mle_restart_fixture) || include("cocoa_benchmark_fixture.jl")

@testset "2022 saved benchmark diagnostic reports" begin
    D = CocoaBenchmarkDiagnostics
    B = D.B
    trial = only(filter(t->t.trial==:ghhp2,B.benchmark_trials(;scope=:all)))
    mktempdir() do temporary
        pilot = mle_restart_fixture(joinpath(temporary,"fixture"),D.CocoaBenchmarkFits,trial)
        output = joinpath(temporary,"diagnostics")
        historical = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits/mle_pilot")
        @test_throws ArgumentError D.diagnose_benchmark(historical;output,trials=[:ghhp2])
        @test !ispath(output)
        for factors in ((),(1.,1.),(1.1,),(-1.,),(NaN,))
            @test_throws ArgumentError D.diagnose_benchmark(pilot;output,tolerance_factors=factors)
            @test !ispath(output)
        end
        @test_throws ArgumentError D.diagnose_benchmark(pilot;output,trials=[:unknown])
        @test !ispath(output)
        rows = D.diagnose_benchmark(pilot;output,trials=[:ghhp2],tolerance_factors=(1.,))
        settings = TOML.parsefile(joinpath(output,"settings.toml"))
        report = TOML.parsefile(joinpath(output,"ghhp2/report_1.toml"))
        candidate = TOML.parsefile(joinpath(pilot,"ghhp2/report.toml"))
        parameter_rows = [split(line,',') for line in
            readlines(joinpath(output,"ghhp2/parameters_1.csv"))[2:end]]
        @test length(rows) == 1 && length(only(rows)) == 13
        @test length(readlines(joinpath(output,"comparison.csv"))) == 2
        @test report["solve_succeeded"] && !report["diagnostics_accept_fit"]
        @test settings["trials"] == ["ghhp2"] && !settings["diagnostics_accept_fit"]
        @test report["observation_count"] == 105
        @test report["nll"] ≈ candidate["candidate_nll"]
        @test report["conditional_sigma_mle"] ≈ candidate["candidate_normalized_rmse"]
        @test length(parameter_rows) == report["information_parameter_count"] == 25
        @test maximum(abs(parse(Float64,r[4])) for r in parameter_rows) ==
            report["projected_gradient_inf"]
        @test length(report["finite_difference_scaled_errors"]) == 2
        @test settings["fit_settings_sha256"] == B.digest(joinpath(pilot,"settings.toml"))
        provenance = TOML.parsefile(joinpath(output,"ghhp2/provenance.toml"))
        @test provenance["files_sha256"]["ghhp2/report.toml"] == B.digest(joinpath(pilot,"ghhp2/report.toml"))
        @test_throws ArgumentError D.diagnose_benchmark(pilot;output)
    end
end
