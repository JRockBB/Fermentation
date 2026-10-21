include(joinpath(pwd(),".julia-depot/p25_brwb3_depth_refinement.jl"))
@testset "Depth refinement comparison preflight" begin
    directory = joinpath(B.ROOT,"data/processed/cocoa_benchmark/fits/bayes_brwb3_depth_pilot1/fit")
    trial = only(filter(t->t.trial==:brwb3,B.benchmark_trials(;scope=:all)))
    report = TOML.parsefile(joinpath(directory,"brwb3/report.toml"))
    row = summarize_run("preflight",directory,trial)
    @test row.total_draws==400 && row.draws_per_chain==100 && row.warmup_per_chain==250 && row.max_depth==9
    @test row.max_rhat==report["max_rhat"] && row.min_ess_bulk==report["min_ess_bulk"] && row.min_ess_tail==report["min_ess_tail"]
    @test row.retained_leapfrog_steps==152048 && row.depth_hits==133 && row.divergences==0
    @test row.min_bulk_ess_per_draw==row.min_ess_bulk/400 && row.min_tail_ess_per_draw==row.min_ess_tail/400
    saved = D.saved_draws(directory,trial)
    block = summarize_block("preflight",saved,1:100)
    @test block.max_rhat==row.max_rhat && block.min_ess_bulk==row.min_ess_bulk && block.min_ess_tail==row.min_ess_tail
    @test block.first_draw==1 && block.last_draw==100 && block.total_draws==400 && block.depth_hits==row.depth_hits
end
