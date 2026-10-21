# Disposable software fixture with synthetic coefficients and freshly evaluated
# scores. Never rewrite a historical fit's implementation hashes to make it load.
function mle_restart_fixture(directory,F,trial)
    B = F.CocoaBenchmarkData
    names = (COCOA_PARAMETERS...,:sigma)
    initial = NamedTuple{names}(ntuple(_->.5,25))
    parameters = merge(initial,(yc6=5.,))
    lower = NamedTuple{names}(ntuple(_->1e-5,25))
    upper = NamedTuple{names}(ntuple(_->5.,25))
    model = cocoa_fit_model(;initial,lower,upper,solver=:rodas5p,
        abstol=1e-9,reltol=1e-8,maxiters=10_000)
    score = F.training_scores(model,trial.data,parameters)
    isnothing(score) && error("Synthetic MLE fixture failed to simulate")
    label = string(trial.trial)
    path = joinpath(directory,label)
    mkpath(path)
    F.write_toml(joinpath(directory,"settings.toml"),Dict(
        "method"=>"mle","scope"=>"all","trials"=>[label],
        "solver"=>"rodas5p","abstol"=>1e-9,"reltol"=>1e-8,
        "solver_maxiters"=>10_000,"mle_bounds"=>[1e-5,5.],
        "fixture"=>"Synthetic software test; no optimization or biological estimate",
        "implementation_sha256"=>Dict(f=>B.digest(joinpath(B.ROOT,f))
            for f in F.LIKELIHOOD_IMPLEMENTATION)))
    F.write_toml(joinpath(path,"input.toml"),Dict(
        "trial"=>label,"source_file"=>trial.source_file,"source_sha256"=>trial.source_sha256,
        "source_branch"=>trial.branch,"hidden_lactate"=>trial.hidden_lactate,
        "initial_normalized"=>trial.data.initial,"times_h"=>trial.data.times,
        "state_order"=>String.(collect(COCOA_STATES))))
    point = Dict(string(k)=>v for (k,v) in pairs(parameters))
    F.write_toml(joinpath(path,"report.toml"),Dict(
        "trial"=>label,"status"=>"synthetic_fixture","accepted"=>false,
        "candidate_start"=>1,"candidate_parameters"=>point,
        "candidate_nll"=>score.nll,"candidate_normalized_rmse"=>score.normalized_rmse))
    F.write_toml(joinpath(path,"start_001.toml"),Dict(
        "parameters"=>point,"accepted"=>false,"status"=>"synthetic_fixture"))
    return directory
end
