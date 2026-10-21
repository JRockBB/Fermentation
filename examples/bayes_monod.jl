# Synthetic P1.6 example: infer growth rate while holding Ks and yields known.
using FermentModels
using OrdinaryDiffEq
using Distributions: LogNormal
using Random
using Statistics
using TOML

function bayes_monod(output=joinpath(@__DIR__, "output", "monod_bayes"))
    spec = ModelSpec(; organisms=(:microbe,), substrates=(:sugar,), products=(:acid,),
        stoichiometry=reshape([1.,-2.,1.],3,1))
    function predict(p, data)
        problem = fermentation_problem(spec, [0.2,4.,0.] .* one(p.mu), (0.,5.),
            (mu_max=[p.mu], Ks=reshape([0.7],1,1)))
        solution = solve(problem, Tsit5(); abstol=1e-10, reltol=1e-9)
        successful_retcode(solution) || return nothing
        return [solution(t)[1] for t in data.times]
    end
    truth, sigma = (mu=0.8,), 0.03
    times = collect(0.25:0.25:5.)
    data = (; times, values=predict(truth,(; times)) + sigma .* randn(MersenneTwister(122),length(times)), sigma)
    model = FitModel(predict; initial=(mu=0.5,), lower=(mu=0.1,), upper=(mu=2.,))
    future_times = collect(0.125:0.25:4.875)
    future = (times=future_times, values=fill(missing,length(future_times)), sigma=sigma)
    posterior = fit_bayes(model, data, (mu=LogNormal(log(0.7),0.3),);
        samples=1000, warmup=700, chains=2, predictive_draws=1000,
        predictive_data=future, rng=MersenneTwister(123))
    heldout = predict(truth,future) + sigma .* randn(MersenneTwister(124),length(future_times))
    bands = posterior.predictive
    coverage = mean(bands.observed.lower .<= heldout .<= bands.observed.upper)
    mkpath(output)
    open(joinpath(output,"observations.csv"),"w") do io
        println(io,"time_h,synthetic_biomass_equiv_mmol_L,sigma_mmol_L")
        for i in eachindex(times)
            println(io,join((times[i],data.values[i],sigma),','))
        end
    end
    open(joinpath(output,"predictive.csv"),"w") do io
        println(io,"time_h,heldout_mmol_L,latent_lower,latent_median,latent_upper,observed_lower,observed_median,observed_upper")
        for i in eachindex(future_times)
            println(io,join((future_times[i],heldout[i],bands.latent.lower[i],bands.latent.median[i],
                bands.latent.upper[i],bands.observed.lower[i],bands.observed.median[i],bands.observed.upper[i]),','))
        end
    end
    open(joinpath(output,"draws.csv"),"w") do io
        println(io,"sample,chain,mu_per_h")
        for chain in axes(posterior.draws.mu,2), sample in axes(posterior.draws.mu,1)
            println(io,join((sample,chain,posterior.draws.mu[sample,chain]),','))
        end
    end
    report = Dict("synthetic"=>true, "script"=>"examples/bayes_monod.jl",
        "dataset"=>"observations.csv", "noise_seed"=>122, "sampling_seed"=>123, "heldout_seed"=>124,
        "mu_truth_per_h"=>truth.mu, "prior_median"=>0.7, "prior_log_sd"=>0.3,
        "fixed_Ks_mmol_L"=>0.7, "heldout_coverage"=>coverage,
        "summary"=>Dict(string(k)=>v for (k,v) in pairs(posterior.summary.mu)),
        "diagnostics"=>Dict(string(k)=>v for (k,v) in pairs(posterior.diagnostics)),
        "settings"=>Dict(string(k)=>(v isa Symbol ? String(v) : v) for (k,v) in pairs(posterior.settings)))
    open(joinpath(output,"report.toml"),"w") do io
        TOML.print(io,report;sorted=true)
    end
    println("Synthetic posterior for mu (1/h): ", posterior.summary.mu)
    println("Diagnostics passed: ", posterior.diagnostics.passed)
    println("Held-out pointwise coverage (20 observations): ", coverage)
    println("Results: ", abspath(output))
    return posterior
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? bayes_monod() : bayes_monod(only(ARGS))
end
