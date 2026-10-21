# Synthetic recovery demonstration for P1.5; no measured food data or parameters.
using FermentModels
using OrdinaryDiffEq
using Random
using TOML

function fit_monod(output=joinpath(@__DIR__, "output", "monod_fit"); starts=10)
    spec = ModelSpec(; organisms=(:microbe,), substrates=(:sugar,), products=(:acid,),
        stoichiometry=reshape([1.,-2.,1.],3,1), conserved_weights=ones(3))
    function predict(parameters, data)
        initial = [0.2,4.,0.] .* one(parameters.mu + parameters.Ks)
        problem = fermentation_problem(spec, initial, (0.,8.),
            (mu_max=[parameters.mu], Ks=reshape([parameters.Ks],1,1)))
        trajectory = solve(problem, Tsit5(); abstol=1e-10, reltol=1e-9)
        successful_retcode(trajectory) || return nothing
        return [trajectory(t)[1] for t in data.times]
    end
    truth = (mu=0.8, Ks=0.7)
    times = collect(0.25:0.25:7.)
    clean = predict(truth, (; times))
    sigma = 0.02
    data = (; times, values=clean + sigma .* randn(MersenneTwister(314), length(times)), sigma)
    model = FitModel(predict; initial=(mu=0.5,Ks=1.), lower=(mu=0.1,Ks=0.05), upper=(mu=2.,Ks=5.))
    result = fit_mle(model, data; starts, rng=MersenneTwister(12))
    result.converged || error("No start converged; inspect the model and optimizer diagnostics")

    mkpath(output)
    open(joinpath(output, "observations.csv"), "w") do io
        println(io, "time_h,synthetic_biomass_equiv_mmol_L,sigma_mmol_L,truth_mmol_L,predicted_mmol_L,residual_mmol_L")
        for i in eachindex(times)
            println(io, join((times[i],data.values[i],sigma,clean[i],result.predictions[i],result.residuals[i]), ','))
        end
    end
    open(joinpath(output, "starts.csv"), "w") do io
        println(io, "start,initial_mu,initial_Ks,final_mu,final_Ks,nll,converged,retcode")
        for (i,run) in enumerate(result.runs)
            final = isnothing(run.parameters) ? (mu=missing,Ks=missing) : run.parameters
            println(io, join((i,run.initial.mu,run.initial.Ks,final.mu,final.Ks,run.objective,run.converged,run.retcode), ','))
        end
    end
    report = Dict{String,Any}(
        "synthetic" => true, "script" => "examples/fit_monod.jl",
        "dataset" => "observations.csv", "noise_seed" => 314, "starts_seed" => 12,
        "starts" => starts, "likelihood" => "independent Gaussian, known sigma",
        "uncertainty_status" => string(result.uncertainty.status), "negative_log_likelihood" => result.objective)
    for (name, unit) in ((:mu,"1/h"), (:Ks,"mmol/L"))
        parameter = Dict{String,Any}("estimate" => result.parameters[name], "truth" => truth[name], "units" => unit)
        if !isnothing(result.uncertainty.intervals)
            parameter["lower_95"] = result.uncertainty.intervals[name].lower
            parameter["upper_95"] = result.uncertainty.intervals[name].upper
        end
        report[string(name)] = parameter
    end
    open(joinpath(output, "report.toml"), "w") do io
        TOML.print(io, report; sorted=true)
    end
    println("Synthetic Monod fit: ", result.parameters)
    println("Local 95% interval status: ", result.uncertainty.status)
    println("Converged starts: ", count(r -> r.converged, result.runs), "/", length(result.runs))
    println("Results: ", abspath(output))
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? fit_monod() : fit_monod(only(ARGS))
end
