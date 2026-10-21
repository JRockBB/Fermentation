using FermentModels
using OrdinaryDiffEq
using TOML
include(joinpath(@__DIR__, "..", "scripts", "prepare_cocoa.jl"))
using .CocoaSourceData: cocoa_trial, TRIALS

# A plug-in trajectory at published, rounded posterior means. This is not a
# posterior predictive distribution; the source archive omits the joint draws.
function run_cocoa_baseline(; output=joinpath(@__DIR__, "output", "cocoa_baseline"))
    mkpath(output)
    root = normpath(joinpath(@__DIR__, ".."))
    report = Dict{String,Any}("interpretation" => "Plug-in simulation at normalized posterior means in supplement Tables S4-S6; no uncertainty bands inferred from marginal summaries")
    for trial in TRIALS
        data = cocoa_trial(trial)
        parameters = load_params(joinpath(root, "params", "cocoa_$(trial)_scaled.toml")).values
        initial = data.initial ./ data.scales
        problem = cocoa_problem(initial, (0.0, last(data.times)), parameters)
        solution = solve(problem, Tsit5(); abstol=1e-11, reltol=1e-10)
        successful_retcode(solution) || error("Cocoa solve failed for $trial: $(solution.retcode)")
        times = collect(0.0:0.5:last(data.times))
        open(joinpath(output, "$trial.csv"), "w") do io
            println(io, "time_h,", join(String.(COCOA_STATES), ','))
            for t in times
                println(io, t, ',', join(solution(t) .* data.scales, ','))
            end
        end
        residuals = reduce(hcat, solution.(data.times))' .- data.observations ./ data.scales'
        report[String(trial)] = Dict("source" => "doi:10.6084/m9.figshare.7152614.v2",
            "normalized_rmse" => sqrt(sum(abs2, residuals)/length(residuals)),
            "sigma_published_mean" => parameters.sigma,
            "minimum_saved_normalized_state" => minimum(reduce(hcat, solution.(times))),
            "final_physical_state" => solution(last(times)) .* data.scales)
        println(trial, ": normalized RMSE = ", report[String(trial)]["normalized_rmse"])
    end
    open(joinpath(output, "report.toml"), "w") do io
        TOML.print(io, report; sorted=true)
    end
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cocoa_baseline()
end
