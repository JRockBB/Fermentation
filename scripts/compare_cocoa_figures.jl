using FermentModels
using OrdinaryDiffEq
using TOML
using SHA
include("prepare_cocoa.jl")
using .CocoaSourceData: cocoa_trial, TRIALS

"""Report the published-curve discrepancy without changing source parameters or tolerances."""
function compare_cocoa_figures(; root=normpath(joinpath(@__DIR__, "..")))
    path = joinpath(root, "data", "raw", "digitized", "moreno2018-figures", "points.csv")
    points = split.(readlines(path)[2:end], ',')
    destination = joinpath(root, "data", "processed", "cocoa", "moreno2018")
    mkpath(destination)
    tolerance = 0.02 # Preselected fraction of the calibrated vertical range.
    report = Dict{String,Any}("tolerance_fraction" => tolerance,
        "reference" => "doi:10.1098/rsos.180964 Figures 3-5; thesis Figures 2.3-2.5",
        "digitized_input" => "data/raw/digitized/moreno2018-figures/points.csv",
        "digitized_sha256" => bytes2hex(open(sha256,path)),
        "script" => "scripts/compare_cocoa_figures.jl",
        "prediction" => "ODE evaluated at printed normalized posterior mean parameters",
        "all_within_tolerance" => true)
    open(joinpath(destination, "figure_comparison.csv"), "w") do io
        println(io, "experiment,variable,time_h,published_curve,simulation,axis_span,absolute_fraction_error,within_tolerance")
        for trial in TRIALS
            data = cocoa_trial(trial; root)
            parameters = load_params(joinpath(root, "params", "cocoa_$(trial)_scaled.toml")).values
            solution = solve(cocoa_problem(data.initial ./ data.scales, (0., 144.), parameters),
                Tsit5(); abstol=1e-11, reltol=1e-10)
            successful_retcode(solution) || error("Failed cocoa comparison solve")
            errors = Float64[]
            details = Dict{String,Any}()
            for (i, state) in enumerate(COCOA_STATES)
                state_errors = Float64[]
                for row in filter(r -> r[1] == String(trial) && r[2] == String(state), points)
                    time, reference, span = parse.(Float64, row[[3,4,8]])
                    prediction = solution(time)[i] * data.scales[i]
                    error_fraction = abs(prediction-reference)/span
                    push!(state_errors, error_fraction)
                    println(io, join((trial, state, time, reference, prediction, span,
                        error_fraction, error_fraction <= tolerance), ','))
                end
                append!(errors, state_errors)
                details[String(state)] = Dict("max_fraction_error" => maximum(state_errors),
                    "points" => length(state_errors), "outside_tolerance" => count(>(tolerance), state_errors))
            end
            details["max_fraction_error"] = maximum(errors)
            details["outside_tolerance"] = count(>(tolerance), errors)
            report[String(trial)] = details
            report["all_within_tolerance"] &= all(<=(tolerance), errors)
        end
    end
    open(joinpath(destination, "figure_comparison.toml"), "w") do io
        TOML.print(io, report; sorted=true)
    end
    println("Published-curve comparison within 2%: ", report["all_within_tolerance"])
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_cocoa_figures()
end
