using FermentModels
using OrdinaryDiffEq
using TOML
using SHA
include("prepare_cocoa.jl")
using .CocoaSourceData: cocoa_trial, published_summaries, TRIALS

"""
Compare trajectories at the printed marginal means and medians with the same
digitized central curves. Neither is a joint-posterior predictive summary.
This diagnostic does not identify the original plotting method or refit data.
"""
function compare_cocoa_summaries(; root=normpath(joinpath(@__DIR__, "..")))
    reference_path = joinpath(root,"data","raw","digitized","moreno2018-figures","points.csv")
    summary_path = joinpath(root,"data","raw","transcribed","cocoa2018","posterior_summaries.csv")
    points = split.(readlines(reference_path)[2:end],',')
    summaries = published_summaries(;root)
    destination = joinpath(root,"data","processed","cocoa","moreno2018")
    mkpath(destination)
    tolerance = 0.02 # Unchanged preselected fraction of the panel's vertical range.
    report = Dict{String,Any}(
        "script" => "scripts/compare_cocoa_summaries.jl",
        "script_sha256" => bytes2hex(open(sha256,@__FILE__)),
        "reference" => "doi:10.1098/rsos.180964 Figures 3-5; supplement Tables S4-S6",
        "reference_sha256" => bytes2hex(open(sha256,reference_path)),
        "parameter_summary_sha256" => bytes2hex(open(sha256,summary_path)),
        "prediction" => "ODE evaluated at printed marginal parameter summaries, not joint-posterior predictions",
        "source_summary_identified" => false,
        "tolerance_fraction" => tolerance, "solver" => "Tsit5",
        "state_units" => "mg/g pulp", "time_units" => "h",
        "abstol" => 1e-11, "reltol" => 1e-10, "julia_version" => string(VERSION),
        "implementation_sha256" => Dict(path => bytes2hex(open(sha256,joinpath(root,path)))
            for path in ("src/cocoa/baseline.jl","scripts/prepare_cocoa.jl")),
        "results" => Dict{String,Any}[])
    open(joinpath(destination,"point_summary_comparison.csv"),"w") do io
        println(io,"experiment,parameter_summary,variable,time_h,published_curve,simulation,axis_span,absolute_fraction_error,within_tolerance")
        for trial in TRIALS
            data = cocoa_trial(trial;root)
            reference = filter(row -> row[1] == String(trial),points)
            length(reference) == 56 || error("Expected 56 digitized points for $trial")
            for summary in (:mean,:q50)
                parameters = NamedTuple{COCOA_PARAMETERS}(Tuple(
                    getproperty(only(filter(row -> row.experiment == trial && row.parameter == name,summaries)),summary)
                    for name in COCOA_PARAMETERS))
                solution = solve(cocoa_problem(data.initial ./ data.scales,(0.,144.),parameters),
                    Tsit5();abstol=1e-11,reltol=1e-10)
                successful_retcode(solution) || error("Failed $trial/$summary comparison solve")
                errors = map(reference) do row
                    index = findfirst(==(Symbol(row[2])),COCOA_STATES)
                    time,observed,span = parse.(Float64,row[[3,4,8]])
                    prediction = solution(time)[index]*data.scales[index]
                    fraction = abs(prediction-observed)/span
                    println(io,join((trial,summary,row[2],time,observed,prediction,span,
                        fraction,fraction <= tolerance),','))
                    fraction
                end
                result = Dict("trial"=>String(trial),"parameter_summary"=>String(summary),
                    "max_fraction_error"=>maximum(errors),"points"=>length(errors),
                    "outside_tolerance"=>count(>(tolerance),errors),
                    "all_within_tolerance"=>all(<=(tolerance),errors),
                    "raw_file"=>replace(data.source_file,'\\'=>'/'),
                    "raw_sha256"=>bytes2hex(open(sha256,joinpath(root,data.source_file))))
                push!(report["results"],result)
                println(trial," ",summary,": maximum=",maximum(errors),
                    "; outside 2%=",count(>(tolerance),errors))
            end
        end
    end
    open(joinpath(destination,"point_summary_comparison.toml"),"w") do io
        TOML.print(io,report;sorted=true)
    end
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_cocoa_summaries()
end
