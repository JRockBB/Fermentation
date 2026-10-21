using TOML
using SHA

"""
Compare saved joint-posterior predictions with digitized published central curves.
No refitting or marginal-parameter resampling occurs. Both posterior mean and
median comparisons are reported because the archived plotting script is absent.
Convergence is required before either comparison can be accepted.
"""
function compare_cocoa_fit(directory; root=normpath(joinpath(@__DIR__, "..")))
    settings = TOML.parsefile(joinpath(directory,"settings.toml"))
    fit = TOML.parsefile(joinpath(directory,"report.toml"))
    settings["method"] == "bayes" || throw(ArgumentError("requires a Bayesian fit"))
    trial = settings["trial"]
    predictions = split.(readlines(joinpath(directory,"predictive.csv"))[2:end], ',')
    reference_path = joinpath(root,"data","raw","digitized","moreno2018-figures","points.csv")
    reference = split.(readlines(reference_path)[2:end], ',')
    reference = filter(row -> row[1] == trial, reference)
    length(reference) == 56 || error("Expected 56 digitized points for $trial")
    tolerance = 0.02 # Same preselected axis-relative criterion as the baseline.
    mean_errors, median_errors = Float64[], Float64[]
    open(joinpath(directory,"curve_comparison.csv"),"w") do io
        println(io,"variable,time_h,published_curve,posterior_mean,posterior_median,axis_span,mean_fraction_error,median_fraction_error")
        for row in reference
            t, value, span = parse.(Float64,row[[3,4,8]])
            matches = filter(p -> p[2] == row[2] && parse(Float64,p[1]) == t, predictions)
            length(matches) == 1 || error("Expected one posterior prediction for $(row[2]) at $t")
            mean_value, median_value = parse.(Float64,only(matches)[[4,5]])
            mean_error, median_error = abs(mean_value-value)/span, abs(median_value-value)/span
            push!(mean_errors,mean_error)
            push!(median_errors,median_error)
            println(io,join((row[2],t,value,mean_value,median_value,span,mean_error,median_error),','))
        end
    end
    report = Dict("trial"=>trial, "tolerance_fraction"=>tolerance,
        "fit_converged"=>fit["converged"], "points"=>length(reference),
        "mean_max_fraction_error"=>maximum(mean_errors),
        "median_max_fraction_error"=>maximum(median_errors),
        "mean_outside_tolerance"=>count(>(tolerance),mean_errors),
        "median_outside_tolerance"=>count(>(tolerance),median_errors),
        "mean_comparison_accepted"=>fit["converged"] && all(<=(tolerance),mean_errors),
        "median_comparison_accepted"=>fit["converged"] && all(<=(tolerance),median_errors),
        "source_summary_identified"=>false,
        "reference_sha256"=>bytes2hex(open(sha256,reference_path)),
        "predictive_sha256"=>bytes2hex(open(sha256,joinpath(directory,"predictive.csv"))),
        "fit_report_sha256"=>bytes2hex(open(sha256,joinpath(directory,"report.toml"))),
        "reference"=>"doi:10.1098/rsos.180964 Figures 3-5",
        "script"=>"scripts/compare_cocoa_fit.jl")
    open(joinpath(directory,"curve_comparison.toml"),"w") do io
        TOML.print(io,report;sorted=true)
    end
    println("$trial: maximum axis-relative mean error ", maximum(mean_errors),
        "; median error ",maximum(median_errors),"; fit converged: ",fit["converged"])
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia --project scripts/compare_cocoa_fit.jl <fit-directory>")
    compare_cocoa_fit(only(ARGS))
end
