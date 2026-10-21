using CairoMakie
using TOML

"""Plot all retained parameter traces from a completed Bayesian cocoa report."""
function cocoa_trace_figures(directory)
    settings = TOML.parsefile(joinpath(directory,"settings.toml"))
    report = TOML.parsefile(joinpath(directory,"report.toml"))
    settings["method"] == "bayes" || error("Requires a Bayesian fit directory")
    lines = readlines(joinpath(directory,"draws.csv"))
    header = split(first(lines),',')
    header[1:2] == ["chain","draw"] || error("Unexpected draw columns")
    parameters = header[3:end]
    rows = split.(lines[2:end],',')
    all(length(row) == length(header) for row in rows) || error("Incomplete draw row")
    chain_ids = parse.(Int,getindex.(rows,1))
    draw_ids = parse.(Int,getindex.(rows,2))
    nchains, nsamples = settings["chains"], settings["samples"]
    length(rows) == nchains*nsamples || error("Draw count differs from fit settings")
    Set(chain_ids) == Set(1:nchains) || error("Unexpected chain identifiers")
    indices = [findall(==(chain),chain_ids) for chain in 1:nchains]
    all(draw_ids[index] == collect(1:nsamples) for index in indices) || error("Incomplete or unordered chain")
    comparisons = split.(readlines(joinpath(directory,"parameter_comparison.csv"))[2:end],',')
    diagnostics = Dict(row[1] => parse.(Float64,row[7:9]) for row in comparisons)
    Set(parameters) == Set(keys(diagnostics)) || error("Parameter columns differ from report")
    colors = Makie.wong_colors()
    ncols = 5
    nrows = cld(length(parameters),ncols)
    figure = Figure(size=(1800,330nrows+150),fontsize=14,figure_padding=30)
    status = report["converged"] ? "convergence diagnostics passed" : "CONVERGENCE CHECKS FAILED"
    Label(figure[1,1:ncols],"Cocoa $(settings["trial"]): retained parameter traces — $status";
        fontsize=24,font=:bold,color=report["converged"] ? :black : :firebrick)
    for (i,parameter) in enumerate(parameters)
        values = parse.(Float64,getindex.(rows,i+2))
        all(isfinite,values) || error("Nonfinite parameter draw")
        rhat,bulk,tail = diagnostics[parameter]
        passed = all(isfinite,(rhat,bulk,tail)) && rhat < report["rhat_threshold"] &&
            min(bulk,tail) >= report["minimum_ess"]
        ess_label(value) = isfinite(value) ? string(round(Int,value)) : "undefined"
        title = "$parameter | R-hat $(round(rhat;digits=5))\nESS bulk $(ess_label(bulk)), tail $(ess_label(tail))"
        axis = Axis(figure[(i-1)÷ncols+2,(i-1)%ncols+1];title,titlecolor=passed ? :black : :firebrick,
            titlesize=16,xlabel="Retained draw",ylabel="Normalized coefficient")
        for chain in 1:nchains
            lines!(axis,draw_ids[indices[chain]],values[indices[chain]];
                color=(colors[mod1(chain,length(colors))],0.75),linewidth=0.8)
        end
    end
    Legend(figure[nrows+2,1:ncols],
        [LineElement(color=colors[mod1(chain,length(colors))],linewidth=2) for chain in 1:nchains],
        ["Chain $chain" for chain in 1:nchains];orientation=:horizontal)
    Label(figure[nrows+3,1:ncols],
        "All retained draws shown; warmup excluded. Red titles fail R-hat < $(report["rhat_threshold"]) or ESS ≥ $(report["minimum_ess"]). Recorded divergences: $(report["divergences"]); depth hits: $(report["depth_hits"]).\nNormalized-coordinate μ and k coefficients have units h⁻¹; ks, yc, and σ are dimensionless. Source: doi:10.1098/rsos.180964.";
        fontsize=13,halign=:left)
    save(joinpath(directory,"traces.png"),figure)
    save(joinpath(directory,"traces.svg"),figure)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia --project=scripts scripts/cocoa_trace_figures.jl <fit-directory>")
    cocoa_trace_figures(only(ARGS))
end
