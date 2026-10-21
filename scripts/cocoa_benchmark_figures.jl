using CairoMakie
using SHA
using TOML

function benchmark_csv(path)
    lines = readlines(path)
    header = split(first(lines),',')
    rows = [split(line,',';keepempty=true) for line in lines[2:end]]
    all(length(row)==length(header) for row in rows) || error("Incomplete CSV row in $path")
    return header,rows
end

"""Render saved benchmark predictions and traces with their diagnostic status visible."""
function cocoa_benchmark_figures(directory)
    directory = normpath(abspath(directory))
    settings = TOML.parsefile(joinpath(dirname(directory),"settings.toml"))
    input = TOML.parsefile(joinpath(directory,"input.toml"))
    report = TOML.parsefile(joinpath(directory,"report.toml"))
    settings["method"] == "bayes" || error("Requires a Bayesian benchmark run")
    report["status"] in ("posterior_diagnostics_passed","posterior_diagnostics_failed") ||
        error("Requires a completed trial report")
    trial = input["trial"]
    status = report["diagnostics_passed"] ? "sampling diagnostics passed" : "SAMPLING DIAGNOSTICS FAILED"
    title_color = report["diagnostics_passed"] ? :black : :firebrick
    header,rows = benchmark_csv(joinpath(directory,"predictive.csv"))
    columns = Dict(name=>index for (index,name) in enumerate(header))
    states = ["G","F","E","L","A","Xy","Xl","Xa"]
    input["state_order"] == states || error("Unexpected source state order")
    labels = ["Glucose","Fructose","Ethanol","Lactic acid","Acetic acid","Yeast","LAB","AAB"]
    figure = Figure(size=(1500,1050),fontsize=16,figure_padding=30)
    Label(figure[1,1:4],"Cocoa benchmark $trial — $status";
        fontsize=23,font=:bold,color=title_color)
    for (i,state) in enumerate(states)
        data = filter(row->row[columns["variable"]]==state,rows)
        values = name->parse.(Float64,getindex.(data,columns[name]))
        times = values("time_h")
        times == input["times_h"] || error("Prediction times differ from source input")
        hidden = state == "L" && input["hidden_lactate"]
        suffix = hidden ? "normalized" : "mg_g"
        axis = Axis(figure[(i-1)÷4+2,(i-1)%4+1];
            title=hidden ? "Lactate (unmeasured)" : labels[i],
            xlabel="Time (h)",ylabel=hidden ? "Normalized latent state" : "mg/g pulp")
        if !hidden
            band!(axis,times,values("observation_lower_$suffix"),values("observation_upper_$suffix");
                color=(:gray,0.25))
        end
        band!(axis,times,values("latent_lower_$suffix"),values("latent_upper_$suffix");
            color=(:firebrick,0.23))
        lines!(axis,times,values("latent_mean_$suffix");color=:firebrick,linewidth=2.5)
        !hidden && scatter!(axis,times,values("observed_$suffix");color=:black,markersize=7)
    end
    Legend(figure[4,1:4],
        [LineElement(color=:firebrick,linewidth=2.5),PolyElement(color=(:firebrick,0.23)),
         PolyElement(color=(:gray,0.25)),MarkerElement(color=:black,marker=:circle)],
        ["Joint-posterior mean trajectory","Latent pointwise 95% interval",
         "Observation pointwise 95% interval","Source observations"];
        orientation=:horizontal,nbanks=2,fontsize=14)
    hidden_note = input["hidden_lactate"] ?
        "Lactate has no measured series or physical scale; its panel shows only the normalized latent state." :
        "All eight source series have physical scales in mg/g pulp."
    Label(figure[5,1:4],
        "Source: Moreno-Zambrano et al. (2022), doi:10.1098/rsos.210274. $hidden_note\n" *
        "Bands use $(settings["predictive_draws"]) retained joint draws. Lines connect observation times; the fixed initial row is excluded.\n" *
        "Gaussian observation bands may extend below zero. Sampling checks do not establish model validity or published-figure agreement.";
        fontsize=13,halign=:left)
    save(joinpath(directory,"predictive.png"),figure)
    save(joinpath(directory,"predictive.svg"),figure)

    _,summaries = benchmark_csv(joinpath(directory,"parameters.csv"))
    parameters = getindex.(summaries,1)
    samples,chains = settings["samples_per_chain"],settings["chains"]
    traces = map(1:chains) do chain
        draw_header,draw_rows = benchmark_csv(joinpath(directory,"chain_$(lpad(chain,3,'0')).csv"))
        draw_header[3:end-2] == parameters || error("Checkpoint parameter order differs")
        length(draw_rows) == samples || error("Incomplete retained chain")
        parse.(Int,getindex.(draw_rows,1)) == fill(chain,samples) || error("Chain identity differs")
        parse.(Int,getindex.(draw_rows,2)) == collect(1:samples) || error("Draw order differs")
        matrix = reduce(hcat,[parse.(Float64,row[3:end-2]) for row in draw_rows])
        all(isfinite,matrix) || error("Nonfinite retained parameter")
        matrix
    end
    colors = Makie.wong_colors()
    ncols,nrows = 5,cld(length(parameters),5)
    traces_figure = Figure(size=(1800,330nrows+200),fontsize=14,figure_padding=30)
    Label(traces_figure[1,1:ncols],"Cocoa benchmark $trial: retained traces — $status";
        fontsize=24,font=:bold,color=title_color)
    for (i,parameter) in enumerate(parameters)
        rhat,bulk,tail = parse.(Float64,summaries[i][6:8])
        passed = all(isfinite,(rhat,bulk,tail)) && rhat < report["rhat_threshold"] &&
            min(bulk,tail) >= report["minimum_ess"]
        ess_label(value) = isfinite(value) ? string(round(Int,value)) : "undefined"
        title = "$parameter | R-hat $(round(rhat;digits=5))\nESS bulk $(ess_label(bulk)), tail $(ess_label(tail))"
        axis = Axis(traces_figure[(i-1)÷ncols+2,(i-1)%ncols+1];
            title,titlecolor=passed ? :black : :firebrick,titlesize=16,
            xlabel="Retained draw",ylabel="Normalized coefficient")
        for chain in 1:chains
            lines!(axis,1:samples,traces[chain][i,:];
                color=(colors[mod1(chain,length(colors))],0.75),linewidth=0.8)
        end
    end
    Legend(traces_figure[nrows+2,1:ncols],
        [LineElement(color=colors[mod1(chain,length(colors))],linewidth=2) for chain in 1:chains],
        ["Chain $chain" for chain in 1:chains];orientation=:horizontal)
    Label(traces_figure[nrows+3,1:ncols],
        "All retained draws; warmup excluded. Red titles fail R-hat < $(report["rhat_threshold"]) or ESS ≥ $(report["minimum_ess"]).\n" *
        "Recorded divergences: $(report["divergences"]); depth hits: $(report["depth_hits"]). " *
        "Normalized-coordinate μ and k have units h⁻¹; ks, yc, and σ are dimensionless.\n" *
        "Source: doi:10.1098/rsos.210274. Priors have positive support with no upper bound.";
        fontsize=13,halign=:left)
    save(joinpath(directory,"traces.png"),traces_figure)
    save(joinpath(directory,"traces.svg"),traces_figure)
    source_files = ["../settings.toml","input.toml","report.toml","parameters.csv",
        "predictive.csv","predictive_selection.csv",
        ["chain_$(lpad(chain,3,'0')).csv" for chain in 1:chains]...]
    file_hash(path) = open(io->bytes2hex(sha256(io)),path)
    open(joinpath(directory,"figure_sources.toml"),"w") do io
        TOML.print(io,Dict("script"=>"scripts/cocoa_benchmark_figures.jl",
            "script_sha256"=>file_hash(@__FILE__),"julia_version"=>string(VERSION),
            "source_sha256"=>Dict(file=>file_hash(joinpath(directory,file)) for file in source_files));
            sorted=true)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia --project=scripts scripts/cocoa_benchmark_figures.jl <trial-directory>")
    cocoa_benchmark_figures(only(ARGS))
end
