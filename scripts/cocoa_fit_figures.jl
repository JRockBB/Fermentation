using CairoMakie
using TOML

"""Plot saved joint-posterior predictions, with their convergence status visible."""
function cocoa_fit_figures(directory; root=normpath(joinpath(@__DIR__,"..")))
    settings = TOML.parsefile(joinpath(directory,"settings.toml"))
    report = TOML.parsefile(joinpath(directory,"report.toml"))
    settings["method"] == "bayes" || error("Requires a Bayesian fit directory")
    trial = settings["trial"]
    rows = split.(readlines(joinpath(directory,"predictive.csv"))[2:end], ',')
    references = split.(readlines(joinpath(root,"data","raw","digitized",
        "moreno2018-figures","points.csv"))[2:end], ',')
    states = ["G","F","E","L","A","Xy","Xl","Xa"]
    labels = ["Glucose","Fructose","Ethanol","Lactic acid","Acetic acid","Yeast","LAB","AAB"]
    title = Dict("camu"=>"Camu: heap fermentation", "box1"=>"Papalexandratou: box 1",
        "box2"=>"Papalexandratou: box 2")[trial]
    status = report["converged"] ? "convergence diagnostics passed" : "CONVERGENCE CHECKS FAILED"
    figure = Figure(size=(1500,1000),fontsize=16,figure_padding=30)
    Label(figure[1,1:4],title * " — " * status;fontsize=23,font=:bold,
        color=report["converged"] ? :black : :firebrick)
    for i in 1:8
        data = filter(row -> row[2] == states[i],rows)
        values = column -> parse.(Float64,getindex.(data,column))
        times = values(1)
        axis = Axis(figure[(i-1) ÷ 4 + 2,(i-1) % 4 + 1];title=labels[i],
            xlabel="Time (h)",ylabel="mg/g pulp",xticks=0:24:144)
        band!(axis,times,values(8),values(9);color=(:gray,0.22))
        band!(axis,times,values(6),values(7);color=(:firebrick,0.23))
        lines!(axis,times,values(4);color=:firebrick,linewidth=2.5)
        scatter!(axis,times,values(3);color=:black,markersize=7)
        points = filter(row -> row[1] == trial && row[2] == states[i],references)
        scatter!(axis,parse.(Float64,getindex.(points,3)),parse.(Float64,getindex.(points,4));
            color=:dodgerblue3,marker=:xcross,markersize=10)
        xlims!(axis,0,144)
    end
    Legend(figure[4,1:4], [LineElement(color=:firebrick,linewidth=2.5),
        PolyElement(color=(:firebrick,0.23)),PolyElement(color=(:gray,0.22)),
        MarkerElement(color=:black,marker=:circle),MarkerElement(color=:dodgerblue3,marker=:xcross)],
        ["Joint-posterior mean trajectory","Latent pointwise 95% interval","Observation pointwise 95% interval",
         "Original observations","Published central curve"];
        orientation=:horizontal,nbanks=2,fontsize=14)
    Label(figure[5,1:4],
        "Source: Moreno-Zambrano et al. (2018), doi:10.1098/rsos.180964. Intervals use retained joint draws, including sampled noise.\nPredictions evaluated at observed times; lines connect those times. Gaussian observation intervals may extend below zero.";
        fontsize=13,halign=:left)
    save(joinpath(directory,"predictive.png"),figure)
    save(joinpath(directory,"predictive.svg"),figure)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia --project=scripts scripts/cocoa_fit_figures.jl <fit-directory>")
    cocoa_fit_figures(only(ARGS))
end
