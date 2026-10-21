using CairoMakie

# Plot existing baseline outputs; no model fitting or parameter changes occur here.
# First run scripts/prepare_cocoa.jl and examples/cocoa_baseline.jl in the main project.
root = normpath(joinpath(@__DIR__, ".."))
states = ["G", "F", "E", "L", "A", "Xy", "Xl", "Xa"]
labels = ["Glucose", "Fructose", "Ethanol", "Lactic acid", "Acetic acid", "Yeast", "LAB", "AAB"]
titles = Dict("camu" => "Camu: heap fermentation", "box1" => "Papalexandratou: box 1", "box2" => "Papalexandratou: box 2")
digitized = split.(readlines(joinpath(root, "data", "raw", "digitized", "moreno2018-figures", "points.csv"))[2:end], ',')
for trial in ("camu", "box1", "box2")
    trajectory_file = joinpath(root, "examples", "output", "cocoa_baseline", "$trial.csv")
    trajectory = reduce(vcat, permutedims.(parse.(Float64, split(line, ',')) for line in readlines(trajectory_file)[2:end]))
    data_file = joinpath(root, "data", "processed", "cocoa", "moreno2018", "$trial.csv")
    observations = split.(readlines(data_file)[2:end], ',')
    figure = Figure(size=(1500, 950), fontsize=16, figure_padding=30)
    Label(figure[1,1:4], titles[trial]; fontsize=24, font=:bold)
    for i in 1:8
        axis = Axis(figure[(i-1) ÷ 4 + 2, (i-1) % 4 + 1]; title=labels[i],
            xlabel="Time (h)", ylabel="mg/g pulp", xticks=0:24:144)
        lines!(axis, trajectory[:,1], trajectory[:,i+1]; color=:firebrick, linewidth=2.5)
        points = filter(row -> row[4] == states[i], observations)
        scatter!(axis, parse.(Float64, getindex.(points, 3)), parse.(Float64, getindex.(points, 5));
            color=:black, markersize=7)
        references = filter(row -> row[1] == trial && row[2] == states[i], digitized)
        scatter!(axis, parse.(Float64, getindex.(references, 3)), parse.(Float64, getindex.(references, 4));
            color=:dodgerblue3, marker=:xcross, markersize=10)
        xlims!(axis, 0, 144)
        ylims!(axis, 0, nothing)
    end
    Legend(figure[4,1:4], [LineElement(color=:firebrick, linewidth=2.5), MarkerElement(color=:black, marker=:circle),
        MarkerElement(color=:dodgerblue3, marker=:xcross)],
        ["Simulation at published parameter means", "Original observations", "Digitized published central curve"]; orientation=:horizontal)
    Label(figure[5,1:4], "Source: Moreno-Zambrano et al. (2018), doi:10.1098/rsos.180964, supplement Tables S4-S6.\nPoint-parameter reproduction; posterior uncertainty is not reconstructed from marginal summaries.";
        fontsize=13, halign=:left)
    destination = joinpath(root, "docs", "figures")
    mkpath(destination)
    save(joinpath(destination, "cocoa_$trial.png"), figure)
    save(joinpath(destination, "cocoa_$trial.svg"), figure)
end
