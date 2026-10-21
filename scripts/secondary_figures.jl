using FermentModels
using CairoMakie

# Synthetic demonstration parameters only; no fitted organism claims.
figure = Figure(size=(1050, 700))
curves = (
    ("Cardinal temperature", "Temperature (°C)", -5:0.1:45, t -> γ_T_ctmi(t; Tmin=0, Topt=30, Tmax=40)),
    ("Gänzle temperature", "Temperature (°C)", -10:0.1:45, t -> γ_T_ganzle(t; a=2, b=2, c=0.2, Tref=40)),
    ("Cardinal pH", "pH", 2:0.01:8, h -> γ_pH_cardinal(h; pHmin=3, pHopt=5, pHmax=7)),
    ("Undissociated acid", "HA (mmol/L)", 0:0.2:120, x -> γ_acid(x; MIC=100, α=2)),
    ("Salt inhibition", "NaCl (% w/w)", 0:0.01:6, x -> γ_salt(x; max=4, α=1)),
)
for (i, (title, xlabel, x, f)) in enumerate(curves)
    axis = Axis(figure[(i-1) ÷ 3 + 1, (i-1) % 3 + 1]; title, xlabel, ylabel="Relative growth factor")
    lines!(axis, x, f.(x); color=:teal, linewidth=2.5)
    ylims!(axis, -0.03, 1.05)
end
Label(figure[3, 1:3], "Synthetic parameters for shape review; not food predictions.")
directory = joinpath(@__DIR__, "..", "docs", "figures")
mkpath(directory)
save(joinpath(directory, "secondary_models.png"), figure)
save(joinpath(directory, "secondary_models.svg"), figure)
