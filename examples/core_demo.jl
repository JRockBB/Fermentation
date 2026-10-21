# P1.1-P1.4 integration example. Every constant below is synthetic.
# These are numerical demonstration units, not CFU or a calibrated food model.
using FermentModels
using OrdinaryDiffEq

function core_demo(output=joinpath(@__DIR__, "output", "core_demo.csv"))
    temperature = TemperatureSchedule([0.,12.,24.], [28.,8.,28.])
    buffer = buffer_from_capacity(; pH0=6., capacity=20., pKa=5., pKw=14.,
                                  acid_pKas=(acid=4.,))
    function environment(u, parameters, t)
        acid = max(u[3], zero(u[3]))
        acidity = pH((acid=acid,), parameters.buffer)
        return γ_T_ctmi(temperature(t); Tmin=0., Topt=30., Tmax=40.) *
               γ_pH_cardinal(acidity; pHmin=3., pHopt=6., pHmax=8.) *
               γ_acid(undissociated_acid(acid, acidity; pKa=4.); MIC=100.)
    end
    # Synthetic biomass equivalents, substrate and acid all use mmol/L here.
    spec = ModelSpec(; organisms=(:microbe,), substrates=(:sugar,), products=(:acid,),
        stoichiometry=reshape([1.,-2.,1.],3,1), modifiers=(environment,),
        log_populations=true, conserved_weights=ones(3))
    parameters = (mu_max=[0.8], Ks=ones(1,1), buffer=buffer)
    problem = fermentation_problem(spec, [0.1,10.,0.], (0.,36.), parameters)
    fraction = feeding_fraction(1)
    feeds = dilution_callback(spec, [12.,24.]; fraction, substrate_addition=[(1-fraction)*10])
    stops = temperature_stops(temperature, problem.tspan)
    solution = solve(problem, Tsit5(); callback=feeds, tstops=stops,
        d_discontinuities=stops, abstol=1e-10, reltol=1e-8, saveat=0.25)
    successful_retcode(solution) || error("Demo integration failed: $(solution.retcode)")
    mkpath(dirname(abspath(output)))
    open(output, "w") do io
        println(io, "time_h,temperature_C,biomass_equiv_mmol_L,sugar_mmol_L,acid_mmol_L,pH")
        for (t, state) in zip(solution.t, solution.u)
            u = physical_state(spec, state)
            row = (t, temperature(t), u..., pH((acid=u[3],), buffer))
            println(io, join(round.(row; digits=8), ','))
        end
    end
    final = physical_state(spec, solution.u[end])
    println("Synthetic demo completed: $(length(solution.t)) saved states over 36 h.")
    println("Final biomass equivalents, sugar, acid (mmol/L): ", round.(final; digits=4))
    println("Final synthetic pH: ", round(pH((acid=final[3],), buffer); digits=4))
    println("CSV: ", abspath(output))
    return solution
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? core_demo() : core_demo(only(ARGS))
end
