@testset "Feeding and temperature schedules" begin
    spec = ModelSpec(; organisms=(:microbe,), substrates=(:sugar,), products=(:acid,),
                     stoichiometry=reshape([1.,-2.,1.],3,1))
    logspec = ModelSpec(; organisms=spec.organisms, substrates=spec.substrates,
                        products=spec.products, stoichiometry=spec.stoichiometry,
                        log_populations=true)
    @test feeding_fraction(1) ≈ 1/3
    @test feeding_fraction(0) == 1
    @test_throws DomainError feeding_fraction(-1)
    @test ForwardDiff.derivative(feeding_fraction, 2.) ≈ -2/25

    physical = [2.,3.,4.]
    @test dilute!(copy(physical), spec; fraction=0.25, substrate_addition=[6.]) == [0.5,6.75,1.]
    logarithmic = initial_state(logspec, physical)
    dilute!(logarithmic, logspec; fraction=0.25, substrate_addition=[6.])
    @test physical_state(logspec, logarithmic) ≈ [0.5,6.75,1.]
    # Feeding remains representable even below the floating-point exp range.
    tiny = [-1000.,3.,4.]
    dilute!(tiny, logspec; fraction=0.5, substrate_addition=[0.])
    @test tiny[1] == -1000. + log(0.5)
    unchanged = copy(physical)
    @test_throws DomainError dilute!(unchanged, spec; fraction=1.1, substrate_addition=[0.])
    @test unchanged == physical
    @test_throws DimensionMismatch dilute!(unchanged, spec; fraction=0.5, substrate_addition=[])
    @test_throws DomainError dilute!(unchanged, spec; fraction=0.5, substrate_addition=[-1.])
    @test_throws DomainError dilute!(logarithmic, logspec; fraction=0., substrate_addition=[1.])
    @test dilute!(copy(physical), spec; fraction=0., substrate_addition=[6.]) == [0.,6.,0.]

    # A zero-growth problem permits large steps; callbacks must still hit every feed.
    frozen = (mu_max=[0.], Ks=ones(1,1))
    times = [0.,1.25,3.,4.]
    callback = dilution_callback(spec, times; fraction=0.5, substrate_addition=[2.])
    times[2] = 100. # callback owns its schedule
    problem = fermentation_problem(spec, physical, (0.,4.), frozen)
    solution = solve(problem, Tsit5(); callback, adaptive=false, dt=4.)
    @test successful_retcode(solution)
    expected = copy(physical)
    for _ in 1:4
        expected = [expected[1]/2, expected[2]/2+2, expected[3]/2]
    end
    @test solution.u[end] ≈ expected
    @test count(==(1.25), solution.t) == 2
    @test count(==(3.), solution.t) == 2
    @test solve(problem, Tsit5(); callback, adaptive=false, dt=4.).u[end] ≈ expected
    @test_throws ArgumentError dilution_callback(spec, [1.,1.]; fraction=0.5, substrate_addition=[2.])
    outside = dilution_callback(spec, [5.]; fraction=0.5, substrate_addition=[2.])
    @test_throws ArgumentError solve(problem, Tsit5(); callback=outside)

    # Synthetic 1:1:1 repeated feeding converges to the same phase of each cycle.
    growth = (mu_max=[1.], Ks=ones(1,1))
    feeds = collect(8.:8.:320.)
    fraction = feeding_fraction(1)
    cycles = solve(fermentation_problem(spec, [0.1,10.,0.], (0.,320.), growth), Tsit5();
        callback=dilution_callback(spec, feeds; fraction, substrate_addition=[(1-fraction)*10]),
        abstol=1e-10, reltol=1e-9)
    @test successful_retcode(cycles)
    @test all(u -> all(>=(0), u), cycles.u)
    postfeed(t) = cycles.u[findlast(==(t), cycles.t)]
    @test maximum(abs.(postfeed(feeds[end])-postfeed(feeds[end-1]))) < 1e-7
    @test maximum(abs.(postfeed(feeds[end-1])-postfeed(feeds[end-2]))) < 1e-7
    @test postfeed(feeds[end])[1] > 1 # distinguish an active culture from washout

    schedule = TemperatureSchedule([0.,2.,4.], [30.,10.,20.])
    @test [schedule(t) for t in [-1.,0.,1.99,2.,3.,4.,5.]] == [30.,30.,30.,10.,10.,20.,20.]
    @test temperature_stops(schedule, (0.,4.)) == [2.]
    linear = TemperatureSchedule([0.,2.,4.], [30.,10.,20.]; interpolation=:linear)
    @test [linear(t) for t in [-1.,1.,2.,3.,5.]] == [30.,20.,10.,15.,20.]
    @test ForwardDiff.derivative(linear, 1.) == -10.
    @test_throws ArgumentError TemperatureSchedule([], [])
    @test_throws ArgumentError TemperatureSchedule([1.,1.], [20.,30.])
    @test_throws ArgumentError TemperatureSchedule([1.,0.], [20.,30.])
    @test_throws DimensionMismatch TemperatureSchedule([0.,1.], [20.])
    @test_throws DomainError TemperatureSchedule([0.], [NaN])
    @test_throws ArgumentError TemperatureSchedule([0.], [20.]; interpolation=:cubic)

    # Resource-unlimited growth has an analytic integral across a temperature step.
    temperature = TemperatureSchedule([0.,2.], [30.,10.])
    response(T) = γ_T_ctmi(T; Tmin=0., Topt=30., Tmax=40.)
    thermal = ModelSpec(; organisms=(:microbe,), substrates=(), products=(),
        stoichiometry=ones(1,1), modifiers=((u,p,t)->response(temperature(t)),))
    p = (mu_max=[0.5], Ks=zeros(0,1))
    stops = temperature_stops(temperature, (0.,4.))
    trajectory = solve(fermentation_problem(thermal, [1.], (0.,4.), p), Tsit5();
        tstops=stops, d_discontinuities=stops, abstol=1e-10, reltol=1e-9)
    @test successful_retcode(trajectory)
    @test trajectory.u[end][1] ≈ exp(0.5*(2response(30.)+2response(10.))) rtol=1e-7
    du = zeros(1)
    build_rhs(thermal)(du, [1.], p, 1.)
    warm = du[1]
    build_rhs(thermal)(du, [1.], p, 3.)
    @test du[1]/warm ≈ response(10.)
end
