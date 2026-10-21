using OrdinaryDiffEq

@testset "Consumer-resource assembly" begin
    matrix = reshape([1., -2., 1.], 3, 1)
    spec = ModelSpec(; organisms=(:microbe,), substrates=(:sugar,), products=(:acid,),
                     stoichiometry=matrix, conserved_weights=ones(3))
    parameters = (mu_max=[0.8], Ks=reshape([1.], 1, 1))
    u0 = [0.1, 10., 0.]
    problem = fermentation_problem(spec, u0, (0., 8.), parameters)
    solution = solve(problem, Tsit5(); abstol=1e-11, reltol=1e-10, saveat=0.25)
    @test successful_retcode(solution)
    @test all(u -> all(>=(0), u), solution.u)
    @test maximum(abs(sum(u)-sum(u0)) for u in solution.u) < 1e-9
    # Analytic integrated Monod law, using C=X+Y*S and K'=Y*Ks.
    C, Kprime = u0[1] + u0[2]/2, parameters.Ks[1]/2
    elapsed(X) = ((C+Kprime)/C * log(X/u0[1]) - Kprime/C * log((C-X)/(C-u0[1]))) / parameters.mu_max[1]
    function analytic_population(t)
        lo, hi = u0[1], C
        for _ in 1:80
            mid = lo + (hi-lo)/2
            (mid == lo || mid == hi) && break
            elapsed(mid) < t ? (lo = mid) : (hi = mid)
        end
        return lo + (hi-lo)/2
    end
    # Compare states: inverse elapsed time is ill-conditioned near exhaustion.
    @test maximum(abs(u[1]-analytic_population(t)) for (u,t) in zip(solution.u, solution.t)) < 1e-8

    logspec = ModelSpec(; organisms=spec.organisms, substrates=spec.substrates, products=spec.products,
                        stoichiometry=matrix, log_populations=true)
    logsolution = solve(fermentation_problem(logspec, u0, (0.,8.), parameters), Tsit5();
                        abstol=1e-11, reltol=1e-10, saveat=0.25)
    @test successful_retcode(logsolution)
    @test maximum(maximum(abs.(physical_state(logspec, v)-u)) for (u,v) in zip(solution.u, logsolution.u)) < 1e-7
    @test physical_state(logspec, initial_state(logspec, u0)) ≈ u0
    @test_throws DomainError initial_state(logspec, [0., 10., 0.])

    rhs! = build_rhs(spec)
    du = zeros(3)
    rhs!(du, [0.,10.,0.], parameters, 0.)
    @test all(iszero, du)
    rhs!(du, [1.,0.,0.], parameters, 0.)
    @test all(iszero, du)
    conto = ModelSpec(; organisms=spec.organisms, substrates=spec.substrates, products=spec.products,
                      stoichiometry=matrix, kinetics=:contois)
    build_rhs(conto)(du, [2.,3.,0.], parameters, 0.)
    @test du ≈ [0.96, -1.92, 0.96]

    # All consumed substrates must be available, and modifier selection is per organism.
    coupled = ModelSpec(; organisms=(:a,:b), substrates=(:s1,:s2), products=(),
        stoichiometry=[1. 0.; 0. 1.; -1. -1.; -1. 0.],
        modifiers=((u,p,t)->0.5, (u,p,t)->0.25))
    p2 = (mu_max=[1.,1.], Ks=ones(2,2))
    du2 = zeros(4)
    build_rhs(coupled)(du2, [1.,1.,1.,0.], p2, 0.)
    @test du2 ≈ [0.,0.125,-0.125,0.]

    endpoint(mu) = solve(fermentation_problem(spec, u0 .* one(mu), (0.,2.),
                              (mu_max=[mu], Ks=parameters.Ks)), Tsit5();
                              abstol=1e-10, reltol=1e-9).u[end][1]
    @test ForwardDiff.derivative(endpoint, 0.8) ≈ (endpoint(0.80001)-endpoint(0.79999))/0.00002 rtol=1e-5
    @test_throws ArgumentError ModelSpec(; organisms=(:x,), substrates=(:s,), products=(),
        stoichiometry=reshape([1.,-2.],2,1), conserved_weights=ones(2))
    @test_throws ArgumentError ModelSpec(; organisms=(:x,), substrates=(:x,), products=(), stoichiometry=zeros(2,1))
    @test_throws DomainError fermentation_problem(spec, u0, (0.,1.), (mu_max=[-1.], Ks=parameters.Ks))
    @test_throws DomainError fermentation_problem(spec, [-1.,1.,0.], (0.,1.), parameters)
    @test_throws DimensionMismatch fermentation_problem(spec, u0, (0.,1.), (mu_max=[1.], Ks=ones(2,1)))
end
