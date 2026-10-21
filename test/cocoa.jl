using OrdinaryDiffEq
using ForwardDiff
using Distributions: logpdf
include(joinpath(@__DIR__, "..", "scripts", "prepare_cocoa.jl"))
using .CocoaSourceData: cocoa_trial, published_summaries, TRIALS

@testset "Cocoa source equations and boundaries" begin
    p = NamedTuple{COCOA_PARAMETERS}(ntuple(_ -> 1.0, 24)) # synthetic algebra check
    du = zeros(8)
    @test cocoa_rhs!(du, ones(8), p, 0.0) === nothing
    @test du == [-1.0, -0.5, 0.5, 0.0, 1.5, 0.0, -0.5, 0.0]
    cocoa_rhs!(du, zeros(8), p, 0.0)
    @test du == zeros(8) # includes the otherwise undefined Contois origin
    for index in 1:8
        state = ones(8)
        state[index] = 0.0
        cocoa_rhs!(du, state, p, 0.0)
        @test du[index] >= 0 # vector field points into the positive orthant
    end
    state = [2., 3., 4., 5., 6., 0., 0., 0.]
    cocoa_rhs!(du, state, p, 0.0)
    @test du == zeros(8) # no spontaneous microbial generation
    cocoa_rhs!(du, [-1., -1., -1., -1., -1., -1., -1., -1.], p, 0.0)
    @test du == zeros(8)
    @test_throws DimensionMismatch cocoa_problem(ones(7), (0., 1.), p)
    @test_throws DomainError cocoa_problem(-ones(8), (0., 1.), p)
    @test_throws DomainError cocoa_problem(ones(8), (0., 1.), merge(p, (ks5=0.,)))
    @test_throws DomainError cocoa_problem(ones(8), (0., 1.), merge(p, (mu1=-1.,)))
    @test_throws ArgumentError cocoa_problem(ones(8), (1., 0.), p)
    @test_throws ArgumentError cocoa_problem(ones(8), (0., 1.), (mu1=1.,))
    problem = cocoa_problem(ones(8), (0., 1.), p)
    @test problem.kwargs[:isoutofdomain]([-1.; ones(7)], p, 0.)
    @test problem.kwargs[:isoutofdomain]([NaN; ones(7)], p, 0.)
end

@testset "Original cocoa inferred-noise likelihood and priors" begin
    source = cocoa_trial(:camu)
    data = cocoa_observations(source.initial, source.times, source.observations)
    @test data.scales == source.scales
    @test length(data.values) == 128
    @test data.values[1:16] == source.observations[:,1] ./ source.scales[1]
    @test data.sigma == :sigma
    initial = load_params(joinpath(@__DIR__, "..", "params", "cocoa_camu_scaled.toml")).values
    lower = NamedTuple{keys(initial)}(ntuple(_ -> 1e-5,25))
    upper = NamedTuple{keys(initial)}(ntuple(_ -> 5.,25))
    model = cocoa_fit_model(; initial, lower, upper)
    predictions = model.predict(initial,data)
    dense_solution = solve(cocoa_problem(data.initial,(data.t0,last(data.times)),initial),
        Tsit5();abstol=1e-6,reltol=1e-6)
    @test predictions ≈ vec([dense_solution(t)[j] for t in data.times,j in 1:8]) atol=1e-14
    response(mu1) = sum(abs2,model.predict(merge(initial,(;mu1)),data)-data.values)
    step = 1e-5
    @test ForwardDiff.derivative(response,initial.mu1) ≈
        (response(initial.mu1+step)-response(initial.mu1-step))/(2step) rtol=2e-4
    squared_error = sum(abs2, predictions-data.values)
    @test gaussian_nll(model,data,initial) ≈ squared_error/(2initial.sigma^2) +
        128log(initial.sigma) + 64log(2pi)
    priors = cocoa_priors()
    @test keys(priors) == (COCOA_PARAMETERS..., :sigma)
    @test all(entry -> minimum(entry.distribution) == 0 && maximum(entry.distribution) == Inf, values(priors))
    @test logpdf(priors.sigma.distribution,1.) ≈ -log(pi)
    @test all(entry -> entry.provenance.source == "doi:10.1098/rsos.180964", values(priors))
    @test_throws DimensionMismatch cocoa_observations(source.initial,source.times,zeros(8,16))
    @test_throws ArgumentError cocoa_observations(source.initial,reverse(source.times),source.observations)
end

@testset "Cocoa stiff-solver agreement and derivatives" begin
    for trial in TRIALS
        source = cocoa_trial(trial)
        data = cocoa_observations(source.initial,source.times,source.observations)
        initial = load_params(joinpath(@__DIR__,"..","params","cocoa_$(trial)_scaled.toml")).values
        lower = NamedTuple{keys(initial)}(ntuple(_ -> 1e-5,25))
        upper = NamedTuple{keys(initial)}(ntuple(_ -> 5.,25))
        explicit = cocoa_fit_model(;initial,lower,upper,abstol=1e-10,reltol=1e-9)
        stiff = cocoa_fit_model(;initial,lower,upper,abstol=1e-10,reltol=1e-9,solver=:rodas5p)
        a,b = explicit.predict(initial,data), stiff.predict(initial,data)
        @test !isnothing(b) && all(>=(0),b)
        @test a ≈ b rtol=1e-6 atol=1e-8
        z = log.(collect(values(initial)))
        decode(z) = NamedTuple{keys(initial)}(Tuple(exp.(z)))
        loss(z,model) = gaussian_nll(model,data,decode(z))
        g = ForwardDiff.gradient(z -> loss(z,stiff),z)
        @test all(isfinite,g)
        @test g ≈ ForwardDiff.gradient(z -> loss(z,explicit),z) rtol=2e-5 atol=1e-6
        direction = cos.(eachindex(z))
        step = 1e-4
        finite_difference = (loss(z+step*direction,stiff)-loss(z-step*direction,stiff))/(2step)
        @test sum(g .* direction) ≈ finite_difference rtol=2e-5 atol=1e-6
        @test_throws ArgumentError cocoa_fit_model(;initial,lower,upper,solver=:unknown)
    end
    # Synthetic, analytically soluble fast mortality with slow remaining states.
    # It exercises a stiff regime independently of any fitted biological values.
    initial = NamedTuple{(COCOA_PARAMETERS...,:sigma)}(ntuple(_ -> 1.,25))
    lower = NamedTuple{keys(initial)}(ntuple(_ -> 1e-5,25))
    upper = NamedTuple{keys(initial)}(ntuple(_ -> 5.,25))
    model = cocoa_fit_model(;initial,lower,upper,solver=:rodas5p,
        abstol=1e-10,reltol=1e-9,maxiters=10_000)
    p = merge(initial,(mu1=0.,mu2=0.,mu3=0.,mu4=0.,mu5=0.,k1=1e5,k2=0.3,k3=0.4))
    data = (initial=ones(8),t0=0.,times=[0.1,1.],values=zeros(16),sigma=:sigma)
    prediction = reshape(model.predict(p,data),2,8)
    @test all(>=(0),prediction)
    @test prediction[:,1:5] ≈ ones(2,5)
    expected = exp.(-data.times .* [p.k1 p.k2 p.k3])
    @test prediction[:,6:8] ≈ expected rtol=1e-6 atol=1e-9
end

@testset "Cocoa analytic mortality and additive growth limits" begin
    base = NamedTuple{COCOA_PARAMETERS}(ntuple(_ -> 1.0, 24))
    p = merge(base, (mu1=0., mu2=0., mu3=0., mu4=0., mu5=0., k1=0.2, k2=0.3, k3=0.4))
    initial = [2., 3., 0.5, 0.7, 0.9, 0.4, 0.3, 0.2]
    sol = solve(cocoa_problem(initial, (0., 4.), p), Tsit5(); abstol=1e-12, reltol=1e-11)
    expected = initial[6:8] .* exp.(-4 .* [p.k1*initial[3], p.k2*initial[4], p.k3*initial[5]^2])
    @test sol(4.)[1:5] ≈ initial[1:5]
    @test sol(4.)[6:8] ≈ expected rtol=1e-9
    # Zero yields keep substrates fixed; both yeast routes contribute to the
    # exact exponential rate, while each zero microbial group stays absent.
    p = merge(base, NamedTuple{Tuple(Symbol("yc",i) for i in 1:11)}(ntuple(_ -> 0., 11)),
        (mu1=0.2, mu2=0.3, k1=0., k2=0., k3=0.))
    initial = [2., 3., 0., 0., 0., 0.1, 0., 0.]
    sol = solve(cocoa_problem(initial, (0., 4.), p), Tsit5(); abstol=1e-12, reltol=1e-11)
    rate = 0.2*2/3 + 0.3*3/4
    @test sol(4.)[6] ≈ 0.1exp(4rate) rtol=1e-9
    @test sol(4.)[[1,2,3,4,5,7,8]] ≈ initial[[1,2,3,4,5,7,8]]
end

@testset "Cocoa published parameters and coordinate invariance" begin
    root = normpath(joinpath(@__DIR__, ".."))
    summaries = published_summaries()
    @test length(summaries) == 75
    for trial in TRIALS
        data = cocoa_trial(trial)
        loaded = load_params(joinpath(root, "params", "cocoa_$(trial)_scaled.toml"))
        p = loaded.values
        @test length(loaded.provenance) == 25
        @test all(row -> getproperty(p, row.parameter) == row.mean,
                  filter(row -> row.experiment == trial, summaries))
        @test all(row -> occursin("10.6084/m9.figshare.7152614.v2", row.source), loaded.provenance)
        physical = cocoa_physical_parameters(p, data.scales)
        roundtrip = cocoa_scaled_parameters(physical, data.scales)
        @test collect(values(roundtrip)) ≈ [getproperty(p, n) for n in COCOA_PARAMETERS]
        @test physical.yc3 ≈ p.yc3 * data.scales[2]/data.scales[6]
        @test physical.ks5 ≈ p.ks5 * data.scales[4]/data.scales[8]
        @test physical.k3 ≈ p.k3 / data.scales[5]^2
        @test !hasproperty(physical, :sigma)
        scaled_sol = solve(cocoa_problem(data.initial ./ data.scales, (0., 144.), p), Tsit5();
            abstol=1e-11, reltol=1e-10, saveat=0.:6.:144.)
        physical_sol = solve(cocoa_problem(data.initial, (0., 144.), physical), Tsit5();
            abstol=data.scales .* 1e-11, reltol=1e-10, saveat=0.:6.:144.)
        @test successful_retcode(scaled_sol) && successful_retcode(physical_sol)
        @test all(>=(0), Array(scaled_sol))
        @test Array(physical_sol) ./ data.scales ≈ Array(scaled_sol) rtol=1e-7 atol=1e-9
        @test all(<=(1e-10), diff(Array(scaled_sol)[1,:]))
        @test all(<=(1e-10), diff(Array(scaled_sol)[2,:]))
        @test all(>=(-1e-10), diff(Array(scaled_sol)[5,:]))
    end
    p = load_params(joinpath(root, "params", "cocoa_camu_scaled.toml")).values
    @test_throws DimensionMismatch cocoa_physical_parameters(p, ones(7))
    @test_throws DomainError cocoa_physical_parameters(p, zeros(8))
end

@testset "Cocoa parameter derivatives through the solver" begin
    data = cocoa_trial(:camu)
    p = load_params(joinpath(@__DIR__, "..", "params", "cocoa_camu_scaled.toml")).values
    function response(mu1)
        sol = solve(cocoa_problem(data.initial ./ data.scales, (0., 12.), merge(p, (; mu1))),
            Tsit5(); abstol=1e-12, reltol=1e-11)
        return sol(12.)[3]
    end
    derivative = ForwardDiff.derivative(response, p.mu1)
    step = 1e-5
    finite_difference = (response(p.mu1 + step) - response(p.mu1 - step))/(2step)
    @test isfinite(derivative)
    @test derivative ≈ finite_difference rtol=1e-5
end

@testset "Published cocoa curves: documented reproduction gap" begin
    points = split.(readlines(joinpath(@__DIR__, "..", "data", "raw", "digitized",
        "moreno2018-figures", "points.csv"))[2:end], ',')
    @test length(points) == 168
    errors = Float64[]
    for trial in TRIALS
        data = cocoa_trial(trial)
        p = load_params(joinpath(@__DIR__, "..", "params", "cocoa_$(trial)_scaled.toml")).values
        solution = solve(cocoa_problem(data.initial ./ data.scales, (0., 144.), p),
            Tsit5(); abstol=1e-11, reltol=1e-10)
        for row in filter(r -> r[1] == String(trial), points)
            i = findfirst(==(Symbol(row[2])), COCOA_STATES)
            time, reference, span = parse.(Float64, row[[3,4,8]])
            push!(errors, abs(solution(time)[i]*data.scales[i]-reference)/span)
        end
    end
    # The published rounded means do not reproduce all plotted central curves.
    # Preserve this acceptance check; do not enlarge its preselected tolerance.
    # Joint source draws and the plotting script are absent. See cocoa_baseline.md.
    @test_broken all(<=(0.02), errors)
end
