@testset "Local sensitivities and unresolved directions" begin
    model = FitModel((p,d)->p.a .+ p.b .* d.times; initial=(a=2.,b=3.),
        lower=(a=0.01,b=0.01), upper=(a=10.,b=10.))
    data = (times=[0.,1.,2.,3.], values=[2.,missing,8.,11.], sigma=[0.5,missing,1.,2.])
    result = local_sensitivities(model,data,(b=3.,a=2.))
    design = [1. 0.; 1. 2.; 1. 3.]
    @test result.parameter_names == (:a,:b)
    @test result.observed_indices == [1,3,4]
    @test result.jacobian ≈ design * Diagonal([2.,3.]) ./ [0.5,1.,2.]
    @test result.status == :full_rank && result.rank == 2
    @test size(result.unresolved_directions) == (2,0)
    @test isfinite(result.condition_number)
    raw = local_sensitivities(model,data;scale=:parameter,standardized=false)
    @test raw.jacobian ≈ design
    @test raw.column_norms.a ≈ sqrt(3.)
    underdetermined = local_sensitivities(model,(times=[1.],values=[5.],sigma=1.))
    @test underdetermined.rank == 1
    @test size(underdetermined.unresolved_directions) == (2,1)
    @test norm(underdetermined.jacobian * underdetermined.unresolved_directions) < 1e-12
    constant = FitModel((p,d)->fill(one(p.a),length(d.values));initial=(a=1.,),lower=(a=0.1,),upper=(a=10.,))
    zero = local_sensitivities(constant,(values=[1.,1.],sigma=1.))
    @test zero.rank == 0 && zero.status == :rank_deficient && isinf(zero.condition_number)
    @test_throws ArgumentError local_sensitivities(model,data;scale=:unsupported)
    @test_throws ArgumentError local_sensitivities(model,data,(c=1.,))
    @test_throws DomainError local_sensitivities(model,data,(a=0.,b=3.))
    @test_throws DomainError local_sensitivities(model,data;rtol=0.)
end

@testset "Sensitivity through a Monod ODE" begin
    spec = ModelSpec(;organisms=(:microbe,),substrates=(:sugar,),products=(:acid,),
        stoichiometry=reshape([1.,-2.,1.],3,1))
    function predict(p,d)
        problem = fermentation_problem(spec,[0.2,4.,0.] .* one(p.mu+p.Ks),(0.,4.),
            (mu_max=[p.mu],Ks=reshape([p.Ks],1,1)))
        trajectory = solve(problem,Tsit5();abstol=1e-11,reltol=1e-10)
        return [trajectory(t)[1] for t in d.times]
    end
    model = FitModel(predict;initial=(mu=0.8,Ks=0.7),lower=(mu=0.1,Ks=0.1),upper=(mu=2.,Ks=5.))
    data = (times=[1.,2.,3.,4.],values=zeros(4),sigma=0.05)
    result = local_sensitivities(model,data)
    h = 1e-5
    reference = hcat([(predict(merge(model.initial,NamedTuple{(name,)}((model.initial[name]*exp(h),))),data)-
                       predict(merge(model.initial,NamedTuple{(name,)}((model.initial[name]*exp(-h),))),data)) /
                       (2h*data.sigma) for name in keys(model.initial)]...)
    @test result.jacobian ≈ reference rtol=1e-5
    @test result.rank == 2
end

@testset "Profiles identify exact parameter confounding" begin
    model = FitModel((p,d)->fill(p.a*p.b,length(d.values));initial=(a=1.,b=1.),
        lower=(a=0.1,b=0.1),upper=(a=10.,b=10.))
    data = (values=fill(2.,12),sigma=0.1)
    sensitivity = local_sensitivities(model,data,(a=1.,b=2.))
    @test sensitivity.status == :rank_deficient && sensitivity.rank == 1
    @test norm(sensitivity.jacobian*sensitivity.unresolved_directions) < 1e-10
    grid = [0.5,1.,2.,4.]
    profile = profile_likelihood(model,data,:a,grid;starts=2,rng=MersenneTwister(67))
    @test profile.reference_valid
    @test profile.status == :flat_on_grid
    @test profile.supported_values == grid
    @test maximum(abs(p.deviance) for p in profile.points) < 1e-6
    @test all(p -> isapprox(p.parameters.b,2/p.value;atol=1e-6), profile.points)
    @test all(p -> length(p.runs) == 2,profile.points)
    @test isnothing(only(profile.regions).lower_bracket) && isnothing(only(profile.regions).upper_bracket)
end

@testset "Profile likelihood agrees with analytical Gaussian profiles" begin
    # At fixed slope, the nuisance intercept has a closed-form optimum.
    model = FitModel((p,d)->p.a .+ p.b .* d.times;initial=(a=1.,b=1.),
        lower=(a=0.1,b=0.1),upper=(a=10.,b=5.))
    times = collect(range(0.,2.;length=20))
    data = (;times,values=2. .+ 1.5 .* times,sigma=0.2)
    grid = collect(1.1:0.1:1.9)
    profile = profile_likelihood(model,data,:b,grid;starts=2,rng=MersenneTwister(4))
    expected = (grid .- 1.5).^2 .* sum(abs2,times .- mean(times)) ./ data.sigma^2
    @test [p.deviance for p in profile.points] ≈ expected atol=1e-6
    @test all(p -> isapprox(p.parameters.a,2+(1.5-p.value)*mean(times);atol=1e-5),profile.points)
    @test profile.threshold ≈ 3.841458820694124 rtol=1e-10
    @test profile.status == :crossings_bracketed
    @test collect(only(profile.regions).lower_bracket) ≈ [1.3,1.4]
    @test collect(only(profile.regions).upper_bracket) ≈ [1.6,1.7]

    level = FitModel((p,d)->fill(p.level,length(d.values));initial=(level=1.,),lower=(level=0.1,),upper=(level=4.,))
    single = profile_likelihood(level,(values=fill(2.,10),sigma=0.2),:level,[1.,2.,3.];starts=1)
    @test [p.deviance for p in single.points] ≈ [250.,0.,250.] atol=1e-6
    @test all(isempty(p.runs) for p in single.points)
    @test single.status == :crossings_bracketed
    limited = profile_likelihood(level,(values=fill(2.,10),sigma=0.2),:level,[2.,2.05];starts=1)
    @test limited.status == :grid_limited
    @test_throws ArgumentError profile_likelihood(level,data,:other,[1.,2.])
    @test_throws ArgumentError profile_likelihood(level,data,:level,[2.,1.])
    @test_throws ArgumentError profile_likelihood(level,data,:level,[1.,1.])
    @test_throws DomainError profile_likelihood(level,data,:level,[0.01,1.])
    @test_throws DomainError profile_likelihood(level,data,:level,[1.,2.];confidence=1.)
end

@testset "Failed profile fits remain explicit gaps" begin
    model = FitModel((p,d)->1.5 < p.a < 2.5 ? nothing : fill(p.a*p.b,length(d.values));
        initial=(a=1.,b=1.),lower=(a=0.1,b=0.1),upper=(a=10.,b=10.))
    profile = profile_likelihood(model,(values=fill(1.,8),sigma=0.1),:a,[0.5,1.,2.,3.,4.];starts=1)
    @test profile.status == :incomplete_grid
    @test !profile.points[3].converged && ismissing(profile.points[3].deviance)
    @test profile.points[3].runs[1].retcode == "PredictionFailure"
    @test length(profile.regions) == 2
    @test isnothing(profile.regions[1].upper_bracket)
    @test isnothing(profile.regions[2].lower_bracket)
    invalid = FitModel((p,d)->nothing;initial=(a=1.,),lower=(a=0.1,),upper=(a=10.,))
    failed = profile_likelihood(invalid,(values=[1.],sigma=1.),:a,[1.,2.];starts=1)
    @test failed.status == :no_reference_fit
    @test !failed.reference.converged && isnothing(failed.supported_values)
end

@testset "A profile can expose an inferior reference optimum" begin
    # A positive response with a global minimum at a=1 and a separate local
    # minimum near a=3. Starting only in the second basin misses the global fit.
    model = FitModel((p,d)->fill((p.a-1)^2*((p.a-3)^2+0.1),length(d.values));
        initial=(a=3.,),lower=(a=0.2,),upper=(a=5.,))
    profile = profile_likelihood(model,(values=[0.],sigma=0.1),:a,[0.8,1.,3.];starts=1)
    @test profile.reference.converged
    @test profile.status == :better_profile_optimum
    @test !profile.reference_valid && isnothing(profile.supported_values)
    @test profile.points[2].deviance < -1
end
