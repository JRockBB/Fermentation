using LinearAlgebra
include("../scripts/cocoa_mle_diagnostics.jl")

@testset "MLE point diagnostics: bounds and noise information" begin
    model = FitModel((p,d)->fill(p.a,3);initial=(a=1.,sigma=1.),
        lower=(a=0.1,sigma=0.1),upper=(a=3.,sigma=3.))
    data = (values=[1.,2.,3.],sigma=:sigma)
    p = (a=2.,sigma=sqrt(2/3))
    d = CocoaMLEDiagnostics.point_diagnostics(model,data,p)
    @test d.projected_gradient_inf < 1e-12
    @test d.information_rank == 2
    @test d.information_factor' * d.information_factor ≈ [18. 0.;0. 6.]
    @test d.conditional_sigma ≈ p.sigma
    @test all(c.max_scaled_error < 1e-6 for c in d.finite_differences)
    @test !any(d.active_lower) && !any(d.active_upper)
    reordered = CocoaMLEDiagnostics.point_diagnostics(model,data,(sigma=p.sigma,a=p.a))
    @test reordered.gradient ≈ d.gradient
    missing_data = (values=[1.,missing,3.],sigma=:sigma)
    missing_result = CocoaMLEDiagnostics.point_diagnostics(model,missing_data,(a=2.,sigma=1.))
    @test missing_result.observation_count == 2
    @test missing_result.information_factor' * missing_result.information_factor ≈ [8. 0.;0. 4.]
    # A true constrained optimum has a nonzero outward gradient at the bound.
    high_data = (values=[3.,4.,5.],sigma=:sigma)
    edge = CocoaMLEDiagnostics.point_diagnostics(model,high_data,(a=3.,sigma=sqrt(5/3)))
    @test edge.gradient[1] < -1
    @test edge.projected_gradient_inf < 1e-12
    @test edge.active_upper == [true,false]
    # Inward descent at the same bound must remain visible.
    wrong = CocoaMLEDiagnostics.point_diagnostics(model,data,(a=3.,sigma=1.))
    @test wrong.gradient[1] > 1
    @test wrong.projected_gradient_inf > 1
    # A redundant pair remains rank deficient after the noise row is added.
    redundant = FitModel((p,d)->fill(p.a*p.b,3);initial=(a=1.,b=1.,sigma=1.),
        lower=(a=.1,b=.1,sigma=.1),upper=(a=3.,b=3.,sigma=3.))
    r = CocoaMLEDiagnostics.point_diagnostics(redundant,data,(a=1.,b=2.,sigma=sqrt(2/3)))
    @test r.information_rank == 2
    @test isinf(r.information_condition)
    @test norm(r.information_factor*r.weakest_direction) < 1e-12
    @test_throws ArgumentError CocoaMLEDiagnostics.point_diagnostics(model,data,(a=4.,sigma=1.))
end
