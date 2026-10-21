using ForwardDiff

@testset "Secondary growth responses" begin
    ctmi(t) = γ_T_ctmi(t; Tmin=0, Topt=30, Tmax=40)
    ganzle(t) = γ_T_ganzle(t; a=2, b=2, c=0.2, Tref=40)
    cardinal(h) = γ_pH_cardinal(h; pHmin=3, pHopt=5, pHmax=7)
    acid(x) = γ_acid(x; MIC=100, α=2)
    salt(x) = γ_salt(x; max=4, α=1)
    for (f, grid, optimum, outside) in (
        (ctmi, -10:0.1:50, 30, [-10, 0, 40, 50]),
        (ganzle, -100:0.1:50, 30, [40, 50]),
        (cardinal, 0:0.01:10, 5, [0, 3, 7, 10]),
        (acid, 0:0.1:150, 0, [100, 150]),
        (salt, 0:0.01:10, 0, [4, 10]),
    )
        @test f(optimum) ≈ 1
        @test all(x -> 0 <= f(x) <= 1, grid)
        @test all(x -> iszero(f(x)), outside)
    end
    # Independent hand calculations for synthetic parameter sets.
    @test ctmi(20) ≈ 2/3
    @test ganzle(20) ≈ 4exp(-2)
    @test cardinal(4) ≈ 3/4
    @test acid(50) ≈ 3/4
    @test salt(2) ≈ 1/2
    @test γ_T_ganzle(20; a=1e250, b=2, c=0.2, Tref=40) ≈ ganzle(20)
    @test ganzle(-1e6) == 0
    @test all(t -> 0 <= ganzle(t) <= 1, 30 .+ (-100:100) .* 1e-9)
    @test γ_T_ctmi(10; Tmin=0, Topt=20, Tmax=40) ≈ 3/4
    @test ctmi(big"20") isa BigFloat

    @testset "Differentiation" begin
        for (f, x) in ((ctmi, 20.0), (ganzle, 20.0), (cardinal, 4.0), (acid, 50.0), (salt, 2.0))
            difference = (f(x + 1e-5) - f(x - 1e-5)) / 2e-5
            @test ForwardDiff.derivative(f, x) ≈ difference rtol=1e-6 atol=1e-10
        end
        @test ForwardDiff.derivative(ganzle, 20.0) ≈ ganzle(20) * (0.2 - 2/20)
        @test ForwardDiff.derivative(ctmi, 30.0) ≈ 0 atol=1e-12
        @test ForwardDiff.derivative(cardinal, 5.0) ≈ 0 atol=1e-12
        @test ForwardDiff.derivative(ctmi, 45.0) == 0
        @test all(isfinite, ForwardDiff.gradient(p -> γ_T_ctmi(20; Tmin=p[1], Topt=p[2], Tmax=p[3]), [0., 30., 40.]))
        @test all(isfinite, ForwardDiff.gradient(p -> γ_T_ganzle(20; a=p[1], b=p[2], c=p[3], Tref=p[4]), [2., 2., 0.2, 40.]))
        @test ForwardDiff.derivative(a -> γ_acid(0; MIC=100, α=a), 2.0) == 0
    end
    @testset "Invalid domains" begin
        @test_throws DomainError γ_T_ctmi(20; Tmin=0, Topt=10, Tmax=40)
        @test_throws DomainError γ_T_ctmi(20; Tmin=30, Topt=10, Tmax=40)
        @test_throws DomainError γ_pH_cardinal(5; pHmin=3, pHopt=7, pHmax=7)
        @test_throws DomainError γ_T_ganzle(20; a=1, b=2, c=-0.2, Tref=40)
        @test_throws DomainError γ_acid(-1; MIC=100)
        @test_throws DomainError γ_acid(1; MIC=0)
        @test_throws DomainError γ_salt(1; max=4, α=0)
        @test_throws DomainError ctmi(NaN)
        @test_throws DomainError cardinal(Inf)
    end
end
