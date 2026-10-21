@testset "Acid speciation and charge balance" begin
    @test undissociated_fraction(5, 5) == 1/2
    @test undissociated_fraction(6, 5) ≈ 1/11
    @test undissociated_fraction(4, 5) ≈ 10/11
    @test undissociated_fraction(1000, 5) == 0
    @test undissociated_fraction(-1000, 5) == 1
    @test undissociated_acid(20, 5; pKa=5) == 10
    @test henderson_hasselbalch(10, 100; pKa=5) ≈ 6
    @test_throws DomainError henderson_hasselbalch(0, 1; pKa=5)
    @test_throws DomainError undissociated_acid(-1, 5; pKa=5)

    # Synthetic buffer: no numerical food-specific assumptions.
    buffer = (total=20.0, pKa=5.0, pH0=6.0, pKw=14.0, acid_pKas=(acid=5.0,))
    @test pH((acid=0.0,), buffer) ≈ 6 atol=1e-10
    @test pH((acid=0,), (total=20, pKa=5, pH0=6, pKw=14, acid_pKas=(acid=5,))) ≈ 6 atol=1e-10
    @test all(<(0), diff([pH((acid=c,), buffer) for c in 0:5:100]))
    capacity = log(10) * buffer.total * (10/11) * (1/11)
    from_capacity = buffer_from_capacity(; pH0=6., capacity, pKa=5., pKw=14., acid_pKas=(acid=5.,))
    @test from_capacity.total ≈ 20
    @test pH((acid=30.,), from_capacity) ≈ pH((acid=30.,), buffer)

    # Invert the mass/charge equations by hand at target pHs.
    H0, OH0 = 10.0^(3-buffer.pH0), 10.0^(3+buffer.pH0-buffer.pKw)
    fixed = buffer.total / (1+10.0^(buffer.pKa-buffer.pH0)) + OH0 - H0
    for target in (3.5, 4.0, 4.5, 5.0, 5.5)
        H, OH = 10.0^(3-target), 10.0^(3+target-buffer.pKw)
        fraction = 1/(1+10.0^(5-target))
        concentration = (H + fixed - OH) / fraction - buffer.total
        @test pH((acid=concentration,), buffer) ≈ target atol=1e-10
    end

    # Unbuffered weak acid: compare with the quadratic approximation, water negligible.
    water = (total=0., pKa=5., pH0=7., pKw=14., acid_pKas=(acid=5.,))
    Ka, concentration = 0.01, 10.0 # mmol/L units
    H = (sqrt(Ka^2 + 4Ka*concentration) - Ka)/2
    # Bound the pH error from neglecting water: delta(H) <= Kw/H.
    water_correction_bound = exp10(6-water.pKw) / (H^2 * log(10))
    @test abs(pH((acid=concentration,), water) - (3-log10(H))) <= water_correction_bound
    mixture = merge(buffer, (acid_pKas=(first=4.0, second=5.0),))
    @test pH((first=10., second=10.), mixture) < pH((first=0., second=10.), mixture)
    @test pH((first=10., second=10.), mixture) ≈ pH((second=10., first=10.), mixture)

    @testset "Implicit differentiation" begin
        f(c) = pH((acid=c,), buffer)
        derivative = ForwardDiff.derivative(f, 30.)
        @test derivative < 0
        @test derivative ≈ (f(30.001)-f(29.999))/0.002 rtol=1e-7
        root = f(30.)
        H, OH = 10^(3-root), 10^(3+root-buffer.pKw)
        fraction = 1/(1+10^(5-root))
        expected = -fraction / (log(10)*(H+OH+50*fraction*(1-fraction)))
        @test derivative ≈ expected rtol=1e-10
        @test isfinite(ForwardDiff.derivative(b -> pH((acid=30.,), merge(buffer, (total=b,))), 20.))
    end
    @test_throws DomainError pH((acid=-1.,), buffer)
    @test_throws ArgumentError pH((other=1.,), buffer)
    @test_throws DomainError pH((acid=30.,), buffer; bracket=(8, 12))
    @test_throws ErrorException pH((acid=30.,), buffer; maxiters=1)
    @test_throws ArgumentError pH((acid=30.,), buffer; atol=0)
    @test_throws DomainError buffer_from_capacity(; pH0=6., capacity=-1., pKa=5., pKw=14., acid_pKas=(acid=5.,))
end
