"""Fraction of a monoprotic acid present as HA: `1/(1 + 10^(pH-pKa))`."""
function undissociated_fraction(pH::Real, pKa::Real)
    _finite_real(pH, pKa)
    difference = pH - pKa
    if difference >= 0
        ratio = exp10(-difference)
        return ratio / (one(ratio) + ratio)
    end
    ratio = exp10(difference)
    return inv(one(ratio) + ratio)
end

"""Undissociated concentration, in the same units as the nonnegative total."""
function undissociated_acid(total::Real, pH::Real; pKa::Real)
    _finite_real(total)
    total >= 0 || throw(DomainError(total, "total acid must be nonnegative"))
    return total * undissociated_fraction(pH, pKa)
end

"""
    henderson_hasselbalch(acid, base; pKa)

Compute `pKa + log10(base/acid)` for positive HA and A⁻ amounts in matching units.
This fallback assumes a known acid/base pair; total acid alone is insufficient.
"""
function henderson_hasselbalch(acid::Real, base::Real; pKa::Real)
    _finite_real(acid, base, pKa)
    acid > 0 && base > 0 || throw(DomainError((acid, base), "both amounts must be positive"))
    return pKa + log10(base) - log10(acid)
end

"""
    buffer_from_capacity(; pH0, capacity, pKa, pKw, acid_pKas)

Construct buffer parameters from its initial pH and local buffer-only capacity
in mmol/L per pH unit. `total = capacity/(ln(10)*f*(1-f))`, where f is the
initial dissociated fraction. `acid_pKas` is a NamedTuple of fermentation-acid
pKas. `pKw` must be supplied for the intended temperature and activity convention.
"""
function buffer_from_capacity(; pH0::Real, capacity::Real, pKa::Real, pKw::Real,
                              acid_pKas::NamedTuple)
    _finite_real(pH0, capacity, pKa, pKw, values(acid_pKas)...)
    capacity >= 0 || throw(DomainError(capacity, "capacity must be nonnegative"))
    fraction = undissociated_fraction(pKa, pH0)
    total = iszero(capacity) ? zero(capacity) : capacity / (log(10) * fraction * (1-fraction))
    isfinite(total) || throw(DomainError(capacity, "buffer capacity is ill-conditioned at this pH"))
    return (; total, pKa, pH0, pKw, acid_pKas)
end

function _buffer_inputs(acids, buffer)
    all(name -> hasproperty(buffer, name), (:total, :pKa, :pH0, :pKw, :acid_pKas)) ||
        throw(ArgumentError("buffer requires total, pKa, pH0, pKw, and acid_pKas"))
    buffer.acid_pKas isa NamedTuple || throw(ArgumentError("acid_pKas must be a NamedTuple"))
    Set(keys(acids)) == Set(keys(buffer.acid_pKas)) ||
        throw(ArgumentError("acid concentrations and acid_pKas must have matching names"))
    _finite_real(buffer.total, buffer.pKa, buffer.pH0, buffer.pKw,
                 values(acids)..., values(buffer.acid_pKas)...)
    buffer.total >= 0 || throw(DomainError(buffer.total, "buffer total must be nonnegative"))
    all(>=(0), values(acids)) || throw(DomainError(acids, "acid totals must be nonnegative"))
end

function _buffer_charge(h, acids, buffer, fixed_charge)
    # pH is defined using mol/L; 10^(3-pH) expresses [H+] directly in mmol/L.
    hydrogen = exp10(3-h)
    hydroxide = exp10(3+h-buffer.pKw)
    fraction = undissociated_fraction(buffer.pKa, h)
    residual = hydrogen + fixed_charge - hydroxide - buffer.total * fraction
    slope = -hydrogen - hydroxide - buffer.total * fraction * (1-fraction)
    for name in keys(acids)
        fraction = undissociated_fraction(buffer.acid_pKas[name], h)
        residual -= acids[name] * fraction
        slope -= acids[name] * fraction * (1-fraction)
    end
    return residual, log(10) * slope
end

"""
    pH(acids::NamedTuple, buffer_params::NamedTuple; bracket=(-2, 16), atol=1e-10, maxiters=200)

Solve electroneutrality for added monoprotic acids (mmol/L) in one effective
monoprotic buffer. Parameters are `total` (mmol/L), `pKa`, `pH0` (initial pH
before adding these acids), `pKw`, and `acid_pKas` with the same names as acids.
Use `buffer_from_capacity` when a local buffer capacity is known instead of total.

The fixed charge is determined from the initial buffer and water equilibrium.
Acids are added as neutral HA; their salts require a different charge balance.
Bisection finds the unique root; two Newton refinements propagate implicit
derivatives through the root for ForwardDiff. Bracketing/convergence failures
raise errors. This is an ideal-solution engine, not a calibrated food model.
"""
function pH(acids::NamedTuple, buffer::NamedTuple; bracket=(-2, 16), atol=1e-10, maxiters=200)
    _buffer_inputs(acids, buffer)
    length(bracket) == 2 || throw(ArgumentError("bracket must contain two pH values"))
    _finite_real(bracket..., atol)
    lo, hi = float.(bracket)
    lo < hi || throw(ArgumentError("bracket must be increasing"))
    atol > 0 || throw(ArgumentError("atol must be positive"))
    maxiters isa Integer && maxiters > 0 || throw(ArgumentError("maxiters must be a positive integer"))
    fixed_charge = buffer.total * undissociated_fraction(buffer.pKa, buffer.pH0) +
                   exp10(3+buffer.pH0-buffer.pKw) - exp10(3-buffer.pH0)
    first(_buffer_charge(lo, acids, buffer, fixed_charge)) >= 0 &&
        first(_buffer_charge(hi, acids, buffer, fixed_charge)) <= 0 ||
        throw(DomainError(bracket, "charge-balance root lies outside the bracket"))

    converged = false
    for _ in 1:maxiters
        midpoint = lo + (hi-lo)/2
        if hi-lo <= atol || midpoint == lo || midpoint == hi
            converged = true
            break
        end
        residual, _ = _buffer_charge(midpoint, acids, buffer, fixed_charge)
        if residual > 0
            lo = midpoint
        else
            hi = midpoint
        end
    end
    converged || error("Charge-balance solve did not converge within $maxiters iterations")
    root = lo + (hi-lo)/2
    # Do not skip these when the primal residual is zero: its tangent may not be.
    for _ in 1:2
        residual, slope = _buffer_charge(root, acids, buffer, fixed_charge)
        root -= residual / slope
    end
    isfinite(root) || error("Nonfinite charge-balance result")
    return root
end
