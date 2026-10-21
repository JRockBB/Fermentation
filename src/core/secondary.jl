function _finite_real(values...)
    all(x -> x isa Real && isfinite(x), values) ||
        throw(DomainError(values, "inputs must be finite real numbers"))
end

"""
    γ_T_ctmi(T; Tmin, Topt, Tmax)

Rosso cardinal temperature response, normalized to one at `Topt` and zero at
and beyond the limits. Implements Rosso et al. (1995), equation 6. Requires
`Tmin < Topt < Tmax` and `Topt - Tmin >= Tmax - Topt`; lower optima produce a
pole rather than a bounded response. Temperatures are in degrees Celsius.
"""
function γ_T_ctmi(T::Real; Tmin::Real, Topt::Real, Tmax::Real)
    _finite_real(T, Tmin, Topt, Tmax)
    Tmin < Topt < Tmax || throw(DomainError((Tmin, Topt, Tmax), "cardinal temperatures must be ordered"))
    Topt - Tmin >= Tmax - Topt || throw(DomainError(Topt, "CTMI optimum must be at or above the midpoint"))
    (T <= Tmin || T >= Tmax) && return zero(T + Tmin + Topt + Tmax)
    numerator = (T - Tmax) * (T - Tmin)^2
    denominator = (Topt - Tmin) * ((Topt - Tmin) * (T - Topt) -
                                  (Topt - Tmax) * (Topt + Tmin - 2T))
    return numerator / denominator
end

"""
    γ_T_ganzle(T; a, b, c, Tref)

Normalized Gänzle et al. (1998) Table 1 model 1a: `a*x^b*exp(-c*x)`, with
`x = Tref - T`, positive `a, b, c`, and optimum `Tref - b/c`. The scale `a`
cancels on normalization. The paper uses a *negative* exponential, correcting
the initial plan. There is no finite lower cardinal limit in this model.
The response is zero at and above `Tref`, and tends to zero as T decreases.
"""
function γ_T_ganzle(T::Real; a::Real, b::Real, c::Real, Tref::Real)
    _finite_real(T, a, b, c, Tref)
    (a > 0 && b > 0 && c > 0) || throw(DomainError((a, b, c), "a, b, c must be positive"))
    T >= Tref && return zero(T + a + b + c + Tref)
    ratio = (Tref - T) * c / b
    # Log-space normalization avoids overflowing the unnormalized peak.
    return exp(b * (log(ratio) - (ratio - one(ratio))))
end

"""
    γ_pH_cardinal(pH; pHmin, pHopt, pHmax)

Rosso et al. (1995), equation 7. Returns one at the optimum and zero at and
outside the ordered cardinal bounds. This is also Gänzle (1998) Table 1 model 2.
"""
function γ_pH_cardinal(pH::Real; pHmin::Real, pHopt::Real, pHmax::Real)
    _finite_real(pH, pHmin, pHopt, pHmax)
    pHmin < pHopt < pHmax || throw(DomainError((pHmin, pHopt, pHmax), "cardinal pHs must be ordered"))
    (pH <= pHmin || pH >= pHmax) && return zero(pH + pHmin + pHopt + pHmax)
    numerator = (pH - pHmin) * (pH - pHmax)
    return numerator / (numerator - (pH - pHopt)^2)
end

function _inhibition(x, limit, α)
    _finite_real(x, limit, α)
    x >= 0 || throw(DomainError(x, "concentration must be nonnegative"))
    limit > 0 || throw(DomainError(limit, "inhibition limit must be positive"))
    α > 0 || throw(DomainError(α, "inhibition exponent must be positive"))
    x >= limit && return zero(x + limit + α)
    iszero(x) && return one(x + limit + α)
    return one(x / limit) - (x / limit)^α
end

"""
    γ_acid(HA; MIC, α=1)

Undissociated-acid inhibition `1 - (HA/MIC)^α`, bounded by zero above the MIC.
Gänzle (1998), Table 1 model 4. HA and MIC must share concentration units.
Negative HA is invalid; the optimum is zero acid. α=1 gives linear inhibition.
"""
γ_acid(HA::Real; MIC::Real, α::Real=1) = _inhibition(HA, MIC, α)

"""
    γ_salt(NaCl; max, α=1)

Power-law inhibition `1 - (NaCl/max)^α`, truncated at zero. NaCl and max are
percent w/w on total vegetable-plus-salt mass. This is an empirical extension
of Gänzle's model 4; salt and ionic strength are not interchangeable units.
The optimum is zero salt; this form does not model salt-stimulated growth.
"""
γ_salt(NaCl::Real; max::Real, α::Real=1) = _inhibition(NaCl, max, α)
