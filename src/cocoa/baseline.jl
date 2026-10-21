"""Source state order: glucose, fructose, ethanol, lactate, acetate, yeast, LAB, AAB."""
const COCOA_STATES = (:G, :F, :E, :L, :A, :Xy, :Xl, :Xa)
const COCOA_PARAMETERS = (:mu1, :mu2, :mu3, :mu4, :mu5,
    :ks1, :ks2, :ks3, :ks4, :ks5, :k1, :k2, :k3,
    :yc1, :yc2, :yc3, :yc4, :yc5, :yc6, :yc7, :yc8, :yc9, :yc10, :yc11)

function _cocoa_parameters(p)
    all(name -> hasproperty(p, name), COCOA_PARAMETERS) ||
        throw(ArgumentError("cocoa parameters require all 24 names in COCOA_PARAMETERS"))
    values = map(name -> getproperty(p, name), COCOA_PARAMETERS)
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x >= 0, values) ||
        throw(DomainError(values, "cocoa kinetic parameters must be finite and nonnegative"))
    all(name -> getproperty(p, name) > 0, (:ks1, :ks2, :ks3, :ks4, :ks5)) ||
        throw(DomainError(p, "cocoa saturation constants must be positive"))
    return NamedTuple{COCOA_PARAMETERS}(values)
end

function _cocoa_scales(scales)
    length(scales) == 8 || throw(DimensionMismatch("one scale per cocoa state is required"))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x > 0, scales) ||
        throw(DomainError(scales, "cocoa scales must be finite and positive"))
    return Tuple(scales)
end

# Supplement S3, pp. 8-9. In particular yc3 uses fructose, not glucose.
function _cocoa_conversion(scales)
    G, F, E, L, A, Xy, Xl, Xa = _cocoa_scales(scales)
    return (one(G), one(G), one(G), one(G), one(G),
        G, F, G, E, L/Xa, inv(E), inv(L), inv(A^2),
        G/Xy, G/Xl, F/Xy, E/Xy, E/Xy, E/Xa, L/Xl, L/Xa, A/Xl, A/Xa, A/Xa)
end

"""
    cocoa_physical_parameters(scaled, scales)

Convert the 24 normalized kinetic coefficients to the source's physical mg/g-pulp
basis, using maxima in `COCOA_STATES` order. Extra fields such as observation
`sigma` are excluded: physical observation SDs are `sigma .* scales`.
No volumetric or CFU conversion is inferred. See supplement S3, pp. 8-9.
"""
function cocoa_physical_parameters(scaled, scales)
    p = _cocoa_parameters(scaled)
    return NamedTuple{COCOA_PARAMETERS}(values(p) .* _cocoa_conversion(scales))
end

"""Inverse of `cocoa_physical_parameters`; returns only the 24 kinetic coefficients."""
function cocoa_scaled_parameters(physical, scales)
    p = _cocoa_parameters(physical)
    return NamedTuple{COCOA_PARAMETERS}(values(p) ./ _cocoa_conversion(scales))
end

"""
    cocoa_rhs!(du, u, p, t)

Eight-state cocoa model, Moreno-Zambrano et al. (2018), equations (2.1)-(2.8),
matching archived `model.stan:26-33`. States follow `COCOA_STATES`; parameters
follow the source names in `COCOA_PARAMETERS`. The same equations apply to
physical or normalized coordinates when state and parameter scales agree.

Negative solver trial states are treated as zero when evaluating rates. The
Contois flux has the continuous value zero at lactate=AAB=0. `cocoa_problem`
validates parameters and rejects negative/nonfinite accepted states. No
temperature, pH, feeding, oxygen, or additional conservation law is introduced.
"""
function cocoa_rhs!(du, u, p, t)
    length(u) == length(du) == 8 || throw(DimensionMismatch("cocoa states have length eight"))
    G, F, E, L, A, Xy, Xl, Xa = ntuple(i -> max(u[i], zero(u[i])), 8)
    v1 = p.mu1 * G/(p.ks1 + G) * Xy
    v2 = p.mu2 * F/(p.ks2 + F) * Xy
    v3 = p.mu3 * G/(p.ks3 + G) * Xl
    v4 = p.mu4 * E/(p.ks4 + E) * Xa
    denominator = p.ks5*Xa + L
    v5 = iszero(denominator) ? zero(p.mu5*L*Xa) : p.mu5 * L/denominator * Xa
    v6 = p.k1 * Xy * E
    v7 = p.k2 * Xl * L
    v8 = p.k3 * Xa * A^2
    du[1] = -p.yc1*v1 - p.yc2*v3
    du[2] = -p.yc3*v2
    du[3] = p.yc4*v1 + p.yc5*v2 - p.yc6*v4
    du[4] = p.yc7*v3 - p.yc8*v5
    du[5] = p.yc9*v3 + p.yc10*v4 + p.yc11*v5
    du[6] = v1 + v2 - v6
    du[7] = v3 - v7
    du[8] = v4 + v5 - v8
    return nothing
end

"""
    cocoa_problem(initial, tspan, parameters)

Construct the source cocoa ODE with finite, nonnegative initial states and
positivity rejection. Time is in hours; use either normalized states and
normalized parameters, or mg/g-pulp states and converted physical parameters.
Initial state storage is promoted to the parameter types for ForwardDiff.
"""
function cocoa_problem(initial, tspan, parameters)
    p = _cocoa_parameters(parameters)
    length(initial) == 8 || throw(DimensionMismatch("cocoa initial state must have eight entries"))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x >= 0, initial) ||
        throw(DomainError(initial, "cocoa initial states must be finite and nonnegative"))
    length(tspan) == 2 && all(x -> x isa Real && isfinite(x), tspan) && tspan[1] < tspan[2] ||
        throw(ArgumentError("tspan must be an increasing finite pair of hours"))
    state_type = promote_type(map(x -> typeof(float(x)), Tuple(initial))..., map(typeof, values(p))...)
    state = state_type.(collect(initial))
    invalid = (u, p, t) -> any(x -> !isfinite(x) || x < 0, u)
    return ODEProblem(cocoa_rhs!, state, tspan, p; isoutofdomain=invalid)
end
