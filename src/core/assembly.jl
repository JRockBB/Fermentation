_unit_modifier(u, p, t) = one(t)

"""
    ModelSpec(; organisms, substrates, products, stoichiometry,
                kinetics=:monod, modifiers=nothing, log_populations=false,
                conserved_weights=nothing)

Small consumer-resource specification with states ordered `[X; S; P]` and one
growth flux per organism. Columns of `stoichiometry` are the changes per unit
population growth: the population block must be the identity, substrate entries
nonpositive, and product entries nonnegative. Units and yields must be consistent
with the chosen population basis; no CFU-to-mass conversion is inferred.

Kinetics are `:monod` or `:contois`, globally or one per organism. Each consumed
substrate contributes a multiplicative limitation term. Per-organism modifiers
are callable `(physical_state, parameters, time) -> factor` in [0,1]. Optional
conservation weights must annihilate every stoichiometric column.
"""
struct ModelSpec{O,S,P,M,K,F,W}
    organisms::O
    substrates::S
    products::P
    stoichiometry::M
    kinetics::K
    modifiers::F
    log_populations::Bool
    conserved_weights::W
end

function ModelSpec(; organisms, substrates, products, stoichiometry,
                   kinetics=:monod, modifiers=nothing, log_populations=false,
                   conserved_weights=nothing)
    organisms, substrates, products = Tuple(organisms), Tuple(substrates), Tuple(products)
    names = (organisms..., substrates..., products...)
    all(x -> x isa Symbol, names) || throw(ArgumentError("state names must be Symbols"))
    length(unique(names)) == length(names) || throw(ArgumentError("state names must be unique"))
    n, ns, np = length(organisms), length(substrates), length(products)
    n > 0 || throw(ArgumentError("at least one organism is required"))
    size(stoichiometry) == (n+ns+np, n) || throw(DimensionMismatch("stoichiometry must have one row per state and one column per organism"))
    _finite_real(stoichiometry...)
    for i in 1:n, j in 1:n
        stoichiometry[i,j] == (i == j ? 1 : 0) || throw(ArgumentError("population stoichiometry must be the identity"))
    end
    all(<=(0), stoichiometry[n+1:n+ns, :]) || throw(ArgumentError("substrate coefficients must be nonpositive"))
    all(>=(0), stoichiometry[n+ns+1:end, :]) || throw(ArgumentError("product coefficients must be nonnegative"))
    rates = kinetics isa Symbol ? ntuple(_ -> kinetics, n) : Tuple(kinetics)
    length(rates) == n && all(in((:monod, :contois)), rates) || throw(ArgumentError("choose Monod or Contois kinetics for each organism"))
    factors = isnothing(modifiers) ? ntuple(_ -> _unit_modifier, n) : Tuple(modifiers)
    length(factors) == n || throw(DimensionMismatch("one modifier per organism is required"))
    weights = isnothing(conserved_weights) ? nothing : collect(conserved_weights)
    if !isnothing(weights)
        length(weights) == length(names) || throw(DimensionMismatch("one conservation weight per state is required"))
        _finite_real(weights...)
        any(!iszero, weights) || throw(ArgumentError("conservation weights must not all be zero"))
        for j in 1:n
            residual = sum(weights[i] * stoichiometry[i,j] for i in eachindex(weights))
            scale = sum(abs(weights[i] * stoichiometry[i,j]) for i in eachindex(weights))
            abs(residual) <= 1e-10 * max(one(scale), scale) || throw(ArgumentError("stoichiometry violates declared conservation weights"))
        end
    end
    return ModelSpec(organisms, substrates, products, copy(stoichiometry), rates,
                     factors, Bool(log_populations), weights)
end

_state_count(spec) = length(spec.organisms) + length(spec.substrates) + length(spec.products)

"""Convert nonnegative physical concentrations/populations into solver coordinates."""
function initial_state(spec::ModelSpec, physical::AbstractVector)
    length(physical) == _state_count(spec) || throw(DimensionMismatch("initial state length differs from specification"))
    _finite_real(physical...)
    all(>=(0), physical) || throw(DomainError(physical, "initial states must be nonnegative"))
    state = float.(physical)
    if spec.log_populations
        n = length(spec.organisms)
        all(>(0), state[1:n]) || throw(DomainError(state[1:n], "log populations require positive initial populations"))
        for i in 1:n
            state[i] = log(state[i])
        end
    end
    return state
end

"""Decode solver coordinates to physical units; population logs are natural logs."""
function physical_state(spec::ModelSpec, state::AbstractVector)
    length(state) == _state_count(spec) || throw(DimensionMismatch("state length differs from specification"))
    physical = copy(state)
    if spec.log_populations
        for i in eachindex(spec.organisms)
            physical[i] = exp(state[i])
        end
    end
    return physical
end

function _growth_parameters(spec, p)
    hasproperty(p, :mu_max) && hasproperty(p, :Ks) || throw(ArgumentError("parameters require mu_max and Ks"))
    n, ns = length(spec.organisms), length(spec.substrates)
    length(p.mu_max) == n || throw(DimensionMismatch("mu_max must have one entry per organism"))
    size(p.Ks) == (ns, n) || throw(DimensionMismatch("Ks must be substrates × organisms"))
    _finite_real(p.mu_max..., p.Ks...)
    all(>=(0), p.mu_max) || throw(DomainError(p.mu_max, "mu_max must be nonnegative"))
    all(>(0), p.Ks) || throw(DomainError(p.Ks, "Ks values must be positive"))
end

"""
    build_rhs(spec) -> f!(du, u, p, t)

Assemble growth fluxes using `p.mu_max` (1/h, vector) and `p.Ks` (substrates ×
organisms). Monod terms are `S/(Ks+S)`; Contois terms are `S/(Ks*X+S)`.
Contois Ks therefore has substrate/population units, rather than substrate units.
Negative trial concentrations are treated as zero for rate evaluation; use
`fermentation_problem` to reject negative accepted solver states.
"""
function build_rhs(spec::ModelSpec)
    n, ns = length(spec.organisms), length(spec.substrates)
    function rhs!(du, u, p, t)
        length(u) == length(du) == _state_count(spec) || throw(DimensionMismatch("state length differs from specification"))
        _growth_parameters(spec, p)
        state = spec.log_populations ? physical_state(spec, u) : u
        fill!(du, zero(eltype(du)))
        for i in 1:n
            population = max(state[i], zero(state[i]))
            factor = spec.modifiers[i](state, p, t)
            _finite_real(factor)
            0 <= factor <= 1 || throw(DomainError(factor, "growth modifiers must lie in [0,1]"))
            rate = p.mu_max[i] * factor
            for j in 1:ns
                if spec.stoichiometry[n+j, i] < 0
                    substrate = max(state[n+j], zero(state[n+j]))
                    if iszero(substrate)
                        rate = zero(rate)
                        break
                    end
                    denominator = spec.kinetics[i] == :monod ? p.Ks[j,i] + substrate : p.Ks[j,i] * population + substrate
                    rate *= substrate / denominator
                end
            end
            flux = rate * population
            du[i] = spec.log_populations ? rate : flux
            for j in n+1:length(du)
                du[j] += spec.stoichiometry[j,i] * flux
            end
        end
        return nothing
    end
    return rhs!
end

"""Build an OrdinaryDiffEq ODEProblem from physical initial values, with positivity rejection."""
function fermentation_problem(spec::ModelSpec, physical, tspan, parameters)
    _growth_parameters(spec, parameters)
    length(tspan) == 2 && tspan[1] < tspan[2] || throw(ArgumentError("tspan must be an increasing pair of hours"))
    _finite_real(tspan...)
    start = spec.log_populations ? length(spec.organisms)+1 : 1
    invalid = (u,p,t) -> any(x -> !isfinite(x), u) || any(<(0), @view u[start:end])
    return ODEProblem(build_rhs(spec), initial_state(spec, physical), tspan, parameters;
                      isoutofdomain=invalid)
end
