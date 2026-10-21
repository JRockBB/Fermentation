"""
    feeding_fraction(N)

Retained fraction `1/(1+2N)` for a 1:N:N starter:flour:water feed by mass.
Applying it to concentrations assumes a consistent mixture basis and hydration.
"""
function feeding_fraction(N::Real)
    _finite_real(N)
    N >= 0 || throw(DomainError(N, "feed ratio must be nonnegative"))
    ratio = float(N)
    return inv(one(ratio) + 2ratio)
end

function _feeding_inputs(spec, fraction, substrate_addition)
    _finite_real(fraction, substrate_addition...)
    0 <= fraction <= 1 || throw(DomainError(fraction, "retained fraction must lie in [0,1]"))
    length(substrate_addition) == length(spec.substrates) ||
        throw(DimensionMismatch("one addition per substrate is required"))
    all(>=(0), substrate_addition) || throw(DomainError(substrate_addition, "substrate additions must be nonnegative"))
    spec.log_populations && iszero(fraction) &&
        throw(DomainError(fraction, "complete dilution cannot be represented by log populations"))
    return nothing
end

"""
    dilute!(state, spec; fraction, substrate_addition)

Retain `fraction` of each physical population/concentration, then add substrate
amounts expressed as concentrations in the final mixture. For fresh-feed
concentrations `Sfresh`, use `substrate_addition = (1-fraction) .* Sfresh`.
Population coordinates follow `spec.log_populations`; no new organisms are added.
Validation and numeric conversion finish before mutating the input.
"""
function dilute!(state::AbstractVector, spec::ModelSpec; fraction::Real, substrate_addition)
    _feeding_inputs(spec, fraction, substrate_addition)
    length(state) == _state_count(spec) || throw(DimensionMismatch("state length differs from specification"))
    _finite_real(state...)
    n = length(spec.organisms)
    firstphysical = spec.log_populations ? n+1 : 1
    all(>=(0), @view state[firstphysical:end]) || throw(DomainError(state, "physical states must be nonnegative"))
    candidate = state .* fraction
    if spec.log_populations
        for i in 1:n
            candidate[i] = state[i] + log(fraction)
        end
    end
    for j in eachindex(spec.substrates)
        candidate[n+j] += substrate_addition[j]
    end
    _finite_real(candidate...)
    converted = convert.(eltype(state), candidate)
    copyto!(state, converted)
    return state
end

"""
    dilution_callback(spec, times; fraction, substrate_addition)

Create a reusable PresetTimeCallback for fixed feeds. Feed times must be unique
and inside the problem's increasing tspan, including either endpoint. The solver
stops at every feed; both sides of each interior jump are saved. A feed at the
initial time is applied during callback initialization.
"""
function dilution_callback(spec::ModelSpec, times; fraction::Real, substrate_addition)
    _feeding_inputs(spec, fraction, substrate_addition)
    stops = sort(collect(times))
    _finite_real(stops...)
    length(unique(stops)) == length(stops) || throw(ArgumentError("feed times must be unique"))
    addition = collect(substrate_addition)
    function initialize(callback, u, t, integrator)
        t0, tf = integrator.sol.prob.tspan
        t0 < tf || throw(ArgumentError("feeding requires an increasing tspan"))
        all(t0 <= stop <= tf for stop in stops) || throw(ArgumentError("feed times must lie inside tspan"))
        derivative_discontinuity!(integrator, true)
        return nothing
    end
    affect! = integrator -> dilute!(integrator.u, spec; fraction, substrate_addition=addition)
    return PresetTimeCallback(stops, affect!; initialize, save_positions=(true, true))
end

"""
    TemperatureSchedule(times, temperatures; interpolation=:step)

Callable temperature input in hours and °C, with strictly increasing knot times.
`:step` is right-continuous; `:linear` interpolates between knots. Both hold the
endpoint temperatures outside the schedule. Pass `temperature_stops(schedule,
tspan)` to the solver's `tstops` and `d_discontinuities` options.
"""
struct TemperatureSchedule{T,V}
    times::T
    temperatures::V
    interpolation::Symbol
    function TemperatureSchedule(times, temperatures; interpolation=:step)
        knots, values = collect(times), collect(temperatures)
        isempty(knots) && throw(ArgumentError("temperature schedule must not be empty"))
        length(knots) == length(values) || throw(DimensionMismatch("one temperature per time is required"))
        _finite_real(knots..., values...)
        all(knots[i] < knots[i+1] for i in 1:length(knots)-1) ||
            throw(ArgumentError("temperature times must be strictly increasing"))
        interpolation in (:step, :linear) || throw(ArgumentError("interpolation must be :step or :linear"))
        return new{typeof(knots),typeof(values)}(knots, values, interpolation)
    end
end

function (schedule::TemperatureSchedule)(t::Real)
    _finite_real(t)
    t <= first(schedule.times) && return first(schedule.temperatures)
    t >= last(schedule.times) && return last(schedule.temperatures)
    i = searchsortedlast(schedule.times, t)
    schedule.interpolation == :step && return schedule.temperatures[i]
    fraction = (t-schedule.times[i]) / (schedule.times[i+1]-schedule.times[i])
    return (1-fraction)*schedule.temperatures[i] + fraction*schedule.temperatures[i+1]
end

"""Return temperature knots strictly inside an increasing solver time span."""
function temperature_stops(schedule::TemperatureSchedule, tspan)
    length(tspan) == 2 || throw(ArgumentError("tspan must be a pair of hours"))
    _finite_real(tspan...)
    tspan[1] < tspan[2] || throw(ArgumentError("tspan must be increasing"))
    return filter(t -> tspan[1] < t < tspan[2], schedule.times)
end
