"""
    cocoa_observations(initial, times, observations; t0=0.0)

Prepare the source likelihood data from physical mg/g-pulp states. `observations`
has times in rows and `COCOA_STATES` in columns. Normalize using column maxima
including the fixed initial state, which is excluded from the likelihood.
The flattened values are state-major (all times for G, then F, etc.).
`sigma=:sigma` requests a shared inferred SD in these normalized coordinates.
"""
function cocoa_observations(initial::AbstractVector, times::AbstractVector,
                            observations::AbstractMatrix; t0=0.0)
    length(initial) == 8 && size(observations) == (length(times),8) ||
        throw(DimensionMismatch("cocoa observations must be times by eight states"))
    t0 isa Real && isfinite(t0) && !isempty(times) &&
        all(x -> x isa Real && isfinite(x), times) && all(>(0), diff([t0;times])) ||
        throw(ArgumentError("observation times must increase strictly after finite t0"))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x >= 0, initial) ||
        throw(DomainError(initial, "initial cocoa concentrations must be finite and nonnegative"))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x), observations) ||
        throw(DomainError(observations, "cocoa observations must be finite"))
    scales = vec(maximum(vcat(initial', observations); dims=1))
    _cocoa_scales(scales)
    return (; t0, times=copy(times), initial=initial ./ scales, scales,
             values=vec(observations ./ scales'), sigma=:sigma)
end

"""
    cocoa_fit_model(; initial, lower, upper, abstol=1e-6, reltol=1e-6,
                     maxiters=100_000_000, solver=:tsit5,
                     maxiters_scope=:solve, max_step_attempts=500)

Build a `FitModel` for the source 24 kinetic coefficients plus inferred `sigma`.
Use `cocoa_observations` for its data. Bounds are explicit numerical choices for
MLE; `fit_bayes` uses the supplied priors without imposing those bounds.
The default Julia integrator is Tsit5. Use `solver=:rodas5p` for stiff parameter
proposals during optimization. This changes the numerical integration method,
not the equations or likelihood. Failed solves return `nothing` to the fit engine.

`:solve` preserves Julia's whole-trajectory limit on attempted steps. Choose
`maxiters_scope=:output_interval` to limit accepted steps between requested output
times, as in Stan's dense-output RK45 interface. This uses one continuous
integrator, without adding time stops or restarting its adaptive controller.
Rejected attempts have a separate `max_step_attempts` guard per accepted step
(500, matching Boost's dense-output retry guard). Integrators and error control
still differ from Stan; equal tolerances do not imply identical trajectories.
"""
function cocoa_fit_model(; initial::NamedTuple, lower::NamedTuple, upper::NamedTuple,
                         abstol=1e-6, reltol=1e-6, maxiters=100_000_000, solver=:tsit5,
                         maxiters_scope=:solve, max_step_attempts=500)
    names = (COCOA_PARAMETERS..., :sigma)
    Set(keys(initial)) == Set(names) || throw(ArgumentError("fit all 24 cocoa kinetics and sigma"))
    all(x -> x isa Real && isfinite(x) && x > 0, (abstol,reltol)) ||
        throw(DomainError((abstol,reltol), "solver tolerances must be positive"))
    maxiters isa Integer && !(maxiters isa Bool) && 0 < maxiters <= typemax(Int) ||
        throw(ArgumentError("maxiters must be a positive machine integer"))
    maxiters_scope in (:solve,:output_interval) ||
        throw(ArgumentError("maxiters_scope must be :solve or :output_interval"))
    max_step_attempts isa Integer && !(max_step_attempts isa Bool) &&
        0 < max_step_attempts <= typemax(Int) ||
        throw(ArgumentError("max_step_attempts must be a positive machine integer"))
    solver in (:tsit5, :rodas5p) || throw(ArgumentError("solver must be :tsit5 or :rodas5p"))
    algorithm = solver == :tsit5 ? Tsit5() : Rodas5P()
    initial = NamedTuple{names}(Tuple(initial[name] for name in names))
    function predict(parameters, data)
        problem = cocoa_problem(data.initial, (data.t0,last(data.times)), parameters)
        solution = if maxiters_scope == :solve
            solved = solve(problem,algorithm;abstol,reltol,maxiters,
                saveat=data.times,save_start=false,dense=false)
            successful_retcode(solved) ? solved : nothing
        else
            _solve_output_intervals(problem,algorithm,data.times;
                abstol,reltol,maxiters,max_step_attempts)
        end
        isnothing(solution) && return nothing
        return vec(permutedims(Array(solution)))
    end
    return FitModel(predict; initial, lower, upper)
end

# Stan RK45's Boost checked_stepper counts successful dense-output steps; its
# checked_observer resets the allowance after each interpolated output. Keep
# OrdinaryDiffEq's controller and cumulative iteration counter intact. Only its
# mutable attempt ceiling changes, to bound retries inside each step! call.
# on_interval is an internal diagnostic hook, also called for the failed interval.
function _solve_output_intervals(problem,algorithm,times;maxiters,max_step_attempts=500,
                                 on_interval=nothing,kwargs...)
    integrator = init(problem,algorithm;maxiters=Int(max_step_attempts),saveat=times,
        save_start=false,dense=false,kwargs...)
    for (interval,time) in enumerate(times)
        accepted = 0
        attempted_start = integrator.iter
        failure = :none
        while integrator.t < time
            if accepted == maxiters
                failure = :output_step_limit
                break
            end
            # Saturation avoids overflow even for an explicitly unbounded guard.
            integrator.opts.maxiters = integrator.iter +
                min(Int(max_step_attempts),typemax(Int)-integrator.iter)
            step!(integrator)
            if !(integrator.sol.retcode in (ReturnCode.Default,ReturnCode.Success))
                failure = integrator.sol.retcode == ReturnCode.MaxIters ?
                    :step_attempt_limit : :solver_failure
                break
            end
            accepted += 1
        end
        isnothing(on_interval) || on_interval((;interval,time,accepted,
            attempted=integrator.stats.naccept+integrator.stats.nreject-attempted_start,failure))
        failure == :none || return nothing
    end
    return integrator.sol
end

"""
    cocoa_priors()

The source's independent positive-truncated Normal(0.5,0.3) priors on normalized
kinetics and half-Cauchy(0,1) prior on shared sigma, with provenance. Support is
unbounded above. These differ from the generic literature log-normal helper.
"""
function cocoa_priors()
    names = (COCOA_PARAMETERS..., :sigma)
    entries = map(names) do name
        distribution = name == :sigma ? Distributions.truncated(Distributions.Cauchy(0.,1.); lower=0.) :
            Distributions.truncated(Distributions.Normal(0.5,0.3); lower=0.)
        provenance = (source="doi:10.1098/rsos.180964",
            location="Equation (2.13); data/raw/cocoa/CocoaModel-master/simulation/model.stan:117-144",
            coordinates="normalized by trial state maxima", parameter=String(name))
        (; distribution, provenance)
    end
    return NamedTuple{names}(entries)
end
