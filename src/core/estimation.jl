"""
    FitModel(predict; initial, lower, upper)

Define an observation model and positive scalar parameters as NamedTuples with
matching names. `predict(parameters, data)` returns one prediction per entry in
`data.values`, in the same units/observation scale. Initial values must lie
strictly inside finite positive bounds. Fixed parameters belong in the closure.

For an expected simulation failure, return `nothing`. Programming errors are
not swallowed. Promote ODE initial states to the parameter scalar type for AD.
"""
struct FitModel{F,I,L,U}
    predict::F
    initial::I
    lower::L
    upper::U
    function FitModel(predict; initial::NamedTuple, lower::NamedTuple, upper::NamedTuple)
        names = keys(initial)
        isempty(names) && throw(ArgumentError("at least one fitted parameter is required"))
        Set(names) == Set(keys(lower)) == Set(keys(upper)) ||
            throw(ArgumentError("initial values and bounds must have matching parameter names"))
        for name in names
            vals = (lower[name], initial[name], upper[name])
            _finite_real(vals...)
            any(x -> x isa Bool, vals) && throw(ArgumentError("parameters must be numerical scalars, not Bool"))
            0 < vals[1] < vals[2] < vals[3] ||
                throw(DomainError(vals, "$name requires 0 < lower < initial < upper"))
        end
        lo = NamedTuple{names}(Tuple(lower[name] for name in names))
        hi = NamedTuple{names}(Tuple(upper[name] for name in names))
        return new{typeof(predict),typeof(initial),typeof(lo),typeof(hi)}(predict, initial, lo, hi)
    end
end

struct _PredictionFailure <: Exception
    message::String
end
Base.showerror(io::IO, err::_PredictionFailure) = print(io, err.message)

_log_parameters(parameters) = log.(collect(promote(float.(values(parameters))...)))
_decode_parameters(model, logarithms) = NamedTuple{keys(model.initial)}(Tuple(exp.(logarithms)))

function _observations(data)
    hasproperty(data, :values) && hasproperty(data, :sigma) ||
        throw(ArgumentError("data requires values and sigma (known SDs or a fitted parameter name)"))
    data.values isa AbstractVector || throw(ArgumentError("data.values must be a vector"))
    indices = findall(!ismissing, data.values)
    isempty(indices) && throw(ArgumentError("at least one nonmissing observation is required"))
    y = collect(data.values[indices])
    _finite_real(y...)
    if data.sigma isa Symbol
        return (; indices, y, sigma=data.sigma, normalizer=nothing)
    end
    sigma = if data.sigma isa Real
        fill(data.sigma, length(indices))
    elseif data.sigma isa AbstractVector
        length(data.sigma) == length(data.values) || throw(DimensionMismatch("sigma must match values"))
        collect(data.sigma[indices])
    else
        throw(ArgumentError("sigma must be a scalar, a matching vector, or a fitted parameter Symbol"))
    end
    _finite_real(sigma...)
    all(>(0), sigma) || throw(DomainError(sigma, "observed standard deviations must be positive"))
    return (; indices, y, sigma, normalizer=sum(log, sigma) + length(indices)*log(2pi)/2)
end

function _noise_sigma(observations, parameters)
    observations.sigma isa Symbol || return observations.sigma
    name = observations.sigma
    hasproperty(parameters, name) || throw(ArgumentError("noise parameter $name is missing from fitted parameters"))
    sigma = getproperty(parameters, name)
    sigma isa Real && isfinite(sigma) && sigma > 0 ||
        throw(_PredictionFailure("invalid fitted noise scale"))
    return fill(sigma, length(observations.indices))
end

function _gaussian_loss(model, data, observations, parameters)
    sigma = _noise_sigma(observations, parameters)
    prediction = _fit_prediction(model, parameters, data, observations.indices)
    residual = (prediction[observations.indices] - observations.y) ./ sigma
    return sum(abs2, residual)/2 + sum(log, sigma) + length(residual)*log(2pi)/2
end

function _fit_prediction(model, parameters, data, indices)
    all(x -> x isa Real && isfinite(x), values(parameters)) ||
        throw(_PredictionFailure("nonfinite parameter proposal"))
    prediction = model.predict(parameters, data)
    isnothing(prediction) && throw(_PredictionFailure("prediction returned nothing (simulation failed)"))
    prediction isa AbstractVector || throw(ArgumentError("predict must return a vector or nothing"))
    length(prediction) == length(data.values) || throw(DimensionMismatch("prediction must match data.values"))
    all(i -> prediction[i] isa Real && isfinite(prediction[i]), indices) ||
        throw(_PredictionFailure("nonfinite or missing prediction at an observed row"))
    return prediction
end

"""
    gaussian_nll(model, data, parameters)

Independent Gaussian negative log likelihood, including normalization. `data`
has `values` and `sigma` (scalar/per-row known SDs, or a Symbol naming a fitted
shared SD parameter), plus any
inputs used by `model.predict`. Missing observations are omitted, never zeroed.
This likelihood does not impose parameter bounds or a prior.
"""
function gaussian_nll(model::FitModel, data, parameters::NamedTuple)
    Set(keys(parameters)) == Set(keys(model.initial)) || throw(ArgumentError("parameter names differ from model"))
    observations = _observations(data)
    return _gaussian_loss(model, data, observations, parameters)
end

function _fit_starts(model, starts, rng)
    names = keys(model.initial)
    lo, hi = _log_parameters(model.lower), _log_parameters(model.upper)
    if starts isa Integer && !(starts isa Bool)
        starts > 0 || throw(ArgumentError("starts must be positive"))
        points = [model.initial]
        # The vector must also accommodate promoted floating-point random starts.
        result = NamedTuple[points...]
        for _ in 2:starts
            push!(result, _decode_parameters(model, lo + rand(rng, length(lo)).*(hi-lo)))
        end
    elseif starts isa AbstractVector && !isempty(starts)
        result = NamedTuple[]
        for point in starts
            point isa NamedTuple && Set(keys(point)) == Set(names) ||
                throw(ArgumentError("explicit starts must be NamedTuples with the model's parameter names"))
            push!(result, NamedTuple{names}(Tuple(point[name] for name in names)))
        end
    else
        throw(ArgumentError("starts must be a positive integer or a nonempty vector of NamedTuples"))
    end
    for point in result, name in names
        _finite_real(point[name])
        model.lower[name] < point[name] < model.upper[name] ||
            throw(DomainError(point[name], "starting $name must lie strictly inside its bounds"))
    end
    return result, lo, hi
end

function _fit_uncertainty(model, data, observations, z, lo, hi)
    # Expected Gaussian information: mean sensitivity plus (when estimated)
    # 2n * grad(log(sigma))*grad(log(sigma))'. Freeze sigma in the mean Jacobian.
    sigma = _noise_sigma(observations, _decode_parameters(model, z))
    J = ForwardDiff.jacobian(z) do point
        prediction = _fit_prediction(model, _decode_parameters(model, point), data, observations.indices)
        return prediction[observations.indices] ./ sigma
    end
    if observations.sigma isa Symbol
        noise_gradient = ForwardDiff.gradient(point -> log(_decode_parameters(model, point)[observations.sigma]), z)
        J = vcat(J, sqrt(2length(observations.indices)) .* noise_gradient')
    end
    active = any(min.(z-lo, hi-z) .<= 1e-5 .* max.(one.(z), hi-lo))
    if !all(isfinite, J)
        return (; status=:nonfinite_derivatives, intervals=nothing, covariance_log=nothing,
                 rank=nothing, singular_values=nothing)
    end
    decomposition = svd(J)
    cutoff = sqrt(eps(eltype(decomposition.S))) * maximum(decomposition.S)
    numerical_rank = count(>(cutoff), decomposition.S)
    status = active ? :active_bounds : numerical_rank < length(z) ? :rank_deficient : :ok
    if status != :ok
        return (; status, intervals=nothing, covariance_log=nothing,
                 rank=numerical_rank, singular_values=decomposition.S)
    end
    covariance = decomposition.V * Diagonal(inv.(decomposition.S.^2)) * decomposition.V'
    # 97.5th percentile of the standard normal: two-sided local 95% intervals.
    halfwidth = 1.959963984540054 .* sqrt.(diag(covariance))
    low, high = z-halfwidth, z+halfwidth
    crosses_bounds = any(low .<= lo) || any(high .>= hi)
    intervals = crosses_bounds ? nothing : NamedTuple{keys(model.initial)}(
        Tuple((lower=exp(a), upper=exp(b)) for (a,b) in zip(low, high)))
    return (; status=crosses_bounds ? :interval_crosses_bounds : :ok, intervals,
             covariance_log=covariance, rank=numerical_rank, singular_values=decomposition.S)
end

"""
    fit_mle(model, data; starts=50, rng=Random.default_rng(), maxiters=1000,
            outer_iterations=20, gtol=1e-7, maxtime=nothing, on_start_complete=nothing)

Multistart bounded Gaussian MLE using Optimization.jl, Fminbox(LBFGS), and
ForwardDiff. Optimize logarithms of positive parameters. The first random start
is `model.initial`; remaining starts are uniform in log bounds. Alternatively,
provide a vector of explicit NamedTuple starts. Pass a seeded RNG for repetition.

`maxiters` limits each inner LBFGS solve; `outer_iterations` limits barrier updates.
`maxtime`, when supplied, is a positive optimizer time budget in seconds per
start. It is checked between optimizer iterations, so a costly model evaluation
or compilation can overrun it. Pair it with model-specific solver step limits.
Each run records elapsed seconds and the optimizer's explicit time-limit flag.
Native optimizer success is recorded separately from acceptance. A fit also
requires an independently recomputed projected log-parameter gradient with
infinity norm at most `gtol`: `z-clamp(z-gradient,lo,hi)`. This handles active
bounds without accepting termination from small steps alone. A failed gradient
evaluation preserves the finite parameter point but cannot pass stationarity.
An optional `on_start_complete(index, run)` callback receives every completed
start, including failed predictions, before the next start begins. Use it to
persist progress; callback errors propagate rather than being treated as fit failures.
Return the best converged parameters, predictions, residuals (observed minus
predicted), standardized residuals, every run, and local 95% uncertainty. Failed
starts remain in `runs`. If none converges, `converged=false`, `best_start=nothing`,
and no estimates or intervals are supplied. Model/programming errors propagate.

Intervals use expected Gaussian information in log parameters, including the
variance contribution when a shared SD is fitted via `data.sigma=:name`.
They are withheld at active bounds, deficient rank, or crossing bounds;
they are not posterior intervals or a guarantee of coverage for nonlinear models.
"""
function fit_mle(model::FitModel, data; starts=50, rng=Random.default_rng(),
                 maxiters=1000, outer_iterations=20, gtol=1e-7, maxtime=nothing,
                 on_start_complete=nothing)
    maxiters isa Integer && maxiters > 0 || throw(ArgumentError("maxiters must be positive"))
    outer_iterations isa Integer && outer_iterations > 0 || throw(ArgumentError("outer_iterations must be positive"))
    _finite_real(gtol)
    gtol > 0 || throw(ArgumentError("gtol must be positive"))
    if !isnothing(maxtime)
        maxtime isa Real && !(maxtime isa Bool) && isfinite(maxtime) && maxtime > 0 ||
            throw(ArgumentError("maxtime must be positive finite seconds or nothing"))
    end
    observations = _observations(data)
    _noise_sigma(observations, model.initial)
    points, lo, hi = _fit_starts(model, starts, rng)
    function loss(z, unused)
        return _gaussian_loss(model, data, observations, _decode_parameters(model, z))
    end
    objective = Optimization.OptimizationFunction(loss, Optimization.AutoForwardDiff())
    runs = NamedTuple[]
    for (index, point) in enumerate(points)
        started = time_ns()
        try
            problem = Optimization.OptimizationProblem(objective, _log_parameters(point); lb=lo, ub=hi)
            solution = Optimization.solve(problem, OptimizationOptimJL.Fminbox(OptimizationOptimJL.LBFGS());
                maxiters=outer_iterations, local_maxiters=maxiters, maxtime,
                g_abstol=gtol, outer_g_abstol=gtol)
            nll = loss(solution.u, nothing)
            time_limit_reached = solution.original.stopped_by.time_limit
            optimizer_converged = successful_retcode(solution)
            gradient_error = nothing
            projected_gradient_inf = try
                gradient = ForwardDiff.gradient(z -> loss(z,nothing),solution.u)
                maximum(abs.(solution.u - clamp.(solution.u-gradient,lo,hi)))
            catch err
                err isa _PredictionFailure || rethrow()
                gradient_error = sprint(showerror,err)
                missing
            end
            stationarity_passed = !ismissing(projected_gradient_inf) &&
                isfinite(projected_gradient_inf) && projected_gradient_inf <= gtol
            converged = optimizer_converged && !time_limit_reached && stationarity_passed &&
                        isfinite(nll) && all(lo .<= solution.u .<= hi)
            push!(runs, (; initial=point, parameters=_decode_parameters(model, solution.u),
                         objective=nll, converged, retcode=string(solution.retcode),
                         optimizer_converged, stationarity_passed, projected_gradient_inf,
                         gradient_tolerance=gtol,
                         iterations=solution.stats.iterations, error=gradient_error,
                         elapsed_seconds=(time_ns()-started)/1e9, time_limit_reached))
        catch err
            err isa _PredictionFailure || rethrow()
            push!(runs, (; initial=point, parameters=nothing, objective=Inf,
                         converged=false, retcode="PredictionFailure", iterations=missing,
                         optimizer_converged=false, stationarity_passed=false,
                         projected_gradient_inf=missing, gradient_tolerance=gtol,
                         error=sprint(showerror, err), elapsed_seconds=(time_ns()-started)/1e9,
                         time_limit_reached=missing))
        end
        isnothing(on_start_complete) || on_start_complete(index, last(runs))
    end
    successful = findall(run -> run.converged, runs)
    if isempty(successful)
        return (; converged=false, parameters=nothing, objective=Inf, best_start=nothing,
                 runs, predictions=nothing, residuals=nothing, standardized_residuals=nothing,
                 observed_indices=observations.indices, uncertainty=(status=:no_converged_fit, intervals=nothing))
    end
    best = successful[argmin([runs[i].objective for i in successful])]
    parameters = runs[best].parameters
    predictions = _fit_prediction(model, parameters, data, observations.indices)
    residuals = [ismissing(y) ? missing : y-prediction for (y,prediction) in zip(data.values, predictions)]
    sigma = data.sigma isa Symbol ? fill(parameters[data.sigma], length(residuals)) :
        data.sigma isa Real ? fill(data.sigma, length(residuals)) : data.sigma
    standardized_residuals = [ismissing(r) ? missing : r/s for (r,s) in zip(residuals, sigma)]
    uncertainty = _fit_uncertainty(model, data, observations, _log_parameters(parameters), lo, hi)
    return (; converged=true, parameters, objective=runs[best].objective, best_start=best,
             runs, predictions, residuals, standardized_residuals,
             observed_indices=observations.indices, uncertainty)
end
