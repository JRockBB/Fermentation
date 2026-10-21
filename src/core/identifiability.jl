"""
    local_sensitivities(model, data, parameters=model.initial;
                        scale=:log, standardized=true, rtol=sqrt(eps(Float64)))

Differentiate observed predictions with ForwardDiff. Columns follow the model's
parameter names; rows follow nonmissing observations. By default derivatives are
with respect to log parameters and divided by measurement sigma. `scale=:parameter`
uses physical parameter units. Return the Jacobian, column norms, singular values,
numerical rank, condition number, and locally unresolved parameter directions.
Rank is a local numerical diagnostic, not a structural identifiability proof.
"""
function local_sensitivities(model::FitModel, data, parameters::NamedTuple=model.initial;
                             scale=:log, standardized::Bool=true, rtol=sqrt(eps(Float64)))
    names = keys(model.initial)
    Set(keys(parameters)) == Set(names) || throw(ArgumentError("parameter names differ from model"))
    scale in (:log, :parameter) || throw(ArgumentError("scale must be :log or :parameter"))
    _finite_real(rtol, values(parameters)...)
    0 < rtol < 1 || throw(DomainError(rtol, "rank tolerance must lie between zero and one"))
    all(>(0), values(parameters)) || throw(DomainError(parameters, "parameters must be positive"))
    observations = _observations(data)
    ordered = NamedTuple{names}(Tuple(parameters[name] for name in names))
    sigma = _noise_sigma(observations, ordered)
    point = scale == :log ? _log_parameters(ordered) : collect(promote(float.(values(ordered))...))
    jacobian = ForwardDiff.jacobian(point) do x
        p = scale == :log ? _decode_parameters(model, x) : NamedTuple{names}(Tuple(x))
        prediction = _fit_prediction(model, p, data, observations.indices)[observations.indices]
        return standardized ? prediction ./ sigma : prediction
    end
    column_norms = NamedTuple{names}(Tuple(sqrt.(vec(sum(abs2, jacobian; dims=1)))))
    common = (; parameters=ordered, parameter_names=names, observed_indices=observations.indices,
                scale, standardized, jacobian, column_norms, rtol)
    if !all(isfinite, jacobian)
        return (; common..., status=:nonfinite_derivatives, rank=nothing,
                 singular_values=nothing, condition_number=Inf, unresolved_directions=nothing)
    end
    decomposition = svd(jacobian; full=size(jacobian,1) < size(jacobian,2))
    cutoff = rtol * maximum(decomposition.S)
    numerical_rank = count(>(cutoff), decomposition.S)
    deficient = numerical_rank < length(names)
    condition_number = deficient ? Inf : maximum(decomposition.S)/minimum(decomposition.S)
    return (; common..., status=deficient ? :rank_deficient : :full_rank,
             rank=numerical_rank, singular_values=decomposition.S, condition_number,
             unresolved_directions=decomposition.V[:,numerical_rank+1:end])
end

function _profile_point(model, data, parameter, value, reference, starts, rng; kwargs...)
    names = keys(model.initial)
    nuisance_names = Tuple(name for name in names if name != parameter)
    assemble = nuisance -> NamedTuple{names}(Tuple(name == parameter ? value : nuisance[name] for name in names))
    if isempty(nuisance_names)
        parameters = assemble(NamedTuple())
        try
            objective = gaussian_nll(model, data, parameters)
            return (; value, parameters, objective, converged=isfinite(objective), runs=NamedTuple[], error=nothing)
        catch err
            err isa _PredictionFailure || rethrow()
            return (; value, parameters=nothing, objective=Inf, converged=false,
                     runs=NamedTuple[], error=sprint(showerror,err))
        end
    end
    subset = p -> NamedTuple{nuisance_names}(Tuple(p[name] for name in nuisance_names))
    initial = NamedTuple{nuisance_names}(Tuple(
        model.lower[name] < reference.parameters[name] < model.upper[name] ?
            reference.parameters[name] : model.initial[name] for name in nuisance_names))
    reduced = FitModel((p,d)->model.predict(assemble(p),d);
        initial, lower=subset(model.lower), upper=subset(model.upper))
    # A profiled noise parameter becomes a fixed known SD in the nuisance fit.
    nuisance_data = data.sigma === parameter ?
        merge((; (name => getproperty(data,name) for name in propertynames(data))...), (sigma=value,)) : data
    fitted = fit_mle(reduced, nuisance_data; starts, rng, kwargs...)
    return (; value, parameters=fitted.converged ? assemble(fitted.parameters) : nothing,
             objective=fitted.objective, converged=fitted.converged, runs=fitted.runs,
             error=fitted.converged ? nothing : "No nuisance-parameter start converged")
end

function _profile_brackets(points, accepted)
    # Only adjacent successful points can bracket a crossing. Never interpolate
    # over failed grid fits or merge disconnected support regions into one interval.
    regions = NamedTuple[]
    i = 1
    while i <= length(points)
        if !accepted[i]
            i += 1
            continue
        end
        first_index = i
        while i < length(points) && accepted[i+1]
            i += 1
        end
        last_index = i
        lower = first_index > 1 && points[first_index-1].converged ?
            (points[first_index-1].value, points[first_index].value) : nothing
        upper = last_index < length(points) && points[last_index+1].converged ?
            (points[last_index].value, points[last_index+1].value) : nothing
        push!(regions, (; values=[points[j].value for j in first_index:last_index], lower_bracket=lower, upper_bracket=upper))
        i += 1
    end
    return regions
end

"""
    profile_likelihood(model, data, parameter, grid; starts=10, confidence=0.95,
                       flat_tol=1e-4, objective_tol=1e-6,
                       rng=Random.default_rng(), maxiters=1000,
                       outer_iterations=20, gtol=1e-7)

Fit a reference MLE, then fix one named parameter at each strictly increasing grid
value and refit all other parameters. Grid values must be inside the model bounds
(endpoints allowed). Every nuisance fit retains its multistart diagnostics.

Return twice the NLL difference from the reference, a nominal chi-square(1)
threshold, supported grid values, and adjacent crossing brackets. No confidence
endpoint is interpolated. Failed fits leave gaps. A grid optimum below the
reference invalidates its support report and is explicitly flagged. Flat or
open-ended profiles concern only the supplied grid and nuisance bounds.
"""
function profile_likelihood(model::FitModel, data, parameter::Symbol, grid;
                            starts=10, confidence=0.95, flat_tol=1e-4, objective_tol=1e-6,
                            rng=Random.default_rng(), maxiters=1000, outer_iterations=20, gtol=1e-7)
    parameter in keys(model.initial) || throw(ArgumentError("unknown profile parameter: $parameter"))
    starts isa Integer && !(starts isa Bool) && starts > 0 || throw(ArgumentError("profile starts must be a positive integer"))
    values = collect(grid)
    length(values) >= 2 || throw(ArgumentError("profile grid requires at least two values"))
    _finite_real(confidence, flat_tol, objective_tol, values...)
    0 < confidence < 1 || throw(DomainError(confidence, "confidence must lie between zero and one"))
    flat_tol >= 0 && objective_tol >= 0 || throw(DomainError((flat_tol,objective_tol), "tolerances must be nonnegative"))
    all(>(0), diff(values)) || throw(ArgumentError("profile grid must be strictly increasing"))
    all(model.lower[parameter] .<= values .<= model.upper[parameter]) ||
        throw(DomainError(values, "profile grid must stay within parameter bounds"))
    threshold = Distributions.quantile(Distributions.Chisq(1), confidence)
    options = (; maxiters, outer_iterations, gtol)
    reference = fit_mle(model, data; starts, rng, options...)
    common = (; parameter, grid=values, reference, confidence, threshold, flat_tol, objective_tol)
    if !reference.converged
        return (; common..., status=:no_reference_fit, reference_valid=false,
                 points=NamedTuple[], supported_values=nothing, regions=NamedTuple[])
    end
    fitted_points = [_profile_point(model,data,parameter,value,reference,starts,rng; options...) for value in values]
    tolerance = objective_tol * max(1,abs(reference.objective))
    reference_valid = all(p -> !p.converged || p.objective >= reference.objective-tolerance, fitted_points)
    points = [(; p..., deviance=p.converged ? 2(p.objective-reference.objective) : missing)
              for p in fitted_points]
    if !reference_valid
        return (; common..., status=:better_profile_optimum, reference_valid,
                 points, supported_values=nothing, regions=NamedTuple[])
    end
    accepted = [p.converged && p.deviance <= threshold for p in points]
    regions = _profile_brackets(points, accepted)
    all_converged = all(p -> p.converged, points)
    # Include the reference (deviance zero) in the flatness check: a constant
    # but poor distant plateau must not be called an unconstrained direction.
    flat = all_converged && maximum(abs(p.deviance) for p in points) <= flat_tol
    status = !all_converged ? :incomplete_grid : flat ? :flat_on_grid :
        isempty(regions) ? :no_supported_grid_points :
        all(r -> !isnothing(r.lower_bracket) && !isnothing(r.upper_bracket), regions) ?
            :crossings_bracketed : :grid_limited
    return (; common..., status, reference_valid, points,
             supported_values=values[accepted], regions)
end
