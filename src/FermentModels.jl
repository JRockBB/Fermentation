"""
    FermentModels

Mechanistic fermentation tools: growth responses, buffer chemistry, ODE assembly,
feeding and environment schedules, and parameter provenance validation.
See `docs/PLAN.md` for food-specific models and calibration work.
"""
module FermentModels

include("core/CoreModels.jl")
include("cocoa/Cocoa.jl")
include("sourdough/Sourdough.jl")
include("vegetables/Vegetables.jl")

using .CoreModels: load_params
using .CoreModels: γ_T_ctmi, γ_T_ganzle, γ_pH_cardinal, γ_acid, γ_salt
using .CoreModels: undissociated_fraction, undissociated_acid, henderson_hasselbalch,
                   buffer_from_capacity, pH
using .CoreModels: ModelSpec, build_rhs, initial_state, physical_state, fermentation_problem
using .CoreModels: feeding_fraction, dilute!, dilution_callback, TemperatureSchedule,
                   temperature_stops
using .CoreModels: FitModel, fit_mle, gaussian_nll
using .CoreModels: literature_prior, fit_bayes, posterior_predictive
using .CoreModels: local_sensitivities, profile_likelihood
using .Cocoa: COCOA_STATES, COCOA_PARAMETERS, cocoa_rhs!, cocoa_problem,
              cocoa_physical_parameters, cocoa_scaled_parameters,
              cocoa_observations, cocoa_fit_model, cocoa_priors

export load_params
export γ_T_ctmi, γ_T_ganzle, γ_pH_cardinal, γ_acid, γ_salt
export undissociated_fraction, undissociated_acid, henderson_hasselbalch,
       buffer_from_capacity, pH
export ModelSpec, build_rhs, initial_state, physical_state, fermentation_problem
export feeding_fraction, dilute!, dilution_callback, TemperatureSchedule, temperature_stops
export FitModel, fit_mle, gaussian_nll
export literature_prior, fit_bayes, posterior_predictive
export local_sensitivities, profile_likelihood
export COCOA_STATES, COCOA_PARAMETERS, cocoa_rhs!, cocoa_problem,
       cocoa_physical_parameters, cocoa_scaled_parameters,
       cocoa_observations, cocoa_fit_model, cocoa_priors

end
