using FermentModels
using Test
using TOML

@testset "FermentModels" begin
    include("parameters.jl")
    include("fetch.jl")
    include("secondary.jl")
    include("chemistry.jl")
    include("assembly.jl")
    include("environment.jl")
    include("estimation.jl")
    include("bayes.jl")
    include("identifiability.jl")
    include("cocoa.jl")
    include("cocoa_solver_budget.jl")
    include("estimated_noise.jl")
    include("fit_budget.jl")
    include("fit_stationarity.jl")
    include("bayes_execution.jl")
    include("cocoa_fit_output.jl")
    include("cocoa_mle_diagnostics.jl")
    include("cocoa_benchmark.jl")
    include("cocoa_benchmark_diagnostics.jl")
    include("cocoa_benchmark_bayes.jl")
    include("cocoa_bayes_diagnostics.jl")
    include("cocoa_bayes_audit.jl")
    include("cocoa_runtime_diagnostics.jl")
end
