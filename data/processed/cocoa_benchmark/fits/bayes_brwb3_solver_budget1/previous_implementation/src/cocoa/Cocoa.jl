module Cocoa

using OrdinaryDiffEq: ODEProblem, solve, Tsit5, Rodas5P, successful_retcode
using ..CoreModels: FitModel
import Distributions

include("baseline.jl")
include("inference.jl")

end
