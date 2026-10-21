module Cocoa

using OrdinaryDiffEq: ODEProblem, solve, init, step!, Tsit5, Rodas5P, successful_retcode,
                      ReturnCode
using ..CoreModels: FitModel
import Distributions

include("baseline.jl")
include("inference.jl")

end
