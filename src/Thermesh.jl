module Thermesh

using CairoMakie
using Gridap
using Printf
using TOML
using Interpolations
using LinearAlgebra: Diagonal




include("constitutivemodels.jl")
include("bcs.jl")
include("simulation.jl")

export Material

end # module Thermesh
