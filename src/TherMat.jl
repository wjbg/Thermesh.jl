# ==============================================================================
# TherMat.jl
# ==============================================================================
#
# A Julia package for thermal constitutive material modeling and calibration.
#
# This package provides:
# - Representations of temperature-dependent material properties.
# - Tools for loading/saving material data from TOML configuration files.
# - (Upcoming) Utilities for fitting constitutive models to experimental data.
#
# Primary focus: Thermal conductivity, density, and specific heat for 
# engineering materials, with specialized support for polymers.
# ==============================================================================

module TherMat

using Interpolations
using LinearAlgebra: Diagonal
using TOML

include("types.jl")
include("models.jl")
include("io.jl")
include("fitting.jl")

export Material, save_material

end # module TherMat
