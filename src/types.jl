# ================================================================
# TherMat.jl - Core Type Definitions
# ================================================================
#
# This file defines the fundamental data structures used in TherMat.jl,
# including the abstract MaterialModel type and the Material struct.
# These types provide the foundation for representing temperature-dependent
# material properties in thermal simulations.
# ================================================================

abstract type MaterialModel end

"""
    Material{K1,K2,K3,D,C}

A constitutive model for thermal simulations.

A `Material` represents temperature-dependent material behavior by
composing three independent constitutive models:

- `k::Tuple{K1,K2,K3}` : conductivity models for the 1, 2, and 3 directions.
- `ρ::D`               : density model
- `cₚ::C`              : specific heat model

Each field must be callable:

    model(T) -> property value
"""
Base.@kwdef struct Material{K1,K2,K3,D,C}
    # Physics
    k::Tuple{K1,K2,K3}
    ρ::D
    cₚ::C

    # Metadata
    name::String = "Unknown Material"
    note::String = ""
end

"""
    Material(k::MaterialModel, ρ::MaterialModel, cp::MaterialModel)

Convenience constructor for generic isotropic materials.
"""
function Material(k::MaterialModel, ρ::MaterialModel, cₚ::MaterialModel)
    return Material(k=(k, k, k), ρ=ρ, cₚ=cₚ)
end

"""
    Material(k::Real, ρ::Real, cp::Real)

Convenience constructor for isothermal isotropic materials.
"""
function Material(k::Real, ρ::Real, cₚ::Real)
    return Material(ConstantModel(k), ConstantModel(ρ), ConstantModel(cₚ))
end

# Accessor
function (m::Material)(T)
    return (
        k = Diagonal((m.k[1](T), m.k[2](T), m.k[3](T))),
        ρ = m.ρ(T),
        cₚ = m.cₚ(T),
    )
end

# ================================================================
# Print functions
# ================================================================

function Base.show(io::IO, m::Material)
    is_iso = (m.k[1] === m.k[2] === m.k[3])
    type_str = is_iso ? "Isotropic" : "Orthotropic"
    print(io, "Material(\"", m.name, "\", ", type_str, ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", m::Material)
    is_iso = (m.k[1] === m.k[2] === m.k[3])

    # Header and Metadata
    println(io, "Material: ", m.name)
    println(io, "  Type:   ", is_iso ? "Isotropic" : "Orthotropic")
    if !isempty(m.note)
        println(io, "  Note:   ", m.note)
    end

    # Density and Specific Heat
    println(io, "  Models:")
    println(io, "    ρ  => ", typeof(m.ρ).name.name)
    println(io, "    cₚ => ", typeof(m.cₚ).name.name)

    # Conductivity
    if is_iso
        println(io, "    k  => ", typeof(m.k[1]).name.name)
    else
        println(io, "    k₁ => ", typeof(m.k[1]).name.name)
        println(io, "    k₂ => ", typeof(m.k[2]).name.name)
        println(io, "    k₃ => ", typeof(m.k[3]).name.name)
    end
    return nothing
end
