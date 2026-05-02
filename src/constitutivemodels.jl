using LinearAlgebra: Diagonal

# ================================================================
# Definitions and interface
# ================================================================

abstract type MaterialModel end

"""
    Material{K1, K2, K3, D, C}

A constitutive model for thermal simulations.

A `Material` represents temperature-dependent material behavior by
composing three independent constitutive models:

- `k::Tuple{K1, K2, K3}` : conductivity models for the 1, 2, and 3 directions.
- `ρ::D`                 : density model
- `cₚ::C`                : specific heat model

Each field must be callable:

    model(T) -> property value

This allows all material properties to vary smoothly with temperature
without introducing additional material types.
"""
Base.@kwdef struct Material{K1, K2, K3, D, C}
    # Physics
    k::Tuple{K1, K2, K3}
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
Material(k::MaterialModel, ρ::MaterialModel, cₚ::MaterialModel) =
    Material(k=(k, k, k), ρ=ρ, cₚ=cₚ)

"""
    Material(k::Real, ρ::Real, cp::Real)

Convenience constructor for isothermal isotropic materials.
"""
Material(k::Real, ρ::Real, cₚ::Real) =
    Material(ConstantModel(k), ConstantModel(ρ), ConstantModel(cₚ))

# Accessor
(m::Material)(T) = (k=k(m, T), ρ=ρ(m, T), cₚ=cₚ(m, T))


# ================================================================
# Print functions
# ================================================================

function Base.show(io::IO, m::Material)
    is_iso = (m.k[1] === m.k[2] === m.k[3])
    type_str = is_iso ? "Isotropic" : "Orthotropic"

    print(io, "Material(\"", m.name, "\", ", type_str, ")")
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
end


# ================================================================
# Constitutive models
# ================================================================

"""
    smooth_heaviside(T, δ)

Smooth approximation of the Heaviside step function.

Used to regularize phase transitions and avoid discontinuities in FEM assembly.
The parameter `δ` must be non-zero and controls the width of the transition region.
"""
smooth_heaviside(T, δ) = 0.5 * (1 + tanh(T / δ))

"""
    ConstantModel(val)

A temperature-independent constitutive model returning a constant value.
"""
struct ConstantModel{S} <: MaterialModel
    val::S
end

(m::ConstantModel)(x) = m.val

"""
    PolynomialModel(coefs)

A temperature-dependent model using an N-th order polynomial.
Coefficients should be provided in increasing order of power (constant first).
"""
struct PolynomialModel{N, S} <: MaterialModel
    coefs::NTuple{N, S}  # Coefficients, lowest order first
end

(m::PolynomialModel)(x) = evalpoly(x, m.coefs)

"""
    PolynomialModel(coefs::Vararg)
    PolynomialModel(coefs::AbstractVector)

Construct a `PolynomialModel` from individual coefficients or a vector.
"""
PolynomialModel(coefs::Vararg{S,N}) where {S,N} = PolynomialModel{N,S}(coefs)
PolynomialModel(coefs::AbstractVector{S}) where {S} =
    PolynomialModel(Tuple(coefs))

"""
    SemiCrystallineCp <: MaterialModel

Specific heat model for semi-crystalline polymers, capturing the glass transition
jump and the slope change at the melting point.
"""
Base.@kwdef struct SemiCrystallineCp{S} <: MaterialModel
    Tg::S   # Glass transition temperature
    Tm::S   # Melt temperature
    c0::S   # Cp at T = 0 [K]
    c1::S   # Slope for T < Tm [K]
    c2::S   # Slope for T > Tm [K]
    Δcₚ::S  # Jump at Tg
    ΔT::S   # Heaviside smoothing window
end

function (m::SemiCrystallineCp)(T)
    m.c0 + m.c1 * T +                                               # T ≤ Tg
    m.Δcₚ * smooth_heaviside(T - m.Tg, m.ΔT) +                     # Tg ≤ T ≤ Tm
    (m.c2 - m.c1) * (T - m.Tm) * smooth_heaviside(T - m.Tm, m.ΔT)  # Tm ≤ T
end


# ================================================================
# Input / Output (TOML Factory)
# ================================================================

"""
    MaterialModel(d::AbstractDict)

Factory function to create a `MaterialModel` from a dictionary block.
Dispatches based on the "type" key.
"""
function MaterialModel(d::AbstractDict)
    m_type = d["type"]

    if m_type == "Constant"
        return ConstantModel(Float64(d["value"]))

    elseif m_type == "Polynomial"
        return PolynomialModel(Float64.(d["coefs"]))

    elseif m_type == "SemiCrystalline"
        return SemiCrystallineCp(
            Tg  = Float64(d["Tg"]),
            Tm  = Float64(d["Tm"]),
            c0  = Float64(d["c0"]),
            c1  = Float64(d["c1"]),
            c2  = Float64(d["c2"]),
            Δcₚ = Float64(d["dCp"]),
            ΔT  = Float64(d["dT"])
        )
    else
        valid_types = ["Constant", "Polynomial", "SemiCrystalline"]
        error("Unknown model type: '$m_type'. Expected one of: $(join(valid_types, ", ")).")
    end
end

"""
    Material(d::AbstractDict)

Construct a `Material` from a dictionary parsed from TOML. Supports metadata and
both isotropic (k) and orthotropic (k1, k2, k3) inputs.
"""
function Material(d::AbstractDict)
    name   = get(d, "name", "Unknown Material")
    note   = get(d, "note", "")

    ρ  = MaterialModel(d["rho"])
    cₚ = MaterialModel(d["cp"])

    if haskey(d, "k1") && haskey(d, "k2") && haskey(d, "k3")
        k = (MaterialModel(d["k1"]),
             MaterialModel(d["k2"]),
             MaterialModel(d["k3"]))
    elseif haskey(d, "k")
        k_iso = MaterialModel(d["k"])
        k = (k_iso, k_iso, k_iso)
    else
        error("Material must define either 'k' (isotropic) or 'k1', 'k2', 'k3' (orthotropic) in TOML.")
    end

    return Material(k=k, ρ=ρ, cₚ=cₚ, name=name, note=note)
end

"""
    Material(filepath::String)

Reads a TOML file and returns a fully initialized `Material` object.
"""
function Material(filepath::AbstractString)
    isfile(filepath) || error("Material file not found: $filepath")
    return Material(TOML.parsefile(filepath))
end


# Local Variables:
# outline-regexp: "# \\([=-]+\\)"
# outline-level: (lambda ()
#                  (let ((s (match-string 1)))
#                    (cond
#                     ((string-match-p "=" s) 1)
#                     ((string-match-p "-" s) 2)
#                     (t 3))))
# eval: (outline-minor-mode 1)
# End:
