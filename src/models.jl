# ================================================================
# TherMat.jl - Scalar Constitutive Models
# ================================================================
#
# This file contains various mathematical models for material properties
# such as thermal conductivity, density, and specific heat.
# These models (Constant, Polynomial, PiecewiseLinear, etc.) are all
# callable functions of temperature.
# ================================================================

"""
    smooth_heaviside(T, δ)

Smooth approximation of the Heaviside step function.
"""
function smooth_heaviside(T, δ)
    return (1 + (tanh ∘ (x-> x/δ))(T))/2
end

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
"""
struct PolynomialModel{N,S} <: MaterialModel
    coefs::NTuple{N,S}  # Coefficients, lowest order first
end

(m::PolynomialModel)(x) = evalpoly(x, m.coefs)

PolynomialModel(coefs::Vararg{S,N}) where {S,N} = PolynomialModel{N,S}(coefs)
function PolynomialModel(coefs::AbstractVector{S}) where {S}
    return PolynomialModel(Tuple(coefs))
end

"""
    PiecewiseLinearModel(x, y)

Linear interpolation between data points with flat extrapolation (clamping).
"""
struct PiecewiseLinearModel{T} <: MaterialModel
    itp::T
end

function PiecewiseLinearModel(x::AbstractVector, y::AbstractVector)
    issorted(x) || error("x-coordinates must be strictly increasing.")
    length(x) == length(y) || error("x and y vectors must be the same length.")
    itp = linear_interpolation(x, y, extrapolation_bc=Flat())
    return PiecewiseLinearModel(itp)
end

(m::PiecewiseLinearModel)(T) = m.itp(T)

"""
    SemiCrystallineCp <: MaterialModel

Specific heat model for semi-crystalline polymers.
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
    return (
        m.c0 + m.c1 * T +                                               # T ≤ Tg
        m.Δcₚ * smooth_heaviside(T - m.Tg, m.ΔT) +                     # Tg ≤ T ≤ Tm
        (m.c2 - m.c1) * (T - m.Tm) * smooth_heaviside(T - m.Tm, m.ΔT)  # Tm ≤ T
    )
end
