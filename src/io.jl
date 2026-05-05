# ================================================================
# TherMat.jl - Input / Output (TOML Factory)
# ================================================================
#
# This file provides utilities for loading and saving material properties
# from and to TOML files. It includes factory methods for constructing
# Material and MaterialModel objects from parsed dictionary data.
# ================================================================

"""
    MaterialModel(d::AbstractDict)

Factory function to create a `MaterialModel` from a dictionary block.
"""
function MaterialModel(d::AbstractDict)
    m_type = d["type"]

    if m_type == "Constant"
        return ConstantModel(Float64(d["value"]))

    elseif m_type == "Polynomial"
        return PolynomialModel(Float64.(d["coefs"]))

    elseif m_type == "PiecewiseLinear"
        return PiecewiseLinearModel(
            Float64.(d["x"]),
            Float64.(d["y"]),
        )

    elseif m_type == "SemiCrystalline"
        return SemiCrystallineCp(
            Tg = Float64(d["Tg"]),
            Tm = Float64(d["Tm"]),
            c0 = Float64(d["c0"]),
            c1 = Float64(d["c1"]),
            c2 = Float64(d["c2"]),
            Δcₚ = Float64(d["dCp"]),
            ΔT = Float64(d["dT"]),
        )
    else
        valid_types = ["Constant", "Polynomial", "PiecewiseLinear", "SemiCrystalline"]
        error("Unknown model type: '$m_type'. Expected one of: $(join(valid_types, ", ")).")
    end
end

# ================================================================
# Serialization (to_dict)
# ================================================================

"""
    to_dict(m::MaterialModel)

Serialize a `MaterialModel` into a dictionary compatible with the TOML format.
"""
to_dict(m::ConstantModel) = Dict("type" => "Constant", "value" => m.val)

to_dict(m::PolynomialModel) = Dict("type" => "Polynomial", "coefs" => collect(m.coefs))

function to_dict(m::PiecewiseLinearModel)
    # Extract data points from interpolation object
    x = m.itp.knots[1]
    y = m.itp.coefs
    return Dict("type" => "PiecewiseLinear", "x" => collect(x), "y" => collect(y))
end

function to_dict(m::SemiCrystallineCp)
    return Dict(
        "type" => "SemiCrystalline",
        "Tg" => m.Tg,
        "Tm" => m.Tm,
        "c0" => m.c0,
        "c1" => m.c1,
        "c2" => m.c2,
        "dCp" => m.Δcₚ,
        "dT" => m.ΔT,
    )
end

"""
    to_dict(m::Material)

Serialize a `Material` into a dictionary compatible with the TOML format.
"""
function to_dict(m::Material)
    d = Dict{String,Any}(
        "name" => m.name,
        "note" => m.note,
        "rho" => to_dict(m.ρ),
        "cp" => to_dict(m.cₚ),
    )

    if m.k[1] === m.k[2] === m.k[3]
        d["k"] = to_dict(m.k[1])
    else
        d["k1"] = to_dict(m.k[1])
        d["k2"] = to_dict(m.k[2])
        d["k3"] = to_dict(m.k[3])
    end
    return d
end

"""
    save_material(filepath::AbstractString, m::Material)

Saves a `Material` object to a TOML file.
"""
function save_material(filepath::AbstractString, m::Material)
    open(filepath, "w") do io
        return TOML.print(io, to_dict(m))
    end
    return nothing
end

"""
    Material(d::AbstractDict)

Construct a `Material` from a dictionary parsed from TOML.
"""
function Material(d::AbstractDict)
    name = get(d, "name", "Unknown Material")
    note = get(d, "note", "")

    ρ = MaterialModel(d["rho"])
    cₚ = MaterialModel(d["cp"])

    if haskey(d, "k1") && haskey(d, "k2") && haskey(d, "k3")
        k = (
            MaterialModel(d["k1"]),
            MaterialModel(d["k2"]),
            MaterialModel(d["k3"]),
        )
    elseif haskey(d, "k")
        k_iso = MaterialModel(d["k"])
        k = (k_iso, k_iso, k_iso)
    else
        msg = "Material must define either 'k' (isotropic)" *
              "or 'k1', 'k2', 'k3' (orthotropic) in TOML."
        error(msg)
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
