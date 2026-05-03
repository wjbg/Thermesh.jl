# ================================================================
# Boundary condition definitions and utilities
# ================================================================
#
# This module defines boundary condition types and associated logic
# for thermal simulations.
#
# It provides:
# - A hierarchy of boundary condition types (Dirichlet, Neumann,
#   Robin, Radiation, Symmetry)
# - Containers for associating boundary conditions with boundary tags
# - Validation rules to ensure physically consistent combinations
# - Utility functions for constructing and managing boundary maps
#
# Boundary conditions are grouped per boundary and validated using
# simple trait-based rules (e.g. "exclusive" conditions such as
# Dirichlet and Symmetry cannot be combined with others).
# ================================================================

# ================================================================
# Boundary condition types
# ================================================================

abstract type BoundaryCondition end

"""Prescribed temperature."""
struct DirichletBC{T} <: BoundaryCondition
    T::T  # Real or Function
end

"""Prescribed heat flux: q = f."""
struct NeumannBC{Q} <: BoundaryCondition
    q::Q  # Real or Function
end

"""Convection: q = h(T∞ - T)."""
struct RobinBC{H, T} <: BoundaryCondition
    h::H
    T∞::T
end

"""Radiation: q = εσ(T∞⁴ - T⁴)."""
struct RadiationBC{E, T} <: BoundaryCondition
    ε::E
    T∞::T
end

"""Symmetry: q = 0."""
struct SymmetryBC <: BoundaryCondition end

"""
    BoundaryConditions

Vector of boundary conditions associated with a single boundary tag.
"""
const BoundaryConditions   = Vector{BoundaryCondition}

"""
    BoundaryConditionMap

Dictionary mapping boundary tags (e.g. IDs or names) to their
corresponding sets of boundary conditions.
"""
const BoundaryConditionMap = Dict{Union{Int, String}, BoundaryConditions}


# ================================================================
# Traits and validition
# ================================================================

is_exclusive(::BoundaryCondition) = false
is_exclusive(::DirichletBC)       = true
is_exclusive(::SymmetryBC)        = true

"""
    is_valid_boundary_conditions(bcs::BoundaryConditions)

Check whether a set of boundary conditions is consistent.

An *exclusive* condition (e.g. Dirichlet or Symmetry) must be the only
condition applied to a boundary.
"""
function is_valid_boundary_conditions(bcs::BoundaryConditions)
    n = length(bcs)
    n_exclusive = count(is_exclusive, bcs)
    return n_exclusive == 0 || (n == 1 && n_exclusive == 1)
end

"""
    is_valid_boundary_map(bcmap::BoundaryConditionMap)

Check whether all boundary condition sets in the map are consistent.
"""
function is_valid_boundary_map(bcmap::BoundaryConditionMap)
    all(is_valid_boundary_conditions, values(bcmap))
end

"""
    validate_boundary_map(bcmap::BoundaryConditionMap)

Validate that all boundary condition sets in `bcmap` are consistent.

Each boundary is checked using `is_valid_boundary_conditions`, and an error
is thrown if any invalid combination of boundary conditions is found.
"""
function validate_boundary_map(bcmap::BoundaryConditionMap)
    for (tag, bcs) in bcmap
        is_valid_boundary_conditions(bcs) ||
            error("Invalid boundary conditions for tag $tag: $bcs")
    end
end


# ================================================================
# Utilities
# ================================================================

"""
    add_bc!(bcmap::BoundaryConditionMap, tag, bc::BoundaryCondition)

Append a single `bc` to the boundary conditions for `tag`.
Initializes the entry if it does not yet exist.
"""
function add_bc!(bcmap::BoundaryConditionMap, tag, bc::BoundaryCondition)
    push!(get!(bcmap, tag, BoundaryConditions()), bc)
end

"""
    add_bc!(bcmap::BoundaryConditionMap, tag, bcs::Vector{<:BoundaryCondition})

Append a list of boundary conditions to the specified `tag`.
"""
function add_bc!(bcmap::BoundaryConditionMap, tag, bcs::Vector{<:BoundaryCondition})
    append!(get!(bcmap, tag, BoundaryConditions()), bcs)
end


"""
    dirichlet_data(bcmap::BoundaryConditionMap)

Extract Dirichlet boundary condition data from a boundary condition map.

Returns a tuple `(tags, values)` where:
- `tags` are the boundary identifiers carrying Dirichlet conditions
- `values` are the corresponding prescribed temperature values or functions

The returned vectors are aligned such that `values[i]` applies to `tags[i]`.
"""
function dirichlet_data(bcmap::BoundaryConditionMap)
    tags, vals = Int[], Any[]
    for (tag, bcs) in bcmap
        for bc in bcs
            if bc isa DirichletBC
                push!(tags, tag)
                push!(vals, bc.T)
            end
        end
    end
    return tags, vals
end
