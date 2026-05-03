using Gridap

"""
    SimulationStage

Defines the parameters for a specific interval of the simulation.
"""
struct SimulationStage
    name::String
    bcmap::BoundaryConditionMap
    dt::Float64          # Timestep size [s]
    duration::Float64    # How long this stage lasts [s]

    function SimulationStage(name, bcmap, dt, duration)
        validate_boundary_map(bcmap)
        new(name, bcmap, dt, duration)
    end
end

# Alias for the full simulation sequence
const SimulationSchedule = Vector{SimulationStage}


struct ThermalDomain
    model::DiscreteModel
    order::Int
    materials::Dict{Union{Int, String}, Material} # Mapping tags to Material properties
    u0::Any                                       # Initial Temperature field
end

# The "Grand Entry" for the solver
struct Simulation
    domain::ThermalDomain
    schedule::SimulationSchedule
end


"""
    build_spaces(stage::SimulationStage, domain::ThermalDomain)

Construct finite element trial and test spaces for a given simulation stage.

This function extracts Dirichlet boundary conditions from `stage.bc_map`
and uses them to impose essential boundary constraints on the trial space.

Returns:
- `U`: Transient trial finite element space with Dirichlet constraints applied
- `V`: Test finite element space

The spatial discretization order is taken from `domain.order`, and the
finite element model from `domain.model`.
"""
function build_spaces(stage::SimulationStage, domain::ThermalDomain)
    reffe = ReferenceFE(lagrangian, Float64, domain.order)
    tags, vals = dirichlet_data(stage.bc_map)

    V = TestFESpace(domain.model, reffe; conformity=:H1, dirichlet_tags=tags)
    U = TransientTrialFESpace(V, vals)
    return U, V
end
