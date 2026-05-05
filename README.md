# TherMat.jl

A Julia package for thermal constitutive material modeling and calibration.

## Overview

**TherMat.jl** provides a structured framework for representing, loading, and (upcoming) calibrating temperature-dependent material properties. It is designed for engineering applications requiring precise thermal characterization, with a particular focus on semi-crystalline polymers and fiber-reinforced composites.

## Features

- **Temperature-Dependent Models**: Support for constant, polynomial, piecewise linear, and specialized semi-crystalline specific heat models.
- **Hierarchical Material Structures**: Combine anisotropic conductivity, density, and specific heat models into a single callable `Material` object.
- **TOML Integration**: Full support for loading and saving material properties from/to TOML configuration files.
- **Blue Style Adherent**: Developed following the [Blue Style](https://github.com/invenia/BlueStyle) guide for clean, idiomatic Julia code.

## Installation

```julia
using Pkg
Pkg.add("TherMat")
```

## Quick Start

```julia
using TherMat

# Load a material from a TOML file
m = Material("path/to/material.toml")

# Evaluate properties at a specific temperature [K]
T = 450.0
props = m(T)

println("Density at $T K: ", props.ρ)
println("Conductivity Matrix: ", props.k)

# Save the material back to disk
save_material("updated_material.toml", m)
```

## Project Structure

- `src/types.jl`: Core data structures and abstract types.
- `src/models.jl`: Mathematical constitutive model implementations.
- `src/io.jl`: TOML loading and serialization logic.
- `src/fitting.jl`: (Placeholder) Tools for experimental data fitting.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details (or choice of license).
