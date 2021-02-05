"""
measurement-scheduling0.jl
Description:
    Contains all of the necessary files.
"""

# List of Necessary Packages

using LinearAlgebra
using JuMP
using Gurobi
using Polyhedra
using HCubature

# Interface Definitions
include("mpc_functions.jl")
include("linear_system.jl")
include("alap_estimation_schedule.jl")
include("fixed_height_planar_lipm.jl")
include("double_integrator_drones.jl")