"""
    constraint_generator.jl
    Description:
        Defines the set of functions that help us work with linear systems.
        Represents a discrete-time linear dynamical system as follows:
         x(t+1) = A x(t) + B u(t) + w(t)
         y(t)   = C x(t) + v(t)
        where w(t) in { w in R^n | ||w||_{infty} <= eta_w } and v(t) in { v in R^p | ||v||_{infty} <= eta_v }
"""

using LinearAlgebra
using JuMP
using Gurobi

struct ConstraintGenerator

end

"""
GetPolytopeInclusionConstraints
Description:

"""
function GetPolytopeInclusionConstraints( Hx , hx , Hy , hy )

    # Get the Dimensions of Each Polytope
    sizeHx = size(Hx)
    if length(sizeHx) ≠ 0
        println("This is a scalar!")
        return NaN
    end

    n_Hx, dim_x = sizeHx
    println(n_Hx)
    println(dim_x)


end