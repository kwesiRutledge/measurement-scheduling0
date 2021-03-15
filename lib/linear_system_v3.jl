"""
    linear_system.jl
    Description:
        Defines the set of functions that help us work with linear systems.
        Represents a discrete-time linear dynamical system as follows:
         x(t+1) = A x(t) + B u(t) + w(t)
         y(t)   = C x(t) + v(t)
         z(t)   = D x(t) + d
        where w(t) in { w in R^n | ||w||_{infty} <= scaleW } and v(t) in { v in R^p | ||v||_{infty} <= scaleV }
"""

using LinearAlgebra
using JuMP
using Gurobi

struct LinearSystem
    A
    B
    C
    D
    k
    d
    polW # The bounded set that process noise can come from
    polV
    polX0
    polU
    list_polS
end


"""
    check
    Description:   
        Checks that the system described by the tuple LinearSystem is consistent.
"""
function check( system_in::LinearSystem )
    # Declare Variables
    n_x = 1
    n_y = 1
    n_u = 1
    n_z = 1
    
    # Algorithm

    sizeA = size(system_in.A)
    if length(sizeA) ≠ 0 #If A is not a scalar, then apply the following checks.
        n_x1 , n_x2 = sizeA
    end
    
    sizeB = size(system_in.B)
    if length(sizeB) ≠ 0
        n_x3, n_u1 = sizeB
    end

    sizeC = size(system_in.C)
    if length(sizeC) ≠ 0 #If C is not a scalar, then apply the following checks.
        n_y1, n_x4 = sizeC
    end
    
    sizeD = size(system_in.D)
    if length(sizeD) ≠ 0
        n_z1, n_x5 = sizeD
    end
    
    n_x6 = length(system_in.k)
    n_z2 = length(system_in.d)
    
    n_x7 = fulldim(system_in.polX0)
    n_x8 = fulldim(system_in.polW)
    
    n_y2 = fulldim(system_in.polV)
    
    n_u2 = fulldim(system_in.polU)
    
    if !(n_x1 == n_x2 == n_x3 == n_x4 == n_x5 == n_x6 == n_x7 == n_x8)
        return false
    end
    
    if !(n_y1 == n_y2)
        return false
    end
    
    if !(n_z1 == n_z2)
        return false
    end
    
    if !(n_u1 == n_u2)
        return false
    end
    
    if any(fulldim.(system_in.list_polS).!=n_z1)
        return false
    end
    
    #All Checks Passed
    return true
end

"""
    x_dim
    Description:
        Computes the dimension of the system.
"""
function x_dim( system_in::LinearSystem )

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    return size(system_in.A,1)
end

"""
    y_dim
    Description:
        Computes the dimension of the system's measurements.
"""
function y_dim( system_in::LinearSystem )

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    return size(system_in.C,1)
end

"""
    z_dim
    Description:
        Computes the dimension of the system's outputs.
"""
function z_dim( system_in::LinearSystem )

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    return size(system_in.D,1)
end

"""
    u_dim
    Description:
        Computes the dimension of the system's outputs.
"""
function u_dim( system_in::LinearSystem )

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    return size(system_in.B,2)
end




