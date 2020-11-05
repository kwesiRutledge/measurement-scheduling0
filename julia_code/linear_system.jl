"""
    linear_system.jl
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

struct LinearSystem
    A
    B
    C
    eta_w # The infinity norm bound of disturbances w
    eta_v # The infinity norm bound of disturbances v
end

"""
    check
    Description:   
        Checks that the system described by the tuple LinearSystem is consistent.
"""
function check( system_in::LinearSystem )
    n_x1 , n_x2 = size(system_in.A)
    if n_x1 ≠ n_x2
        #This indicates that system_in.A is not square
        return false
    end

    n_y, n_x3 = size(system_in.C)
    if n_x1 ≠ n_x3
        #This indicates that system_in.A is of different dimension than Checks
        return false
    end

    #All Checks Passed
    return true
end

"""
    n_x
    Description:
        Computes the dimension of the system.
"""
function n_x( system_in::LinearSystem )

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    return size(system_in.A,1)

end

"""
    n_y
    Description:
        Computes the dimension of the system's outputs.
"""
function n_y( system_in::LinearSystem )

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    return size(system_in.C,1)

end

"""
    find_ALAP_time_from_X0_to_Bound
    Description:
        Identifies the amount of time it takes for the 
"""
function find_ALAP_time_from_X0_to_bound( system_in::LinearSystem, eta_x0 , bound_in ; iterLimit=10 )
    #Constants
    dim_x = n_x(system_in)
    dim_y = n_y(system_in)

    model = Model(Gurobi.Optimizer)

    #Start the Loop on time horizon T 
    for T = 1:iterLimit

        #Clear the model?
        empty!(model)

        #Create the Optimization Variables
        @variable(model,Q[1:n_x*T,1:n_x*T])
        @variable(model,r[1:n_x*T,1])
        @variable(model,Lambda[1:2*n_x*(T+1),1:T*(n_w+n_v)+n_x0])
        @variable(model,alpha0[1:T+1,1])

        #Constraints
        @constraint(model,alpha0.>=0)
        @constraint(model,Lambda.>=zeros(size(Lambda)))

        S = compute_Sw( system_in.A , T )
        C_bar = compute_C_M( system_in.C , [T-1] , T )
        J = compute_Sx0( system_in.A , T)

        P_xi_w = S + S*Q*C_bar*S
        P_xi_v = S * Q 
        P_xi_xi0 = ( I + S*Q*C_bar)*J

        R_T_mT = [ I ; -Matrix(1I, n_x*(T+1), n_x*(T+1)) ]

        #Optimize!
        optimize!(model)

        print("termination_status = $(termination_status(model))\n")
    end

end