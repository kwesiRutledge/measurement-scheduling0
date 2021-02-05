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
    W # The bounded set that process noise can come from
    V # The bounded set that measurement noise can come from
end

"""
    check
    Description:   
        Checks that the system described by the tuple LinearSystem is consistent.
"""
function check( system_in::LinearSystem )
    # Declare Variables
    n = 1
    m = 1
    p = 1
    
    # Algorithm

    sizeA = size(system_in.A)
    if length(sizeA) ≠ 0 #If A is not a scalar, then apply the following checks.
        n_x1 , n_x2 = sizeA
        if n_x1 ≠ n_x2
            #This indicates that system_in.A is not square
            return false
        end
        n = n_x1 # Assign dimension.
    end

    sizeC = size(system_in.C)
    if length(sizeC) ≠ 0 #If C is not a scalar, then apply the following checks.
        n_y, n_x3 = sizeC
        if n_x1 ≠ n_x3
            #This indicates that system_in.A is of different dimension than Checks
            return false
        end
        p = n_y # Assign dimension of the output.
    end

    # Check the Type of Disturbance Sets
    A_symbol = Symbol("A")
    b_symbol = Symbol("b")

    W_type = typeof(system_in.W)
    if !( A_symbol in fieldnames(W_type) && b_symbol in fieldnames(W_type) )
        # If this is true, then the disturbance set is not a Polyhedron
        # in H representation.
        return false
    end

    if fulldim(system_in.W) ≠ n
        # If the W Polyhedron is not defined in the proper dimension,
        # then return an error.
        return false
    end

    V_type = typeof(system_in.V)
    if !( A_symbol in fieldnames(V_type) && b_symbol in fieldnames(V_type) )
        # If this is true, then the disturbance set is not a Polyhedron
        # in H representation.
        return false
    end

    if fulldim(system_in.V) ≠ p
        # If the V Polyhedron is not defined in the proper dimension,
        # then return an error.
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
        Computes the dimension of the system's outputs.
"""
function y_dim( system_in::LinearSystem )

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    return size(system_in.C,1)

end

"""
    define_simple_eta_HPolyhtope
    Description:

"""
function define_simple_eta_HPolytope( system_in::LinearSystem , eta_x0 , T_Horizon_in )
    #Constants
    n_x = x_dim(system_in)
    n_y = y_dim(system_in)

    n_w = n_x
    n_v = n_y

    T = T_Horizon_in

    #Create the Polytopes defined by system_in.eta_w and system_in.eta_v
    H_w = system_in.W.A
    h_w = system_in.W.b
    dim_H_w = length(h_w)

    H_v = system_in.V.A
    h_v = system_in.V.b
    dim_H_v = length(h_v)

    H_x0 = [I(n_x) ; -I(n_x)]
    h_x0 = eta_x0*ones(2*n_x,1)
    dim_H_x0 = 2*n_x

    #Algorithm
    H_eta = [kron(Matrix(1I,T,T), H_w) zeros(T*dim_H_w,T*n_y + n_x) ;
             zeros(T*dim_H_v, T*n_x) kron(Matrix(1I,T,T), H_v) zeros(T*dim_H_v,n_x);
             zeros(dim_H_x0,T*(n_x+n_y)) H_x0]

    h_eta = [kron(ones(T,1),h_w); kron(ones(T,1),h_v); h_x0];

    return H_eta, h_eta

end

"""
    find_est_error_bound_from_time_0_to_T
    Description:
        Finds the bound B which bounds all estimation error values between time t=0
        and time T.
"""
function find_est_error_bound_from_time_0_to_T( system_in::LinearSystem , T::Int , eta_x0 )

    #Define Constants
    model = Model(Gurobi.Optimizer) # Create Optimization model
    H_eta, h_eta = define_simple_eta_HPolytope( system_in , eta_x0 , T ) # Create Eta Polyhedron (Polyhedron of all disturbances over time horizon)
    n_Heta = size(H_eta,1)

    n_x = x_dim(system_in)
    n_y = y_dim(system_in)

    S = compute_Sw( system_in.A , T )
    # C_bar = compute_C_M( system_in.C , [T-1] , T )
    J = compute_Sx0( system_in.A , T)

    R_T = Matrix{Float64}([zeros(n_x,T*n_x) I(n_x)])

    #Create the Optimization Variables
    @variable(model,alpha0)
    @variable(model,Lambda[1:2*n_x,1:n_Heta])

    # Create Constraints
    @constraint(model,alpha0.>=0)
    @constraint(model,Lambda.>=zeros(size(Lambda)))

    LHS1 = Lambda * H_eta
    # println(typeof(S))
    # println(typeof(J))
    # println(typeof([ S  zeros(size(S,1),n_y*T) J ]))
    # println(S)
    # println(typeof(R_T))
    # println( R_T*[ S  zeros(size(S,1),n_y*T) J ] )
    RHS1 = [I(n_x); -I(n_x)] * R_T * [ S  zeros(size(S,1),n_y*T) J ]
    @constraint(model, LHS1 .== RHS1 )

    LHS2 = Lambda * h_eta
    RHS2 = alpha0*ones(2*n_x,1)
    @constraint(model, LHS2 .<= RHS2)

    # Create objective
    @objective(model, Min, alpha0)

    #Optimize!
    set_silent(model)
    optimize!(model)

    return termination_status(model) , objective_value(model)

end

"""
    evaluate_schedule_wrt
    Description:
        Evaluates the schedule with respect to (i) a desired linear system, (ii) the schedule
        (iii) the time horizon, and (iv) one of several objectives including:
            1 - Maximum Infinity Norm of Estimation error 
            2 - Sum of the Estimation Error bounds at each time
    Usage:
        obj_val0 , feas_flag0 = evaluate_schedule_wrt( ls0, schedule0 , T0 , 1 , eta_x0 )
        
        
"""
function evaluate_schedule_wrt( system_in::LinearSystem , schedule_in , time_horizon_in , objective_flag::Int , x0_description )

    # Input Processing

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end

    for schedule_val in schedule_in
        if schedule_val > time_horizon_in
            throw(ArgumentError("The schedule given contains times which are outside of the time horizon! (" + str(schedule_val) + ")" ))
        end
    end

    expected_objective_flag_range = (1,2)

    if (objective_flag < expected_objective_flag_range[1]) | (objective_flag > expected_objective_flag_range[2])
        throw(ArgumentError("The input objective_flag is not in the desired range $(expected_objective_flag_range).\n"))
    end 

    # Constants

    eta_x0 = x0_description #Fix this later when eta_description becomes a polytope?
    H_eta, h_eta = define_simple_eta_HPolytope( system_in , eta_x0 , time_horizon_in ) # Create Eta Polyhedron (Polyhedron of all disturbances over time horizon)
    n_Heta = size(H_eta,1)

    n_x = x_dim(system_in)
    n_y = y_dim(system_in)

    #Clear the model?
    empty!(model)

    #Create the Optimization Variables
    @variable(model,Q[1:n_x*time_horizon_in,1:n_y*time_horizon_in])
    @variable(model,alpha0[1:time_horizon_in+1,1])
    @variable(model,alpha1)

    #Constraints
    @constraint(model,alpha0.>=0)
    @constraint(model,alpha1.>=0)

    S = compute_Sw( system_in.A , time_horizon_in )
    C_bar = compute_C_M( system_in.C , schedule_in , time_horizon_in )
    J = compute_Sx0( system_in.A , time_horizon_in)

    P_xi_w = S + S*Q*C_bar*S
    P_xi_v = S * Q 
    P_xi_xi0 = ( I + S*Q*C_bar)*J
    P_xiTilde = 0

    R_T_mT = [ I ; -Matrix(1I, n_x*(time_horizon_in+1), n_x*(time_horizon_in+1)) ]
    
    R_0_T = I( n_x*(time_horizon_in+1) )
    R_T = [ zeros(n_x,n_x*time_horizon_in) I(n_x) ]

    if objective_flag == 1 #Maximum Estimation Error Bound Objective
        #
        @variable(model,Lambda[1:2*n_x*(time_horizon_in+1),1:n_Heta])
        @constraint(model,Lambda.>=zeros(size(Lambda)))

        LHS1 = Lambda * H_eta
        RHS1 = [I(n_x*(time_horizon_in+1)) ; -I(n_x*(time_horizon_in+1))] * R_0_T * [ P_xi_w P_xi_v P_xi_xi0 ]
        @constraint(model, LHS1 .== RHS1 )

        LHS2 = Lambda * h_eta
        RHS2 = alpha1 * ones(2*n_x*(time_horizon_in+1),1)
        @constraint(model, LHS2 .<= RHS2)

    elseif objective_flag == 2 #Sum of estimation error bounds at each time steps

        # Variable Definitions
        @variable(model,Lambda[1:2*n_x*(time_horizon_in+1),1:n_Heta])
        @constraint(model,Lambda.>=zeros(size(Lambda)))

        LHS1 = Lambda * H_eta
        RHS1 = [I(n_x*(time_horizon_in+1)) ; -I(n_x*(time_horizon_in+1))] * R_0_T * [ P_xi_w P_xi_v P_xi_xi0 ]
        @constraint(model, LHS1 .== RHS1 )

        LHS2 = Lambda * h_eta
        RHS2 = kron(alpha0,ones(2*n_x))
        @constraint(model, LHS2 .<= RHS2)

    end

    # Define Objective

    if objective_flag == 1
        @objective(model,Min,alpha1)
    elseif objective_flag == 2
        @objective(model,Min,sum(alpha0))
    end

    # Optimize!
    optimize!(model)

    print("termination_status = $(termination_status(model))\n")

    return objective_value(model) , termination_status(model)

end