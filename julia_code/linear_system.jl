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
    sizeA = size(system_in.A)
    if length(sizeA) ≠ 0 #If A is not a scalar, then apply the following checks.
        n_x1 , n_x2 = sizeA
        if n_x1 ≠ n_x2
            #This indicates that system_in.A is not square
            return false
        end
    end

    sizeC = size(system_in.C)
    if length(sizeC) ≠ 0 #If C is not a scalar, then apply the following checks.
        n_y, n_x3 = sizeC
        if n_x1 ≠ n_x3
            #This indicates that system_in.A is of different dimension than Checks
            return false
        end
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
    H_w = [ I(n_w) ; -I(n_w) ]
    h_w = system_in.eta_w*ones(2*n_w,1)
    dim_H_w = 2*n_w

    H_v = [ I(n_v) ; -I(n_v) ]
    h_v = system_in.eta_v*ones(2*n_v,1)
    dim_H_v = 2*n_v

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

"""
    get_alap_estimation_schedule
    Description:
        Uses one of a variety of methods to find an ALAP schedule.
"""
function get_alap_estimation_schedule( system_in , TimeHorizon , MeasurementBudget , method_number::Int )



end

"""
    alap_estimation_schedule_alg1
    Description:
        This algorithm takes the time horizon, divides it into |M|+1 intervals.
        It is agnostic of the system and is strictly based on the numbers provided.
    Inputs:
        TimeHorizon - An integer representing the final time instant of the problem.
        MeasurementBudget - An integer representing the number of measurements that the estimator
                            can take during the time [0,TimeHorizon-1].
"""
function alap_estimation_schedule_alg1( TimeHorizon , MeasurementBudget )

    # Equally Divide the Time Horizon Up Into M+1 blocks of length TimeChunk (if possible)
    FloatTimeChunk = TimeHorizon / (MeasurementBudget+1)
    IntTimeChunk = floor(FloatTimeChunk);

    M = []

    if IntTimeChunk == FloatTimeChunk
        # If there are exactly MeasurementBudget+1 windows of length IntTimeChunk,
        # then the measurement times are straight forward:
        for M_idx = 1:MeasurementBudget
            append!(M,Int(M_idx*IntTimeChunk) )
        end

    else
        # If there aren't exactly MeasurementBudget+1 windows of length IntTimeChunk,
        # then use a small window initially before doing periodic chunks.
        initialChunk = rem(TimeHorizon,MeasurementBudget+1)
        append!(M,Int( initialChunk ))
        for M_idx = 2:MeasurementBudget
            append!(M,Int(initialChunk + (M_idx-1)*IntTimeChunk) )
        end
    end

    return M

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