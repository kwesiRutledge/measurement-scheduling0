"""
    alap_estimation_schedule.jl
    Description:
        This file contains functions related to As-Late-As-Possible Measurement Scheduling
        for the linear systems defined as LinearSystem objects.
"""

"""
    norm_of_est_error_after
    Description:
        Returns the estimation error for the system (in terms of infinity norm) after 
        T time steps for the input system system_in.
    Usage:
        objective_val , termination_status = norm_of_est_error_after( T , ls0 , eta_x0 )
"""
function norm_of_est_error_after( time_horizon_in::Int , system_in::LinearSystem , x0_description )
    
    # Input Checking

    if !check( system_in )
        throw(ArgumentError("The input LinearSystem is not valid. Please check its dimensions."))
    end
    
    if length(x0_description) == 1
        # This means that x0_description is eta_x0
        eta_x0 = x0_description 
    end

    # Constants

    n_x = x_dim( system_in )
    n_y = y_dim( system_in )

    model = Model(Gurobi.Optimizer)

    # Define Disturbance Polytope

    H_eta, h_eta = define_simple_eta_HPolytope( system_in , eta_x0 , time_horizon_in ) # Create Eta Polyhedron (Polyhedron of all disturbances over time horizon)
    n_Heta = size(H_eta,1)

    #Clear the model?
    empty!(model)

    #Create the Optimization Variables
    @variable(model,alpha0)
    @variable(model,Lambda[1:2*n_x,1:n_Heta])

    #Create Optimization Constraints
    @constraint(model,alpha0.>=0)
    @constraint(model,Lambda.>=0)

    S = compute_Sw( system_in.A , time_horizon_in )
    # C_bar = compute_C_M( system_in.C , [time_horizon_in-1] , time_horizon_in )
    J = compute_Sx0( system_in.A , time_horizon_in)

    R_T = [ zeros(n_x,n_x*time_horizon_in) I(n_x) ]

    LHS1 = Lambda * H_eta
    RHS1 = [I(n_x) ; -I(n_x)] * R_T * [ S  zeros(size(S,1),n_y*time_horizon_in) J ]
    @constraint(model, LHS1 .== RHS1 )

    LHS2 = Lambda * h_eta
    RHS2 = alpha0 * ones(2*n_x,1)
    @constraint(model, LHS2 .<= RHS2)

    # Create Objective
    @objective(model,Min,alpha0)

    #Optimize!
    set_silent(model)
    optimize!(model)

    print("termination_status = $(termination_status(model))\n")

    return objective_value(model) , termination_status(model)

end

"""
    find_time_from_X0_to_bound
    Description:
        Identifies the amount of time it takes for the system to reach a bound bound_in.
    Usage:
        T_star = find_time_from_X0_to_bound( ls0 , eta_x0 , bound_in )
"""
function find_time_from_X0_to_bound( system_in::LinearSystem, eta_x0 , bound_in ; iterLimit=10 )
    #Constants
    dim_x = x_dim(system_in)
    dim_y = y_dim(system_in)

    # model = Model(Gurobi.Optimizer)

    # Define Disturbance Polytope
    H_eta, h_eta = define_simple_eta_HPolytope( system_in , eta_x0 , T ) # Create Eta Polyhedron (Polyhedron of all disturbances over time horizon)
    n_Heta = size(H_eta,1)

    #Start the Loop on time horizon T 
    for T = 1:iterLimit

        obj_val0 , term_status0 = norm_of_est_error_after( T , system_in , eta_x0 )

        if obj_val0 > bound_in
            return T-1
        end
        
    end

    return -1

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