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
    IntTimeChunk = ceil(FloatTimeChunk);

    M = []

    for M_idx = 1:MeasurementBudget
        append!(M,TimeHorizon - Int(M_idx*IntTimeChunk) )
    end

    return reverse(M)

end

"""
sample_polytope
Description:
    This function samples a polytope as specified by polytope_in. 
"""
function sample_polytope( polytope_in::Polyhedron )
    
    # Constants

    vrep_of_polytope_in = vrep(polytope_in)
    num_vertices = size(vrep_of_polytope_in.V,1)

    # Algorithm
    theta = rand(Float64,num_vertices)
    theta = theta/sum(theta) #This vector sums to one

    # Return our sample a random and convex combination of the vertices
    return transpose(tempV2.V)*theta 

end

"""
ConcreteScheduleParameters
Description:
    Defines a simple object which can store some of the 
"""
struct ConcreteScheduleParameters
    ControlSchedule::Array{Float64}
    MeasurementSchedule::Array{Float64}
end

function instantiateASAP( system , N_m , N_c , T )
    """
    Construct a JuMP model of the MILP problem.
    If enforceRegular, sigma[t]=1 for Round( t=k*(T-1)/(N-1) ), t=0,...,N-1. Or t=Round( (T-1)/2 ) if N=1.
    If minimize scaling, minimize alpha such that z_t in alpha*polS for all t
    """
    
    n_x=x_dim(system)
    n_y=y_dim(system)
    n_z=z_dim(system)
    n_u=u_dim(system)
    
    # avoid errors when enforceRegular
    N_m=min(T,N_m)
    N_c=min(T,N_c)
    
    polW=system.polW
    polV=system.polV
    polX0=system.polX0
    
    polEta=etaPolyhedra(polW,polV,polX0,T)
    (H_eta,h_eta)=getHrep(polEta)
    
    polU=system.polU
    (H_U,h_U)=getHrep(polU)
    H_all_U=kron(Matrix(I,T,T),H_U)
    h_all_U=kron(ones(T),h_U)
    
    pol_all_S=prod(system.list_polS) # cartesian product polS_0 X polS_1 X ... X polS_T
    (H_all_S,h_all_S)=getHrep(pol_all_S)
    
    n_w=nhalfspaces(polW)
    n_v=nhalfspaces(polV)
    n_x0=nhalfspaces(polX0)
    n_S=nhalfspaces(pol_all_S)
    n_U=nhalfspaces(polU) # different than n_u
    n_eta=nhalfspaces(polEta)
    
    model = Model()
    if minimizeScaling
        @variable(model,alpha)
    else
        alpha=1
    end
    @variable(model,Q[1:n_u*T,1:n_y*T])
    @variable(model,r[1:n_u*T])
    @variable(model,Lambda_rob[1:n_S,1:n_eta])
    @variable(model,Gamma[1:T*n_U,1:n_eta])
    
    if enforceRegular
        sigma_meas=regularSampling(T,N_m) # not a JuMP variable
        sigma_control=regularSampling(T,N_c)
    else
        @variable(model,sigma_meas[1:T],Bin) # sigma_meas[i] is sigma^m_(i-1)
        @variable(model,sigma_control[1:T],Bin) # sigma_control[i] is sigma^c_(i-1)
        
        # Budget constraints
        @constraint(model,sum(sigma_meas)<=N_m)
        @constraint(model,sum(sigma_control)<=N_c)
    end
    
    P_z_all,z_tilde, P_u_all,u_tilde, C_bar,S_bar=computePmatrix(system,T, Q,r)
    
    # Safety (robustness) constraint
    @constraint(model,Lambda_rob*H_eta .== H_all_S*P_z_all)
    @constraint(model,Lambda_rob*h_eta .<= alpha*h_all_S-H_all_S*z_tilde)
    
    # Bounded inputs constraint
    @constraint(model,Gamma*H_eta .== H_all_U*P_u_all)
    @constraint(model,Gamma*h_eta .<= h_all_U-H_all_U*u_tilde)

    # Measurement constraint (indicator constraint)
    for t=0:T-1
        I_t=1+t*n_y:(t+1)*n_y
        @constraint(model,  Q[:,I_t] .<= ones(n_u*T,n_y)*sigma_meas[t+1]*BigM )
        @constraint(model, -Q[:,I_t] .<= ones(n_u*T,n_y)*sigma_meas[t+1]*BigM )
    end

    # Control constraint (indicator constraint)
    # special case t=0
    J_t_1=1:n_u # J_{t-1}
    @constraint(model,  Q[J_t_1,:] .<= ones(n_u,n_y*T)*sigma_control[1]*BigM)
    @constraint(model, -Q[J_t_1,:] .<= ones(n_u,n_y*T)*sigma_control[1]*BigM)
    @constraint(model,  r[J_t_1] .<= ones(n_u)*sigma_control[1]*BigM)
    @constraint(model, -r[J_t_1] .<= ones(n_u)*sigma_control[1]*BigM)
    for t=1:T-1
        J_t=1+t*n_u:(t+1)*n_u
        @constraint(model,  Q[J_t,:]-Q[J_t_1,:] .<= ones(n_u,n_y*T)*sigma_control[t+1]*BigM)
        @constraint(model, -Q[J_t,:]+Q[J_t_1,:] .<= ones(n_u,n_y*T)*sigma_control[t+1]*BigM)
        @constraint(model,  r[J_t]-r[J_t_1] .<= ones(n_u)*sigma_control[t+1]*BigM)
        @constraint(model, -r[J_t]+r[J_t_1] .<= ones(n_u)*sigma_control[t+1]*BigM)
        J_t_1=J_t
    end

    # positivity constraints
    if minimizeScaling
        @constraint(model,alpha >= 0)
    end
    @constraint(model,Lambda_rob .>= zeros(size(Lambda_rob)))
    @constraint(model,Gamma .>= zeros(size(Gamma)))
    
    
    for blk_row_idx = 1:div(size(Q,1), n_u)-1 # should it be n_u instead of n_x
        LHS_bloc = Q[(blk_row_idx-1)*n_u.+(1:n_u),blk_row_idx*n_y+1:end]
        @constraint(model,LHS_bloc .== zeros(size(LHS_bloc)))
    end
    
    if minimizeScaling
        @objective(model,Min,alpha)
    end
    
    return model, C_bar, S_bar, P_z_all, z_tilde, P_u_all, u_tilde, polEta # 5 last quantities are to compute polZ = P_z_all*polEta + z_tilde and polU = P_u_all*polEta + u_tilde
end