"""
Contains all functions used to solve the As Safe As Possible (ASAP) Problem and the As Late As Possible (ALAP) Problem.
"""

using Polyhedra # Polyhedra are represented using Polyhedra.jl
using JuMP
using LinearAlgebra
using Gurobi #, Cbc, GLPK
using SparseArrays
using Random

BigM=1e6 # for indicator constraints

function hypercube(n_dims)
# Return the H-representation of a n_dim dimensional hypercube.
    halfSpaces=Array{HalfSpace{Float64,Array{Float64,1}}}(undef, 2*n_dims)
    for i=1:n_dims
        normalVector1=zeros(n_dims)
        normalVector1[i]=1
        halfSpaces[2*i-1]=HalfSpace(normalVector1,1)

        normalVector2=zeros(n_dims)
        normalVector2[i]=-1
        halfSpaces[2*i]=HalfSpace(normalVector2,1)
    end
    return hrep(halfSpaces)
end

function hyperbox(scalingVector)
# Return the H-representation of a n_dim dimensional hyperbox centered around the origin.
    
    n_dims=length(scalingVector)
    
    halfSpaces=Array{HalfSpace{Float64,Array{Float64,1}}}(undef, 2*n_dims)
    for i=1:n_dims
        normalVector1=zeros(n_dims)
        normalVector1[i]=1
        halfSpaces[2*i-1]=HalfSpace(normalVector1,scalingVector[i])

        normalVector2=zeros(n_dims)
        normalVector2[i]=-1
        halfSpaces[2*i]=HalfSpace(normalVector2,scalingVector[i])
    end
    return hrep(halfSpaces)
end

# To improve
function getHrep(pol)
# Return the H and h of a polyhedra under H-representation.
    z=collect(allhalfspaces(pol))
    dim=fulldim(pol)
    nHalfSpaces=nhalfspaces(pol)
    H=zeros(nHalfSpaces,dim)
    h=zeros(nHalfSpaces)
    for i=1:length(z)
        a=z[i].a'
        beta=z[i].β
        H[i,:]=a
        h[i]=beta
    end
    return H,h
end


"""
sample_polytope
Description:
    This function samples a polytope as specified by polytope_in. 
"""
function sample_polytope( v_polytope_in::VRepresentation, fromVertices::Bool=false )
    # To have field .V in v_polytope_in
    v_polytope_in=convert(MixedMatVRep{Float64,Array{Float64,2}},v_polytope_in)
    
    # Constants
    num_vertices = size(v_polytope_in.V,1)
    
    if fromVertices # return a random vertex
        return v_polytope_in.V[rand(1:num_vertices),:]
    else
        # Algorithm
        theta = randexp(Float64,num_vertices)
        theta = theta/sum(theta) # This vector sums to one. Normalized exponential random variables are uniformly distributed over the symplex.

        # Return our sample a random and convex combination of the vertices
        return transpose(v_polytope_in.V)*theta
    end
end

# To improve
function etaPolyhedra(polW,polV,polX0,T)
# return polW^T x polV^T x polX0
    polEta=polX0
    for t=1:T
        polEta=polV*polEta
    end
    for t=1:T
        polEta=polW*polEta
    end
    return polEta
end

function regularSampling(T,N)
    # enforce regular measurements
    if N==0
        t=[]
    elseif N==1
        t=[round(Int8,(T-1)/2)]
    else
        t=round.(Int8,0:(T-1)/(N-1):T-1)
    end
    sigma=zeros(Bool,T)
    sigma[t.+1].=1
    return sigma
end

function simulateTrajectories(system,F,f;numberSamples=1)
# simulate 'numberSamples' trajectories. Noises are chosen randomly in polytopes
# Works only if all polyhedra (except S) are hypercubes
    n_x=x_dim(system)
    n_y=y_dim(system)
    n_z=z_dim(system)
    n_u=u_dim(system)
    
    T=length(system.list_polS)-1
    
    print("Compute V-representation of W... ")
    vRep_polW=vrep(polyhedron(system.polW))
    println("Done.")
    
    print("Compute V-representation of V... ")
    vRep_polV=vrep(polyhedron(system.polV))
    println("Done.")
    
    print("Compute V-representation of X0... ")
    vRep_polX0=vrep(polyhedron(system.polX0))
    println("Done.")
    
    x=zeros(numberSamples,n_x,T+1) # we don't need to stock the history of x
    z=zeros(numberSamples,n_z,T+1)
    u=zeros(numberSamples,n_u,T)
    for i=1:numberSamples
        y=zeros(n_y,T) # preallocation
        
        # compute x_0
        x[i,:,0+1]= sample_polytope(vRep_polX0,mod(i,2)==0)
        
        # compute z_0
        z[i,:,0+1]=system.D*x[i,:,0+1] + system.d
        
        for t=0:T-1
            # compute y_t
            v=sample_polytope(vRep_polV,mod(i,2)==0) # v_t
            y[:,t+1]=system.C*x[i,:,t+1] + v

            # compute u_t
            u[i,:,t+1]=f[:,t+1]
            for tau=0:t
                u[i,:,t+1]+=F[:,:,t+1,tau+1]*y[:,tau+1]
            end

            # compute x_t+1
            w=sample_polytope(vRep_polW,mod(i,2)==0) # w_t
            x[i,:,t+2]=system.A*x[i,:,t+1] + system.B*u[i,:,t+1] + system.k + w
            
            # compute z_t+1
            z[i,:,t+2]=system.D*x[i,:,t+2] + system.d
        end
    end
    return z,u
end


function computeBounds(system,P_z_all_val, z_tilde_val, P_u_all_val, u_tilde_val, polEta)
# compute theoretical bounds on z and u using [z; u]=P*eta+[z_tilde; u_tilde]
    dim_eta=fulldim(polEta)

    mod=Model(Gurobi.Optimizer)
    @variable(mod, eta[1:dim_eta] in polEta)
    
    # bounds on z
    z=P_z_all_val*eta+z_tilde_val
    len=length(z)
    
    lower_bounds_z=zeros(len)
    upper_bounds_z=zeros(len)

    for i=1:len
        print("\n\n === Compute bounds on z for entry $i ===\n\n")
        @objective(mod,Min,z[i])
        optimize!(mod)
        lower_bounds_z[i]=objective_value(mod)

        @objective(mod,Max,z[i])
        optimize!(mod)
        upper_bounds_z[i]=objective_value(mod)
    end
    
    # bounds on u
    u=P_u_all_val*eta+u_tilde_val
    len=length(u)
    
    lower_bounds_u=zeros(len)
    upper_bounds_u=zeros(len)

    for i=1:len
        print("\n\n === Compute bounds on u for entry $i ===\n\n")
        @objective(mod,Min,u[i])
        optimize!(mod)
        lower_bounds_u[i]=objective_value(mod)

        @objective(mod,Max,u[i])
        optimize!(mod)
        upper_bounds_u[i]=objective_value(mod)
    end
    
    # bounds of S
    n_z=z_dim(system)
    list_polS=system.list_polS
    T=length(list_polS)-1
    
    lower_bounds_S=zeros(n_z,T+1)
    upper_bounds_S=zeros(n_z,T+1)
    
    for t=0:T
        mod=Model(Gurobi.Optimizer)
        @variable(mod, s[1:n_z] in list_polS[t+1])
        for i=1:n_z
            print("\n\n === Compute bounds of S_$t for component $i ===\n\n")
            @objective(mod,Min,s[i])
            optimize!(mod)
            lower_bounds_S[i,t+1]=objective_value(mod)
            
            @objective(mod,Max,s[i])
            optimize!(mod)
            upper_bounds_S[i,t+1]=objective_value(mod)
        end
    end
    
    # bounds of U (independent of t)
    n_u=u_dim(system)
    
    lower_bounds_U=zeros(n_u)
    upper_bounds_U=zeros(n_u)
    
    mod=Model(Gurobi.Optimizer)
    @variable(mod, input[1:n_u] in system.polU)
    for i=1:n_u
        print("\n\n === Compute bounds of U for component $i ===\n\n")
        @objective(mod,Min,input[i])
        optimize!(mod)
        lower_bounds_U[i]=objective_value(mod)
        
        @objective(mod,Max,input[i])
        optimize!(mod)
        upper_bounds_U[i]=objective_value(mod)
    end
            
    # reshaping
    lower_bounds_z=reshape(lower_bounds_z,(n_z,T+1))
    upper_bounds_z=reshape(upper_bounds_z,(n_z,T+1))
    lower_bounds_u=reshape(lower_bounds_u,(n_u,T))
    upper_bounds_u=reshape(upper_bounds_u,(n_u,T))
    
    return lower_bounds_z, upper_bounds_z, lower_bounds_u, upper_bounds_u, lower_bounds_S, upper_bounds_S, lower_bounds_U, upper_bounds_U
end

function computePmatrix(system,T, Q,r)
# Return P_z_all=[P_z_w P_z_v P_z_x0], z_tilde, P_u_all=[P_u_w P_u_v P_u_x0] and u_tilde
    
    A=system.A
    B=system.B
    C=system.C
    D=system.D
    k=system.k
    d=system.d
    
    n_x=x_dim(system)
    n_y=y_dim(system)
    n_z=z_dim(system)
    n_u=u_dim(system)
    
    # compute the H matrix
    H=zeros(n_x*(T+1),n_x*T)
    l=one(A) 
    H[n_x+1:2*n_x,1:n_x]=l # A^0
    for i=1:T-1 # loop over the lines
        l=[A*l[1:n_x,1:n_x] l] # A*l[1:nx,1:nx] is A^i
        H[(i+1)*n_x+1:(i+2)*n_x, 1:(i+1)*n_x]=l
    end

    # compute the S matrix (not related to the safety set S)
    S=zeros(n_x*(T+1),n_u*T)
    l=B
    S[n_x+1:2*n_x,1:n_u]=l # A^0*B
    for i=1:T-1 # loop over the lines
        l=[A*l[1:n_x,1:n_u] l] # A*l[1:nx,1:nx] is A^i
        S[(i+1)*n_x+1:(i+2)*n_x, 1:(i+1)*n_u]=l
    end

    # compute C_bar
    C_bar=[kron(Matrix(I,T,T),C) zeros(T*n_y,n_x)]

    # compute the J matrix
    J=zeros((T+1)*n_x,n_x)
    l=one(A) # Identity matrix
    J[1:n_x,1:n_x]=l
    for i = 1:T
        l=A*l # l=A^i
        J[i*n_x+1:(i+1)*n_x, 1:n_x] = l 
    end

    # compute k_tilde
    k_tilde=kron(ones(T),k)

    # compute N
    N=Matrix(I,T*n_y,T*n_y)
    
    # compute d_tilde
    d_tilde=kron(ones(T+1),d)
    
    # compute D_bar
    D_bar=kron(Matrix(I,T+1,T+1),D)
    
    
    # Compute matrix P
    P_z_w=D_bar*(H+S*Q*C_bar*H)
    P_z_v=D_bar*S*Q*N
    P_z_x0=D_bar*( Matrix(I,n_x*(T+1),n_x*(T+1)) + S*Q*C_bar )*J
    P_z_all=[P_z_w P_z_v P_z_x0]
    
    z_tilde=D_bar*( Matrix(I,n_x*(T+1),n_x*(T+1)) + S*Q*C_bar )*H*k_tilde + D_bar*S*r + d_tilde
    
    P_u_w=Q*C_bar*H
    P_u_v=Q*N
    P_u_x0=Q*C_bar*J
    P_u_all=[P_u_w P_u_v P_u_x0]
    
    u_tilde=Q*C_bar*H*k_tilde+r
    
    return P_z_all, z_tilde, P_u_all, u_tilde, C_bar, S
end

# ASAP functions

function instantiateASAP(system,N_m,N_c,T;minimizeScaling=false,enforceRegular=false)
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

function getResultsASAP(model,C_bar,S_bar,N_m,N_c,T;minimizeScaling=false,enforceRegular=false)
# Returns the solution of the As Safe As Possible problem (Control gains, optimal alpha,  measurement and control signals)
    if minimizeScaling
        alpha_opt=value(model[:alpha])
    else
        alpha_opt=1
    end
    
    if enforceRegular
        sigma_meas_opt=regularSampling(T,N_m)
        sigma_control_opt=regularSampling(T,N_c)
    else
        sigma_meas_opt=value.(model[:sigma_meas])
        sigma_meas_opt=sigma_meas_opt.>=0.9 # to have binary

        sigma_control_opt=value.(model[:sigma_control])
        sigma_control_opt=sigma_control_opt.>=0.9 # to have binary
    end
    
    Q_opt=value.(model[:Q]);
    r_opt=value.(model[:r]);
    
    l1,l2=size(Q_opt)
    n_u=convert(Int64,l1/T)
    n_y=convert(Int64,l2/T)
    
    #F_bigMatrix = inv(I + Q_opt*C_bar*S_bar ) * Q_opt # slow and unstable
    F_bigMatrix = sparse(I + Q_opt*C_bar*S_bar) \ Q_opt
    
    F=zeros(n_u,n_y,T,T) # possible with a reshape ? Or with only one loop ?
    for t=0:T-1
        for tau=0:t
            F[:,:,t+1,tau+1]=F_bigMatrix[1+t*n_u:(t+1)*n_u , 1+tau*n_y:(tau+1)*n_y]
        end
    end

    #f=inv(I + Q_opt*C_bar*S_bar ) * r_opt; # slow and unstable
    f=sparse(I + Q_opt*C_bar*S_bar ) \ r_opt
    f=reshape(f,(n_u,T))
    
    return F, f, alpha_opt, sigma_meas_opt, sigma_control_opt
end

# ALAP functions

function binarySearchALAP(system,N_m,N_c,T_max;enforceRegular=false)
# Solve the As Late As Possible Problem by doing a binary search. In the binary search, each evaluation requires to solve one instance of ASAP.
# We assume that polX0 is included in polS, then T >= 1
    
    original_list_polS=copy(system.list_polS) # keep a copy because the list will be modified when solving ASAP for T < T_max.
    
    
    # solve for T=1
    current_T=1
    
    append!(empty!(system.list_polS), original_list_polS[1:current_T+1])
    
    ASAPmodel_lb, C_bar_lb, S_lb, _, _, _= instantiateASAP(system,N_m,N_c,current_T; enforceRegular=enforceRegular,minimizeScaling=false)
    
    set_optimizer(ASAPmodel_lb, Gurobi.Optimizer)
    optimize!(ASAPmodel_lb)
    print("\nBINARY SEARCH: current_T=$(current_T).\n\n")
    
    if termination_status(ASAPmodel_lb) != MOI.TerminationStatusCode(1) # not feasible
        print("\nCAUTION: Maybe T_opt=0: We have assumed that X_0 is in S.\n")
        
        # restore the original system.list_polS
        append!(empty!(system.list_polS), original_list_polS)
        
        return ASAPmodel_lb, 0, C_bar_lb, S_lb
    end
    
    # solve for T=T_max
    current_T=T_max
    
    append!(empty!(system.list_polS), original_list_polS[1:current_T+1])
    
    ASAPmodel, C_bar, S, _,_,_=instantiateASAP(system,N_m,N_c,current_T; enforceRegular=enforceRegular, minimizeScaling=false)
    set_optimizer(ASAPmodel, Gurobi.Optimizer)
    optimize!(ASAPmodel)
    print("\nBINARY SEARCH: current_T=$(current_T).\n\n")
    
    if termination_status(ASAPmodel) == MOI.TerminationStatusCode(1) # feasible
        
        # restore the original system.list_polS
        append!(empty!(system.list_polS), original_list_polS)
        
        return ASAPmodel, current_T, C_bar, S
    end
    
    # lb <= T_opt < ub holds at each iteration
    lb=1
    ub=T_max
    while ub-lb>1
        current_T=round(Int8,(ub+lb)/2)
        
        append!(empty!(system.list_polS), original_list_polS[1:current_T+1])
        
        ASAPmodel, C_bar, S, _, _, _=instantiateASAP(system,N_m,N_c,current_T; enforceRegular=enforceRegular,minimizeScaling=false)
        set_optimizer(ASAPmodel, Gurobi.Optimizer)
        optimize!(ASAPmodel)
        if termination_status(ASAPmodel) != MOI.TerminationStatusCode(1) # not feasible
            ub=current_T
        else
            # lb is the new best guess for the optimal T
            lb=current_T
            ASAPmodel_lb=ASAPmodel
            C_bar_lb=C_bar
            S_lb=S
        end
        print("\nBINARY SEARCH: current_T=$(current_T). T_opt is in [$(lb) $(ub)[\n\n")
    end
    
    # restore the original system.list_polS
    append!(empty!(system.list_polS), original_list_polS)
    
    # lb is the optimal T
    return ASAPmodel_lb, lb, C_bar_lb, S_lb
end


function getResultsALAP(system, ASAPmodel,T_max, C_bar_ASAP,S_ASAP, N_m,N_c,T_opt; enforceRegular=false)
# Returns the solution of the As Long As Possible problem (Control gains, optimal alpha,  measurement and control signals)
    
    F_ASAP, f_ASAP, alpha_ASAP, sigma_meas_ASAP, sigma_control_ASAP = getResultsASAP(ASAPmodel,C_bar_ASAP,S_ASAP,N_m,N_c,T_opt;minimizeScaling=false,enforceRegular=enforceRegular)
    
    n_x=x_dim(system)
    n_y=y_dim(system)
    #n_z=z_dim(system)
    n_u=u_dim(system)
    
    sigma_meas_opt=[sigma_meas_ASAP; zeros(Int8,T_max-T_opt)]
    sigma_control_opt=[sigma_control_ASAP; zeros(Int8,T_max-T_opt)]
    
    F=zeros(n_u, n_y, T_max, T_max)
    F[:,:,1:T_opt,1:T_opt]=F_ASAP
    
    f=zeros(n_u,T_max)
    f[:,1:T_opt]=f_ASAP
    
    
    # Impose zero-order old
    for t=T_opt:T_max-1 # index of t in matrices is t+1
        F[:,:,t+1,:]=F[:,:,(t-1)+1,:]
        f[:,t+1]=f[:,(t-1)+1]
    end
    
    
    Q_ASAP=value.(ASAPmodel[:Q])
    r_ASAP=value.(ASAPmodel[:r])
    
    Q=zeros(T_max*n_u,T_max*n_y)
    Q[1:T_opt*n_u,1:T_opt*n_y]=Q_ASAP
    
    r=zeros(T_max*n_u)
    r[1:T_opt*n_u]=r_ASAP
    
    # Impose zero-order hold
    for t=T_opt:T_max-1 # index of t in matrices is t+1
        Q[1+t*n_u:(t+1)*n_u,:]=Q[1+(t-1)*n_u:t*n_u,:]
        r[1+t*n_u:(t+1)*n_u]=r[1+(t-1)*n_u:t*n_u]
    end
    
    P_z_all, z_tilde, P_u_all, u_tilde, _, _ = computePmatrix(system,T_max, Q,r)
    
    polEta=etaPolyhedra(system.polW,system.polV,system.polX0,T_max)
    
    return F, f, sigma_meas_opt, sigma_control_opt, P_z_all, z_tilde, P_u_all, u_tilde, polEta # 5 last quantities are to compute polz = P_z_all*polEta + z_tilde and polu = P_u_all*polEta + u_tilde
end




""" 
======== OLD FUNCTIONS ==========
"""

"""
function instantiateALAP(A,B,C,k,N_m,N_c,T_max,alpha)

    n_y,n_x=size(C)
    _,n_u=size(B)

    polW=hypercube(n_x)
    polV=hypercube(n_y)
    polX0=hypercube(n_x)
    polS=alpha*hypercube(n_x)

    n_w=nhalfspaces(polW)
    n_v=nhalfspaces(polV)
    n_x0=nhalfspaces(polX0)
    n_s=nhalfspaces(polS)

    polEta=etaPolyhedra(polW,polV,polX0,T_max)
    n_eta=nhalfspaces(polEta)
    
    model = Model()
    @variable(model,Q[1:n_u*T_max,1:n_y*T_max])
    @variable(model,r[1:n_u*T_max])
    @variable(model,Lambda_rob[1:T_max+1,1:n_s,1:n_eta])
    @variable(model,sigma_meas[1:T_max],Bin) # simga_meas[i] is sigma^m_(i-1)
    @variable(model,sigma_control[1:T_max],Bin) # simga_control[i] is sigma^c_(i-1)
    @variable(model,delta[1:T_max+1],Bin)
    @variable(model,mu[1:T_max+1],Bin)

    # compute the H matrix
    H=zeros(n_x*(T_max+1),n_x*T_max)
    l=one(A) 
    H[n_x+1:2*n_x,1:n_x]=l # A^0
    for i=1:T_max-1 # loop over the lines
        l=[A*l[1:n_x,1:n_x] l] # A*l[1:nx,1:nx] is A^i
        H[(i+1)*n_x+1:(i+2)*n_x, 1:(i+1)*n_x]=l
    end

    # compute the S matrix
    S=zeros(n_x*(T_max+1),n_u*T_max)
    l=B
    S[n_x+1:2*n_x,1:n_u]=l # A^0*B
    for i=1:T_max-1 # loop over the lines
        l=[A*l[1:n_x,1:n_u] l] # A*l[1:nx,1:nx] is A^i
        S[(i+1)*n_x+1:(i+2)*n_x, 1:(i+1)*n_u]=l
    end

    # compute C_bar
    C_bar=[kron(Matrix(I,T_max,T_max),C) zeros(T_max*n_y,n_x)]

    # compute the J matrix
    J=zeros((T_max+1)*n_x,n_x)
    l=one(A) # Identity matrix
    J[1:n_x,1:n_x]=l
    for i = 1:T_max
        l=A*l # l=A^i
        J[i*n_x+1:(i+1)*n_x, 1:n_x] = l 
    end

    # compute k_tilde
    k_tilde=kron(ones(T_max),k)

    # compute N
    N=Matrix(I,T_max*n_y,T_max*n_y)

    # Compute matrix P
    P_x_w=H+S*Q*C_bar*H
    P_x_v=S*Q*N
    P_x_x0=( Matrix(I,n_x*(T_max+1),n_x*(T_max+1)) + S*Q*C_bar )*J
    x_tilde=(Matrix(I,n_x*(T_max+1),n_x*(T_max+1))+S*Q*C_bar)*H*k_tilde+S*r

    P_u_w=Q*C_bar*H
    P_u_v=Q*N
    P_u_x0=Q*C_bar*J
    u_tilde=Q*C_bar*H*k_tilde+r

    # find H-representations
    (H_eta,h_eta)=getHrep(polEta)
    (H_S,h_S)=getHrep(polS)
    #H_all_S=kron(Matrix(I,T_max+1,T_max+1),H_S)
    #h_all_S=kron(ones(T_max+1),h_S)

    # Budget constraints
    @constraint(model,sum(sigma_meas)<=N_m)
    @constraint(model,sum(sigma_control)<=N_c)

    # Constraints on delta and mu
    @constraint(model,delta.<=mu)
    @constraint(model,delta[1:T_max].>=delta[2:T_max+1]) # delta is decreasing

    # Robustness constraint (indicator constraint)
    for t=0:T_max
        R_t=[zeros(n_x,n_x*t) I zeros(n_x,n_x*(T_max-t))]
        @constraint(model,Lambda_rob[t+1,:,:]*H_eta .<= H_S*R_t*[P_x_w P_x_v P_x_x0] + ones(n_s,(T_max+1)*n_x+T_max*n_y)*(1-mu[t+1])*BigM)

        LHS=Lambda_rob[t+1,:,:]*H_eta + ones(n_s,(T_max+1)*n_x+T_max*n_y)*(1-mu[t+1])*BigM # to avoid a bug
        @constraint(model,LHS .>= H_S*R_t*[P_x_w P_x_v P_x_x0])

        @constraint(model,Lambda_rob[t+1,:,:]*h_eta .<= h_S-H_S*R_t*x_tilde + ones(n_s)*(1-mu[t+1])*BigM)
    end

    # Measurement constraint (indicator constraint)
    for t=0:T_max-1
        I_t=1+t*n_y:(t+1)*n_y
        @constraint(model,  Q[:,I_t] .<= ones(n_u*T_max,n_y)*sigma_meas[t+1]*BigM )
        @constraint(model, -Q[:,I_t] .<= ones(n_u*T_max,n_y)*sigma_meas[t+1]*BigM )
    end

    # Control constraint (indicator constraint)
    # special case t=0
    J_t_1=1:n_u # J_{t-1}
    @constraint(model,  Q[J_t_1,:] .<= ones(n_u,n_y*T_max)*sigma_control[1]*BigM)
    @constraint(model, -Q[J_t_1,:] .<= ones(n_u,n_y*T_max)*sigma_control[1]*BigM)
    @constraint(model,  r[J_t_1] .<= ones(n_u)*sigma_control[1]*BigM)
    @constraint(model, -r[J_t_1] .<= ones(n_u)*sigma_control[1]*BigM)
    for t=1:T_max-1
        J_t=1+t*n_u:(t+1)*n_u
        @constraint(model,  Q[J_t,:]-Q[J_t_1,:] .<= ones(n_u,n_y*T_max)*sigma_control[t+1]*BigM)
        @constraint(model, -Q[J_t,:]+Q[J_t_1,:] .<= ones(n_u,n_y*T_max)*sigma_control[t+1]*BigM)
        @constraint(model,  r[J_t]-r[J_t_1] .<= ones(n_u)*sigma_control[t+1]*BigM)
        @constraint(model, -r[J_t]+r[J_t_1] .<= ones(n_u)*sigma_control[t+1]*BigM)
        J_t_1=J_t
    end

    # positivity constraint
    @constraint(model,Lambda_rob .>= zeros(size(Lambda_rob)))

    # low diagonal constraint
    for blk_row_idx = 1:div(size(Q,1), n_x)-1
        LHS_bloc = Q[(blk_row_idx-1)*n_x.+(1:n_x),blk_row_idx*n_y+1:end]
        @constraint(model,LHS_bloc .== zeros(size(LHS_bloc)))
    end

    @objective(model,Max,sum(delta))
    
    return model, C_bar, S
end

function greedyALAP(A,B,C,k,N_m,N_c,T_max,alpha)
    model,_,_=instantiateALAP(A,B,C,k,1,1,T_max,alpha)
    set_optimizer(model, Gurobi.Optimizer)
    optimize!(model)
    
    ind_meas=findall(value.(model[:sigma_meas]).>=0.9)
    ind_control=findall(value.(model[:sigma_control]).>=0.9)
    
    N_min=min(N_m,N_c)
    for N_loc=2:N_min
        model,_,_=instantiateALAP(A,B,C,k,N_loc,N_loc,T_max,alpha)
        @constraint(model,model[:sigma_meas][ind_meas].==1)
        @constraint(model,model[:sigma_control][ind_control].==1)
        set_optimizer(model, Gurobi.Optimizer)
        optimize!(model)
        ind_meas=findall(value.(model[:sigma_meas]).>=0.9)
        ind_control=findall(value.(model[:sigma_control]).>=0.9)
    end
    
    if N_m>N_c
        for N_m_loc=N_min+1:N_m
            model,_,_=instantiateALAP(A,B,C,k,N_m_loc,N_c,T_max,alpha)
            @constraint(model,model[:sigma_meas][ind_meas].==1)
            @constraint(model,model[:sigma_control][ind_control].==1)
            set_optimizer(model, Gurobi.Optimizer)
            optimize!(model)
            ind_meas=findall(value.(model[:sigma_meas]).>=0.9)
        end
    elseif N_c>N_m
        for N_c_loc=N_min+1:N_c
            model,_,_=instantiateALAP(A,B,C,k,N_m,N_c_loc,T_max,alpha)
            @constraint(model,model[:sigma_meas][ind_meas].==1)
            @constraint(model,model[:sigma_control][ind_control].==1)
            set_optimizer(model, Gurobi.Optimizer)
            optimize!(model)
            ind_control=findall(value.(model[:sigma_control]).>=0.9)
        end
    end
    
    return model
end



function getResultsALAP_old(model,C_bar,S)
    sigma_meas=model[:sigma_meas]
    sigma_control=model[:sigma_control]
    Q=model[:Q]
    r=model[:r]
    
    T_opt=round(Int8,objective_value(model)-1)
    sigma_meas_opt=value.(sigma_meas)
    sigma_meas_opt=sigma_meas_opt.>=0.9 # to have binary
    sigma_control_opt=value.(sigma_control)
    sigma_control_opt=sigma_control_opt.>=0.9 # to have binary
    Q_opt=value.(Q);
    r_opt=value.(r);
    
    T_max=length(sigma_meas_opt)
    l1,l2=size(Q)
    n_u=convert(Int64,l1/T_max)
    n_y=convert(Int64,l2/T_max)
    
    F_bigMatrix = inv(I + Q_opt*C_bar*S ) * Q_opt; # slow
    F=zeros(n_u,n_y,T_max,T_max) # possible with a reshape ? Or with only one loop ?
    for t=0:T_opt-1 # values for t>=T_opt are arbitrary and are set to 0
        for tau=0:t
            F[:,:,t+1,tau+1]=F_bigMatrix[1+t*n_u:(t+1)*n_u , 1+tau*n_y:(tau+1)*n_y]
        end
    end

    f=inv(I + Q_opt*C_bar*S ) * r_opt; # slow
    f=reshape(f,(n_u,T_max))
    f[:,T_opt+1:end].=0 # values for t>=T_opt are arbitrary and are set to 0
    
    return F, f, T_opt, sigma_meas_opt, sigma_control_opt
end


function instantiateASAP_timeVaryingS(system,N_m,N_c,T;minimizeScaling=false,enforceRegular=false)
    #WRONG union safety set approach
    #This version of instatiateASAP can deal with time-varying S and with S=Union S_i.
    #
    #Construct a JuMP model of the MILP problem.
    #If enforceRegular, sigma[t]=1 for Round( t=k*(T-1)/(N-1) ), t=0,...,N-1. Or t=Round( (T-1)/2 ) if N=1.
    #If minimize scaling, minimize alpha such that z_t in alpha*polS for all t
    
    
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
    polU=system.polU

    n_w=nhalfspaces(polW)
    n_v=nhalfspaces(polV)
    n_x0=nhalfspaces(polX0)
    n_U=nhalfspaces(polU) # different than n_u
    
    polEta=etaPolyhedra(polW,polV,polX0,T)
    n_eta=nhalfspaces(polEta)

    model = Model()
    if minimizeScaling
        @variable(model,alpha)
    else
        alpha=1
    end
    @variable(model,Q[1:n_u*T,1:n_y*T])
    @variable(model,r[1:n_u*T])
    @variable(model,Gamma[1:T*n_U,1:n_eta])
    
    if enforceRegular
        sigma_meas=regularSampling(T,N_m) # not a JuMP variable
        sigma_control=regularSampling(T,N_c)
    else
        @variable(model,sigma_meas[1:T],Bin) # simga_meas[i] is sigma^m_(i-1)
        @variable(model,sigma_control[1:T],Bin) # simga_control[i] is sigma^c_(i-1)
        
        # Budget constraints
        @constraint(model,sum(sigma_meas)<=N_m)
        @constraint(model,sum(sigma_control)<=N_c)
    end
    
    P_z_all,z_tilde,P_u_all,u_tilde, C_bar,S=computePmatrix(system,T, Q,r)

    # find H-representations
    (H_eta,h_eta)=getHrep(polEta)
    
    (H_U,h_U)=getHrep(polU)
    H_all_U=kron(Matrix(I,T,T),H_U)
    h_all_U=kron(ones(T),h_U)
    
    # Safety (robustness) constraint
    list_S=system.list_polS
    
    for t=0:T
        N_S_t=length(list_S[t+1]) # Number of safety sets at time t
        mu = @variable(model,[1:N_S_t],Bin)
        @constraint(model,sum(mu)>=1) # sum_i mu_t^i >= 1
        I_t=1+t*n_z:(t+1)*n_z
        for i=1:N_S_t # add indicator constraints
            polS_t_i=list_S[t+1][i]
            H_S_t_i,h_S_t_i=getHrep(polS_t_i)
            n_S_t_i=nhalfspaces(polS_t_i)
            
            
            Lambda_t_i=@variable(model,[1:n_S_t_i,1:n_eta])
            @constraint(model,Lambda_t_i.>=zeros(size(Lambda_t_i))) # Lambda>=0
            
            # constraints on H
            LHS1 = Lambda_t_i*H_eta
            RHS1 = H_S_t_i*P_z_all[I_t,:] + ones(size(LHS1))*(1-mu[i])*BigM # ERROR HERE
            
            @constraint(model,LHS1 .<= RHS1)
            
            LHS2 = Lambda_t_i*H_eta + ones(size(LHS1))*(1-mu[i])*BigM
            
            @constraint(model, LHS2 .>= H_S_t_i*P_z_all[I_t,:])

            # constraint on h
            LHS3 = Lambda_t_i*h_eta

            @constraint(model, LHS3 .<= alpha*h_S_t_i - H_S_t_i*z_tilde[I_t] + ones(size(LHS3))*(1-mu[i])*BigM)
        end
    end
    
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
    @constraint(model,Gamma .>= zeros(size(Gamma)))

    # low diagonal constraint
    for blk_row_idx = 1:div(size(Q,1), n_x)-1
        LHS_bloc = Q[(blk_row_idx-1)*n_x.+(1:n_x),blk_row_idx*n_y+1:end]
        @constraint(model,LHS_bloc .== zeros(size(LHS_bloc)))
    end

    if minimizeScaling
        @objective(model,Min,alpha)
    end
    
    return model, C_bar, S, P_z_all, z_tilde, P_u_all, u_tilde, polEta # 5 last quantities are to compute polZ = P_z_all*polEta + z_tilde and polU = P_u_all*polEta + u_tilde
end

"""
