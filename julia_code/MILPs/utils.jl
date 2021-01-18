using Polyhedra
using JuMP
using LinearAlgebra

function hypercube(n_dims)
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

# Very unefficient (and probably stupid)
function getHrep(pol)
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


function etaPolyhedra(polW,polV,polX0,T) # could be faster
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


function instantiateASAP(A,B,C,k,N_m,N_c,T)
    BigM=1e5
    
    n_y,n_x=size(C)
    _,n_u=size(B)

    polW=hypercube(n_x)
    polV=hypercube(n_y)
    polX0=hypercube(n_x)
    polS=hypercube(n_x)

    n_w=nhalfspaces(polW)
    n_v=nhalfspaces(polV)
    n_x0=nhalfspaces(polX0)
    n_s=nhalfspaces(polS)

    polEta=etaPolyhedra(polW,polV,polX0,T)
    n_eta=nhalfspaces(polEta)

    model = Model()
    @variable(model,alpha)
    @variable(model,Q[1:n_u*T,1:n_y*T])
    @variable(model,r[1:n_u*T])
    @variable(model,Lambda_rob[1:(T+1)*n_s,1:n_eta])
    @variable(model,sigma_meas[1:T],Bin) # simga_meas[i] is sigma^m_(i-1)
    @variable(model,sigma_control[1:T],Bin) # simga_control[i] is sigma^c_(i-1)

    # compute the H matrix
    H=zeros(n_x*(T+1),n_x*T)
    l=one(A) 
    H[n_x+1:2*n_x,1:n_x]=l # A^0
    for i=1:T-1 # loop over the lines
        l=[A*l[1:n_x,1:n_x] l] # A*l[1:nx,1:nx] is A^i
        H[(i+1)*n_x+1:(i+2)*n_x, 1:(i+1)*n_x]=l
    end

    # compute the S matrix
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

    # Compute matrix P
    P_x_w=H+S*Q*C_bar*H
    P_x_v=S*Q*N
    P_x_x0=( Matrix(I,n_x*(T+1),n_x*(T+1)) + S*Q*C_bar )*J
    x_tilde=(Matrix(I,n_x*(T+1),n_x*(T+1))+S*Q*C_bar)*H*k_tilde+S*r

    P_u_w=Q*C_bar*H
    P_u_v=Q*N
    P_u_x0=Q*C_bar*J
    u_tilde=Q*C_bar*H*k_tilde+r

    # find H-representations
    (H_eta,h_eta)=getHrep(polEta)
    (H_S,h_S)=getHrep(polS)
    H_all_S=kron(Matrix(I,T+1,T+1),H_S)
    h_all_S=kron(ones(T+1),h_S)

    # Budget constraints
    @constraint(model,sum(sigma_meas)==N_m)
    @constraint(model,sum(sigma_control)==N_c)

    # Robustness constraint
    @constraint(model,Lambda_rob*H_eta .== H_all_S*[P_x_w P_x_v P_x_x0])
    @constraint(model,Lambda_rob*h_eta .<= alpha*h_all_S-H_all_S*x_tilde)

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

    # positivity constraint
    @constraint(model,alpha >= 0)
    @constraint(model,Lambda_rob .>= zeros(size(Lambda_rob)))

    # low diagonal constraint
    for blk_row_idx = 1:div(size(Q,1), n_x)-1
        LHS_bloc = Q[(blk_row_idx-1)*n_x.+(1:n_x),blk_row_idx*n_y+1:end]
        @constraint(model,LHS_bloc .== zeros(size(LHS_bloc)))
    end

    @objective(model,Min,alpha)
    
    return model, C_bar, S
end

function instantiateALAP(A,B,C,k,N_m,N_c,T_max,alpha)
    BigM=1e5

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
    @constraint(model,sum(sigma_meas)==N_m)
    @constraint(model,sum(sigma_control)==N_c)

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


function getResultsASAP(model,C_bar,S)
    alpha=model[:alpha]
    sigma_meas=model[:sigma_meas]
    sigma_control=model[:sigma_control]
    Q=model[:Q]
    r=model[:r]
    
    alpha_opt=value(alpha)
    sigma_meas_opt=value.(sigma_meas)
    sigma_meas_opt=sigma_meas_opt.>=0.9 # to have binary
    sigma_control_opt=value.(sigma_control)
    sigma_control_opt=sigma_control_opt.>=0.9 # to have binary
    Q_opt=value.(Q);
    r_opt=value.(r);
    
    T=length(sigma_meas_opt)
    l1,l2=size(Q)
    n_u=convert(Int64,l1/T)
    n_y=convert(Int64,l2/T)
    

    F_bigMatrix = inv(I + Q_opt*C_bar*S ) * Q_opt; # slow
    F=zeros(n_u,n_y,T,T) # possible with a reshape ? Or with only one loop ?
    for t=0:T-1
        for tau=0:t
            F[:,:,t+1,tau+1]=F_bigMatrix[1+t*n_u:(t+1)*n_u , 1+tau*n_y:(tau+1)*n_y]
        end
    end

    f=inv(I + Q_opt*C_bar*S ) * r_opt; # slow
    f=reshape(f,(n_u,T))
    
    return F, f, alpha_opt, sigma_meas_opt, sigma_control_opt
end

function getResultsALAP(model,C_bar,S)
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


function simulateTrajectories(A,B,C,k,F,f;numberSamples=1)
    # Works only if all polyhedra are hypercubes
    n_y,n_x=size(C)
    T=size(f)[2]
    
    x=zeros(numberSamples,n_x,T+1)
    for i=1:numberSamples
        y=zeros(n_y,T)
        x[i,:,0+1]=rand(Float64,n_x)*2 .-1 # x_0
        for t=0:T-1
            # compute y
            v=rand(Float64,n_y)*2 .-1
            y[:,t+1]=C*x[i,:,t+1]+v

            # compute u
            u=f[:,t+1]
            for tau=0:t
                u+=F[:,:,t+1,tau+1]*y[:,tau+1]
            end

            # compute x
            w=rand(Float64,n_x)*2 .-1
            x[i,:,t+2]=A*x[i,:,t+1]+B*u+k+w
        end
    end
    return x
end











