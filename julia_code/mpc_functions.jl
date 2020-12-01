"""
    mpc_functions.jl
    Description:
        File contains some of the commonly used shortcuts for creating the large block matrices for
        MPC instances.
"""

using LinearAlgebra
using Test

"""
    compute_Sx0
    Description:
        Construct the matrix which describes how the initial condition (x0) propagates 
        to all of the states of a time horizon T trajectory problem.
    Inputs:
        A: A square matrix of some dimension.
    Assumptions:
        This assumes that the initial time in the trajectory is t=0 and the input T 
        is a scalar defining the final time in the time horizon.
"""
function compute_Sx0( A::Matrix{<:Number} , T )
    J = I(size(A,1))
    for i = 1:T
        J = [J; (A^i)]
    end
    return J
end

function compute_Sx0( A::Number , T::Int )
    J = 1.0
    for i = 1:T
        J = [J; (A^i)]
    end
    return J
end

"""
    compute_Sw
    Description:
        Construct the matrix which describes how the disturbances to a linear system (w)
        propagates to all of the states of a time horizon T trajectory problem.
    Inputs:
        A: A square matrix of some dimension.
    Assumptions:
        This assumes that the initial time in the trajectory is t=0 and the input T 
        is a scalar defining the final time in the time horizon.
"""
function compute_Sw( A::Matrix{<:Number} , TimeHorizon )
    #Compute dimension constants
    n_x,n_test=size(A)
    if n_x ≠ n_test
        throw(ArgumentError("Expected argument A to be a square matrix."))
    end

    #Algorithm
    S = zeros(n_x*(TimeHorizon+1),n_x*TimeHorizon)
    current_row = I( size(A,1) )
    S[n_x+1:2*n_x,1:n_x] = current_row
    for i = 1:TimeHorizon-1
        current_row = [A*current_row[1:n_x,1:n_x] current_row]
        S[(i+1)*n_x+1:(i+2)*n_x, 1:(i+1)*n_x] = current_row
    end
    return S
end

function compute_Sw( A::Number , TimeHorizon )
    #This should be a scalar/

    # Constants
    floatA = A + 0.0

    #Algorithm
    S = zeros(TimeHorizon+1,TimeHorizon)
    S[2,1] = 1
    current_row = [1]
    for i = 1:TimeHorizon-1
        current_row = [floatA*current_row[1] current_row]
        S[(i+1)+1:(i+2), 1:(i+1)] = current_row
    end
    # println(S)
    return S
end

"""
    compute_C_M(C,M,T)
    Description:
        Constructs the measurement matrix C which is dependent on the provided schedule (M)
        as well as the total time horizon (T).
    Inputs:

    Assumptions:
"""
function compute_C_M( C , M , T)
    #Compute Dimension constants
    n_y,n_x=size(C)

    # construct C_sched (=C_bar)
    M=M.+1 # to have indices
    diag=zeros(T)
    diag[M].=1
    dMat=Diagonal(diag)
    C_sched=[kron(dMat,C) zeros(n_y*T,n_x)]

    #Return C matrix
    return C_sched
end