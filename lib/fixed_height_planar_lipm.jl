"""
fixed_height_planar_lipm.jl
Description:
    This file contains the definition of our simple linear system for the
    Fixed Height Planar Linear Inverted Pendulum from `Balancing and Step Recovery Capturability
    via Sums-of-Squares Optimization.” Robotics: Science and Systems XIII, July, 2017, Cambridge, 
    Massachusetts, Robotics: Science and Systems Foundation, 2017. Version: Author's final manuscript`
    [ https://hdl.handle.net/1721.1/124432 ].
"""

using HCubature

"""
get_lipm
Description:
    This function retrieves a standard Fixed Height LIPM defined with the parameters discussed in the
    paper referenced above.
Usage:
    (continuous_system, discrete_system) = get_lipm()
"""
function get_lipm(T)
    # Constants
    g = 9.81; #Gravitational Constant m/s^2
    
    zbar_cm = 1;      #1 meter
    r_foot = 0.5;    #0.05 meters

    dt = 0.1

    # Create Continuous Time System
    A = [   0.0  1.0 ;
            (g/(zbar_cm))    0 ]
    B = [ 0.0 ; (g/(zbar_cm))*r_foot ]
    C = [1.0 0.0; 0.0 1.0]
    D = [1.0 0.0; 0.0 1.0]
    
    n_x=size(A,1)
    n_u=size(B,2)
    n_y=size(C,1)
    n_z=size(D,1)
    
    k = zeros(n_x)
    d = zeros(n_z)
    
    
    eta_w = 0.05 #r_foot/5.0
    W = eta_w*hypercube(n_x)
    
    #W = MixedMatHRep(HalfSpace([1.0, 0.0], 0.0) ∩ HalfSpace([0.0, 1.0], eta_w) ∩ HalfSpace([-1.0, 0.0], 0.0) ∩ HalfSpace([0.0, -1.0], eta_w) )
    
    
    eta_v = 0.01
    V = eta_v*hypercube(n_y)
    
    #H_v = [ 1.0 0 ; 0 1.0 ; -1.0 0 ; 0.0 -1.0 ]
    #h_v = eta_v * [1.0;1.0;1.0;1.0]
    #V = hrep(H_v,h_v)
    
    
    eta_x0 = 0.1
    X0 = eta_x0*hypercube(n_x)
    
    
    #H_u = [-1.0 ; 1.0]
    #H_u=reshape(H_u, length(H_u), 1) # convert array into matrix
    #h_u = [0.0 ; r_foot]
    #U = hrep(H_u,h_u)
    U = hypercube(n_u) # [-1,1]
    
    
    eta_x_cm = 0.75
    eta_x_cm_dot = 5
    S = (eta_x_cm*hypercube(1)) * (eta_x_cm_dot*hypercube(1))
    
    list_S=[]
    for t=0:T
        push!(list_S,S)
    end
    
    # continuous time system
    c_system = LinearSystem(A,B,C,D,k,d,W,V,X0,U,list_S)

    # Create Discrete Time System

    Ad = exp(A*dt)
    Bd = zeros((2,1))
    est_err = 0.0

    let f = x -> transpose([1; 0.0])*exp(A.*x)*B
        (Bd[1],est_err) = hquadrature( f , 0.0 , dt )
    end
    print("Estimated error in calculating first entry of discretized B: ")
    println(est_err)

    let f = x-> transpose([0.0; 1])*exp(A.*x)*B
        (Bd[2],est_err) = hquadrature(f, 0.0, dt)
    end
    print("Estimated error in calculating second entry of discretized B: ")
    println(est_err)

    d_system = LinearSystem(Ad,Bd,C,D,k,d,W,V,X0,U,list_S)

    return  d_system, c_system
end