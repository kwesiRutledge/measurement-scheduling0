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
function get_lipm()
    # Constants
    g = 10; #Gravitational Constant
    
    zbar_cm = 1;      #1 meter
    r_foot = 0.05;    #0.05 meters

    dt = 0.1

    # Create Continuous Time System
    A = [   0  1 ;
            (g/(zbar_cm))    0 ]
    B = [ 0 ; (g/(zbar_cm))*r_foot ]
    C = [1.0 0; 0.0 1.0]

    eta_w = r_foot/5.0
    W = MixedMatHRep(HalfSpace([1.0, 0.0], 0.0) ∩ HalfSpace([0.0, 1.0], eta_w) ∩ HalfSpace([-1.0, 0.0], 0.0) ∩ HalfSpace([0.0, -1.0], eta_w) )
    eta_v = 0.1
    H_v = [ 1.0 0 ; 0 1.0 ; -1.0 0 ; 0.0 -1.0 ]
    h_v = eta_v * [1.0;1.0;1.0;1.0]
    V = hrep(H_v,h_v)

    c_system = LinearSystem(A,B,C,W,V)

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

    d_system = LinearSystem(Ad,Bd,C,W,V)

    return (c_system, d_system)


end