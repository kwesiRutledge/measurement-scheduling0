"""
gouble_integrator_drones.jl
Description:
    This file contains the definition of our simple linear system for a fleet of 
    drones where each drone is represented by a double integrator.
"""

"""
get_drones
Description:
    Creates a system representing the number of drones that the user desires.
"""
function get_drones( num_drones )
    # Constants
    dt = 0.1
    
    # Create A Matrix
    num_states = num_drones*4
    Ad = zeros((num_states,num_states))
    Ad = Ad + I
    # There is likely a faster way to do this. I don't see it yet...
    for drone_index = 1:num_drones
        Ad[(drone_index-1)*4+1,(drone_index-1)*4+3] = dt
        Ad[(drone_index-1)*4+2,(drone_index-1)*4+4] = dt
    end
    
    # Create B Matrix
    Bd = zeros((num_states,num_drones*2))
    for drone_index = 1:num_drones
        Bd[(drone_index-1)*4+3,(drone_index-1)*2+1] = dt
        Bd[(drone_index-1)*4+4,(drone_index-1)*2+2] = dt
    end

    # Create C Matrix
    C = zeros((num_drones*2,num_states))
    for drone_index = 1:num_drones
        C[(drone_index-1)*2+1,(drone_index-1)*4+3] = 1
        C[(drone_index-1)*2+2,(drone_index-1)*4+4] = 1
    end

    # Create the Performance Matrices D_x, D_u
    D = zeros((num_drones*2,num_states))
    for drone_index = 1:num_drones
        D[(drone_index-1)*2+1,(drone_index-1)*2+1] = 1
        D[(drone_index-1)*2+2,(drone_index-1)*2+2] = 1
    end
    
    d=zeros(size(D,1))

    # Creating Disturbance Set

    eta_w = 0.5
    H_w = [ 1.0 0 ; 0 1.0 ; -1.0 0 ; 0.0 -1.0 ]
    h_w = eta_w * [1.0;1.0;1.0;1.0]
    W = hrep(H_w,h_w)

    eta_v = 0.2
    H_v = [ 1.0 0 ; 0 1.0 ; -1.0 0 ; 0.0 -1.0 ]
    h_v = eta_v * [1.0;1.0;1.0;1.0]
    V = hrep(H_v,h_v)

    # Create Input Set U
    eta_u = 50
    U = hRepRectangle(-eta_u,-eta_u,eta_u,eta_u)

    # Create Continuous Time Linear System

    d_system = LinearSystem(Ad,Bd,C,D,d,W,V,U)

    return d_system

end