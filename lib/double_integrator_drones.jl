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
function get_drones( num_drones, T ) # T is the time horizon
    
    if mod(T,2)!=0
        println("ERROR, T must be even.")
        return
    end
    
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
        C[(drone_index-1)*2+1,(drone_index-1)*4+1] = 1
        C[(drone_index-1)*2+2,(drone_index-1)*4+2] = 1
    end
    
    # Create D Matric
    D = zeros((num_drones*2,num_states))
    for drone_index = 1:num_drones
        D[(drone_index-1)*2+1,(drone_index-1)*4+1] = 1
        D[(drone_index-1)*2+2,(drone_index-1)*4+2] = 1
    end
    
    k=zeros(num_states)
    d=zeros(size(D,1))
    
    # scaling w
    eta_w = 0.05
    scalingW = [eta_w; eta_w; eta_w; eta_w]
    polW = hyperbox(scalingW)
    polW = cartesianPower(polW,num_drones)
    
    # scaling v
    eta_v = 0.05
    scalingV = [eta_v; eta_v]
    polV = hyperbox(scalingV)
    polV = cartesianPower(polV,num_drones)
    
    # scaling x0
    eta_x0 = 1.0
    scalingX0 = [eta_x0; eta_x0; 0; 0]
    polX0 = hyperbox(scalingX0)
    polX0 = cartesianPower(polX0,num_drones)
    
    # Creating the set constraining inputs: U
    eta_u = 20.0
    polU = hRepRectangle(-eta_u,-eta_u,eta_u,eta_u)
    polU = cartesianPower(polU,num_drones)
    
    # creating a list of Safety set: S_t = list_S[t+1] with the constraint that z_t mus be in S_t.
    horizontalRectangle=hRepRectangle(-1.5,-1.5,6.5,1.5)
    
    list_S=[]
    for t=0:T/2-1
        push!(list_S, cartesianPower(horizontalRectangle,num_drones) )
    end
    
    verticalRectangle=hRepRectangle(3.5,-1.5,6.5,6.5)
    
    for t=T/2:T-1
        push!(list_S, cartesianPower(verticalRectangle,num_drones) )
    end
    
    finalRectangle=hRepRectangle(4.5,4.5,5.5,5.5) # 4.0,4.0,6.0,6.0
    push!(list_S, cartesianPower(finalRectangle,num_drones) )
    
    # Create Continuous Time Linear System
    d_system = LinearSystem(Ad,Bd,C,D,k,d,polW,polV,polX0,polU,list_S)
    
    polZ0=hRepRectangle(-eta_x0,-eta_x0,eta_x0,eta_x0) # used for plotting
    
    return d_system, polZ0
end

function hRepRectangle(x_min,y_min,x_max,y_max)
    H=[-1 0; 1 0; 0 -1; 0 1]
    h=[-x_min; x_max; -y_min; y_max]

    return hrep(H,h)
end

function cartesianPower(pol,n)
    # return the n-th cartesian power of the polyhedron "pol". Return pol if n<2
    power=pol
    for i=1:n-1
       power*=pol 
    end
    return power
end