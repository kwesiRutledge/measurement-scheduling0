using LinearAlgebra

struct KalmanFilter
    transition_matrix        # A
    transition_covariance    # Q
    observation_matrix       # C
    observation_covariance   # R
    initial_state_covariance # P0
end



function loadKF()
    # Define parameters of the problem
    transition_matrix=[0.2 0; 0 1.3]
    transition_covariance=[1 0; 0 1]
    observation_matrix=[0.5 1.5; 0 1]
    observation_covariance=[1 0; 0 1]
    initial_state_covariance=[1 0; 0 1]
    
    return KalmanFilter(transition_matrix, transition_covariance, observation_matrix, observation_covariance, initial_state_covariance)
end

function timeUpdate(transition_matrix,transition_covariance,state_covariance)
    return transition_matrix*state_covariance*transition_matrix'+transition_covariance
end

function observationUpdate(observation_matrix,observation_covariance,state_covariance)
    return inv( inv(state_covariance) + observation_matrix'*inv(observation_covariance)*observation_matrix )
end

function scheduleCost(kalman_filter,sched,T)
    costs=zeros(T)
    P=kalman_filter.initial_state_covariance
    for t=1:T
        P=timeUpdate(kalman_filter.transition_matrix,kalman_filter.transition_covariance,P)
        
        if t in sched
            P=observationUpdate(kf.observation_matrix,kf.observation_covariance,P)
        end
        costs[t] = tr(P)
    end
    cost=sum(costs)
    return cost, costs
end