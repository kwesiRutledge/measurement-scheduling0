classdef IntermittentKalmanFilter
	%
	%
	%
	%Description:
	%	A dynamical system as defined by Antoine Aspeel and Kwesi Rutledge in the measurement scheduling
	%	problem. See document for details.

	properties
		DynamicalSystem;

		t=0; %Current time

		%Current Estimates
		xh_t_t; %x-hat(t|t)
		P_t_t; %Covariance matrix P(t|t)

		%Estimate Ingredients
		xh_t_tm1;
		P_t_tm1;
	end

	methods
		%%%%%%%%%%%%%%%%%
		%% Constructor %%
		%%%%%%%%%%%%%%%%%
		function ikf = IntermittentKalmanFilter( DynamicalSystem_in )
			%Description:
			%	Initializes a linear dynamical system with no input, an offset term, and a noise term

			allowed_nargins = [ 1 ];

			%% Input processing %%
			if ~any(nargin == allowed_nargins)
				error([ num2str(nargin) ' number of inputs is not supported.']);
            end

            %% Algorithm %%

            ikf.DynamicalSystem = DynamicalSystem_in;

            %Initialize Estimate Variables
            ikf.xh_t_tm1 = DynamicalSystem_in.x0_bar;
            ikf.P_t_tm1 = DynamicalSystem_in.P0_bar;

            ikf.xh_t_t = NaN( DynamicalSystem_in.n_x, 1 );
            ikf.P_t_t = NaN( DynamicalSystem_in.n_x, DynamicalSystem_in.n_x );

		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Simulation Functions %%
		%%%%%%%%%%%%%%%%%%%%%%%%%%

		function [ P_history ] = simulate_schedule( varargin )
			%Description:
			%	Simulates the filter and estimator working together over a given
			%	time horizon (T) or until after the schedule has elapsed and a
			%	covariance bound (P_bound) is violated.
			%	
			%Usage:
			%	[ P_history ] = ikf.simulate_schedule( M , 'TimeHorizon' , T )
			%	[ P_history ] = ikf.simulate_schedule( M , 'CovarianceBound' , P_bound )

			%% Input Processing %%

			obj = varargin{1};
			M = varargin{2};

			switch varargin{3}
				case 'TimeHorizon'
					T = varargin{4};
				case 'CovarianceBound'
					P_bound = varargin{4};
				otherwise
					error(['Unrecognized option: ' varargin{3} ])
			end

			%% Algorithm %%

			switch varargin{3}
				case 'TimeHorizon'
					%Do step 0
					[ xhat_t_t , P_t_t ] = obj.correction_update( M );
					obj.xh_t_t = xhat_t_t;
					obj.P_t_t = P_t_t;

					P_history = [P_t_t];

					%Update State
					obj.DynamicalSystem.update_state();

					%Iterate through every time step.
					for tau = 1:T
						%Update internal filter time.
						obj.t = obj.t+1;

						%Perform Correction Step
						[ xh_t_tm1 , P_t_tm1 ] = obj.time_update();
						obj.xh_t_tm1 = xh_t_tm1;
						obj.P_t_tm1 = P_t_tm1;

						%Perform Time Update.
						[ xhat_t_t , P_t_t ] = obj.correction_update(M);
						obj.xh_t_t = xhat_t_t;
						obj.P_t_t = P_t_t;

						%Update State
						obj.DynamicalSystem.update_state();

						%Save P matrix
						P_history(:,:,end+1) = P_t_t; 
					end

				case 'CovarianceBound'
					error('This part is not done yet.')
			end


		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Kalman Estimator Update Functions %%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		function [ xh_t_t , P_t_t ] = correction_update( obj , M )
			%Description:
			%	Computes the correction step of the intermittent kalman filter when the 
			%	measurement schedule M is given.

			%% Constants %%

			C = obj.DynamicalSystem.C;
			R = obj.DynamicalSystem.R;

			n_x = obj.DynamicalSystem.n_x;

			%% Algorithm %%

			if any( obj.t == M )
				%If the current time is in the measurement schedule.
				K_t = obj.P_t_tm1 * C'* ( C*obj.P_t_tm1*C' + R )^(-1);
			else
				%If the current time is not in the measurement schedule.
				K_t = zeros( n_x );
			end

			%Final Results
			z_t = obj.DynamicalSystem.get_output();

			xh_t_t = obj.xh_t_tm1 + K_t * ( z_t - C * obj.xh_t_tm1 );
			P_t_t = (eye(n_x) - K_t * C )*obj.P_t_tm1;

		end


		function [ xh_tp1_t , P_tp1_t ] = time_update( obj )
			%Description:
			%	Computes the time update from the intermittent kalman filter estimator.

			%% Constants %%

			A = obj.DynamicalSystem.A;
			b = obj.DynamicalSystem.b;
			Q = obj.DynamicalSystem.Q;

			%% Algorithm %%
			xh_tp1_t = A * obj.xh_t_t + b;
			P_tp1_t = A * obj.P_t_t * A' + Q;

		end

	end

end