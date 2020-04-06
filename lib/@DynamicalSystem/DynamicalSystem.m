classdef DynamicalSystem
	%ds = DynamicalSystem( A , b , C, Q, R , x0_bar , P0_bar )
	%
	%
	%Description:
	%	A dynamical system as defined by Antoine Aspeel and Kwesi Rutledge in the measurement scheduling
	%	problem. See document for details.

	properties
		A;
		b;
		C;
		Q;
		R;
		x0_bar;
		P0_bar;

		n_x; %Dimension of the state
		n_y; %Dimension of the output

		%State
		x;
	end

	methods
		%%%%%%%%%%%%%%%%%
		%% Constructor %%
		%%%%%%%%%%%%%%%%%
		function ds = DynamicalSystem( A , b , C, Q, R , x0_bar , P0_bar )
			%Description:
			%	Initializes a linear dynamical system with no input, an offset term, and a noise term

			allowed_nargins = [ 7 ];

			%% Input processing %%
			if ~any(nargin == allowed_nargins)
				error([ num2str(nargin) ' number of inputs is not supported.']);
            end
            
            ds.n_x = size(A,1);
            ds.n_y = size(C,1);

            if any( size(A) ~= [ds.n_x,ds.n_x] )
            	error('A matrix must be a square matrix.')
            end

            if any( size(b) ~= [ds.n_x,1] )
            	error('b matrix must be a vector that has n rows and 1 column.')
            end

            if any( size(C) ~= [ds.n_y,ds.n_x] )
            	error('C matrix must be a vector with n_y rows and n_x columns.')
            end

            if any( size(Q) ~= [ds.n_x,ds.n_x] )
            	error('Q must be a square matrix with n_x rows and columns.')
            end

            if any( size(R) ~= [ds.n_y,ds.n_y] )
            	error('R must be a square matrix with n_y rows and columns.')
            end

            if any( size(x0_bar) ~= [ds.n_x,1] )
            	error('x0_bar must be a vector with n_x elements.')
            end

            if any( size(P0_bar) ~= [ds.n_x,ds.n_x] )
            	error('P0_bar must be a square matrix with n_x rows and columns.')
            end

            %% Algorithm %%

            ds.A = A;
            ds.b = b;
            ds.C = C;
            ds.Q = Q;
            ds.R = R;
            ds.x0_bar = x0_bar;
            ds.P0_bar = P0_bar;

            %Initialize State
            ds.x = mvnrnd(x0_bar,P0_bar)';

		end

		function [ y ] = get_output( obj )
			%Description:
			%	Creates an output of the system's current state using the equation
			%		y = Cx + v, where v ~ N(0,R)

			%% Algorithm %%
			y = obj.C * obj.x + mvnrnd(zeros(obj.n_y,1),obj.R)';

		end

		function [ x_tp1 ] = update_state( obj )
			%Description:
			%	Creates the state of the system at the next time instant according to the equation
			%		x^+ = Ax + b + w, where w ~ N(0,Q)
			%	Also updates the value of the current state.

			%% Algorithm %%
			x_tp1 = obj.A * obj.x + obj.b + mvnrnd( zeros(obj.n_x,1) , obj.Q );

			obj.x = x_tp1;

		end
	end

end