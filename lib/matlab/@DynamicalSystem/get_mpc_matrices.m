function [varargout] = get_mpc_matrices(varargin)
	% get_mpc_matrices
	%	Description:
	%		Creates MPC matrices that are needed for the dynamical system in our work.
	%		This is currently designed to handle:
	%		- A missing observation system being part of the switched system
	%		- A system with unique disturbance matrices as being part of the switched system
	%
	%	Usage:
	%		[J,S,C_bar] = get_mpc_matrices(ds,T,M)
	%		[J,S,C_bar] = ds.get_mpc_matrices(T,M)
	%
	%	Inputs:
	%		sys_arr - 	An array of structs, each containing system matrices and an initial 
	%					condition for the desired system.
	%					Required fields: .A,.B,.C,.x0
	%
	%		T - 		Time horizon for the MPC matrices.
	%					When T is given (i.e. a single integer is given as the second input),
	%					it is assumed that the switching sequence is: 
	%		
	%		sigma -		An array of integers.
    %                   This should be a single sequence of integer values
    %                   which defines the value of the discrete state at
    %                   each time state (and thus which dynamics object in
    %                   the sys_arr should be used).
    %
    %		lcsas -		A Language Constrained, Switched Affine System object representing the
    %					the switched system.
	%
	%	Outputs:
	%		H - 	This is the matrix which defines how the process disturbances (w(t)) affect the
	%				system's state trajectory.
	%
	%		S - 	This is the matrix which defines how the inputs ( u(t) ) affect the system's
	%				state trajectory.
	%
	%		C_bar - ... defines how the state of the system is transmitted to the measurement trajectory.
	%
	%		J - 	... defines how the initial state (x(t_0)) affects the system's state trajectory.
	%
	%	Assumptions:
	%		We assume that the piecewise affine systems that are given as input maintain constant dimension of
	%		their disturbance (w), input (u), etc. That means that if u(t) is a 2 dimensional input at time 1
	%		then it always is.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if (nargin ~= [2,3])
		error(['Inappropriate number of arguments. (Received ' num2str(nargin) ')'])
	end

	ds = varargin{1};
	T  = varargin{2};
	if nargin > 2
		M  = varargin{3};
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Defining the Matrices %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Find S, and J Matrices
	S = ds.calc_w_effect_mat(T);
	J = ds.calc_J_mat(T);
	
	varargout{1} = J;
	varargout{2} = S;


	%%%%%%%%%%%%%%%%%%%%%%
	%% Optional Outputs %%
	%%%%%%%%%%%%%%%%%%%%%%

	%Output Matrix	
	if nargin > 2
		%Calculate Big C Matrix
		C_at_each_n = {};
		for t = 0:T-1
			C_at_each_n{t+1} = ds.C*( any(t == M) ); 
		end

		C_bar = [ blkdiag(C_at_each_n{:}) zeros( [ds.n_y , ds.n_x] * [ T 0 ; 0 1 ] ) ];

		varargout{3} = C_bar;
	end

end