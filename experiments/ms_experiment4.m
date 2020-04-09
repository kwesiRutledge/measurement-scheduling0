function [results] = ms_experiment4(varargin)
	%ms_experiment4.m
	%Description:
	%	Create a dynamics object and play around with it.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	A1 = diag([0.2;1.3]);
	b = rand([2,1]);
	C = [ 0.5, 1.5; 0 , 1 ];

	P = 1.7*eye(dim);
	Q1 = diag([1.5;0.7]);
	R1 = diag([0.1;0.2]);

	x0_mean = zeros(dim,1);

	results.constants.A1 = A1;
	results.constants.b = b;
	results.constants.C = C;
	results.constants.P = P;
	results.constants.Q1 = Q1;
	results.constants.R1 = R1;
	results.constants.x0_mean = x0_mean;

	experiment_name = 'ms_experiment4';

	plot_conv_sets = false;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	disp(['Beginning ' experiment_name '.'])
	disp(' ')
	
	dyn0 = DynamicalSystem( A1 , b , C , Q1 , R1 , x0_mean , P );	
	
	eta_x0 = 1.7; eta_w = 1.5; eta_v = 0.3;
	P_x0 = eta_x0 * Polyhedron('lb',-ones(1,dyn0.n_x),'ub',ones(1,dyn0.n_x));
	P_w = eta_w * Polyhedron('lb',-ones(1,dyn0.n_x),'ub',ones(1,dyn0.n_x));
	P_v = eta_v * Polyhedron('lb',-ones(1,dyn0.n_y),'ub',ones(1,dyn0.n_y));

	results.exp1.P_x0 = P_x0;
	results.exp1.P_w = P_w;
	results.exp1.P_v = P_v;

	T0 = 10;
	M = [2,5,8];

	[J,S,C_bar] = dyn0.get_mpc_matrices(T0,M);

	results.exp1.J = J;
	results.exp1.S = S;
	results.exp1.C_bar = C_bar;

	results.dyn0 = dyn0;

end