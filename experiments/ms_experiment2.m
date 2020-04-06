function [results] = ms_experiment2(varargin)
	%ms_experiment2.m
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

	experiment_name = 'ms_experiment2';

	plot_conv_sets = false;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	disp(['Beginning ' experiment_name '.'])
	disp(' ')
	
	dyn0 = DynamicalSystem( A1 , b , C , Q1 , R1 , x0_mean , P );
	
	ifk0 = IntermittentKalmanFilter( dyn0 );

	rand_schedule0 = [ 0, 3, 7];
	rand_schedule1 = [ 2, 4, 6];

	P_traj = ifk0.simulate_schedule(rand_schedule0,'TimeHorizon',10);
	results.exp1.schedule = rand_schedule0;
	results.exp1.P_traj = P_traj;

	P_traj = ifk0.simulate_schedule(rand_schedule1,'TimeHorizon',10);
	results.exp2.schedule = rand_schedule1;
	results.exp2.P_traj = P_traj;


	results.dyn0 = dyn0;
	results.ifk0 = ifk0;


end