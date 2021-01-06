function [results] = ms_experiment7(varargin)
	%ms_experiment7.m
	%Description:
	%	Create a dynamics object and tries to find the best schedule
	%	of all schedules with a given length.

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

	dyn0 = DynamicalSystem( A1 , b , C , Q1 , R1 , x0_mean , P );	
	results.constants.dyn0 = dyn0;

	eta_x0 = 1.7; eta_w = 1.5; eta_v = 0.3;
	P_x0 = eta_x0 * Polyhedron('lb',-ones(1,dyn0.n_x),'ub',ones(1,dyn0.n_x));
	P_w = eta_w * Polyhedron('lb',-ones(1,dyn0.n_x),'ub',ones(1,dyn0.n_x));
	P_v = eta_v * Polyhedron('lb',-ones(1,dyn0.n_y),'ub',ones(1,dyn0.n_y));

	results.constants.P_x0 = P_x0;
	results.constants.P_w = P_w;
	results.constants.P_v = P_v;

	T0 = 10;
	card_M = 3;

	experiment_name = 'ms_experiment7';

	plot_conv_sets = false;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%% Test synthesize_feedback_and_schedule

	script_start = tic;

	disp(['Beginning ' experiment_name '.'])
	disp(' ')

	disp('1. Call MILP Function.')

	schedule_synthesis_start = tic;

	[ ~ , ~ , opt_schedule , opt_cost ] = dyn0.synthesize_feedback_and_schedule( T0 , card_M , P_x0 , P_w , P_v , 1 );
		
	schedule_synthesis_time = toc(schedule_synthesis_start);

	%Print results
	disp(['  + Cost: ' num2str(opt_cost) ])
	disp(['  + Time: ' num2str(schedule_synthesis_time) ])

	results.Test1.Cost = opt_cost;
	results.Test1.OptimalSchedule = opt_schedule;

end