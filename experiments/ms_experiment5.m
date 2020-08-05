function [results] = ms_experiment5(varargin)
	%ms_experiment5.m
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

	experiment_name = 'ms_experiment5';

	plot_conv_sets = false;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	script_start = tic;

	disp(['Beginning ' experiment_name '.'])
	disp(' ')

	disp('1. Creating All Possible Schedules.')
	time_instants = [0:T0-1];
	all_schedules = nchoosek( time_instants , card_M );

	disp(['- There are ' num2str(size(all_schedules,1)) ' schedules in the list.' ])

	disp('2. Identifying the costs for each schedule.')
	schedule_costs = NaN( size(all_schedules,1) , 1 );
	schedule_synthesis_times = zeros( size(all_schedules,1) , 1 );

	for schedule_idx = 1:size(all_schedules,1)
		disp(['- Schedule #' num2str(schedule_idx) ])
		schedule_synthesis_start = tic;

		temp_schedule = all_schedules(schedule_idx,:); %Grab a schedule.

		[ ~ , ~ , schedule_costs(schedule_idx) ] = synthesize_fhae( dyn0 , T0 , temp_schedule , P_x0 , P_w , P_v , 0 );
		
		schedule_synthesis_times( schedule_idx ) = toc(schedule_synthesis_start);

		%Print results
		disp(['  + Cost: ' num2str(schedule_costs(schedule_idx)) ])
		disp(['  + Time: ' num2str(schedule_synthesis_times(schedule_idx)) ])

	end

	disp('3. Determining the Optimal Schedule Based on Exhaustive Search')

	[~,opt_schedule_idcs] = min( schedule_costs );
	disp('- Optimal schedule is:')
	all_schedules( opt_schedule_idcs(1) , : )
	disp(['- Optimal cost = ' num2str( schedule_costs( opt_schedule_idcs(1) ) ) ])


	script_time = toc(script_start);

	results.exp2.schedule_costs = schedule_costs;
	results.exp2.schedule_synthesis_times = schedule_synthesis_times;
	results.script_time = script_time;

end