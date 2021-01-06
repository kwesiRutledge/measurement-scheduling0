function [results] = ms_experiment6(varargin)
	%ms_experiment6.m
	%Description:
	%	Demonstrates the MILP solver available in Gurobi.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	Px = Polyhedron('lb',-ones(1,dim),'ub',ones(1,dim));
	Py = Polyhedron('lb',-ones(1,dim),'ub',[1,2]);

	results.Parameters.dim = dim;
	results.Parameters.Px = Px;

	experiment_name = 'ms_experiment6';

	plot_conv_sets = false;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	disp(['Beginning ' experiment_name '.'])
	disp(' ')
	

	disp('This is a sanity check to verify that Gurobi understands how to solve MILP problems.')
	disp(' ')

	%% Create Simple MILP
	
	%Create Optimization Variables
	b = binvar(2,1);
	x = sdpvar(2,1,'full');
	y = sdpvar(2,1,'full');

	%Create constraints

	c = [0;1];

	objective = [ c'*(y+x) ];

	constraints = [ Px.A * x <= Px.b*b(1) ] + [ Py.A * y <= Py.b*b(2) ] + [ sum(b) == 1 ];

	%Solve Problem

	disp('The solution to this optimization should be 2.')
	ops0 = sdpsettings('verbose',1);

	diagnostics = optimize( constraints , -objective,  ops0)

	results.Experiment1.diagnostics = diagnostics;
	results.Experiment1.OptimalObjective = value(objective);

	%% Test synthesize_feedback_and_schedule

	disp(' ')
	disp('Experiment2')

	

end