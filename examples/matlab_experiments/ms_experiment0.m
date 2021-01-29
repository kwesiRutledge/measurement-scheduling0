function [results] = ms_experiment0(varargin)
	%ms_experiment0.m
	%Description:
	%	Tests the function L_AQ().

	%% Constants %%

	dim = 2;

	A1 = diag([0.2;0.5]);
	A2 = eye(dim);

	Q = 0.3*eye(dim);
	P = 1.7*eye(dim);

	results.constants.A1 = A1;
	results.constants.A2 = A2;
	results.constants.Q = Q;
	results.constants.P = P;

	experiment_name = 'ms_experiment0';

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	disp(['Beginning ' experiment_name '.'])
	disp(' ')
	disp('L_{A1,Q}(P) = ')
	L_AQ1 = L_AQ(A1,Q,P)

	disp('L_{A2,Q}(P) = ')
	L_AQ2 = L_AQ(A2,Q,P)

	results.L_AQ1 = L_AQ1;
	results.L_AQ2 = L_AQ2;

	%% Determine How Each Compares to their Predecessor
	disp('Is P \succeq L_{A1,Q}(P)?')
	if all(eig( L_AQ1 - P ) >= 0)
		disp('- Yes!')
	else
		disp('- No!')
	end
	disp(' ')

	disp('Is P \succeq L_{A2,Q}(P)?')
	if all(eig( L_AQ2 - P ) >= 0)
		disp('- Yes!')
	else
		disp('- No!')
	end

	%% Use Optimization to Find Out When/If We Can Reconstruct A Previous Covariance

	P_bar = diag([1.7;3.7]);
	P_prev = sdpvar(dim,dim,'symmetric');

	objective = [ trace(P_prev) ]

	update_constraint = [ A1*P_prev*A1' + Q == P_bar ];

	ops0 = sdpsettings('verbose',1);

	diagnostics = optimize(update_constraint , objective , ops0);

	value(P_prev)

	results.experim2.diagnostics = diagnostics;
	results.experim2.P_prev = value(P_prev);

end