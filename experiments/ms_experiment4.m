function [results] = ms_experiment4(varargin)
	%ms_experiment4.m
	%Description:
	%	Create a dynamics object and solving the problem from Section 4 of the Bounded_Error document.

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
	results.dyn0 = dyn0;

	T0 = 10;
	M = [2,5,8];

	[J,S,C_bar] = dyn0.get_mpc_matrices(T0,M);

	results.exp1.J = J;
	results.exp1.S = S;
	results.exp1.C_bar = C_bar;

	disp('1. Obtained MPC matrices.')

	%% Define Optimization %%
	disp('2. Starting to construct optimization.')

	eta_x0 = 1.7; eta_w = 1.5; eta_v = 0.3;
	P_x0 = eta_x0 * Polyhedron('lb',-ones(1,dyn0.n_x),'ub',ones(1,dyn0.n_x));
	P_w = eta_w * Polyhedron('lb',-ones(1,dyn0.n_x),'ub',ones(1,dyn0.n_x));
	P_v = eta_v * Polyhedron('lb',-ones(1,dyn0.n_y),'ub',ones(1,dyn0.n_y));

	disp('- Created the disturbance sets.')

	results.exp1.P_x0 = P_x0;
	results.exp1.P_w = P_w;
	results.exp1.P_v = P_v;

	%% Create Optimization Variables %%

	T = 10;
	[J,S,C_bar] = dyn0.get_mpc_matrices(T,M);

	Q = sdpvar(dyn0.n_x*T,dyn0.n_x*T,'full');
	r = sdpvar(dyn0.n_x*T,1);

	P_wT = 1; P_vT = 1;
	for t = 1:T
		%Construct the sets by products
		P_wT = P_wT * P_w;
		P_vT = P_vT * P_v;
	end

	P_eta = P_wT * P_vT * P_x0;

	Lambda = sdpvar( 2*dyn0.n_x*(T+1) , size(P_eta.A,1) , 'full' );

	%Calculate Big C Matrix
	b = binvar(T,1);
	% C_at_each_n = {};
	% for t = 0:T-1
	% 	C_at_each_n{t+1} = dyn0.C*b(t+1); 
	% end

	% C_bar = [ blkdiag(C_at_each_n{:}) zeros( [dyn0.n_y , dyn0.n_x] * [ T 0 ; 0 1 ] ) ];

	%Alpha-bar
	alpha0 = sdpvar(T+1,1,'full');
	alpha_bar = kron(alpha0,ones(dyn0.n_x,1));

	alpha_bar = [ alpha_bar ; alpha_bar ];

	disp('- Created optimization variables.')

	%% Write Optimization Constraints %%

	cg = ConstraintGenerator();

	P_xi_w = S + S*Q*C_bar*S;
	P_xi_v = S*Q;
	P_xi_xi0 = ( eye(size(S,1)) + S*Q*C_bar )*J;

	R_T_mT = [ eye(dyn0.n_x*(T+1)) ; -eye(dyn0.n_x*(T+1)) ];

	dual_constraints = 	[ Lambda * P_eta.A == R_T_mT * [ P_xi_w , P_xi_v , P_xi_xi0 ] ] + ...
						[ Lambda * P_eta.b <= alpha_bar - R_T_mT * S*r ];

	low_diag_constr = cg.get_block_lower_diagonal_constraint_on(Q,[dyn0.n_x,dyn0.n_y]);

	positivity_constr = [ Lambda >= 0 , alpha0 >= 0 ];

	k0 = 2;
	schedule_constr = [ sum( b ) <= k0 ];

	%% Write Optimization Problem %%

	ops0 = sdpsettings('verbose',1);

	diagnostics = optimize(	dual_constraints+low_diag_constr+positivity_constr+schedule_constr , ...
							sum(alpha0), ...
							ops0)

	results.exp1.Q = value(Q);
	results.exp1.r = value(r);
	results.exp1.F = inv(eye( size(Q ,1)) + value(Q)*C_bar*S ) * value(Q);
	results.exp1.f = inv(eye( size(Q ,1)) + value(Q)*C_bar*S ) * value(r);
	results.exp1.alpha0 = value(alpha0);
	results.exp1.objective = sum(value(alpha0));

	[ F2 , f2 , cost2 ] = synthesize_fhae( dyn0 , T , M , P_x0 , P_w , P_v );

	results.exp2.F = F2;
	results.exp2.f = f2;
	results.exp2.objective = cost2;

end