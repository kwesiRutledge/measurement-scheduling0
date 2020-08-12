function [ F , f , M_opt , opt_cost ] = synthesize_feedback_and_schedule( varargin )
%synthesize_feedback_and_schedule.m
%Description:
%	Synthesizes a Finite Horizon Affine Estimator (uses the structure of Skaf and Boyd's 2010 work), i.e.
%	for a system
%		x(t+1) = Ax(t) + b + w(t)
%		y(t)   = Cx(t) + v(t)
%	where w(t) \in P_w, v(t) \in P_v, and all other matrices are known. One can derive an estimator of the following form:
%		xh(t+1) = Axh(t) + b, 	if t \notin M
%				= Axh(t) + b - \sum_{tau=0}^{t} F_{t,\tau} (y(\tau) - yh(\tau)) + f_t
%		yh(t)   = Cxh(t)
%	Then the estimation error system is defined as:
%		e(t+1) = Ae(t), 											if t \notin M
%			   = Ae(t) + \sum_{tau=0}^{t} F_{t,\tau} e(\tau) + f_t, if t \in M
%
%Example Usage:
%	[F,f,M_out,cost_out] = synthesize_fhae( ds , T , MCardinality , P_x0 , P_w , P_v )
%
%Inputs:
%	

	%% Input Processing %%
	ds = varargin{1};
	T = varargin{2};
	MCardinality = varargin{3};

	P_x0 = varargin{4};
	P_w = varargin{5};
	P_v = varargin{6};

	if ~isa(P_x0,'Polyhedron')
		error('P_x0 must be a Polyhedron() object from MPT3.')
	end

	if ~isa(P_w,'Polyhedron')
		error('P_w must be a Polyhedron() object from MPT3.')
	end

	if ~isa(P_v,'Polyhedron')
		error('P_v must be a Polyhedron() object from MPT3.')
	end

	if nargin > 6
		verbosity = varargin{7};
	end		

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	if ~exist('verbosity')
		verbosity = 0;
	end

	time_instants = [0:T-1];
	all_schedules = nchoosek( time_instants , MCardinality );

	bigM = 10^5;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%% Create Disturbance Sets %%

	P_wT = 1; P_vT = 1;
	for t = 1:T
		%Construct the sets by products
		P_wT = P_wT * P_w;
		P_vT = P_vT * P_v;
	end
	P_eta = P_wT * P_vT * P_x0;

	%% Create Optimization Variables %%

	[J,S,~] = ds.get_mpc_matrices(T,all_schedules(1,:) ); %Define MPC Matrice

	for schedule_idx = 1:size(all_schedules,1)
		%Get the schedule at index
		M_i = all_schedules(schedule_idx,:);

		%Get the C-bar matrix for this schedule
		[~,~,C_bar{schedule_idx}] = ds.get_mpc_matrices(T,M_i);

		%
		Q{schedule_idx} = sdpvar(ds.n_x*T,ds.n_x*T,'full');
		r{schedule_idx} = sdpvar(ds.n_x*T,1);

		Lambda{schedule_idx} = sdpvar( 2*ds.n_x*(T+1) , size(P_eta.A,1) , 'full' );

	end

	alpha0 = sdpvar(T+1,1,'full');
	alpha_bar = kron(alpha0,ones(ds.n_x,1));

	alpha_bar = [ alpha_bar ; alpha_bar ];

	s = binvar(size(all_schedules,1),1);

	if verbosity > 0
		disp('- Created optimization variables.')
	end

	%% Write Optimization Constraints %%

	cg = ConstraintGenerator();

	dual_constraints = [];
	low_diag_constr = [];
	positivity_constr = [];

	for schedule_idx = 1:size(all_schedules,1)

		P_xi_w = S + S*Q{schedule_idx}*C_bar{schedule_idx}*S;
		P_xi_v = S*Q{schedule_idx};
		P_xi_xi0 = ( eye(size(S,1)) + S*Q{schedule_idx}*C_bar{schedule_idx} )*J;

		R_T_mT = [ eye(ds.n_x*(T+1)) ; -eye(ds.n_x*(T+1)) ];

		dual_constraints = 	dual_constraints + ...
							[ Lambda{schedule_idx} * P_eta.A == R_T_mT * [ P_xi_w , P_xi_v , P_xi_xi0 ] ] + ...
							[ Lambda{schedule_idx} * P_eta.b <= alpha_bar + bigM*(1-s(schedule_idx))*ones(size(alpha_bar)) - R_T_mT * S*r{schedule_idx} ];

		low_diag_constr = low_diag_constr + cg.get_block_lower_diagonal_constraint_on(Q{schedule_idx},[ds.n_x,ds.n_y]);

		positivity_constr = positivity_constr + [ Lambda{schedule_idx} >= 0 , alpha0 >= 0 ];

	end

	binary_constraint = [ sum(s) == 1 ];

	%% Write Optimization Problem %%

	ops0 = sdpsettings('verbose',verbosity);
	diagnostics = optimize(	dual_constraints+low_diag_constr+positivity_constr+binary_constraint , ...
							sum(alpha0), ...
							ops0);

	%% Define Outputs %%

	%Find the selected gains
	optimal_schedule_idx = find( value(s) , 1 );

	Q0 = value(Q{optimal_schedule_idx});
	r0 = value(r{optimal_schedule_idx});

	F = inv(eye(size(Q0,1)) + Q0*C_bar{optimal_schedule_idx}*S ) * Q0;
	f = inv(eye(size(Q0,1)) + Q0*C_bar{optimal_schedule_idx}*S ) * r0;

	opt_cost = sum( value(alpha0) );

	M_opt = all_schedules(optimal_schedule_idx,:);


end