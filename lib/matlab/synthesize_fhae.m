function [ F , f , opt_cost ] = synthesize_fhae( varargin )
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
%Example Usage:
%	[F,f,cost_out] = synthesize_fhae( ds , T , M , P_x0 , P_w , P_v )
%
%Inputs:
%	

	%% Input Processing %%
	ds = varargin{1};
	T = varargin{2};
	M = varargin{3};

	P_x0 = varargin{4};
	P_w = varargin{5};
	P_v = varargin{6};

	if ~all( T >= M )
		warning('All measurement times should be less than or equal to the time horizon.')
	end

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

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%% Create Optimization Variables %%

	[J,S,C_bar] = ds.get_mpc_matrices(T,M); %Define MPC Matrices

	Q = sdpvar(ds.n_x*T,ds.n_x*T,'full');
	r = sdpvar(ds.n_x*T,1);

	P_wT = 1; P_vT = 1;
	for t = 1:T
		%Construct the sets by products
		P_wT = P_wT * P_w;
		P_vT = P_vT * P_v;
	end
	P_eta = P_wT * P_vT * P_x0;

	Lambda = sdpvar( 2*ds.n_x*(T+1) , size(P_eta.A,1) , 'full' );

	alpha0 = sdpvar(T+1,1,'full');
	alpha_bar = kron(alpha0,ones(ds.n_x,1));

	alpha_bar = [ alpha_bar ; alpha_bar ];

	if verbosity > 0
		disp('- Created optimization variables.')
	end

	%% Write Optimization Constraints %%

	cg = ConstraintGenerator();

	P_xi_w = S + S*Q*C_bar*S;
	P_xi_v = S*Q;
	P_xi_xi0 = ( eye(size(S,1)) + S*Q*C_bar )*J;

	R_T_mT = [ eye(ds.n_x*(T+1)) ; -eye(ds.n_x*(T+1)) ];

	dual_constraints = 	[ Lambda * P_eta.A == R_T_mT * [ P_xi_w , P_xi_v , P_xi_xi0 ] ] + ...
						[ Lambda * P_eta.b <= alpha_bar - R_T_mT * S*r ];

	low_diag_constr = cg.get_block_lower_diagonal_constraint_on(Q,[ds.n_x,ds.n_y]);

	positivity_constr = [ Lambda >= 0 , alpha0 >= 0 ];

	%% Write Optimization Problem %%

	ops0 = sdpsettings('verbose',verbosity);
	diagnostics = optimize(	dual_constraints+low_diag_constr+positivity_constr , ...
							sum(alpha0), ...
							ops0);

	%% Define Outputs %%

	Q0 = value(Q);
	r0 = value(r);

	F = inv(eye(size(Q0,1)) + Q0*C_bar*S ) * Q0;
	f = inv(eye(size(Q0,1)) + Q0*C_bar*S ) * r0;

	opt_cost = sum( value(alpha0) );


end