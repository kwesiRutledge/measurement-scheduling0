function [results] = ms_experiment1(varargin)
	%ms_experiment1.m
	%Description:
	%	Create a linearization of the covariance dynamics.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	A1 = diag([0.2;1.3]);

	P = 1.7*eye(dim);
	Q1 = diag([1.5;0.7]);

	results.constants.A1 = A1;
	results.constants.P = P;
	results.constants.Q1 = Q1;

	experiment_name = 'ms_experiment1';

	plot_conv_sets = false;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	disp(['Beginning ' experiment_name '.'])
	disp(' ')
	
	[Ap,q] = get_vector_cov_dynamics(A1,Q1)

	results.exp1.Ap = Ap;
	results.exp1.q = q;

	%% Draw the PSD Cone in R3
	disp('2. Using YALMIP to plot the PSD cone in R3.')

	p_bar = sdpvar(3,1,'full');
	constr = [ [p_bar(1),p_bar(2);p_bar(2),p_bar(3)] >= 0 ];
	Y = YSet(p_bar,constr);

	if plot_conv_sets
		figure;
		plot(Y)

		saveas(gcf,'results/experiment1/transformed_psd_cone','epsc')
		disp('- Successfully plotted the full cone.')
	end

	clear p_bar constr

	p_bar = sdpvar(3,1,'full');
	constr = [ [p_bar(1),p_bar(2);p_bar(2),p_bar(3)] >= 0 ] + [ [p_bar(1),p_bar(2);p_bar(2),p_bar(3)] <= 3*eye(dim) ];
	Y2 = YSet(p_bar,constr);

	if plot_conv_sets
		figure;
		plot(Y2)
		saveas(gcf,'results/experiment1/bounded_psd_cone','epsc')
		disp('- Successfully plotted a bounded version of the full cone.')
	end

	results.exp2.Y = Y;
	results.exp2.Y2 = Y2;

	%% Plot the path of the Covariance in this 3-D space
	disp('3. Plotting the trajectory of one matrix.')

	th = deg2rad(10);
	A2 = 0.95*[cos(th),-sin(th); sin(th), cos(th)];

	[Ap2,q] = get_vector_cov_dynamics(A2,Q1);

	%Create the Path of the covariance
	p0 = [1.7;0;1.7];
	T = 500;
	p = [p0];
	for t = [1:T-1]
		p = [p, Ap2*p(:,t)+q ];
	end 

	%Create the set of PSD Matrices
	diag_val = -10;
	bounding_diag_val = 12;
	P_offset = bounding_diag_val*eye(dim);
	p_bar = sdpvar(3,1,'full');
	constr = [ [p_bar(1),p_bar(2);p_bar(2),p_bar(3)] <= P_offset ] + [ [p_bar(1),p_bar(2);p_bar(2),p_bar(3)] >= diag_val*eye(dim) ];
	Y3 = YSet(p_bar,constr);

	grid_n = 10;
	lw0 = 2; %Default LineWidth = 0.5
	figure;
	hold;
	scatter3(p(1,:),p(2,:),p(3,:),'rx','LineWidth',lw0)
	Y3.plot('Alpha',0.5,'Color','lightblue','linestyle','--','Grid',grid_n)
	axis([ 	min([0,p(1,:)]) , max([bounding_diag_val+1,p(1,:)]) , ...
	 		min([0,p(2,:)]) , max([bounding_diag_val+1,p(2,:)]) , ...
	 		min([0,p(3,:)]) , max([bounding_diag_val+1,p(3,:)]) ])

	set(gcf,'units','Normalized','Position',[0 0 1 1])
	saveas(gcf,['results/experiment1/covariance_matrix_traj' ],'epsc')

	savefig('results/experiment1/covariance_matrix_traj_fig')

	view(0,90)
	saveas(gcf,['results/experiment1/covariance_matrix_traj_view2' ],'epsc')

	results.exp3.p = p;
	results.exp3.T = T;
	results.exp3.p0 = p0;
	results.exp3.Y3 = Y3;


end