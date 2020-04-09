function [S] = calc_w_effect_mat(ds,T)
	%calc_w_effect_mat
	%	Description:
	%		Calculates the matrix that describes how the trajectory of noise
	%		[w(t_0),w(t_0+1),...,w(T-1)] affects the state trajectory
	%		[x(t_0);x(t_0+1); ... ;x(T-1);x(T)]
	%
	%	Usage:
	%		S = calc_u_effect_mat(ds,T);

	%% Constants %%

	A = ds.A;

	%% Algorithm %%		

	%Constants
	S = zeros(ds.n_x*(T+1),ds.n_x*T);

	for i = 1: T

		temp = [];
		for k = 1:i
			temp = [temp A^(i-k)];
		end
		
		S([i*ds.n_x+1:(i+1)*ds.n_x],[1:ds.n_x*i]) = temp;
	
	end

end