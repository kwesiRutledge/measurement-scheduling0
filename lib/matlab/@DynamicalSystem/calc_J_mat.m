function [J] = calc_J_mat(obj,T)
	%calc_J_mat
	%	Description:
	%		Calculates the matrix that describes how the initial state x(t_0)
	%		affects the state trajectory [x(t_0);x(t_0+1); ... ;x(T-1);x(T)]
	%
	%	Usage:
	%		J = calc_x0_mat(ds,T)

	%Constants
	A = obj.A;

	%%Input Checking

	if nargin ~= 2
		error('Algorithm was made to accept only 2 arguments.')
	end

	%% Algorithm

	%Process Inputs

	%Creating matrix.
	J = eye(obj.n_x);

	for i = 1 : T
		J = [J; (A^i)];
	end
	
end