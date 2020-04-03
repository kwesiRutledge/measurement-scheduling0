function [Abar,q] = get_vector_cov_dynamics( A , Q )
	%Description:
	%	The dynamics for the estimation error covariance matrix over time is written in the following form in the system matrices
	%		P(t+1) = A P(t) A' + Q
	%	at open loop.
	%	This is linear and thus can be written in terms of linear equations like this:
	%		p(t+1) = Abar p(t) + q
	%

	%% Constants %%

	dim = size(A,1);

	%% Algorithm %%

	if dim == 2
		a1 = A(1,1);
		a2 = A(1,2);
		a3 = A(2,1);
		a4 = A(2,2);

		Abar = [ 	a1^2,	2*a1*a2,		a2^2	;
					a1*a3,  a2*a3+a1*a4,	a2*a4	;
					a3^2,	2*a3*a4,		a4^2];

		q = [ Q(1,1) ; Q(1,2) ; Q(2,2) ];

	else
		error('This function was designed for 2-dimensional systems!')
	end

end