function [ res ] = L_AQ(A,Q,P)
	%Description:
	%	Assume that the initial estimation error covariance is P.
	%	Update the estimation error covariance according to the linear equation:
	%		L_AQ(P) = A P A' + Q

	res = A*P*A'+Q;

end