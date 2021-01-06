%test_L_AQ.m
%Description:
%	Tests the function L_AQ().

function tests = test_L_AQ
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test_L_AQ1(testCase)
	%Description:
	%

	%% Include Libraries
	addpath(genpath('../../lib/matlab/'))

	%% Constants
	P1 = [1,2;3,4];
	Q1 = [5,6;7,8];
	A1 = eye(2);

	%% Algorithm
	assert( all(all( [6,8;10,12] == L_AQ( A1 , Q1 , P1 ))) )
    
function test_L_AQ2(testCase)
	%Description:
	%   Providing a test case where L_AQ(A,Q,P) of a matrix P is not >=
    %   P.

	%% Include Libraries
	addpath(genpath('../../lib/matlab/'))

	%% Constants
    dim = 2;
	Q = 0.3*eye(dim);
	P = 1.7*eye(dim);
	A1 = diag([0.2;0.5]);

	%% Algorithm
    % The eigenvalues
	assert( any(eig( L_AQ(A1,Q,P) - P ) < 0 ) )
    
function test_L_AQ3(testCase)
	%Description:
	%   Providing a test case where L_AQ(A,Q,P) of a matrix P IS >=
    %   P.

	%% Include Libraries
	addpath(genpath('../../lib/matlab/'))

	%% Constants
    dim = 2;
	Q = 0.3*eye(dim);
	P = 1.7*eye(dim);
	A2 = eye(dim);

	%% Algorithm
    % The eigenvalues
	assert( all(eig( L_AQ(A2,Q,P) - P ) >= 0 ) )
