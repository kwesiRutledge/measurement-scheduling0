%test_DynamicalSystem.m
%Description:
%	Tests the function DynamicalSystem().

function tests = test_DynamicalSystem
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test_DynamicalSystem1(testCase)
	%Description:
	%

	%% Include Libraries
	addpath(genpath('../../lib/matlab/'))

	%% Constants
    
    dim = 2;

	A1 = diag([0.2;1.3]);
	b = rand([2,1]);
	C = [ 0.5, 1.5; 0 , 1 ];

	P = 1.7*eye(dim);
	Q1 = diag([1.5;0.7]);
	R1 = diag([0.1;0.2]);
    
    x0_mean = zeros(dim,1);
    
    dyn0 = DynamicalSystem( A1 , b , C , Q1 , R1 , x0_mean , P );

	%% Algorithm
	assert( all(all(  dyn0.A == A1 )) & all(all(  dyn0.b == b )) & ...
            all(all(  dyn0.C == C )) & all(all(  dyn0.Q == Q1 )) & ...
            all(all(  dyn0.R == R1 )) & all(all(  dyn0.x0_bar == x0_mean )) & ...
            all(all(  dyn0.P0_bar == P )) )
    
