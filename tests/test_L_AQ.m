%test_L_AQ.m
%Description:
%	Tests the function L_AQ().

function tests = test_L_AQ
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test1(testCase)
	%Description:
	%

	%% Include Libraries
	addpath(genpath('../lib/'))

	%% Constants
	P1 = [1,2;3,4];
	Q1 = [5,6;7,8];
	A1 = eye(2);

	%% Algorithm
	assert( all(all( [6,8;10,12] == L_AQ( A1 , Q1 , P1 ))) )
