function [results] = run_experiments( varargin )
	% 	run_experiments
	%		How to run the ms_experimentX.m scripts programmatically
	%

	%% Add Library Folder to the Path
	if isempty(strfind(path,'./lib/'))
		addpath('./lib/')
	end

	if isempty(strfind(path,'./experiments/'))
		addpath('./experiments/')
	end

	add_libraries('tbxmanager')

	%% Constants
	base_name = 'ms_experiment';

	test_nums = varargin{1};
	if nargin < 2
		test_inputs = cell(size(varargin{1}));
	else
		test_inputs = varargin{2};
	end

	%%Run tests from 
	for k = 1 : length(test_nums)
		results{k} = eval([base_name num2str(test_nums(k)) '(' test_inputs{k} ')' ]);
	end

end