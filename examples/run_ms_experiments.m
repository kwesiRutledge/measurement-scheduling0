function [results] = run_experiments( varargin )
	% 	run_experiments
	%		How to run the ms_experimentX.m scripts programmatically
	%

	%% Add Library Folder to the Path
	matlab_library_rel_path = '../lib/matlab/';
	if isempty(strfind(path,matlab_library_rel_path))
		addpath(matlab_library_rel_path)
	end

	matlab_experiments_folder = 'matlab_experiments'
	if isempty(strfind(path,['./' matlab_experiments_folder ]))
		addpath(['./' matlab_experiments_folder ])
	end

	add_libraries('tbxmanager','gurobi')

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