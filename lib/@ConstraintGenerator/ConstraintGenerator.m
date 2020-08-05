classdef ConstraintGenerator
	%Description:
	%	This class abstracts away how we make constraints for some optimization problems.
	properties

	end

	methods
		%Constructor
		function [generator] = constr_gen()
			%Description:
			%	The constructor for this class does nothing.

			disp('ConstraintGenerator() does not initialize any member variables.')
			% generator = [];
		
		end

		function [ blk_low_diag_constr ] = get_block_lower_diagonal_constraint_on( obj , Q_in , block_dims )
			%Description:
			%	Assumes that the variable Q_in is a YALMIP variable.
			%	If Q_in is a matrix with block_dims(1) x block_dims(2) blocks, then the function returns a constraint
			%	that is lower block diagonal according to the blocks.

			%% Input Processing %%

			n_u = block_dims(1);
			n_y = block_dims(2);

			if (rem( size(Q_in,1) , n_u ) ~= 0) || (rem( size(Q_in,2) , n_y ) ~= 0)
				error(['The input matrix cannot be decomposed into blocks of dimension ' num2str(n_u) ' x ' num2str(n_y) '.' ])
			end

			%% Algorithm %%
			blk_low_diag_constr = [];
			for blk_row_idx = 1:(size(Q_in,1) / n_u)-1
				%disp(['blk_row_idx = ' num2str(blk_row_idx) ])
				blk_low_diag_constr = blk_low_diag_constr + [Q_in((blk_row_idx-1)*n_u+[1:n_u],[blk_row_idx*n_y+1:end]) == 0];
			end

		end
		
	end
end