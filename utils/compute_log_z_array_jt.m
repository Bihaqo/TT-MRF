function [logZ] = compute_log_z_array_jt(modelCellArray)
	% Compute logZ for a cell array of models with the Junction tree method
	% from the libDAI library (this method is exact). Shows progress bar.
	% 

	counter = cell(size(modelCellArray));
	counter(:) = num2cell((1:numel(modelCellArray)) / numel(modelCellArray));

	function [logZ] = jtreeFunc(Model, counter)
		persistent waitbarH;
		if isempty(waitbarH) || ~ishandle(waitbarH)
		  waitbarH = waitbar(0, 'Computing logZ with the Junction Tree method.'); 
		end

		[logZ, ~, ~] = dai_jtree(Model.libdaiFactors, { 1 }, '[updates=HUGIN]');
		
		waitbar(counter, waitbarH);
		if counter == 1
			% Last call.
		    close(waitbarH);
		end
	end
	[logZ] = cellfun(@jtreeFunc, modelCellArray, counter);
end
