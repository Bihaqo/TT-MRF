function [logZ, time, logAbsErrorBound] = compute_log_z_array(modelCellArray, roundingPrecision)
	% Compute logZ for a cell array of models with the Tensor Train based
	% method. If the Tensor Train method fails for some model, show warning,
	% save the model into the errors.mat file, return 0 as logZ and continue.
	% Shows progress bar.
	% 

	counter = cell(size(modelCellArray));
	counter(:) = num2cell((1:numel(modelCellArray)) / numel(modelCellArray));

	function [logZ, logAbsErrorBound, time] = ttFunc(Model, counter)
		persistent waitbarH;
		if isempty(waitbarH) || ~ishandle(waitbarH)
		  waitbarH = waitbar(0, ['Computing logZ with the Tensor Train method ', ...
			  	'(precision is ', num2str(roundingPrecision), ').']); 
		end

		try
			[logZ, logAbsErrorBound, time] = compute_log_z(Model, roundingPrecision);
		catch err
			warning('TT method didn''t succeed on the current model, proably out', ...
					'of memory error. Model saved to the errors.mat file.');
			logZ = 0;
			logAbsErrorBound = inf;
			time = 0;
			if exist('errors.mat', 'file')
				load('errors.mat');
				badModels{end + 1} = Model;
			else
				badModels = {Model};
			end
			save('errors.mat', 'badModels');
		end
			
		disp(['Computed loZ with TT method in ', num2str(time), ' seconds.']);
		
		waitbar(counter, waitbarH);
		if counter == 1
			% Last call.
		    close(waitbarH);
		end
	end
	[logZ, logAbsErrorBound, time] = cellfun(@ttFunc, modelCellArray, counter);
end
