function [logZ_TT, logAbsErrorBound_TT, time_TT, logZ_JT, ...
		  logZ_BP, time_BP, logZ_MF, time_MF, ...
		  logZ_TREEEP, time_TREEEP, logZ_AIS, time_AIS] = compute_answers(modelCellArray, ttRoundingPrecision)

	% Compute log partition function with different methods
	% for all models from modelCellArray.

	% Tensor Train.
	ttFunc = @(Model) compute_log_z(Model, ttRoundingPrecision);
	[logZ_TT, logAbsErrorBound_TT, time_TT] = cellfun(ttFunc, modelCellArray);

	% Junction tree (exact method).
	function [logZ_JT] = jtreeFunc(Model)
		[logZ_JT, ~, ~] = dai_jtree(Model.libdaiFactors, { 1 }, '[updates=HUGIN]');
	end
	[logZ_JT] = cellfun(@jtreeFunc, modelCellArray);

	% Belief Propagation.
	function [logZ_BP, time_BP] = bpFunc(Model)
		tic;
		[logZ_BP, ~, ~] = dai(Model.libdaiFactors, 'BP', '[maxiter=1000,tol=1e-8,verbose=0,logdomain=1,updates=SEQMAX]');
		time_BP = toc;
	end
	[logZ_BP, time_BP] = cellfun(@bpFunc, modelCellArray);

	% Mean Field.
	function [logZ_MF, time_MF] = mfFunc(Model)
		tic;
		[logZ_MF, ~, ~] = dai(Model.libdaiFactors, 'MF', '[maxiter=1000,tol=1e-8,verbose=0,updates=NAIVE]');
		time_MF = toc;
	end
	[logZ_MF, time_MF] = cellfun(@mfFunc, modelCellArray);

	% Tree Expectation Propagation.
	function [logZ_TREEEP, time_TREEEP] = treeepFunc(Model)
		tic;
		[logZ_TREEEP, ~, ~] = dai(Model.libdaiFactors, 'TREEEP', '[maxiter=1000,tol=1e-8,type=ORG,verbose=0]');
		time_TREEEP = toc;
	end
	[logZ_TREEEP, time_TREEEP] = cellfun(@treeepFunc, modelCellArray);

	% Annealing important sampling.
	function [logZ_AIS, time_AIS] = aisFunc(Model)
		if strcmp(Model.edgesType, 'number')
			J = Model.edgeWeights(1);
			tic;
			logZ_AIS = ZgibbsIsing4(Model.unaryWeights, J, 1/Model.temperature, 1000, 70);
			time_AIS = toc;
		else
			warning('AIS does not support multiple pairwise weights');
			logZ_AIS = 0;
			time_AIS = 0;
		end
	end
	[logZ_AIS, time_AIS] = cellfun(@aisFunc, modelCellArray);
end
