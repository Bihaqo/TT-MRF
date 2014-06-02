function [logZ, logAbsErrorBound, time] = compute_log_z(Model, roundingPrecision, varargin)
	% Find logarithm of the partition function Z for the graphical model using Tensor Train decomposition.
	% 
	
	tic;    
    factors = tt_factors(Model);
    if nargout == 1
	    logZ = tt_log_sum_prod(factors, roundingPrecision, varargin{:});
	else
	    [logZ, logAbsErrorBound] = tt_log_sum_prod(factors, roundingPrecision, varargin{:});
	end
    time = toc;
end