function [unaryMarginals, logZ, time] = compute_marginals(Model, roundingPrecision, varargin)
    % Approximate the marginal distributions and the partition function of
    % an MRF defined by a using the Tensor Train approach.
    %
    % Optional arguments:
    %       o mv     -- [exact] or amen. That method to use for matrix by vector multiplication.
    %       o rmax   -- [inf] maximal rank for intermediate tensors.
    %       o verb   -- [1] verbosity level, 0-silent, 1-full info.
    %
    % Output:
    %       o unaryMarginals -- d x n matrix, marginal distribution for each variable.
    %       o logZ           -- d x 1 vector, log(Z) computed by summing unnormalized marginal distribution.
    %

    tic;
    factors = tt_factors(Model);
    [unaryMarginals, logZ] = tt_marginals_from_factors(factors, roundingPrecision, varargin{:});
    time = toc;
end
