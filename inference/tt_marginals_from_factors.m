function [unaryMarginals, logZ] = tt_marginals_from_factors(factorArray, roundingPrecision, varargin)
    % Approximate the marginal distributions and the partition function of
    % an MRF defined by a product of factors from the factorArray.
    % The factors should be given in the TT-format.
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

    % Default values
    isAmenMatrixByVector = false;
    rmax = inf;
    verboseLevel = 1;

    % Read parameters from input sequence
    for i=1:2:length(varargin)-1
        switch lower(varargin{i})
            case 'rmax'
                rmax = varargin{i+1};
            case 'mv'
                method = varargin{i+1};
                if strcmp(method, 'exact')
                    isAmenMatrixByVector = false;
                elseif strcmp(method, 'amen')
                    isAmenMatrixByVector = true;
                else
                    error('Unsupported matrix-by-vector multiplication method.');
                end
            case 'verb'
                verboseLevel = varargin{i+1};

            otherwise
                warning('Unrecognized option: %s\n',varargin{i});
        end
    end

    if isAmenMatrixByVector && rmax ~= inf
        error('You cannot specify rmax with amen matrix-by-vector multiplication method.')
    end

    numVariables = factorArray{1}.d;
    numFactors = length(factorArray);

    factorRankArray = zeros(numVariables + 1, numFactors);
    factorModeSizeArray = zeros(numVariables, 1);
    factorCoreArray = cell(numVariables, numFactors);
    for iFactor = 1 : numFactors
        for iVariable = 1 : numVariables
            currCore = factorArray{iFactor}{iVariable};
            [r1, n1, r2] = size(currCore);
            factorRankArray(iVariable, iFactor) = r1;
            factorModeSizeArray(iVariable) = n1;
            factorCoreArray{iVariable, iFactor} = currCore;
        end
    end
    factorRankArray(numVariables+1, :) = 1;

    if any(factorModeSizeArray ~= factorModeSizeArray(1))
        error('MRFs with different mode sizes are unsupported.');
    end


    matrixArray = cell(numVariables, 1);
    for iVariable = 1:numVariables
        currMatrix = [];
        currModeSize = factorModeSizeArray(iVariable);
        rowRankArray = factorRankArray(iVariable, :);
        colRankArray = factorRankArray(iVariable + 1, :);
        tp1 = tt_zeros(rowRankArray .* colRankArray);
        tp1 = tt_matrix(tp1, rowRankArray(:), colRankArray(:));
        for variableValue = 1 : currModeSize
            for iFactor = 1 : numFactors
                currCore = factorCoreArray{iVariable, iFactor};
                [r1, ~, r2] = size(currCore);
                cr1 = reshape(currCore(:, variableValue, :), [1, r1, r2, 1]);
                tp1{iFactor} = cr1;
            end
            currMatrix = currMatrix + tp1;
        end
        matrixArray{iVariable} = currMatrix;
    end

    isTransposeLeftPref = false;

    leftPrefixArray = cell(numVariables, 1);
    leftPrefixLogNorm = zeros(numVariables, 1);
    if isAmenMatrixByVector
        f = matrixArray{1}.tt;
    else
        if isTransposeLeftPref
            f = matrixArray{1}';
        else
            f = matrixArray{1};
        end
    end
    if isTransposeLeftPref
        leftPrefixArray{1} = f';
    else
        leftPrefixArray{1} = f;
    end
    for iVariable = 2 : numVariables
        if isAmenMatrixByVector
            f = amen_mv(transpose(matrixArray{iVariable}), f, roundingPrecision, 'verb', verboseLevel);
        else
            if isTransposeLeftPref
                f = matrixArray{iVariable}' * f;
            else
                f = f * matrixArray{iVariable};
            end
            f = round(f, roundingPrecision, rmax);
        end
        currNorm = norm(f);
        leftPrefixLogNorm(iVariable) = leftPrefixLogNorm(iVariable - 1) + log(currNorm);
        f = f / currNorm;
        if isTransposeLeftPref
            leftPrefixArray{iVariable} = f';
        else
            leftPrefixArray{iVariable} = f;
        end
        if verboseLevel >= 1
            fprintf('i = %d, erank = %3.1f, max_rank = %d \n', iVariable, erank(f), max(rank(f)));
        end
    end

    rightPostfixArray = cell(numVariables, 1);
    rightPostfixLogNorm = zeros(numVariables, 1);
    if isAmenMatrixByVector
        f = matrixArray{numVariables}.tt;
    else
        f = matrixArray{numVariables};
    end
    rightPostfixArray{numVariables} = f;
    for iVariable = (numVariables - 1) : -1 : 1
        if isAmenMatrixByVector
            f = amen_mv(matrixArray{iVariable}, f, roundingPrecision, 'verb', verboseLevel);
        else
            f = matrixArray{iVariable} * f;
            f = round(f, roundingPrecision, rmax);
        end
        currNorm = norm(f);
        rightPostfixLogNorm(iVariable) = rightPostfixLogNorm(iVariable + 1) + log(currNorm);
        f = f / currNorm;
        rightPostfixArray{iVariable} = f;
        if verboseLevel >= 1
            fprintf('i = %d, erank = %3.1f, max_rank = %d \n', iVariable, erank(f), max(rank(f)));
        end
    end

    numValues = factorModeSizeArray(iVariable);
    unaryMarginals = zeros(numVariables, numValues);
    logZ = zeros(numVariables, 1);
    for iVariable = 1 : numVariables
        rowRankArray = factorRankArray(iVariable, :);
        colRankArray = factorRankArray(iVariable + 1, :);
        tp1 = tt_zeros(rowRankArray .* colRankArray);
        tp1 = tt_matrix(tp1, rowRankArray(:), colRankArray(:));
        if iVariable > 1
            currLeft = leftPrefixArray{iVariable - 1};
            currLeftLogNorm = leftPrefixLogNorm(iVariable - 1);
        else
            currLeft = 1;
            currLeftLogNorm = 0;
        end
        if iVariable < numVariables
            currRight = rightPostfixArray{iVariable + 1};
            currRightLogNorm = rightPostfixLogNorm(iVariable + 1);
        else
            currRight = 1;
            currRightLogNorm = 0;
        end

        for variableValue = 1 : numValues
            for iFactor = 1 : numFactors
                currCore = factorCoreArray{iVariable, iFactor};
                [r1, n1, r2] = size(currCore);
                cr1 = reshape(currCore(:, variableValue, :), [1, r1, r2, 1]);
                tp1{iFactor} = cr1;
            end
            if isAmenMatrixByVector
                if iVariable > 1
                    leftProd = amen_mv(transpose(tp1), currLeft, roundingPrecision, 'verb', verboseLevel);
                else
                    leftProd = tp1.tt;
                end
                if iVariable < numVariables
                    prob = dot(leftProd, currRight);
                else
                    prob = full(leftProd);
                end
                unaryMarginals(iVariable, variableValue) = prob;
            else
                unaryMarginals(iVariable, variableValue) = full(currLeft * tp1 * currRight);
            end
        end
        logZ(iVariable) = currLeftLogNorm + currRightLogNorm;
        logZ(iVariable) = logZ(iVariable) + log(sum(unaryMarginals(iVariable, :)));
    end

    % Normalize marginal distributions.
    unaryMarginals = bsxfun(@times, unaryMarginals, 1 ./ sum(unaryMarginals, 2));
end
