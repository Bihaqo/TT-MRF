function [logSum, logAbsDiffBound] = tt_log_sum_prod(factorArray, roundingPrecision, varargin)
    % Find approximate log sum of elementwise product of tensors:
    % tt_log_sum_prod(factorArray, roundingPrecision) is an approximation to:
    %   log(sum(factorArray{1} .* factorArray{2} .* ... .* factorArray{end}))
    %
    % Optional arguments:
    %       o order  -- ['right_to_left'] in that order we would multiply vectors and matrices.
    %                                   Could be 'left_to_right' or 'right_to_left'.
    %       o rmax   -- [inf] maximal rank for intermediate tensors.
    %       o mv     -- ['exact'] or 'amen'. Matrix by vector multiplication method.
    %                                'exact' -- multiply exactly and then apply tt_rounding
    %                                'amen' -- use amen_mv method
    %       o verb   -- [1] verbosity level, 0-silent, 1-full info.
    %

    % Default values
    multiplicationOrder = 'right_to_left';
    rmax = inf;
    isAmenMatrixByVector = false;
    verboseLevel = 1;

    isMvSpecified = false;


    % Read parameters from input sequence
    for i=1:2:length(varargin)-1
        switch lower(varargin{i})
            case 'order'
                multiplicationOrder = varargin{i+1};
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

                isMvSpecified = true;
            case 'verb'
                verboseLevel = varargin{i+1};

            otherwise
                warning('Unrecognized option: %s\n',varargin{i});
        end
    end

    if isAmenMatrixByVector && rmax ~= inf
        error('You cannot specify rmax with amen matrix-by-vector multiplication method.')
    end

    if nargout >= 2 && isAmenMatrixByVector
        error('Bounds on absolute error of the answer is not supported with amen matrix-by-vector multiplication method.');
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

    logNormEstimateArray = zeros(numVariables, 1);
    matrixArray = cell(numVariables, 1);
    for iVariable = 1:numVariables
        currMatrix = [];
        currModeSize = factorModeSizeArray(iVariable);
        rowRankArray = factorRankArray(iVariable, :);
        colRankArray = factorRankArray(iVariable + 1, :);
        tp1 = tt_zeros(rowRankArray .* colRankArray);
        tp1 = tt_matrix(tp1, rowRankArray(:), colRankArray(:));
        currNormEstimate = 0;
        for variableValue = 1 : currModeSize
            tp1LogNormEstimate = 0;
            for iFactor = 1 : numFactors
                currCore = factorCoreArray{iVariable, iFactor};
                [r1, ~, r2] = size(currCore);
                cr1 = reshape(currCore(:, variableValue, :), [1, r1, r2, 1]);
                tp1{iFactor} = cr1;

                tp1LogNormEstimate = tp1LogNormEstimate + log(norm(reshape(cr1, r1, r2)));
            end
            currMatrix = currMatrix + tp1;
            currNormEstimate = currNormEstimate + exp(tp1LogNormEstimate);
        end
        matrixArray{iVariable} = currMatrix;
        logNormEstimateArray(iVariable) = log(currNormEstimate);
    end

    if verboseLevel >= 1
        sufficientPrecision = exp(log(roundingPrecision)  - sum(logNormEstimateArray) - log(numVariables - 1));
        disp(['To guarantee ', num2str(roundingPrecision), ' abs error use ', num2str(sufficientPrecision), ' rounding precision.']);
    end


    logAbsDiffBound = 0;
    fLogNormArray = zeros(numVariables, 1);
    if (strcmp(multiplicationOrder, 'left_to_right'))
        absDiffBound = 0;
        if isAmenMatrixByVector
            f = matrixArray{1}.tt;
        else
            f = matrixArray{1};
        end
        for iVariable = 2 : numVariables
            fNorm = norm(f);
            f = f / fNorm;
            fLogNormArray(iVariable) = log(fNorm);
            if isAmenMatrixByVector
                f = amen_mv(transpose(matrixArray{iVariable}), f, roundingPrecision, 'verb', verboseLevel);
            else
                f = f * matrixArray{iVariable};
                fRounded = round(f, roundingPrecision, rmax);
                absDiffBound = absDiffBound + exp(log(norm(f - fRounded)) + sum(logNormEstimateArray(iVariable+1 : end)) - sum(fLogNormArray(1 : iVariable-1)));
                f = fRounded;
            end
            if verboseLevel >= 1
                fprintf('i = %d, erank = %3.1f, max_rank = %d \n', iVariable, erank(f), max(rank(f)));
            end
        end
        logAbsDiffBound = log(absDiffBound);
    elseif (strcmp(multiplicationOrder, 'right_to_left'))
        absDiffBound = 0;
        if isAmenMatrixByVector
            f = matrixArray{numVariables}.tt;
        else
            f = matrixArray{numVariables};
        end
        for iVariable = (numVariables-1):-1:1
            fNorm = norm(f);
            f = f / fNorm;
            fLogNormArray(iVariable) = log(fNorm);
            if isAmenMatrixByVector
                f = amen_mv(matrixArray{iVariable}, f, roundingPrecision, 'verb', verboseLevel);
            else
                f = matrixArray{iVariable} * f;
                fRounded = round(f, roundingPrecision, rmax);
                absDiffBound = absDiffBound + exp(log(norm(f - fRounded)) + sum(logNormEstimateArray(1 : iVariable-1)) - sum(fLogNormArray(iVariable+1 : end)));
                f = fRounded;
            end
            if verboseLevel >= 1
                fprintf('i = %d, erank = %3.1f, max_rank = %d \n', iVariable, erank(f), max(rank(f)));
            end
        end
        logAbsDiffBound = log(absDiffBound);
    else
        error('Unsupported multiplication order.')
    end

    logSum = sum(fLogNormArray) + log(full(f));
end
