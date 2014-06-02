function Model = generate_spin_glass_model(n, m, temperature, varargin)
    % Build problem instance for spin glass model.
    % 
    % Input:
    %   n, m        -- positive integers, size of the grid
    %   temperature -- number, temperature of the model
    %   J           -- [1] pairwise weights, can be number (all weights equal this number)
    %                                        can be 'rand' (all weights are generated from uniform distribution)
    %   J_distr     -- [-1, 1] vector with two numbers, support of the J uniform distribution.
    %   H           -- [0] unary weights, can be number (all weights equal this number)
    %                                     can be 'rand' (all weights are generated from uniform distribution)
    %                                     can be matrix of size n by m with desired weights.
    %   

    [edgesType, edgesDistr, unaryType, unaryDistr] = process_options(varargin, 'J', 1, 'J_distr', [-1, 1],...
                                                            'H', 0, 'H_distr', [-1, 1]);


    if strcmp(unaryType, 'rand')
        unaryWeights = unifrnd(unaryDistr(1), unaryDistr(2), n, m);
        unaryTextType = 'rand';
    elseif isscalar(unaryType)
        unaryWeights = unaryType * ones(n, m);
        unaryTextType = 'number';
    elseif ismatrix(unaryType)
        if size(unaryType) == [n, m]
            unaryWeights = unaryType;
            unaryTextType = 'matrix';
        else
            error('H has inappropriate dimensions.')
        end
    else
        error('Unknown H type.');
    end

    [edgesFrom, edgesTo] = get_neighborhood4(n, m);
    numEdges = length(edgesFrom);
    if isscalar(edgesType)
        edgeWeights = edgesType * ones(numEdges, 1);
        edgesTextType = 'number';
    elseif strcmp(edgesType, 'rand')
        edgeWeights = unifrnd(edgesDistr(1), edgesDistr(2), numEdges, 1);
        edgesTextType = 'rand';
    else
        error('Unknown J type.');
    end

    numNodes = n * m;
    % There are unary and pairwise factors.
    numFactors = numNodes + numEdges;
    libdaiFactors = cell(numFactors, 1);

    unaryConfigurationMatrix = [-1; 1];
    iFactor = 1;
    for iNode = 1 : numNodes
        if abs(unaryWeights(iNode)) > 1e-20
            libdaiFactors{iFactor}.Member = iNode - 1;
            libdaiFactors{iFactor}.P = -unaryConfigurationMatrix * unaryWeights(iNode);
            iFactor = iFactor + 1;
        end
    end

    pairwiseConfigurationMatrix = [ 1    -1;
                           -1     1];
    for iEdge = 1 : numEdges
        if abs(edgeWeights(iEdge)) > 1e-20
            libdaiFactors{iFactor}.Member = [edgesFrom(iEdge) - 1, edgesTo(iEdge) - 1];
            libdaiFactors{iFactor}.P = -pairwiseConfigurationMatrix * edgeWeights(iEdge);
            if (edgesFrom(iEdge) > edgesTo(iEdge))
                % LibDAI format assumes that nodes in Member field are in sorted order.
                libdaiFactors{iFactor}.Member = sort(libdaiFactors{iFactor}.Member);
                libdaiFactors{iFactor}.P = libdaiFactors{iFactor}.P';
            end
            iFactor = iFactor + 1;
        end            
    end
    % There can be less than numFactors factors, since some weights are zero.
    % Remove unused factors.
    libdaiFactors(iFactor : end) = [];
    numFactors = iFactor - 1;

    % Convert potentials to factors.
    for iFactor = 1 : numFactors
        libdaiFactors{iFactor}.P = exp( -libdaiFactors{iFactor}.P / temperature);
    end

    Model.libdaiFactors = libdaiFactors;
    Model.numNodes = numNodes;
    Model.modeSizes = 2 * ones(1, numNodes);
    Model.numFactors = numFactors;

    Model.type = 'Spin glass';
    Model.grid_n = n;
    Model.grid_m = m;
    Model.temperature = temperature;
    Model.unaryWeights = unaryWeights;
    Model.unaryType = unaryTextType;
    if strcmp(unaryType, 'rand')
        Model.unaryDistr = unaryDistr;
    end
    Model.edgeWeights = edgeWeights;
    Model.edgesType = edgesTextType;
    if strcmp(edgesType, 'rand')
        Model.edgesDistr = edgesDistr;
    end

    Model.description = [int2str(n), 'x', int2str(m), ' spin glass grid, temp = ',...
                    num2str(temperature), ', J '];
    if isscalar(edgesType)
        Model.description = [Model.description, '= ', num2str(edgesType)];
    elseif strcmp(edgesType, 'rand')
        Model.description = [Model.description, 'is uniform on [', num2str(edgesDistr(1)), ', ', num2str(edgesDistr(2)), ']'];
    end
    
    Model.description = [Model.description, ', H '];
    if strcmp(unaryType, 'rand')
        Model.description = [Model.description, 'is uniform on [', num2str(unaryDistr(1)), ', ', num2str(unaryDistr(2)), '].'];
    elseif isscalar(unaryType)
        Model.description = [Model.description, '= ', num2str(unaryType), '.'];
    elseif ismatrix(unaryType)
        Model.description = [Model.description, 'is provided by user.'];
    end
end
