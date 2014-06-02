function Model = load_opengm_model(opengmFilename, dataset)

    if exist(opengmFilename, 'file') == 0
        error(['File ', opengmFilename, ' does not exist!'])
    end

    if nargin == 1
        warning('Dataset is unspecified, using "gm".')
        dataset = 'gm';
    end

    gm = openGMModel;
    gm.load(opengmFilename, dataset);

    numNodes = gm.numberOfVariables();
    modeSizes = zeros(1, numNodes);
    for iNode = 1:numNodes
        modeSizes(iNode) = gm.numberOfLabels(iNode - 1);
    end

    numFactors = gm.numberOfFactors();
    libdaiFactors = cell(numFactors, 1);
    for iFactor = 1:numFactors
        % In OpenGM potentials (-logarithms of factors) called factors.
        % Thats why we are getting currPotentialTable with getFactorTable method.
        [currPotentialTable, vars] = gm.getFactorTable(iFactor - 1);
        assert(issorted(vars), 'Error: we assumed that each factor in OpenGM model ',...
                               'has sorted list of variables, but it is not the case!');
        libdaiFactors{iFactor}.Member = vars;
        libdaiFactors{iFactor}.P = exp(-currPotentialTable);
    end

    Model.libdaiFactors = libdaiFactors;
    Model.numNodes = numNodes;
    Model.modeSizes = modeSizes;
    Model.numFactors = numFactors;

    Model.type = 'OpenGM';
    Model.description = ['OpenGM model loaded from file "', opengmFilename, '" (dataset "', dataset, '")'];
end