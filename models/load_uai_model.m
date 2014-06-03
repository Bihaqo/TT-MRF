function Model = load_uai_model(uaiFilename)

    if exist(uaiFilename, 'file') == 0
        error(['File ', uaiFilename, ' does not exist!'])
    end


    fid = fopen(uaiFilename);
    mrfType = fscanf(fid, '%s', 1);
    if ~strcmp(mrfType, 'MARKOV')
        error(['Unsupported MRF type:', mrfType])
    end
    numNodes = fscanf(fid, '%d', 1);
    modeSizes = fscanf(fid, '%d', numNodes)';

    numFactors = fscanf(fid, '%d', 1);
    factorVars = cell(numFactors, 1);
    for iFactor = 1:numFactors
        numCurrFactorVars = fscanf(fid, '%d', 1);
        factorVars{iFactor} = fscanf(fid, '%d', numCurrFactorVars) + 1;
    end

    libdaiFactors = cell(numFactors, 1);
    for iFactor = 1:numFactors
        currFactorSize = fscanf(fid, '%d', 1);
        currFactorTable = fscanf(fid, '%f', currFactorSize);
        currModeSizes = modeSizes(factorVars{iFactor});
        if length(currModeSizes) > 1
            % Memory order in uai format is different from the matlab one, so we need to rearrange dimensions.
            currFactorTable = reshape(currFactorTable, currModeSizes(end:-1:1));
            new_order = length(currModeSizes):-1:1;
            currFactorTable = permute(currFactorTable, new_order);
        end
        libdaiFactors{iFactor}.Member = factorVars{iFactor}' - 1;
        libdaiFactors{iFactor}.P = currFactorTable;
    end
    fclose(fid);

    Model.libdaiFactors = libdaiFactors;
    Model.numNodes = numNodes;
    Model.modeSizes = modeSizes;
    Model.numFactors = numFactors;

    Model.type = 'UAI';
    Model.description = ['UAI model loaded from file "', uaiFilename, '"'];
end
