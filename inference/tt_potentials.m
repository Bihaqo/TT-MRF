function potentials = tt_potentials(Model)
    % Find TT representation of model potentials (minus logarithm of factors).
    % 

    potentials = cell(Model.numFactors, 1);
    for iFactor = 1 : Model.numFactors
        vars = Model.libdaiFactors{iFactor}.Member;
        % TT related routines assume that variable numbers are starting from 1.
        vars = vars + 1;
        factorTable = Model.libdaiFactors{iFactor}.P;
        factorTable(factorTable < 1e-10) = 1e-10;
        potentialTable = -log(factorTable);
        currPotential = tt_tensor(potentialTable);
        % Build enlarged tensor with nonessential dimensions.
        potentials{iFactor} = add_non_essential_dims(currPotential, Model.modeSizes, vars);
    end
end