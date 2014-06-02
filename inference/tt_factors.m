function factors = tt_factors(Model)
    % Find TT representation of model factors.
	% 

    factors = cell(Model.numFactors, 1);
    for iFactor = 1 : Model.numFactors
        vars = Model.libdaiFactors{iFactor}.Member;
        % TT related routines assume that variable numbers are starting from 1.
        vars = vars + 1;
        table = Model.libdaiFactors{iFactor}.P;
        currFactor = tt_tensor(table);
        
        % Build enlarged tensor with nonessential dimensions.
        factors{iFactor} = add_non_essential_dims(currFactor, Model.modeSizes, vars);
    end
end