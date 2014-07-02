function modelCellArray = generate_data(n, m, temperatureArray, J, numAverages)
	% Generate cell array of spin glass models with grid size equals to n x m
	% and all pairwise weights fixed to J.
	% For each temperature numAverages models with different random
	% unary weights are generted (each time from the uniform distribution on [-1, 1]).
	% 
	% modelCellArray is a cell array of size numAverages by length(temperatureArray).
	% modelCellArray{i, j} has temperature equals to temperatureArray(j).
	% 

	numTemerature = length(temperatureArray);
	modelCellArray = cell(numAverages, numTemerature);
	for iAverage = 1 : numAverages
		unaryType = 'rand';
		for iTemerature = 1 : numTemerature
			currTemperature = temperatureArray(iTemerature);
			currModel = generate_spin_glass_model(n, m, currTemperature, 'J', J, 'H', unaryType);
			modelCellArray{iAverage, iTemerature} = currModel;
		end
	end
end
