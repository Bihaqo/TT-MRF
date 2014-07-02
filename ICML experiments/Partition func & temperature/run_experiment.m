% Reproduce the partition function experiment from the ICML paper.
% 
% This script generates an Ising models dataset where temperature varies from 10^-1 to 10^3,
% computes the partition function with different approaches and plot the comparison.
% 
% You can access original dataset used in the ICML paper here:
% https://www.dropbox.com/s/6onvwjz4avva72m/data.mat
% 

ttRoundingPrecision = 1e-10;
n = 10;
m = 10;
temperatureArray = 10.^[-1:0.1:3];
J = 1;
numAverages = 50;

% Fix seed for reproducibility.
rng(0);

modelCellArray = generate_data(n, m, temperatureArray, J, numAverages);
save('data.mat', 'modelCellArray', 'n', 'm', 'temperatureArray', 'J', 'numAverages');

[logZ_TT, logAbsErrorBound_TT, time_TT, logZ_JT, ...
		  logZ_BP, time_BP, logZ_MF, time_MF, ...
		  logZ_TREEEP, time_TREEEP, logZ_AIS, time_AIS] = compute_answers(modelCellArray, ttRoundingPrecision);

save('results.mat', 'temperatureArray', 'logZ_TT', 'logAbsErrorBound_TT', 'time_TT', 'logZ_JT', 'logZ_BP', 'time_BP', 'logZ_MF', 'time_MF', 'logZ_TREEEP', 'time_TREEEP', 'logZ_AIS', 'time_AIS');


plot_comparison(temperatureArray, logZ_TT, logZ_JT, logZ_BP, logZ_MF, logZ_TREEEP, logZ_AIS);
