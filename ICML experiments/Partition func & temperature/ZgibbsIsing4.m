function logZ = ZgibbsIsing4(H, J, beta, num_pdfs, num_iter)
% Estimate the partition function with Gibbs sampling.
%
% logZ = ZgibbsIsing4(H, J, beta, num_pdfs, num_iter)
%
% Input
%   H -- N x M matrix, external magnetic field
%   J -- number, global pairwise wight;
%   beta -- 1 x num_temp vector, 1/(kT) values;
%   num_pdfs -- integer, number of intermediate distributions;
%   num_iter -- integer, number of iterations;
%
% Output
%   logZ -- 1 x num_temp vector, logarithms of the partition function for all provided beta values;
% 
% 
% Author: Dmitry Kropotov (kropotov -at- bayesgroup.ru)
% 

if nargin ~= 5
    error('Wrong number of input arguments.');
end

if ~isnumeric(H) || (length(size(H))>2)
    error('External magnetic field parameter H should be a matrix.');
end

if ~isnumeric(J) || (length(size(J))>2) || (norm(size(J) - [1 1])~=0)
    error('Pairwise weight J should be a number.');
end

if ~isnumeric(beta) || (length(size(beta))>2) || (size(beta,1)>1) || any(beta<=0)
    error('Beta parameter should be a vector of real values.');
end

if ~isnumeric(num_pdfs) || (length(size(num_pdfs))>2) || (norm(size(num_pdfs)-[1 1])~=0) || (fix(num_pdfs)~=num_pdfs) || (num_pdfs < 1)
    error('Number of intermediate distributions is not a natural number.');
end

if ~isnumeric(num_iter) || (length(size(num_iter))>2) || (norm(size(num_iter)-[1 1])~=0) || (fix(num_iter)~=num_iter) || (num_iter<1)
    error('Number of iterations is not a natural number.');
end


rng('Shuffle');

t_start = tic;

[N,M] = size(H);
NM = N*M;
num_beta = length(beta);

% 1 - outside element
index = reshape(2 : (NM + 1), [N, M]);

numN = 4;
neighbors = ones(N, M, numN);
neighbors(1 : end - 1, :, 1) = index(2 : end, :);
neighbors(2 : end, :, 2) = index(1 : end - 1, :);
neighbors(:, 1 : end - 1, 3) = index(:, 2 : end);
neighbors(:, 2 : end, 4) = index(:, 1 : end - 1);

nN = ones(NM + 1, numN);
for c = 1 : numN
    nN(2 : end, c) = reshape(neighbors(:, :, c), [NM, 1]);
end

hN = [0; H(:)];

alpha = 0:1/(num_pdfs-1):1;

w = zeros(num_iter, num_beta);
fprintf('Total %d iterations:', num_iter);
for iter = 1:num_iter

    if mod(iter, 30) == 0
        fprintf('\n');
    end
    fprintf(' %d', iter);
    
    % Generate points from uniform distrition.
    currX = round(rand(NM + 1, num_beta)) * 2 - 1;
    currX(1, :) = 0;

    energy = zeros(num_pdfs-1, num_beta);
    
    % Compute energies of generated configurations.
    for c = 1:numN
        energy(1,:) = energy(1,:) - sum(currX .* currX(nN(:,c),:), 1);
    end
    energy(1,:) = J*energy(1,:)/2 - sum(currX .* repmat(hN, [1, num_beta]), 1);

    for i_alpha = 2:length(alpha)-1

        % Do one step of Gibbs sampling for p_alpha distribution.
        for i = 2 : NM + 1
            p1 = 1 ./ (1 + exp(-2 * alpha(i_alpha)*beta .*(J * sum(currX(nN(i, :), :), 1) + hN(i))));
            currX(i, :) = 1 - 2 * (rand(1, num_beta) > p1);
        end
        % Compute energy of current configuration.
        for c = 1 : numN
            energy(i_alpha, :) = energy(i_alpha, :) - sum(currX .* currX(nN(:, c), :), 1);
        end
        energy(i_alpha, :) = energy(i_alpha, :) * J / 2 - sum(currX .* repmat(hN, [1, num_beta]), 1);
    end
    w(iter,:) = -beta.*sum(energy, 1)/(num_pdfs-1);
end
fprintf('\n');

wmax = max(w, [], 1);
logZ = log(sum(exp(bsxfun(@minus, w, wmax)),1)) + wmax - log(num_iter) + NM*log(2);

t_finish = toc(t_start);
fprintf('Total time = %.4f\n', t_finish);
