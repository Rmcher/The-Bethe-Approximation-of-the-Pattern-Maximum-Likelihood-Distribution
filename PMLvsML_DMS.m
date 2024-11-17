% Author: WU Binghong
% Date: 2024.Nov.17

% Sec 1: Define the Alphabet of DMS and message.
% Alphabet
k = 500;
% Message length
n = 750;

% Sec 2: Random distribution p ~ exp(1), sorting and normalization. And random pr, sorting and normalization.
% True distribution
lambda = 1;
p = -log(rand(k, 1)) / lambda;  % using rand() to produce exp(1) distribution.
p = sort(p, 'descend');
p = p / sum(p);

% Self define special distribution
R = 250; % Number of guess p distribution segment, no more than k.
pr = randomGenerateP(k, R);
prchange = [];
prchange = [prchange, pr];

% Sec 3: Generate X with length n with p.
X = datasample(1:k, n, 'Weights', p);
X = sort(X, 'descend');

% Sec 4: Generate mu vector, mu_i representing the amount of x_i.
mu = zeros(k, 1);
m = 0;
for i = 1:k
    mu(i) = sum(X == i);
    if mu(i) ~= 0
        m = m+1;
    end
end
mu_len = length(mu);

mu = sort(mu, 'descend');

% Sec 5：Generate Vandermonde Matrix, real one and guess one.
groundtruthmatrix = zeros(k, k);
for i = 1:k
    for j = 1:mu_len
        groundtruthmatrix(i, j) = p(i) ^ mu(j);
    end
end

% Initial the 1st random guess of matrix with distribution pr
PMLmatrix = zeros(k, k);
for i = 1:k
    for j = 1:mu_len
        PMLmatrix(i, j) = pr(i) ^ mu(j);
    end
end

% get the non-repeated cols' index
[~, uniqueCols, ~] = unique(PMLmatrix', 'rows');
uniqueCols = sort(uniqueCols, 'ascend');
colSegmentSizes = [diff(uniqueCols); k - uniqueCols(end) + 1];

[~, uniqueRows, ~] = unique(PMLmatrix, 'rows');
uniqueRows = sort(uniqueRows, 'ascend');
rowSegmentSizes = [diff(uniqueRows); k - uniqueRows(end) + 1];

% Sec 6: Compute the perm (permB).
gammamatrix = zeros(k,k);
tolerance = 1e-4;

for t = 1:1000

    % For every iteration, update new guess matrix with updated pr
    % PMLmatrixnew = zeros(k, k);
    for i = 1:k
        for j = 1:mu_len
            PMLmatrixnew(i, j) = pr(i) ^ mu(j);
        end
    end

    [permB, Lambdal, Lambdar, Zf_left, Zf_right, MessageCellLeft, MessageCellRight, Ze] = computeBethePermanent_double_block_log_all(PMLmatrixnew, k);

    LambdaLeft = reshape(Lambdal(end, :), [k, k]);
    LambdaRight = reshape(Lambdar(end, :), [k, k]);

    for uRowIdx = 1:length(uniqueRows)
        for uColIdx = 1:length(uniqueCols)
            i = uniqueRows(uRowIdx);
            j = uniqueCols(uColIdx);

            % gammamatrix(i, j) = (VLeft(i, j) * VRight(i, j)) / (VLeft(i, j) * VRight(i, j) + 1);
            gammamatrix(i, j) = (1) / (LambdaLeft(i, j) * LambdaRight(i, j) + 1);

            rowRange = i:(i + rowSegmentSizes(uRowIdx) - 1);
            colRange = j:(j + colSegmentSizes(uColIdx) - 1);
            gammamatrix(rowRange, colRange) = gammamatrix(i, j);
        end
    end

    % Step 2: Update pr
    pr = gammamatrix * mu / n;
    prchange = [prchange, pr];

    % step 3: Check iteration stop condition
    difference = max(max(abs(prchange(:, end) - prchange(:, end-1))));  % Difference between VLeft at t and t-2
    if difference < tolerance
        fprintf('Requirement met at Iteration t = %d\n', t);
        break;
    end

    fprintf('Iteration processed at t = %d\n', t);

end

% Sec 7: ML estimation for source distribution
% Estimate the source distribution directly using Maximum Likelihood ML estimation based on symbol frequency
pr_ml = mu / sum(mu);

zero_index = find(pr_ml == 0, 1);

% Sex 8: Combine PML and ML results in a single plot for comparison
figure;

% Plot Ground Truth, BPML, and ML
semilogy(1:k, p, '-o', 'LineWidth', 1.5);
hold on;
semilogy(1:k, pr, '-x', 'LineWidth', 1.5);
h_ml = semilogy(1:k, pr_ml, '-s', 'LineWidth', 1.5);  % Store handle for ML plot

% Add vertical line at the transition point to zero for ML
if ~isempty(zero_index)
    plot([zero_index, zero_index], [pr_ml(zero_index - 1), min(p)], '-s', 'Color', get(h_ml, 'Color'), 'LineWidth', 1.5);
end

% Set title, labels, and legend
title('Comparison of Ground truth, BPML, and ML Distributions');
xlabel('Index');
ylabel('Probability (log scale)');

legend('Ground truth', 'BPML', 'ML');
grid on;

% Set y-axis lower limit to the minimum value of p
ylim([min(p), 1]);

hold off;

filename_combined = sprintf('./results/PML_ML_k_%d_n_%d_t_%d_R_%d.fig', k, n, t, R);
saveas(gcf, filename_combined);
