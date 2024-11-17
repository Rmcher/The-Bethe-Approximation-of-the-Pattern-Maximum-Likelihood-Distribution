function pr = randomGenerateP(k, R)
    uniqueRows = [1; sort(randperm(k-1, R-1)'+1)];
    rowSegmentSizes = diff([uniqueRows; k + 1]);

    random_probs = rand(R, 1);
    random_probs = sort(random_probs, 'descend');
    random_probs = random_probs / sum(random_probs);

    pr = zeros(k, 1);

    for i = 1:length(rowSegmentSizes)
        segment_size = rowSegmentSizes(i);
        segment_prob = random_probs(i);
        pr(uniqueRows(i):uniqueRows(i) + segment_size - 1) = segment_prob;
    end

    pr = pr / sum(pr);
end
