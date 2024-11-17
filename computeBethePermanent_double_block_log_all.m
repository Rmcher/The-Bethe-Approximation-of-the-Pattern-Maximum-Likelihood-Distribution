function [permB, Lambdal, Lambdar, Zf_left, Zf_right, MessageCellLeft, MessageCellRight, Ze] = computeBethePermanent_double_block_log_all(A, n)

    [~, uniqueCols, ~] = unique(A', 'rows');
    uniqueCols = sort(uniqueCols, 'ascend');
    colSegmentSizes = [diff(uniqueCols); n - uniqueCols(end) + 1];
    
    [~, uniqueRows, ~] = unique(A, 'rows');
    uniqueRows = sort(uniqueRows, 'ascend');
    rowSegmentSizes = [diff(uniqueRows); n - uniqueRows(end) + 1];

    MessageCellRight = cell(n, n);
    MessageCellLeft = cell(n, n);
    for i = 1:n
        for j = 1:n
            MessageCellLeft{i, j} = [-Inf; -Inf];
            MessageCellRight{i, j} = [-Inf; -Inf];
        end
    end

    LambdaRight = zeros(n, n);
    LambdaLeft = zeros(n, n);
    Lambdal = [];
    Lambdar = [];

    r = log(1/3); 
    for i = uniqueRows'
        for j = uniqueCols'
            vector = [r; log(1 - exp(r))];
            MessageCellRight{i, j} = vector;
        end
    end 

    t = 0;
    tolerance = 1e-4;

    while true
        t = t + 1;

        if mod(t, 2) == 1
            for u = 1:length(uniqueRows)
                i = uniqueRows(u);
                segSize = rowSegmentSizes(u);

                for j = uniqueCols'
                    sumTerm = -Inf;
                    for v = 1:length(uniqueRows)
                        i_prime = uniqueRows(v);
                        factor = rowSegmentSizes(v);

                        if i_prime == i && segSize == 1
                            continue;
                        end

                        productTerm = MessageCellRight{i_prime, j}(2) + 0.5 * log(A(i_prime, j));

                        for w = 1:length(rowSegmentSizes)
                            refrow = uniqueRows(w);
                            refSize = rowSegmentSizes(w);

                            effectiveSize = refSize - (refrow == i) - (refrow == i_prime);

                            if effectiveSize > 0
                                productTerm = productTerm + effectiveSize * MessageCellRight{refrow, j}(1);
                            end
                        end

                        temp_sum = log(factor - (i == i_prime)) + productTerm;
                        sumTerm = logsumexp(sumTerm, temp_sum);
                    end

                    MessageCellLeft{i, j}(1) = sumTerm;

                    productTerm = 0.5 * log(A(i, j));
                    for q = 1:length(uniqueRows)
                        i_prime = uniqueRows(q);
                        factor = rowSegmentSizes(q);

                        if i_prime == i && factor == 1
                            continue;
                        end

                        effectiveFactor = factor - (i_prime == i);

                        if effectiveFactor > 0 && A(i_prime, j) ~= 0
                            productTerm = productTerm + effectiveFactor * MessageCellRight{i_prime, j}(1);
                        end
                    end

                    MessageCellLeft{i, j}(2) = productTerm;

                    total = logsumexp(MessageCellLeft{i, j}(1), MessageCellLeft{i, j}(2));
                    MessageCellLeft{i, j} = MessageCellLeft{i, j} - total;

                    LambdaLeft(i, j) = exp(MessageCellLeft{i, j}(1) - MessageCellLeft{i, j}(2));
                end
            end

            Lambdal = [Lambdal; LambdaLeft(:)'];

        else
            for u = 1:length(uniqueCols)
                j = uniqueCols(u);
                segSize = colSegmentSizes(u);

                for i = uniqueRows'
                    sumTerm = -Inf;

                    for v = 1:length(uniqueCols)
                        j_prime = uniqueCols(v);
                        factor = colSegmentSizes(v);

                        if j_prime == j && segSize == 1
                            continue;
                        end

                        productTerm = MessageCellLeft{i, j_prime}(2) + 0.5 * log(A(i, j_prime));

                        for w = 1:length(colSegmentSizes)
                            refCol = uniqueCols(w);
                            refSize = colSegmentSizes(w);

                            effectiveSize = refSize - (refCol == j) - (refCol == j_prime);

                            if effectiveSize > 0
                                productTerm = productTerm + effectiveSize * MessageCellLeft{i, refCol}(1);
                            end
                        end

                        temp_sum = log(factor - (j == j_prime)) + productTerm;
                        sumTerm = logsumexp(sumTerm, temp_sum);
                    end

                    MessageCellRight{i, j}(1) = sumTerm;

                    productTerm = 0.5 * log(A(i, j));
                    for q = 1:length(uniqueCols)
                        j_prime = uniqueCols(q);
                        factor = colSegmentSizes(q);

                        if j_prime == j && factor == 1
                            continue;
                        end

                        effectiveFactor = factor - (j_prime == j);

                        if effectiveFactor > 0 && A(i, j_prime) ~= 0
                            productTerm = productTerm + effectiveFactor * MessageCellLeft{i, j_prime}(1);
                        end
                    end

                    MessageCellRight{i, j}(2) = productTerm;

                    total = logsumexp(MessageCellRight{i, j}(1), MessageCellRight{i, j}(2));
                    MessageCellRight{i, j} = MessageCellRight{i, j} - total;

                    LambdaRight(i, j) = exp(MessageCellRight{i, j}(1) - MessageCellRight{i, j}(2));
                end
            end

            Lambdar = [Lambdar; LambdaRight(:)'];
        end

        if t > 2 && mod(t, 2) == 1
            difference = max(max(abs(Lambdal(end, :) - Lambdal(end-1, :))));
            if difference < tolerance || t > 300000
                break;
            end
        end
    end

    log_Zf_left = 0;
    for x = 1:length(uniqueRows)
        i = uniqueRows(x);
        rowsegSize = rowSegmentSizes(x);
        sum_log_i = -Inf;

        for u = 1:length(uniqueCols)
            j = uniqueCols(u);
            segSize = colSegmentSizes(u);

            log_term = 0.5 * log(A(i, j)) + MessageCellLeft{i, j}(2);

            log_product_j_prime = 0;
            for v = 1:length(uniqueCols)
                j_prime = uniqueCols(v);
                factor = colSegmentSizes(v) - (j == j_prime);

                if j == j_prime && segSize == 1
                    continue;
                end

                if factor > 0
                    log_product_j_prime = log_product_j_prime + factor * MessageCellLeft{i, j_prime}(1);
                end
            end

            log_term = log_term + log_product_j_prime;

            if sum_log_i == -Inf
                sum_log_i = log(segSize) + log_term;
            else
                max_log = max(sum_log_i, log(segSize) + log_term);
                sum_log_i = max_log + log(1 + exp(-abs(sum_log_i - (log(segSize) + log_term))));
            end
        end

        log_Zf_left = log_Zf_left + rowsegSize * sum_log_i;
    end
    Zf_left = exp(log_Zf_left);

    log_Zf_right = 0;
    for x = 1:length(uniqueCols)
        j = uniqueCols(x);
        colsegSize = colSegmentSizes(x);
        sum_log_j = -Inf;

        for u = 1:length(uniqueRows)
            i = uniqueRows(u);
            segSize = rowSegmentSizes(u);

            log_term = 0.5 * log(A(i, j)) + MessageCellRight{i, j}(2);
            
            log_product_i_prime = 0;
            for v = 1:length(uniqueRows)
                i_prime = uniqueRows(v);
                factor = rowSegmentSizes(v) - (i == i_prime);

                if i == i_prime && segSize == 1
                    continue;
                end

                if factor > 0
                    log_product_i_prime = log_product_i_prime + factor * MessageCellRight{i_prime, j}(1);
                end
            end

            log_term = log_term + log_product_i_prime;

            if sum_log_j == -Inf
                sum_log_j = log(segSize) + log_term;
            else
                max_log = max(sum_log_j, log(segSize) + log_term);
                sum_log_j = max_log + log(1 + exp(-abs(sum_log_j - (log(segSize) + log_term))));
            end
        end

        log_Zf_right = log_Zf_right + colsegSize * sum_log_j;
    end
        Zf_right = exp(log_Zf_right);

    log_Ze = 0;
    for u = 1:length(uniqueCols)
        j = uniqueCols(u);
        colsegSize = colSegmentSizes(u);

        for v = 1:length(uniqueRows)
            i = uniqueRows(v);
            rowsegSize = rowSegmentSizes(v);

            innerTerm = logsumexp(MessageCellRight{i, j}(1) + MessageCellLeft{i, j}(1), ...
                                  MessageCellRight{i, j}(2) + MessageCellLeft{i, j}(2));

            if innerTerm > -Inf
                log_Ze = log_Ze + rowsegSize * colsegSize * innerTerm;
            else
                log_Ze = -Inf;
                break;
            end
        end
    end
    Ze = exp(log_Ze);

    log_Zf = log_Zf_left + log_Zf_right;
    permB = exp(log_Zf - log_Ze);

end

function log_sum = logsumexp(log_a, log_b)
    if log_a == -Inf
        log_sum = log_b;
    elseif log_b == -Inf
        log_sum = log_a;
    else
        max_log = max(log_a, log_b);
        log_sum = max_log + log(1 + exp(-abs(log_a - log_b)));
    end
end

