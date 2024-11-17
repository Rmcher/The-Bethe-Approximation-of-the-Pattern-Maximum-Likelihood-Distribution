function [permB, Vl, Vr, Zf_left, Zf_right, MessageCellLeft, MessageCellRight, Ze]  = computeBethePermanent_block_log(A, n)
    % Author: WU Binghong
    % Date: 2024.Oct.31
    
    %--------------------------------------------------------------------

    % Get the non-repeated cols' index
    [~, uniqueCols, ~] = unique(A', 'rows');
    uniqueCols = sort(uniqueCols, 'ascend');
    segmentSizes = [builtin('diff', uniqueCols); n - uniqueCols(end) + 1];

    % Get the non-repeated rows' index
    % [~, uniqueRows, ~] = unique(A, 'rows');
    % uniqueRows = sort(uniqueRows, 'ascend');
    % rowSegmentSizes = [diff(uniqueRows); n - uniqueRows(end) + 1];

    % Initialize the n*n message cell and V symbol matrix, cell{i, j} 
    % contain the message on x(i, j), in the form of vector = [u(0), u(1)]T
    MessageCellRight = cell(n, n);
    MessageCellLeft = cell(n, n);
    for i = 1:n
        for j = 1:n
            MessageCellLeft{i, j} = [0; 0];
            MessageCellRight{i, j} = [0; 0];
        end
    end

    VRight = zeros(n, n);
    VLeft = zeros(n, n);
    
    % Store values for plotting
    Vl = [];
    Vr = [];

    

    % Random Initialize the message on each edge for t=0, sending right,
    % i.e. x, ensuring sum=1. For u_{ij}(x_ij) = 1, MessageCell{i, j}(2)
    % In the block SPA, we need randomly generate the message, but all the
    % messages on different edges should be the same.

    % Fix seed for debugging
    rng(42);

    % r = rand();
    r = 1/3;
    for i = 1:n
        for u = 1:length(uniqueCols)
            j = uniqueCols(u);
            vector = [r; 1 - r];
            MessageCellRight{i, j} = vector;
            % VRight{1,1}(i, j) = r/(1-r);
        end
    end 
    
    %--------------------------------------------------------------------
    % Message sending algorithm in iteration
    
    % Iteration setup
    t = 0;  % Initial time step
    tolerance = 1e-4;  % Stopping condition tolerance for VLeft comparison

    while true
        t = t + 1;

        if mod(t, 2) == 1  % t is odd: update MessageCellLeft and VLeft
            % Update the MessageCellLeft and VLeft based on current MessageCellRight
            

            % modify the following, let calculation happens only when j
            % equals uniqueCols and do repeatition operation for the
            % gap.

            % Step 1: Only update the special cols in uniqueCols.
            for i = 1:n
                for u = 1:length(uniqueCols)
                    j = uniqueCols(u);

                    sumTerm = 0;
                    for i_prime = 1:n
                        if i_prime ~= i && A(i_prime, j) ~= 0
                            productTerm = 1;
                            for i_double_prime = 1:n
                                if i_double_prime ~= i && i_double_prime ~= i_prime && A(i_double_prime, j) ~= 0
                                    productTerm = productTerm * MessageCellRight{i_double_prime, j}(1);
                                end
                            end
                            productTerm = productTerm * MessageCellRight{i_prime, j}(2) * sqrt(A(i_prime, j));
                            sumTerm = sumTerm + productTerm;
                        end
                    end
                    MessageCellLeft{i, j}(1) = sumTerm;

                    % Optimized Calculation for MessageCellRight{i, j}(2)
                    productTerm = sqrt(A(i, j));
                    for i_prime = 1:n
                        if i_prime ~= i && A(i_prime, j) ~= 0
                            productTerm = productTerm * MessageCellRight{i_prime, j}(1);
                        end
                    end

                    MessageCellLeft{i, j}(2) = productTerm;
                    
                    total = MessageCellLeft{i, j}(1) + MessageCellLeft{i, j}(2);
                    MessageCellLeft{i, j} = MessageCellLeft{i, j} / total;
                    
                    VLeft(i, j) = MessageCellLeft{i, j}(1) / MessageCellLeft{i, j}(2);

                    fprintf('Computed normal t=%d, i=%d, j=%d.\n', t, i, j);
                end
            end

            % Step 2: Duplicate the special cols to the following cols until the next special col.
            % refCol = uniqueCols(1);
            % 
            % for j = 1:n
            %     if ~ismember(j, uniqueCols)
            %         % If j is not in uniqueCols, make duplication from refCol
            %         MessageCellLeft(:, j) = MessageCellLeft(:, refCol);
            %         VLeft{(t + 1) / 2, 1}(:, j) = VLeft{(t + 1) / 2, 1}(:, refCol);
            %         fprintf('Copied entire column t=%d, j=%d from refCol=%d.\n', t, j, refCol);
            %     else
            %         % If j is in uniqueCols, update anchor refCol
            %         refCol = j;
            %         fprintf('Set refCol to j=%d.\n', j);
            %     end
            % end

            % Store V for plot
            Vl = [Vl; VLeft(:)'];
            
        else  % t is even: update MessageCellRight and VRight
            % Update the MessageCellRight and VRight based on current MessageCellLeft
            for u = 1:length(uniqueCols)
                j = uniqueCols(u);
                segSize = segmentSizes(u);  % 当前列块的大小

                for i = 1:n
                    % Initialize MessageCellRight{i, j}(1)
                    sumTerm = 0;

                    % 遍历 j_prime
                    for v = 1:length(uniqueCols)
                        j_prime = uniqueCols(v);
                        factor = segmentSizes(v);

                        % 跳过 j_prime = j 的情况，如果 j 是唯一列
                        if j_prime == j && segSize == 1
                            continue;  % 唯一列情况，跳过 j_prime = j
                        end

                        % 计算 productTerm，包含 MessageCellLeft{i, j_prime}(2) 和 sqrt(A(i, j_prime))
                        productTerm = MessageCellLeft{i, j_prime}(2) * sqrt(A(i, j_prime));

                        % 计算分段乘积，根据分段数量优化乘法操作
                        for w = 1:length(segmentSizes)
                            refCol = uniqueCols(w);
                            refSize = segmentSizes(w);

                            % 如果 refCol 是 j 或 j_prime，需要减去占位列
                            effectiveSize = refSize - (refCol == j) - (refCol == j_prime);

                            % 更新 productTerm，通过有效列数乘以对应的 MessageCellLeft 值
                            if effectiveSize > 0
                                productTerm = productTerm * (MessageCellLeft{i, refCol}(1) ^ effectiveSize);
                            end
                        end

                        % 将 productTerm 加到 sumTerm
                        sumTerm = sumTerm + (factor - (j == j_prime)) * productTerm;
                    end

                    % 更新 MessageCellRight{i, j}(1)
                    MessageCellRight{i, j}(1) = sumTerm;

                    % Calculate MessageCellRight{i,j}(2)
                    productTerm = sqrt(A(i, j));
                    for q = 1:length(uniqueCols)
                        j_prime = uniqueCols(q);
                        factor = segmentSizes(q);

                        % 跳过 j_prime = j 的情况，如果 j 是唯一列
                        if j_prime == j && factor == 1
                            continue;  % 唯一列情况，跳过 j_prime = j
                        end

                        effectiveFactor = factor - (j_prime == j); % 如果 j_prime 和 j 相等，则减去 1

                        % 如果 effectiveFactor 大于 0，则计算并累乘
                        if effectiveFactor > 0 && A(i, j_prime) ~= 0
                            productTerm = productTerm * (MessageCellLeft{i, j_prime}(1) ^ effectiveFactor);
                        end
                    end
                    
                    MessageCellRight{i, j}(2) = productTerm;

                    % Normalize MessageCellRight{i, j}
                    total = MessageCellRight{i, j}(1) + MessageCellRight{i, j}(2);
                    MessageCellRight{i, j} = MessageCellRight{i, j} / total;

                    % Update VRight
                    VRight(i, j) = MessageCellRight{i, j}(1) / MessageCellRight{i, j}(2);

                    fprintf('Optimized calculation done for t=%d, i=%d, j=%d.\n', t, i, j);
                end
            end

            % Store V for plot
            Vr = [Vr; VRight(:)'];

        end

        %----------------------------------------------------------------
        % Check stopping condition
        if t > 2 && mod(t, 2) == 1  % Only check VLeft in odd iterations
            difference = max(max(abs(Vl(end, :) - Vl(end-1, :))));  % Difference between VLeft at t and t-2
            if difference < tolerance || t > 300000
                % disp(['Iteration converged at t = ', num2str(t)]);
                break;  % Stop if the change is smaller than tolerance
            end
        end
    end

    %--------------------------------------------------------------------
    % Calculate permB = Zf/Ze using optimized segment calculations with adjusted factors
    log_Zf_left = 0;  % Initialize Zf_left
    log_Zf_right = 0; % Initialize Zf_right

    % Calculate log_Zf_left
    for i = 1:n
        sum_log_i = -Inf;  % Initialize for log-sum-exp accumulation

        for u = 1:length(uniqueCols)
            j = uniqueCols(u);
            segSize = segmentSizes(u);  % Size of the repeated block for column j

            % Calculate the primary term for this unique column in log-space
            log_term = log(sqrt(A(i, j))) + log(MessageCellLeft{i, j}(2));

            % Calculate the product for other columns in this row in log-space
            log_product_j_prime = 0;
            for v = 1:length(uniqueCols)
                j_prime = uniqueCols(v);
                factor = segmentSizes(v) - (j == j_prime);  % Adjust factor if j == j_prime

                % Skip if j_prime matches j and it's a single column
                if j == j_prime && segSize == 1
                    continue;
                end

                % Skip if factor becomes zero (i.e., if j == j_prime and not a block)
                if factor > 0
                    log_product_j_prime = log_product_j_prime + factor * log(MessageCellLeft{i, j_prime}(1));
                end
            end

            % Add the product term to the main term in log-space
            log_term = log_term + log_product_j_prime;

            % Accumulate in sum_log_i using manual log-sum-exp
            if sum_log_i == -Inf
                sum_log_i = log(segSize) + log_term;  % Initialization case
            else
                max_log = max(sum_log_i, log(segSize) + log_term);
                sum_log_i = max_log + log(1 + exp(-abs(sum_log_i - (log(segSize) + log_term))));
            end
        end

        % Add the row log-sum to log_Zf_left
        log_Zf_left = log_Zf_left + sum_log_i;
    end
    Zf_left = exp(log_Zf_left);  % Convert from log-space if needed

    % Zf_right
    for u = 1:length(uniqueCols)
        j = uniqueCols(u);
        segSize = segmentSizes(u);

        sum_log_j = -Inf;
        for i = 1:n
            log_term = log(sqrt(A(i, j))) + log(MessageCellRight{i, j}(2));
            log_product_i_prime = 0;
            for i_prime = 1:n
                if i_prime ~= i
                    log_product_i_prime = log_product_i_prime + log(MessageCellRight{i_prime, j}(1));
                end
            end
            log_term = log_term + log_product_i_prime;
            if sum_log_j == -Inf
                sum_log_j = log_term;  % Initialization case
            else
                max_log = max(sum_log_j, log_term);
                sum_log_j = max_log + log(1 + exp(-abs(sum_log_j - log_term)));
            end
        end
        log_Zf_right = log_Zf_right + segSize * sum_log_j;
    end
    Zf_right = exp(log_Zf_right);

    % Calculate Ze by iterating only over uniqueCols and segments
    log_Ze = 0;
    for u = 1:length(uniqueCols)
        j = uniqueCols(u);
        segSize = segmentSizes(u);

        % Calculate the product for the current unique column j
        for i = 1:n
            % Compute the inner term for each (i, j)
            innerTerm = MessageCellRight{i, j}(1) * MessageCellLeft{i, j}(1) + MessageCellRight{i, j}(2) * MessageCellLeft{i, j}(2);

            % Ensure innerTerm is positive before taking log
            if innerTerm > 0
                log_innerTerm = log(innerTerm);
                % Accumulate the log value for Ze with segment size
                log_Ze = log_Ze + segSize * log_innerTerm;
            else
                log_Ze = -Inf;  % If innerTerm is 0, set log_Ze to -Inf to handle log of zero
                break;
            end
        end
    end

    Ze = exp(log_Ze);
    % Final calculation of Zf
    log_Zf = log_Zf_left + log_Zf_right;
    % permB = Zf / Ze;
    permB = exp(log_Zf - log_Ze);
    return;
end