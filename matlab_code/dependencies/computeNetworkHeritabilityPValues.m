function [p_values, actual_diffs,actual_diffs_perc] = computeNetworkHeritabilityPValues(movieNetHerit, restNetHerit, numPermutations)

numNetworks = 17;
numDays = 2;
p_values = zeros(numNetworks, numDays);
actual_diffs = zeros(numNetworks, numDays);

for day = 1:numDays
    actual_diff = mean(movieNetHerit(:,:,day), 2) - mean(restNetHerit(:,:,day), 2);
    actual_diffs(:, day) = actual_diff;
    actual_diffs_perc(:, day) = actual_diffs(:, day)./mean(restNetHerit(:,:,day), 2);
    perm_diffs = zeros(numNetworks, numPermutations);

    for perm = 1:numPermutations
        rng(perm);
        for net = 1:numNetworks
            % Randomly shuffle heritability values for each network
            combined = [movieNetHerit(net,:,day), restNetHerit(net,:,day)];
            shuffled = combined(randperm(numel(combined)));

            % Split shuffled data back into 'movie' and 'rest'
            shuffled_movie = shuffled(1:numNetworks);
            shuffled_rest = shuffled(numNetworks+1:end);

            % Compute difference in average heritability after permutation
            perm_diffs(net, perm) = mean(shuffled_movie) - mean(shuffled_rest);
        end
    end

    % Calculate p-values
    for net = 1:numNetworks
        p_values(net, day) = mean(abs(perm_diffs(net, :)) >= abs(actual_diff(net)));
    end
end
end
