function [fc_herit, fc_herit_perm,fc_herit_jack,fc_similarity] = computeNetworkHeritability(fc_square, net_ids, subjects,kinship,covari,num_perm,fam_ids)

num_nets            = max(net_ids);
n2                  = num_nets + ((num_nets * (num_nets - 1)) / 2);
numsubjects         = length(subjects);
fc_similarity       = zeros(numsubjects, numsubjects, n2, 2);
combination_counter = 0;
hvec                = @(x) reshape(x(triu(true(size(x(:,:,1))), 1)), [], 1);
unique_networks     = unique(net_ids);

% Preallocate fc_vectors
max_elements = numsubjects * max(cellfun(@numel, arrayfun(@(x) find(net_ids == x), unique_networks, 'UniformOutput', false)));
fc_vectors   = NaN(max_elements, numsubjects); 

% Loop over each scan
for scan = 1:size(fc_square, 4)
    % Loop over each network combination
    for net1 = 1:length(unique_networks)
        for net2 = net1:length(unique_networks)
            combination_counter = combination_counter + 1;

            indices_net1 = find(net_ids == unique_networks(net1));
            indices_net2 = find(net_ids == unique_networks(net2));
            fc_submatrix = reshape(fc_square(indices_net1, indices_net2, :, scan), [], numsubjects);
            
            % Store in fc_vectors
            num_elements                  = size(fc_submatrix, 1);
            fc_vectors(1:num_elements, :) = fc_submatrix;
            valid_indices                 = ~any(isnan(fc_vectors), 2);
            fc_vectors_truncated          = fc_vectors(valid_indices, :);
            corr_matrix                   = corr(fc_vectors_truncated, 'rows', 'pairwise');

            % Assign the correlation matrix to the corresponding slice of the fc_similarity matrix
            fc_similarity(:, :, combination_counter, scan) = corr_matrix;
        end
    end
    combination_counter = 0;
end

fc_herit      = zeros(num_nets, size(fc_square,4));
fc_herit_perm = zeros(num_nets, size(fc_square,4));
fc_herit_jack = zeros(num_nets,length(unique(fam_ids)), size(fc_square,4));

for scan_id = 1:size(fc_square,4)
    parfor net = 1:n2
        [fc_herit(net,scan_id), fc_herit_perm(net,scan_id),~, fc_herit_jack(net,:,scan_id)] = h2_mat(squeeze(fc_similarity(:,:,net,scan_id)), kinship(subjects,subjects), covari(subjects,[1,2,scan_id+2]), num_perm,fam_ids);
    end
end

end
