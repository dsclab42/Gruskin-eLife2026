function fcMagnitude = computeNetworkMagnitude(fc_square, net_ids, subjects, task_type, align_type, curr_parc)
% computeNetworkMagnitude  Per-subject FC "strength" (magnitude) for each network combination.
%
%   RECONSTRUCTED FILE — the original was missing from the release. Please
%   review against your archived outputs before relying on it. The two modeling
%   choices most worth checking are flagged with "ASSUMPTION" comments below.
%
%   For each day/scan and each of the 153 unique 17-network combinations
%   (within- and between-network), this averages the parcel-level FC values
%   belonging to that combination, per subject, to yield a single FC strength.
%   These per-subject strengths are the phenotypes fed to SOLAR: the function
%   writes one CSV per scan, named
%       fc_net_<task_type>_<align_type>_parc_<curr_parc>_scan_<scan_id>.csv
%   as an [nCombinations x nSubjects] headerless matrix. r_code/calc_h2.R reads
%   these with header = FALSE and transposes them (so its rows become subjects
%   and its columns become the network combinations "net1".."net153").
%
%   Inputs:
%     fc_square  - [nParcels x nParcels x nSubjects x nDays] square FC matrices
%                  (as produced by computeFunctionalConnectivity)
%     net_ids    - [nParcels x 1] network assignment per parcel (1..17), from
%                  processNetworkData
%     subjects   - subject index vector (used only for count/order; length must
%                  equal size(fc_square,3))
%     task_type  - 'movie' | 'rest'   (filename field)
%     align_type - 'anatomical' | 'piecewise' | 'connectivity' (filename field)
%     curr_parc  - parcellation resolution as a char row vector, e.g. '400'
%
%   Output:
%     fcMagnitude - [nCombinations x nSubjects x nDays] FC strengths (Pearson r),
%                   in the same combination order used by computeNetworkHeritability
%                   and vector_to_symmetric_grid.
%
%   See also: computeNetworkHeritability, computeFunctionalConnectivity, get_mag_se

% --- USER CONFIGURATION (edit for your own system) ---
desktop_dir = '/path/to/desktop';            % originally /Users/davidgruskin/Desktop
r_input_dir = [desktop_dir '/isc_heritability_r'];  % SOLAR input tree = calc_h2.R base_dir
% -----------------------------------------------------
out_dir = fullfile(r_input_dir, 'fc');        % calc_h2.R reads paste0(base_dir,"/fc/fc_net_...")
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Fisher-z helpers (match the rest of the pipeline). atanh/tanh are used so this
% function does not depend on conv_r2z/conv_z2r being on the path.
r2z = @(r) atanh(min(max(r, -0.999999), 0.999999));
z2r = @(z) tanh(z);

unique_networks = unique(net_ids);
num_nets        = numel(unique_networks);
n_combos        = num_nets + (num_nets * (num_nets - 1)) / 2;   % 153 for 17 networks
numsubjects     = size(fc_square, 3);
num_days        = size(fc_square, 4);

if numsubjects ~= numel(subjects)
    warning('computeNetworkMagnitude:subjectMismatch', ...
        'size(fc_square,3) = %d but numel(subjects) = %d; using size(fc_square,3).', ...
        numsubjects, numel(subjects));
end

fcMagnitude = NaN(n_combos, numsubjects, num_days);

for scan = 1:num_days
    combo = 0;
    for net1 = 1:num_nets
        for net2 = net1:num_nets
            combo = combo + 1;

            idx1  = find(net_ids == unique_networks(net1));
            idx2  = find(net_ids == unique_networks(net2));
            block = fc_square(idx1, idx2, :, scan);            % [k1 x k2 x nSubj]
            k1    = numel(idx1);
            k2    = numel(idx2);

            % ASSUMPTION 1 (connection set): for between-network combos use every
            % parcel-to-parcel connection; for within-network combos use each
            % unique pair once and exclude the diagonal (self-FC = 1). This matches
            % the paper's "all unique pairs of one parcel from each network".
            if net1 == net2
                mask = triu(true(k1), 1);
            else
                mask = true(k1, k2);
            end

            vals = reshape(block, k1 * k2, numsubjects);       % column-major, matches mask(:)
            vals = vals(mask(:), :);                           % [nPairs x nSubj]

            % ASSUMPTION 2 (averaging): Fisher-z average across parcel pairs, then
            % convert back to r, consistent with how FC values are averaged
            % elsewhere in this project.
            fcMagnitude(combo, :, scan) = z2r(nanmean(r2z(vals), 1));
        end
    end

    % Write the SOLAR input phenotype file for this scan: rows = combinations,
    % columns = subjects, no header (calc_h2.R transposes it).
    out_file = fullfile(out_dir, sprintf('fc_net_%s_%s_parc_%s_scan_%d.csv', ...
        task_type, align_type, curr_parc, scan));
    writematrix(fcMagnitude(:, :, scan), out_file);
end

end
