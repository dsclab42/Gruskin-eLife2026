%%%%%%%%
% Author: David Gruskin
% Contact: dcg2153@cumc.columbia.edu
% Project: Heritability of movie-evoked brain activity and connectivity
% Description: This is the main analysis and visualization script for the
% FC component of the HCP 7T ISC heritability project.
%%%%%%%%


% Add paths
% ------------------------------------------------------------------------
desktop_dir = '/path/to/desktop';
data_dir    = '/path/to/data';    
% ------------------------------------------------------------------------
addpath(fullfile(fileparts(mfilename('fullpath')), 'dependencies'))  % helper functions shipped in matlab_code/dependencies
addpath([desktop_dir '/isc_heritability_r'])

% Set constants
num_nets  = 17;
num_edges = 153;
num_trs   = [1432, 1409];
num_parcs = 10;
baseDir   = [data_dir '/isc_heritability/data/'];
base_dir  = baseDir;  
load([data_dir '/isc_heritability/covari_motion.mat']) % Note: covari_motion columns correspond to Age, Gender, Movie Day 1 Mean FD, Movie Day 2 Mean FD, Rest Day 1 Mean FD, Rest Day 2 Mean FD

network_names_kong = {'Default A';'Default B';'Default C';'Language';'Control A';'Control B';'Control C';'Ventral Attention A';'Ventral Attention B';'Dorsal Attention A';'Dorsal Attention B';'Auditory';'Somatomotor A';'Somatomotor B';'Visual A';'Visual B';'Visual C'};
network_names_orig = {'DefaultA';'DefaultB';'DefaultC';'Language';'ContA';'ContB';'ContC';'SalVenAttnA';'SalVenAttnB';'DorsAttnA';'DorsAttnB';'Aud';'SomMotA';'SomMotB';'VisualA';'VisualB';'VisualC'};

% Read and process network data
filePath = [desktop_dir '/Schaefer2018_400Parcels_Kong2022_17Networks_order_info.txt'];

[net_ids, parcelColors] = processNetworkData(filePath, 400, network_names_orig); % Add all network names

set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');

%% Anatomical FC Profile analyses

% FC profile heritability
if exist(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit.mat'),'file') == 2
    load(strcat(base_dir,'anatomical/outputs/parc/rest_fc_herit.mat'))
    load(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit.mat'))

    load(strcat(base_dir,'anatomical/outputs/parc/movie_fc_similarity.mat'))
    load(strcat(base_dir,'anatomical/outputs/parc/rest_fc_similarity.mat'))

    load(strcat(base_dir,'anatomical/outputs/parc/rest_fc_herit_jack.mat'))
    load(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit_jack.mat'))

    load(strcat(base_dir,'anatomical/outputs/parc/rest_fc_herit_perm.mat'))
    load(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit_perm.mat'))
else
    rest_subj_data  = cell(num_days,1);
    movie_subj_data = cell(num_days,1);
    for scan_id = 1:num_days
        data                     = loadData(length(subjs_rest), baseDir, 'anatomical','rest','parc','400',scan_id,num_trs);
        rest_subj_data{scan_id}  = data(kong_transform,:,:);
        data                     = loadData(length(subjs_movie), baseDir, 'anatomical','movie','parc','400',scan_id,num_trs);
        movie_subj_data{scan_id} = data(kong_transform,:,:);
    end

    % Compute functional connectivity
    [rest_fc, rest_fc_square]   = computeFunctionalConnectivity(rest_subj_data, length(subjs_rest));
    [movie_fc, movie_fc_square] = computeFunctionalConnectivity(movie_subj_data, length(subjs_movie));

    % Compute network magnitude
    [rest_fcMagnitude]  = computeNetworkMagnitude(rest_fc_square, net_ids, subjs_rest,'rest','anatomical','400');
    [movie_fcMagnitude] = computeNetworkMagnitude(movie_fc_square, net_ids, subjs_movie,'movie','anatomical','400');

    % Compute network heritability (change last argument from [] to fam_ids to get
    % jackknife results)
    [rest_fc_herit, rest_fc_herit_perm,rest_fc_herit_jack,rest_fc_similarity]     = computeNetworkHeritability(rest_fc_square, net_ids, subjs_rest,kinship,covari_motion(:,[1, 2, 5, 6]),num_perm,[]);
    [movie_fc_herit, movie_fc_herit_perm,movie_fc_herit_jack,movie_fc_similarity] = computeNetworkHeritability(movie_fc_square, net_ids, subjs_movie,kinship,covari_motion(:,[1,2,3,4]),num_perm,[]);

    save(strcat(base_dir,'anatomical/outputs/parc/rest_fc_herit.mat'),'rest_fc_herit')
    save(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit.mat'),'movie_fc_herit')
    save(strcat(base_dir,'anatomical/outputs/parc/movie_fc_similarity.mat'),'movie_fc_similarity')
    save(strcat(base_dir,'anatomical/outputs/parc/rest_fc_similarity.mat'),'rest_fc_similarity')
    save(strcat(base_dir,'anatomical/outputs/parc/rest_fc_herit_jack.mat'),'rest_fc_herit_jack')
    save(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit_jack.mat'),'movie_fc_herit_jack')

    if num_perm == 10000
        save(strcat(base_dir,'anatomical/outputs/parc/rest_fc_herit_perm.mat'),'rest_fc_herit_perm')
        save(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit_perm.mat'),'movie_fc_herit_perm')
    end
end

% Difference in FC heritability for Movie vs. Rest data
fc_herit_task_diff = movie_fc_herit - rest_fc_herit;
if num_perm  == 10000
    fc_herit_task_diff_perm   = zeros(num_edges, num_days, num_perm);
    movie_rest_subjs_indices  = find(ismember(subjs_movie, subjs_rest));
    fc_herit_task_perm_square = zeros(num_nets, num_nets, num_perm, num_days);
    fc_herit_task_diff_square = zeros(num_nets, num_nets, num_days);

    fc_herit_task_diff_perm_square_net_mean = zeros(num_nets, num_perm, num_days);
    fc_herit_task_diff_square_net_mean_pval = zeros(num_nets, num_days);
    fc_herit_task_square_net_mean_pval_fdr  = zeros(num_nets, num_days);

    rest_covari  = covari_motion(subjs_rest,[1, 2, 5, 6]);
    movie_covari = covari_motion(subjs_rest,[1, 2, 3, 4]);

    parfor iter = 1:num_perm
        for scan_id = 1:num_days
            [shuffled_rest_fc, shuffled_movie_fc,shuffled_rest_covari,shuffled_movie_covari,swap_indices] = shuffleFCData(rest_fc_square(:,:,:,scan_id), movie_fc_square(:, :, movie_rest_subjs_indices, scan_id),rest_covari(:,[1,2,scan_id+2]),movie_covari(:,[1,2,scan_id+2]), iter);

            [rest_fc_herit_shuffle, ~,~,~]   = computeNetworkHeritability(shuffled_rest_fc, net_ids, 1:length(subjs_rest),kinship(subjs_rest,subjs_rest),shuffled_rest_covari,0,[]);
            [movie_fc_herit_shuffle, ~,~,~]  = computeNetworkHeritability(shuffled_movie_fc, net_ids, 1:length(subjs_rest),kinship(subjs_rest,subjs_rest),shuffled_movie_covari,0,[]);

            fc_herit_task_diff_perm(:, scan_id, iter) = movie_fc_herit_shuffle - rest_fc_herit_shuffle;
        end

    end

    for scan_id = 1:num_days
        for iter = 1:num_perm
            fc_herit_task_perm_square(:, :, iter, scan_id) = vector_to_symmetric_grid(fc_herit_task_diff_perm(:, scan_id, iter), num_nets);
        end

        fc_herit_task_diff_square(:, :, scan_id)               = vector_to_symmetric_grid(fc_herit_task_diff(:, scan_id), num_nets);
        fc_herit_task_diff_perm_square_net_mean(:, :, scan_id) = squeeze(mean(fc_herit_task_perm_square(:, :, :, scan_id), 2));

        for net = 1:num_nets
            fc_herit_task_diff_square_net_mean_pval(net, scan_id) = sum(abs(fc_herit_task_diff_perm_square_net_mean(net, :, scan_id)) > abs(nanmean(fc_herit_task_diff_square(net,:, scan_id),2))) / num_perm;
        end
    end

    for scan_id = 1:num_days
        fc_herit_task_square_net_mean_pval_fdr(:,scan_id) = fdr_bh(fc_herit_task_diff_square_net_mean_pval(:,scan_id));
    end
    fc_herit_task_square_net_mean_pval_fdr_sig = fc_herit_task_square_net_mean_pval_fdr(:,1).*fc_herit_task_square_net_mean_pval_fdr(:,2);
    save(strcat(base_dir,'anatomical/outputs/parc/fc_herit_task_diff_perm.mat'),'fc_herit_task_diff_perm')
    save(strcat(base_dir,'anatomical/outputs/parc/fc_herit_task_diff_square_net_mean_pval.mat'),'fc_herit_task_diff_square_net_mean_pval')
end

% Further analyses for paper
movie_fc_herit = load(strcat(base_dir,'anatomical/outputs/parc/movie_fc_herit.mat'));
movie_fc_herit = movie_fc_herit.movie_fc_herit;

% Percent increase in heritability from Rest to Movie data
fc_herit_task_diff_perc = zeros(num_nets,2);
fc_herit_task_diff      = zeros(num_edges,2);

for scan_id = 1:num_days
    fc_herit_task_diff_perc(:,scan_id) = nanmean(vector_to_symmetric_grid((movie_fc_herit(:,scan_id)./rest_fc_herit(:,scan_id))-1,max(net_ids)));
    fc_herit_task_diff(:,scan_id)      = movie_fc_herit(:,scan_id) - rest_fc_herit(:,scan_id);
end

fc_herit_task_diff_perc_sig = fc_herit_task_diff_perc(fc_herit_task_square_net_mean_pval_fdr_sig==1,:);

mean_fc_herit_task_diff_perc_sig = mean(fc_herit_task_diff_perc_sig);
min_fc_herit_task_diff_perc_sig  = min(fc_herit_task_diff_perc_sig);
max_fc_herit_task_diff_perc_sig  = max(fc_herit_task_diff_perc_sig);

mean_movie_fc_herit = mean(movie_fc_herit);
std_movie_fc_herit  = std(movie_fc_herit);

% Actual test-retest reliabtility between FC heritability patterns from Day 1 and
% Day 2
[movie_fc_herit_trt,~] = corr(movie_fc_herit,'rows','complete','type','spearman');

% Significance testing for FC heritability TRT
movie_net_herit_trt_perm = zeros(num_perm,1);
for perm = 1:10000
    rng(perm)

    % Shuffle the heritability values of the second day
    shuffled_herit = movie_fc_herit(randperm(size(movie_fc_herit, 1)), 2);

    % Compute the correlation for the shuffled data
    movie_net_herit_trt_perm(perm) = corr(movie_fc_herit(:, 1), shuffled_herit,'type','spearman');
end

movie_trt_perm_p = mean(abs(movie_net_herit_trt_perm) >= abs(movie_fc_herit_trt(1,2)));


movie_fc_herit_pvalues          = movie_fc_herit_perm;
movie_fc_herit_pvalues_fdr      = fdr_bh(movie_fc_herit_perm);
movie_fc_herit_pvalues_fdr_perc = sum(movie_fc_herit_pvalues_fdr(:,1).*movie_fc_herit_pvalues_fdr(:,2))/num_edges;

movie_fc_herit_square = zeros(num_nets,num_nets,num_days);
rest_fc_herit_square  = zeros(num_nets,num_nets,num_days);

movie_fc_herit_jack_square = zeros(num_nets,num_nets,length(unique(movie_fam_ids)),num_days);
rest_fc_herit_jack_square  = zeros(num_nets,num_nets,length(unique(movie_fam_ids)),num_days);

for scan_id = 1:num_days
    movie_fc_herit_square(:,:,scan_id) = vector_to_symmetric_grid(movie_fc_herit(:,scan_id),max(net_ids));
    rest_fc_herit_square(:,:,scan_id)  = vector_to_symmetric_grid(rest_fc_herit(:,scan_id),max(net_ids));
    for fam_id = 1:length(unique(movie_fam_ids))
        movie_fc_herit_jack_square(:,:,fam_id,scan_id) = vector_to_symmetric_grid(movie_fc_herit_jack(:,fam_id,scan_id),max(net_ids));
        rest_fc_herit_jack_square(:,:,fam_id,scan_id)  = vector_to_symmetric_grid(rest_fc_herit_jack(:,fam_id,scan_id),max(net_ids));
    end
end

movie_fc_se = zeros(num_parcs,num_days);

for scan_id = 1:num_days
    [~,movie_fc_se(1,scan_id)] = compareHeritability(nanmean(movie_fc_herit(:,scan_id)),nanmean(movie_fc_herit_jack(:,:,scan_id),1),nanmean(rest_fc_herit(:,scan_id)),nanmean(movie_fc_herit_jack(:,:,scan_id),1));
end
%% Anatomical FC Magnitude analysis

fc_movie_magnitude_anatomical = cell(1,num_days);
fc_rest_magnitude_anatomical  = cell(1,num_days);

for scan_id = 1:num_days
    fc_movie_magnitude_anatomical{1,scan_id} = readtable(strcat([data_dir "/isc_heritability/data/solar/fc/fc_herit_net_movie_anatomical_parc_400_scan_"], num2str(scan_id), ".csv"));
    fc_rest_magnitude_anatomical{1,scan_id}  = readtable(strcat([data_dir "/isc_heritability/data/solar/fc/fc_herit_net_rest_anatomical_parc_400_scan_"], num2str(scan_id), ".csv"));
end

movie_fc_herit_mag     = zeros(num_edges,num_days);
rest_fc_herit_mag      = zeros(num_edges,num_days);
task_diff_fc_herit_mag = zeros(num_edges,num_days);

for scan_id = 1:num_days
    task_diff_fc_herit_mag(:,scan_id) = fc_movie_magnitude_anatomical{1,scan_id}.Var - fc_rest_magnitude_anatomical{1,scan_id}.Var;
    movie_fc_herit_mag(:,scan_id)     = fc_movie_magnitude_anatomical{1,scan_id}.Var;
    rest_fc_herit_mag(:,scan_id)      = fc_rest_magnitude_anatomical{1,scan_id}.Var;
end


% Actual test-retest reliabtility between FC heritability patterns from Day 1 and
% Day 2
[fc_movie_magnitude_anatomical_trt_r,~] = corr(fc_movie_magnitude_anatomical{1,1}.Var,fc_movie_magnitude_anatomical{1,2}.Var,'type','spearman');
[fc_rest_magnitude_anatomical_trt_r,~]  = corr(fc_rest_magnitude_anatomical{1,1}.Var,fc_rest_magnitude_anatomical{1,2}.Var,'type','spearman');

% Significance testing for FC heritability TRT
movie_mag_trt_perm = zeros(num_perm,1);
rest_mag_trt_perm  = zeros(num_perm,1);

for perm = 1:num_perm
    rng(perm)
    shuffled_herit_movie_mag = movie_fc_herit(randperm(size(fc_movie_magnitude_anatomical{1,2}.Var, 1)), 2);
    shuffled_herit_rest_mag  = movie_fc_herit(randperm(size(fc_rest_magnitude_anatomical{1,2}.Var, 1)), 2);

    movie_mag_trt_perm(perm) = corr(fc_movie_magnitude_anatomical{1,1}.Var, shuffled_herit_movie_mag,'type','spearman','rows','complete');
    rest_mag_trt_perm(perm)  = corr(fc_rest_magnitude_anatomical{1,1}.Var, shuffled_herit_rest_mag,'type','spearman','rows','complete');
end

pvalue_movie_mag = sum(abs(movie_mag_trt_perm) >= abs(fc_movie_magnitude_anatomical_trt_r))/num_perm;
pvalue_rest_mag  = sum(abs(rest_mag_trt_perm) >= abs(fc_rest_magnitude_anatomical_trt_r))/num_perm;

movie_fc_herit_mag_pval = zeros(num_edges,num_days);
movie_fc_herit_mag      = zeros(num_edges,num_days);

% Get useful stats for paper
for scan_id = 1:num_days
    movie_fc_herit_mag_pval(:,scan_id) = fc_movie_magnitude_anatomical{1,scan_id}.pval;
    movie_fc_herit_mag(:,scan_id)      = fc_movie_magnitude_anatomical{1,scan_id}.Var;
end

fc_movie_magnitude_pval_fdr      = fdr_bh(movie_fc_herit_mag_pval);
fc_movie_magnitude_pval_fdr_perc = sum(fc_movie_magnitude_pval_fdr(:,1).*fc_movie_magnitude_pval_fdr(:,2))/num_edges;

fc_movie_mag_herit_mean = mean(movie_fc_herit_mag);
fc_movie_mag_herit_std  = std(movie_fc_herit_mag);

movie_fc_herit_mag_square = zeros(num_nets,num_nets,num_days);
rest_fc_herit_mag_square  = zeros(num_nets,num_nets,num_days);

for scan_id = 1:num_days
    movie_fc_herit_mag_square(:,:,scan_id) = vector_to_symmetric_grid(movie_fc_herit_mag(:,scan_id),max(net_ids));
    rest_fc_herit_mag_square(:,:,scan_id)  = vector_to_symmetric_grid(rest_fc_herit_mag(:,scan_id),max(net_ids));
end

[anatomical_mag_pval, anatomical_mag_actual_diffs,anatomical_actual_diffs_perc] = computeNetworkHeritabilityPValues(movie_fc_herit_mag_square, rest_fc_herit_mag_square, 10000);

anatomical_pval_magnitude_fdr      = fdr_bh(anatomical_mag_pval);
anatomical_pval_magnitude_fdr_perc = sum(anatomical_pval_magnitude_fdr(:,1).*anatomical_pval_magnitude_fdr(:,2))/num_nets;

mag_dir = [data_dir '/isc_heritability/data/solar/fc'];

movie_mag_se = zeros(1,2);
for scan_id = 1:2
    movie_mag_se(1,scan_id) = get_mag_se(mag_dir,'movie', 'anatomical','400',scan_id);
end

rest_mag_se = zeros(1,2);
for scan_id = 1:2
    rest_mag_se(1,scan_id) = get_mag_se(mag_dir,'rest', 'anatomical','400',scan_id);
end

%% Piecewise FC Profile analysis

if exist(strcat(base_dir,'piecewise/outputs/parc/piecewise_herit.mat'),'file') == 2
    load(strcat(base_dir,'piecewise/outputs/parc/piecewise_similarity.mat'))
    load(strcat(base_dir,'piecewise/outputs/parc/piecewise_herit.mat'))
    load(strcat(base_dir,'piecewise/outputs/parc/piecewise_herit_jack.mat'))
    load(strcat(base_dir,'piecewise/outputs/parc/piecewise_herit_perm.mat'))
else
    % Load and process data
    piecewise_subj_data = cell(num_parcs,1);
    for parc_id = 1:num_parcs
        curr_parc = parc_reses{parc_id};
        for scan_id = 1:num_days
            piecewise_subj_data{parc_id,1}{scan_id} = loadData(num_subj_movie, baseDir, 'piecewise','movie','parc',curr_parc,scan_id,num_trs);
        end
    end

    % Compute functional connectivity
    piecewise_fc        = cell(num_parcs,1);
    piecewise_fc_square = cell(num_parcs,1);

    for parc_id = 1:length(parc_reses)
        [piecewise_fc{parc_id,1}, piecewise_fc_square{parc_id,1}] = computeFunctionalConnectivity(piecewise_subj_data{parc_id,1}, num_subj_movie);
    end

    % Compute network heritability
    piecewise_herit      = cell(num_parcs,1);
    piecewise_herit_perm = cell(num_parcs,1);
    piecewise_herit_jack = cell(num_parcs,1);
    piecewise_similarity = cell(num_parcs,1);
    piecewise_herit_avg  = zeros(num_parcs,num_days);

    for parc_id = 1:length(parc_reses)
        [piecewise_herit{parc_id,1}, piecewise_herit_perm{parc_id,1},piecewise_herit_jack{parc_id,1},piecewise_similarity{parc_id,1}] = computeNetworkHeritability(piecewise_fc_square{parc_id,1}, net_ids, subjs_movie,kinship,covari_motion,num_perm,[]);
         piecewise_herit_avg(parc_id,:) = mean(piecewise_herit{parc_id,1});
    end

    if num_perm == 10000
        save(strcat(base_dir,'piecewise/outputs/parc/piecewise_similarity.mat'),'piecewise_similarity')
        save(strcat(base_dir,'piecewise/outputs/parc/piecewise_herit.mat'),'piecewise_herit')
        save(strcat(base_dir,'piecewise/outputs/parc/piecewise_herit_jack.mat'),'piecewise_herit_jack')
        save(strcat(base_dir,'piecewise/outputs/parc/piecewise_herit_perm.mat'),'piecewise_herit_perm')
    end
end

piecewise_net_herit_square    = cell(num_parcs,1);
piecewiseNetHerit_jack_square = cell(num_parcs,1);
for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        piecewise_net_herit_square{parc_id,1}(:,:,scan_id) = vector_to_symmetric_grid(piecewise_herit{parc_id,1}(:,scan_id),max(net_ids));
        for fam_id = 1:length(unique(movie_fam_ids))
            piecewiseNetHerit_jack_square{parc_id}(:,:,fam_id,scan_id) = vector_to_symmetric_grid(piecewise_herit_jack{parc_id,1}(:,fam_id,scan_id),max(net_ids));
        end
    end
end

% Compute network magnitude
for parc_id = 1:length(parc_reses)
    movie_fcMagnitude{parc_id,1} = computeNetworkMagnitude(piecewise_fc_square{parc_id,1}, net_ids, subjs_movie,'movie','piecewise',parc_reses{parc_id});
end

% Get standard errors
piecewise_fc_se = zeros(num_parcs,num_days);
for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        [~,~,piecewise_fc_se(parc_id,scan_id)] = compareHeritability(nanmean(movie_fc_herit(:,scan_id)),nanmean(movie_fc_herit_jack(:,:,scan_id),1),piecewise_herit_avg(parc_id,scan_id),nanmean(piecewise_herit_jack{parc_id,1}(:,:,scan_id),1));
    end
end

%% Piecewise FC Magnitude analysis

fc_movie_magnitude_piecewise     = cell(num_parcs,num_days);
fc_movie_magnitude_piecewise_avg = zeros(num_parcs,num_days);

for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        fc_movie_magnitude_piecewise{parc_id,scan_id}     = readtable(strcat([data_dir '/isc_heritability/data/solar/fc/fc_herit_net_movie_piecewise_parc_'],parc_reses{parc_id},'_scan_', num2str(scan_id), '.csv'));
        fc_movie_magnitude_piecewise_avg(parc_id,scan_id) = nanmean(fc_movie_magnitude_piecewise{parc_id,scan_id}.Var);
    end
end

piecewise_fc_magnitude                 = zeros(num_edges,num_parcs,num_days);
piecewise_anatomical_fc_magnitude_diff = zeros(num_edges,num_parcs,num_days);

for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        piecewise_fc_magnitude(:,parc_id,scan_id)                 = fc_movie_magnitude_piecewise{parc_id,scan_id}.Var;
        piecewise_anatomical_fc_magnitude_diff(:,parc_id,scan_id) = fc_movie_magnitude_piecewise{parc_id,scan_id}.Var - movie_fc_herit_mag(:,scan_id);
    end
end


piecewise_net_herit_mag_square = zeros(num_nets,num_nets,num_days);
for scan_id = 1:num_days
    piecewise_net_herit_mag_square(:,:,scan_id) = vector_to_symmetric_grid(fc_movie_magnitude_piecewise{1,scan_id}.Var,max(net_ids));
end

piecewise_mag_corrected_pval = zeros(num_nets,2);
for scan_id = 1:num_days
    piecewise_mag_corrected_pval(:,scan_id) = fdr_bh(piecewise_mag_pval(:,scan_id));
end

piecewise_mag_se = zeros(num_parcs,num_days);
for parc_id = 1:num_parcs
    curr_parc = parc_reses{parc_id};
    for scan_id = 1:num_days
        piecewise_mag_se(parc_id,scan_id) = get_mag_se(mag_dir,'movie', 'piecewise',curr_parc,scan_id);
    end
end


%% Connectivity FC Profile analysis

if exist(strcat(base_dir,'connectivity/outputs/parc/connectivity_herit.mat'),'file') == 2
    load(strcat(base_dir,'connectivity/outputs/parc/connectivity_similarity.mat'))
    load(strcat(base_dir,'connectivity/outputs/parc/connectivity_herit.mat'))
    load(strcat(base_dir,'connectivity/outputs/parc/connectivity_herit_jack.mat'))
    load(strcat(base_dir,'connectivity/outputs/parc/connectivity_herit_perm.mat'))
else
    % Load and process data
    connectivity_subj_data = cell(num_parcs,1);
    for parc_id = 1:num_parcs
        curr_parc = parc_reses{parc_id};
        for scan_id = 1:num_days
            connectivity_subj_data{parc_id,1}{scan_id} = loadData(num_subj_movie, baseDir, 'connectivity','movie','parc',curr_parc,scan_id,num_trs);
        end
    end

    % Compute functional connectivity
    connectivity_fc        = cell(num_parcs,1);
    connectivity_fc_square = cell(num_parcs,1);

    for parc_id = 1:length(parc_reses)
        [connectivity_fc{parc_id,1}, connectivity_fc_square{parc_id,1}] = computeFunctionalConnectivity(connectivity_subj_data{parc_id,1}, num_subj_movie);
    end

    % Compute network heritability
    connectivity_herit      = cell(num_parcs,1);
    connectivity_herit_perm = cell(num_parcs,1);
    connectivity_herit_jack = cell(num_parcs,1);
    connectivity_similarity = cell(num_parcs,1);
    connectivity_herit_avg  = zeros(num_parcs,num_days);

    for parc_id = 1:length(parc_reses)
        [connectivity_herit{parc_id,1}, connectivity_herit_perm{parc_id,1},connectivity_herit_jack{parc_id,1},connectivity_similarity{parc_id,1}] = computeNetworkHeritability(connectivity_fc_square{parc_id,1}, net_ids, subjs_movie,kinship,covari_motion,num_perm,[]);
         connectivity_herit_avg(parc_id,:) = mean(connectivity_herit{parc_id,1});
    end

    if num_perm == 10000
        save(strcat(base_dir,'connectivity/outputs/parc/connectivity_similarity.mat'),'connectivity_similarity')
        save(strcat(base_dir,'connectivity/outputs/parc/connectivity_herit.mat'),'connectivity_herit')
        save(strcat(base_dir,'connectivity/outputs/parc/connectivity_herit_jack.mat'),'connectivity_herit_jack')
        save(strcat(base_dir,'connectivity/outputs/parc/connectivity_herit_perm.mat'),'connectivity_herit_perm')
    end
end

connectivity_net_herit_square    = cell(num_parcs,1);
connectivityNetHerit_jack_square = cell(num_parcs,1);
for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        connectivity_net_herit_square{parc_id,1}(:,:,scan_id) = vector_to_symmetric_grid(connectivity_herit{parc_id,1}(:,scan_id),max(net_ids));
        for fam_id = 1:length(unique(movie_fam_ids))
            connectivityNetHerit_jack_square{parc_id}(:,:,fam_id,scan_id) = vector_to_symmetric_grid(connectivity_herit_jack{parc_id,1}(:,fam_id,scan_id),max(net_ids));
        end
    end
end


% Compute network magnitude
for parc_id = 1:length(parc_reses)
    movie_fcMagnitude{parc_id,1} = computeNetworkMagnitude(connectivity_fc_square{parc_id,1}, net_ids, subjs_movie,'movie','connectivity',parc_reses{parc_id});
end

% Get standard errors
connectivity_fc_se = zeros(num_parcs,num_days);
for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        [~,~,connectivity_fc_se(parc_id,scan_id)] = compareHeritability(nanmean(movie_fc_herit(:,scan_id)),nanmean(movie_fc_herit_jack(:,:,scan_id),1),connectivity_herit_avg(parc_id,scan_id),nanmean(connectivity_herit_jack{parc_id,1}(:,:,scan_id),1));
    end
end


%% Connectivity FC Magnitude analysis

fc_movie_magnitude_connectivity     = cell(num_parcs,num_days);
fc_movie_magnitude_connectivity_avg = zeros(num_parcs,num_days);

for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        fc_movie_magnitude_connectivity{parc_id,scan_id}     = readtable(strcat([data_dir '/isc_heritability/data/solar/fc/fc_herit_net_movie_connectivity_parc_'],parc_reses{parc_id},'_scan_', num2str(scan_id), '.csv'));
        fc_movie_magnitude_connectivity_avg(parc_id,scan_id) = nanmean(fc_movie_magnitude_connectivity{parc_id,scan_id}.Var);
    end
end

connectivity_fc_magnitude                 = zeros(num_edges,num_parcs,num_days);
connectivity_anatomical_fc_magnitude_diff = zeros(num_edges,num_parcs,num_days);

for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        connectivity_fc_magnitude(:,parc_id,scan_id)                 = fc_movie_magnitude_connectivity{parc_id,scan_id}.Var;
        connectivity_anatomical_fc_magnitude_diff(:,parc_id,scan_id) = fc_movie_magnitude_connectivity{parc_id,scan_id}.Var - movie_fc_herit_mag(:,scan_id);
    end
end


connectivity_net_herit_mag_square = zeros(num_nets,num_nets,num_days);
for scan_id = 1:num_days
    connectivity_net_herit_mag_square(:,:,scan_id) = vector_to_symmetric_grid(fc_movie_magnitude_connectivity{1,scan_id}.Var,max(net_ids));
end

connectivity_mag_corrected_pval = zeros(num_nets,2);
for scan_id = 1:num_days
    connectivity_mag_corrected_pval(:,scan_id) = fdr_bh(connectivity_mag_pval(:,scan_id));
end

connectivity_mag_se = zeros(num_parcs,num_days);
for parc_id = 1:num_parcs
    curr_parc = parc_reses{parc_id};
    for scan_id = 1:num_days
        connectivity_mag_se(parc_id,scan_id) = get_mag_se(mag_dir,'movie', 'connectivity',curr_parc,scan_id);
    end
end

%% Anatomical FC Profile dyadic similarity
fc_similarity_data_path  = strcat(base_dir, 'anatomical/outputs/parc/fc_similarity.mat');
fc_similarity            = load(fc_similarity_data_path);
fc_similarity_fieldnames = fieldnames(fc_similarity);
fc_similarity            = fc_similarity.(fc_similarity_fieldnames{1});

% Initialize anatomical similarity matrices
anatomical_fc_similarity_mz_net_all = zeros(length(pairwise_mz_id), num_edges, num_days);
anatomical_fc_similarity_dz_net_all = zeros(length(pairwise_dz_id), num_edges, num_days);
anatomical_fc_similarity_nz_net_all = zeros(length(pairwise_nz_id), num_edges, num_days);

% Process fc_similarity data for each scan and network
for scan_id = 1:num_days
    for net = 1:size(fc_similarity, 3)
        fc_similarity_temp                   = squeeze(fc_similarity(:,:,net,scan_id));
        [m, n]                               = size(fc_similarity_temp);
        fc_similarity_temp(triu(true(m, n))) = NaN; % Set upper triangle to NaN

        % Store data for each group
        anatomical_fc_similarity_mz_net_all(:,net,scan_id) = fc_similarity_temp(pairwise_mz_id);
        anatomical_fc_similarity_dz_net_all(:,net,scan_id) = fc_similarity_temp(pairwise_dz_id);
        anatomical_fc_similarity_nz_net_all(:,net,scan_id) = fc_similarity_temp(pairwise_nz_id);
    end

    % Replace diagonal values with NaN and calculate gray net all similarity
    fc_similarity(fc_similarity == 1) = NaN;
    anatomical_fc_similarity_gray_net_all(:,scan_id) = conv_z2r(squeeze(nanmean(nanmean(conv_r2z(fc_similarity(:,:,:,scan_id))))));
end

% Concatenate and process data for plotting
all_xz_fc_similarity = cat(2, ...
    conv_z2r(squeeze(nanmean(conv_r2z(anatomical_fc_similarity_mz_net_all), 1))), ...
    conv_z2r(squeeze(nanmean(conv_r2z(anatomical_fc_similarity_dz_net_all), 1))), ...
    conv_z2r(squeeze(nanmean(conv_r2z(anatomical_fc_similarity_nz_net_all), 1))));

% Initialize matrices for standard deviations
mz_std_mat = zeros(size(anatomical_fc_similarity_mz_net_all, 2), num_days);
dz_std_mat = mz_std_mat;
nz_std_mat = mz_std_mat;
sort_index = mz_std_mat;

% Compute standard deviations and prepare data for plotting
for scan_id = 1:num_days
    mz_std_mat(:, scan_id) = conv_z2r(nanstd(conv_r2z(anatomical_fc_similarity_mz_net_all(:,:,scan_id)))' / sqrt(sum(~isnan(anatomical_fc_similarity_mz_net_all(:,1,scan_id)))));
    dz_std_mat(:, scan_id) = conv_z2r(nanstd(conv_r2z(anatomical_fc_similarity_dz_net_all(:,:,scan_id)))' / sqrt(sum(~isnan(anatomical_fc_similarity_dz_net_all(:,1,scan_id)))));
    nz_std_mat(:, scan_id) = conv_z2r(nanstd(conv_r2z(anatomical_fc_similarity_nz_net_all(:,:,scan_id)))' / sqrt(sum(~isnan(anatomical_fc_similarity_nz_net_all(:,1,scan_id)))));
end

% Loop through each scan for sorting and plotting
for scan_id = 1:num_days
    [~, sort_index(:, scan_id)] = sort(anatomical_fc_similarity_gray_net_all(:, scan_id), 'ascend');
    sorted_means                = all_xz_fc_similarity(sort_index(:, scan_id), [scan_id, scan_id + 2, scan_id + 4]);

    mz_std_mat_sorted = mz_std_mat(sort_index(:, scan_id), scan_id);
    dz_std_mat_sorted = dz_std_mat(sort_index(:, scan_id), scan_id);
    nz_std_mat_sorted = nz_std_mat(sort_index(:, scan_id), scan_id);
    sorted_std_errors = cat(2, mz_std_mat_sorted, dz_std_mat_sorted, nz_std_mat_sorted);

    % Plotting routine for each scan
    subplot(1, 2, scan_id);
    x         = 1:size(anatomical_fc_similarity_gray_net_all, 1); % New x-axis vector
    color_vec = [[100/255, 143/255, 255/255]; [220/255, 38/255, 127/255]; [255/255, 176/255, 0/255]]; % Colors for each line

    % Plot each line with shaded error
    for i = 1:size(sorted_means, 2)
        boundedline(x, sorted_means(:, i), sorted_std_errors(:, i), 'cmap', color_vec(i, :), 'alpha', 'linewidth', 2);
        hold on;
    end

    % Axes and plot customization
    ylim([0, 1]); xlim([0, size(anatomical_fc_similarity_gray_net_all, 1)]);
    xlabel('netel Rank'); ylabel('fc_similarity (r)');
    yticks([0, 0.25, .5, .75, 1]);
    ax = gca;
    set(ax, 'TickDir', 'out', 'FontSize', 20, 'LineWidth', 3, 'Layer', 'bottom');
    grid off;
    legend({'', 'MZ', '', 'DZ', '', 'UR', ''});
    pbaspect([1, 1, 1]);
end

% Save the figure
set(gcf, 'position', [100, 100, 1200, 1200]);
saveas(gcf, [data_dir '/isc_heritability/figs/fig1_fc_similaritylines_avgsort_net.fig']);
close all;
% Start MATLAB's parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

% Data arrays
dataArrays = {anatomical_fc_similarity_mz_net_all, anatomical_fc_similarity_dz_net_all, anatomical_fc_similarity_nz_net_all};

% Preallocate arrays for t-scores and p-values
fc_dyad_pval = zeros(num_nets, 6);

% Parallel computation over nets
parfor net = 1:num_nets
    temp_differences = zeros(1, 6); % Temporary array for actual differences
    temp_pvals       = zeros(1, 6); % Temporary array for p-values

    for day = 1:num_days
        % Actual data for each group
        dataMZ = dataArrays{1}(:, net, day);
        dataDZ = dataArrays{2}(:, net, day);
        dataNZ = dataArrays{3}(:, net, day);

        % Calculate actual differences
        diffMZDZ = nanmean(dataMZ) - nanmean(dataDZ);
        diffMZNZ = nanmean(dataMZ) - nanmean(dataNZ);
        diffDZNZ = nanmean(dataDZ) - nanmean(dataNZ);

        actualDifferences = [diffMZDZ, diffMZNZ, diffDZNZ];

        % Initialize permutation differences
        permDifferences = zeros(num_perm, 3);

        for perm = 1:num_perm
            % Shuffle group identities while maintaining proportions
            allData = [dataMZ; dataDZ; dataNZ];
            shuffledData = allData(randperm(length(allData)));
            numMZ = numel(dataMZ);
            numDZ = numel(dataDZ);

            shuffledMZ = shuffledData(1:numMZ);
            shuffledDZ = shuffledData(numMZ + 1:numMZ + numDZ);
            shuffledNZ = shuffledData(numMZ + numDZ + 1:end);

            % Calculate differences for shuffled data
            permDifferences(perm, 1) = mean(shuffledMZ) - mean(shuffledDZ);
            permDifferences(perm, 2) = mean(shuffledMZ) - mean(shuffledNZ);
            permDifferences(perm, 3) = mean(shuffledDZ) - mean(shuffledNZ);
        end

        % Calculate p-values
        for comp = 1:3
            temp_pvals((day - 1) * 3 + comp) = mean(abs(permDifferences(:, comp)) >= abs(actualDifferences(comp)));
        end

        % Store actual differences for this day
        temp_differences((day - 1) * 3 + 1) = diffMZDZ;
        temp_differences((day - 1) * 3 + 2) = diffMZNZ;
        temp_differences((day - 1) * 3 + 3) = diffDZNZ;
    end

    % Update differences and p-values for each net
    fc_dyad_diffs(net, :) = temp_differences;
    fc_dyad_pval(net, :)  = temp_pvals;
end

% Apply FDR correction
fc_dyad_pval_fdr = fdr_bh(fc_dyad_pval);

% Extract specific p-values
fc_dyad_pval_mzdz = fc_dyad_pval(:,[1,4]);
fc_dyad_pval_mznz = fc_dyad_pval(:,[2,5]);
fc_dyad_pval_dznz = fc_dyad_pval(:,[3,6]);

fc_dyad_pval_fdr_mzdz = fc_dyad_pval_fdr(:,[1,4]);
fc_dyad_pval_fdr_mznz = fc_dyad_pval_fdr(:,[2,5]);
fc_dyad_pval_fdr_dznz = fc_dyad_pval_fdr(:,[3,6]);

% Calculate percentage of significant values
perc_sig_mzdz_fc = sum(fc_dyad_pval_fdr_mzdz(:,1) .* fc_dyad_pval_fdr_mzdz(:,2)) / num_edges;
perc_sig_mznz_fc = sum(fc_dyad_pval_fdr_mznz(:,1) .* fc_dyad_pval_fdr_mznz(:,2)) / num_edges;
perc_sig_dznz_fc = sum(fc_dyad_pval_fdr_dznz(:,1) .* fc_dyad_pval_fdr_dznz(:,2)) / num_edges;

% Additional statistics for the paper
mznz_perc_fc = nanmean(nanmean((nanmean(conv_r2z(anatomical_fc_similarity_mz_net_all),1) - nanmean(conv_r2z(anatomical_fc_similarity_nz_net_all),1)) ./ nanmean(conv_r2z(anatomical_fc_similarity_nz_net_all),1)));
mzdz_perc_fc = nanmean(nanmean((nanmean(conv_r2z(anatomical_fc_similarity_mz_net_all),1) - nanmean(conv_r2z(anatomical_fc_similarity_dz_net_all),1)) ./ nanmean(conv_r2z(anatomical_fc_similarity_dz_net_all),1)));
dznz_perc_fc = nanmean(nanmean((nanmean(conv_r2z(anatomical_fc_similarity_dz_net_all),1) - nanmean(conv_r2z(anatomical_fc_similarity_nz_net_all),1)) ./ nanmean(conv_r2z(anatomical_fc_similarity_nz_net_all),1)));

% Preparing data for pcolor plots
mz_fc_z_avg = squeeze(nanmean(conv_r2z(anatomical_fc_similarity_mz_net_all), 1));
dz_fc_z_avg = squeeze(nanmean(conv_r2z(anatomical_fc_similarity_dz_net_all), 1));
nz_fc_z_avg = squeeze(nanmean(conv_r2z(anatomical_fc_similarity_nz_net_all), 1));

mzdz_fc_avg = conv_z2r(mz_fc_z_avg - dz_fc_z_avg);
mznz_fc_avg = conv_z2r(mz_fc_z_avg - nz_fc_z_avg);
dznz_fc_avg = conv_z2r(dz_fc_z_avg - nz_fc_z_avg);

mz_fc_avg = conv_z2r(mz_fc_z_avg);
dz_fc_avg = conv_z2r(dz_fc_z_avg);
nz_fc_avg = conv_z2r(nz_fc_z_avg);


% Define plot labels for group differences
plotLabels = {'DZ-NZ', 'MZ-DZ', 'MZ-NZ'};

% Creating pcolor plots
fig_counter = 1;
for scan_id = 1:num_days
    % Plot DZ-NZ, MZ-DZ, and MZ-NZ group differences
    plot_groups = {dznz_fc_avg, mzdz_fc_avg, mznz_fc_avg};
    for i = 1:3
        grid_out = vector_to_symmetric_grid(plot_groups{i}(:,scan_id), max(net_ids));
        ax       = subplot(2, 3, fig_counter);
        mask     = triu(true(size(grid_out)), 1);

        grid_out(mask) = NaN; % Set upper triangle to NaN
        grid_out       = flipud(grid_out); % Flip grid
        grid_out(num_nets+1,:) = NaN; % Set specific rows and columns to NaN

        pcolor(grid_out); % Create pcolor plot
        caxis([-.4, 0.4]); % Set color axis limits
        colormap(ax, bluewhitered(256)); % Set colormap
        axis square; % Set axis to square

        % Set the tick labels
        set(ax, 'YTick', 1:max(net_ids), 'YTickLabel', flip(network_names_kong));
        set(ax, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
        if scan_id == 1
            set(ax, 'XTickLabel', {}); % Clear x-axis and y-axis labels for specific plots
        end
        if i>1
            set(ax, 'YTickLabel', {});
        end
        if i == 1
            ylabel('Networks');
        elseif i == 3
            c = colorbar;
            c.Ticks = [-0.4,  0,  0.4];
            caxis([-0.4, 0.4]);
            tickLabels = {'-0.4',  '0',  '0.4'};
            c.TickLabels = tickLabels;
        end
        title(['Day ' num2str(scan_id) ' - ' plotLabels{i}]);

        fig_counter = fig_counter + 1;
    end
end

% Customize axes and colorbars properties
allAxes      = findall(gcf, 'Type', 'axes');
allColorbars = findall(gcf, 'Type', 'colorbar');
for i = 1:length(allAxes)
    set(allAxes(i), 'LineWidth', 2, 'TickLength', [0 0]); % Set properties for axes
end
for i = 1:length(allColorbars)
    set(allColorbars(i), 'LineWidth', 2); % Set properties for colorbars
end

% Save the figure
saveas(gcf, [data_dir '/isc_heritability/figs/fig1_heatmap_anatomical_fc_profile.fig']);
close all;


plotLabels = {'MZ', 'DZ', 'NZ'};
% Creating pcolor plots
fig_counter = 1;
for scan_id = 1:num_days
    plot_groups = {mz_fc_avg, dz_fc_avg, nz_fc_avg};
    for i = 1:3
        grid_out = vector_to_symmetric_grid(plot_groups{i}(:,scan_id), max(net_ids));
        ax       = subplot(2, 3, fig_counter);
        mask     = triu(true(size(grid_out)), 1);

        grid_out(mask) = NaN; % Set upper triangle to NaN
        grid_out       = flipud(grid_out); % Flip grid
        grid_out(num_nets+1,:) = NaN; % Set specific rows and columns to NaN

        pcolor(grid_out); % Create pcolor plot
        caxis([.25, 1]); % Set color axis limits
        colormap(ax, whitered(256)); % Set colormap
        axis square; % Set axis to square

        % Set the tick labels
        set(ax, 'YTick', 1:max(net_ids), 'YTickLabel', flip(network_names_kong));
        set(ax, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
        if scan_id == 1
            set(ax, 'XTickLabel', {}); % Clear x-axis and y-axis labels for specific plots
        end
        if i>1
            set(ax, 'YTickLabel', {});
        end
        if i == 1
            ylabel('Networks');
        elseif i == 3
            c = colorbar;
            c.Ticks = [.25,.5,.75,1];
            caxis([.25, 1]);
            tickLabels = {'.25','0.5','0.75','1'};
            c.TickLabels = tickLabels;
        end
        title(['Day ' num2str(scan_id) ' - ' plotLabels{i}]);

        fig_counter = fig_counter + 1;
    end
end

% Customize axes and colorbars properties
allAxes      = findall(gcf, 'Type', 'axes');
allColorbars = findall(gcf, 'Type', 'colorbar');
for i = 1:length(allAxes)
    set(allAxes(i), 'LineWidth', 2, 'TickLength', [0 0]); % Set properties for axes
end
for i = 1:length(allColorbars)
    set(allColorbars(i), 'LineWidth', 2); % Set properties for colorbars
end

% Save the figure
saveas(gcf, [data_dir '/isc_heritability/figs/figS1_heatmap_anatomical_fc_profile.fig']);
close all;

%% Other figures

% FIGURE 3

fig_counter = 1;
for scan_id = 1:num_days

    % First figure
    grid_out = vector_to_symmetric_grid(movie_fc_herit(:,scan_id),max(net_ids));

    ax1 = subplot(2,3,fig_counter);
    % Mask upper triangle
    mask           = triu(true(size(grid_out)), 1);
    grid_out(mask) = NaN; % Set upper triangle to NaN
    grid_out       = flipud(grid_out);
    grid_out(num_nets+1,:) = NaN;
    grid_out(:,num_nets+1) = NaN;

    pcolor(grid_out);

    caxis([.15, 0.45]);

    set(gca,'YTick', 1:max(net_ids), 'YTickLabel', flip(network_names_kong));

    set(gca, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
    if scan_id == 1
        set(ax1, 'XTickLabel', {});
        xlabel(ax1,'')
    end
    ylabel('Networks');
    title(['Day ' num2str(scan_id)])
    colormap(ax1,whitered(256));
    axis square;

    fig_counter = fig_counter+1;

    % Next figure
    grid_out = vector_to_symmetric_grid(rest_fc_herit(:,scan_id),max(net_ids));

    ax2 = subplot(2,3,fig_counter);
    % Mask upper triangle
    mask           = triu(true(size(grid_out)), 1);
    grid_out(mask) = NaN; % Set upper triangle to NaN
    grid_out       = flipud(grid_out);
    grid_out(num_nets+1,:) = NaN;
    grid_out(:,num_nets+1) = NaN;
    % pcolor plot
    pcolor(grid_out);

    c       = colorbar;
    c.Ticks = [.15 .3 0.45];
    caxis([.15, 0.45]);

    tickLabels   = {'.15', '0.3', '0.45'};
    c.TickLabels = tickLabels;

    set(gca, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
    if scan_id == 1
        set(ax2, 'XTickLabel', {});
        xlabel(ax2,'')
    end
    set(ax2, 'YTickLabel', {});
    title(['Day ' num2str(scan_id)])
    colormap(ax2,whitered(256));

    axis square;

    fig_counter = fig_counter + 1;

    % Next figure
    grid_out = vector_to_symmetric_grid(fc_herit_task_diff(:,scan_id),max(net_ids));

    ax3 = subplot(2,3,fig_counter);
    % Mask upper triangle
    mask           = triu(true(size(grid_out)), 1);
    grid_out(mask) = NaN; % Set upper triangle to NaN
    grid_out       = flipud(grid_out);
    grid_out(num_nets+1,:) = NaN;
    grid_out(:,num_nets+1) = NaN;

    pcolor(grid_out);

    c = colorbar;
    c.Ticks = [-0.3,  0,  0.3];
    caxis([-0.3, 0.3]);

    % Optionally, set tick labels if desired
    tickLabels   = {'-0.3',  '0',  '0.3'};
    c.TickLabels = tickLabels;

    set(gca, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
    if scan_id == 1
        set(ax3, 'XTickLabel', {});
        xlabel(ax3,'')
    end
    set(ax3, 'YTickLabel', {});
    title(['Day ' num2str(scan_id)])
    colormap(ax3,bluewhitered(256));

    axis square;
    fig_counter = fig_counter + 1;

end
sgtitle('Movie FC h^2');

% Find all axes and colorbars in the figure
allAxes      = findall(gcf, 'Type', 'axes');
allColorbars = findall(gcf, 'Type', 'colorbar');

% Set line width and tick length for all axes
for i = 1:length(allAxes)
    set(allAxes(i), 'LineWidth', 2);      % Set line width to 2 for subplots
    set(allAxes(i), 'TickLength', [0 0]); % Set tick length to [0 0] for subplots
end

% Set line width for all colorbars
for i = 1:length(allColorbars)
    set(allColorbars(i), 'LineWidth', 2); % Set line width to 2 for colorbars
end

saveas(gcf,[data_dir '/isc_heritability/figs/fig3_heatmap_anatomical_fc_profile.fig'])
close all

% FIGURE S2


fig_counter = 1;
for scan_id = 1:num_days

    % First figure
    grid_out = vector_to_symmetric_grid(fc_movie_magnitude_anatomical{1,scan_id}.Var,max(net_ids));

    ax1 = subplot(2,3,fig_counter);
    % Mask upper triangle
    mask           = triu(true(size(grid_out)), 1);
    grid_out(mask) = NaN; % Set upper triangle to NaN
    grid_out       = flipud(grid_out);
    grid_out(num_nets+1,:) = NaN;
    grid_out(:,num_nets+1) = NaN;

    pcolor(grid_out);


    caxis([0, 0.75]);

    set(gca, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong, 'YTick', 1:max(net_ids), 'YTickLabel', flip(network_names_kong));
    if scan_id == 1
        set(ax1, 'XTickLabel', {});
        xlabel(ax1,'')
    end
    ylabel('Networks');
    title(['Day ' num2str(scan_id)])
    colormap(ax1,whitered(256));
    axis square;

    fig_counter = fig_counter+1;

    % Next figure
    grid_out = vector_to_symmetric_grid(fc_rest_magnitude_anatomical{1,scan_id}.Var,max(net_ids));

    ax2 = subplot(2,3,fig_counter);
    % Mask upper triangle
    mask           = triu(true(size(grid_out)), 1);
    grid_out(mask) = NaN; % Set upper triangle to NaN
    grid_out       = flipud(grid_out);
    grid_out(num_nets+1,:) = NaN;
    grid_out(:,num_nets+1) = NaN;

    pcolor(grid_out);

    c       = colorbar;
    c.Ticks = [0 .25 0.5 0.75];
    caxis([0, 0.75]);
    set(ax2, 'YTickLabel', {});

    tickLabels   = {'0', '0.25', '0.5','0.75'};
    c.TickLabels = tickLabels;
    set(gca, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
    if scan_id == 1
        set(ax2, 'XTickLabel', {});
        xlabel(ax2,'')
    end
    title(['Day ' num2str(scan_id)])
    colormap(ax2,whitered(256));

    axis square;

    fig_counter = fig_counter + 1;

    % Next figure
    grid_out = vector_to_symmetric_grid(fc_movie_magnitude_anatomical{1,scan_id}.Var-fc_rest_magnitude_anatomical{1,scan_id}.Var,max(net_ids));

    ax3 = subplot(2,3,fig_counter);
    % Mask upper triangle
    mask           = triu(true(size(grid_out)), 1);
    grid_out(mask) = NaN; % Set upper triangle to NaN
    grid_out       = flipud(grid_out);
    grid_out(num_nets+1,:) = NaN;
    grid_out(:,num_nets+1) = NaN;

    pcolor(grid_out);


    c       = colorbar;
    c.Ticks = [-0.6 -0.3,  0,  0.3 0.6];
    caxis([-0.6, 0.6]);

    % Optionally, set tick labels if desired
    tickLabels = {'-0.6','-0.3',  '0',  '0.3','0.6'};
    c.TickLabels = tickLabels;

    set(gca, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
    set(ax3, 'YTickLabel', {});
    if scan_id == 1
        set(ax3, 'XTickLabel', {});
        xlabel(ax3,'')
    end
    title(['Day ' num2str(scan_id)])
    colormap(ax3,bluewhitered(256));

    axis square;
    fig_counter = fig_counter + 1;

end
sgtitle('Movie FC h^2');

% Find all axes and colorbars in the figure
allAxes      = findall(gcf, 'Type', 'axes');
allColorbars = findall(gcf, 'Type', 'colorbar');

% Set line width and tick length for all axes
for i = 1:length(allAxes)
    set(allAxes(i), 'LineWidth', 2);      % Set line width to 2 for subplots
    set(allAxes(i), 'TickLength', [0 0]); % Set tick length to [0 0] for subplots
end

% Set line width for all colorbars
for i = 1:length(allColorbars)
    set(allColorbars(i), 'LineWidth', 2); % Set line width to 2 for colorbars
end
saveas(gcf,[data_dir '/isc_heritability/figs/figS2_heatmap_anatomical_fc_magnitude.fig'])
close all

%% Comparison figures

% FIGURE 5 SCATTER
x_extended           = linspace(-100, 1200, 300);
piecewise_fc_fits    = cell(num_days, 1);
connectivity_fc_fits = cell(num_days, 1);

for scan_id = 1:num_days
    subplot(2,1,scan_id)
    hold on

    % Calculate means and SEMs
    mean_piecewise = piecewise_herit_avg(:,scan_id);
    sem_piecewise  = piecewise_fc_se(:,scan_id);

    mean_connectivity = connectivity_herit_avg(:,scan_id);
    sem_connectivity  = connectivity_fc_se(:,scan_id);

    mean_movie_fc_herit = nanmean(movie_fc_herit(:,scan_id));
    sem_movie_fc_herit  = nanmean(movie_fc_se(:,scan_id));
    mean_piecewise      = [mean_piecewise;mean_movie_fc_herit];
    mean_connectivity   = [mean_connectivity;mean_movie_fc_herit];


    customModel        = fittype('a*x^b + c', 'independent', 'x', 'dependent', 'y');
    options            = fitoptions('Method', 'NonlinearLeastSquares');
    options.Lower      = [-Inf, -Inf, -Inf]; % Adjust lower bounds as necessary
    options.Upper      = [Inf, Inf, Inf]; % Adjust upper bounds as necessary
    options.StartPoint = [-.1, .1, .5]; % Adjust start points as necessary

    % Fit custom model for piecewise data
    [fitresult_piecewise, gof_piecewise] = fit([avg_areas;0], mean_piecewise, customModel, options);
    piecewise_fc_fits{scan_id}           = struct('FitResult', fitresult_piecewise, 'GOF', gof_piecewise);

    % Fit custom model for connectivity data
    [fitresult_connectivity, gof_connectivity] = fit([avg_areas;0], mean_connectivity, customModel, options);
    connectivity_fc_fits{scan_id}              = struct('FitResult', fitresult_connectivity, 'GOF', gof_connectivity);

    % Additional plotting and analysis code remains the same
    % Ensure to replace any instance of plotting or evaluating the fit with the new model

    % For example, calculating new y values for the extended x range using the fit results
    y_extended_piecewise    = feval(fitresult_piecewise, x_extended);
    y_extended_connectivity = feval(fitresult_connectivity, x_extended);

    % Plotting the extended power law lines
    plot(x_extended, y_extended_piecewise, 'LineWidth', 4, 'Color', [254/255 97/255 0/255]);
    plot(x_extended, y_extended_connectivity, 'LineWidth', 4, 'Color', [120/255 94/255 240/255]);

    % boundedline plots
    hl1 = boundedline(avg_areas, mean_piecewise(1:10), sem_piecewise, 'alpha', 'cmap', [254/255 97/255 0/255],'LineWidth',.00001);
    hl2 = boundedline(avg_areas, mean_connectivity(1:10), sem_connectivity, 'alpha', 'cmap', [120/255 94/255 240/255],'LineWidth',.00001);

    % Scatter plots for mean values
    scatter(avg_areas, mean_piecewise(1:10), 100, 'filled', 'MarkerFaceColor', [254/255 97/255 0/255], 'MarkerEdgeColor', 'k');
    scatter(avg_areas, mean_connectivity(1:10), 100, 'filled', 'MarkerFaceColor', [120/255 94/255 240/255], 'MarkerEdgeColor', 'k');
    scatter(0, mean_movie_fc_herit, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    % SEM plot for movie_fc_herit_mag (as a vertical line)
    errorbar(0, mean_movie_fc_herit, sem_movie_fc_herit, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4);

    % Axes settings
    ylim([.1 .4]); xlim([-100 1200]);
    xlabel('Average Parcel Area (mm^2)'); ylabel('Average Heritability (h^2)');
    yticks([.1 .2 .3 .4 ]); ax = gca; set(ax, 'TickDir', 'out'); set(ax, 'FontSize', 20); set(ax, 'LineWidth', 3); set(ax, 'Layer', 'bottom');
    xticks([0 400 800 1200]);
    grid off; legend({'Response Hyperalignment', 'Connectivity Hyperalignment'}, 'Location', 'best'); pbaspect([2 1 1]); set(gcf, 'position', [100,100,600,1200]);

    piecewise_perc(1,scan_id)    = 1-(mean_piecewise(1)/mean_movie_fc_herit);
    connectivity_perc(1,scan_id) = 1-(mean_connectivity(1)/mean_movie_fc_herit);

    piecewise_conf(1,scan_id)    = 1-(mean_piecewise(1)/mean_movie_fc_herit - (sem_piecewise(1)/mean_movie_fc_herit*1.96));
    piecewise_conf(2,scan_id)    = 1-(mean_piecewise(1)/mean_movie_fc_herit +(sem_piecewise(1)/mean_movie_fc_herit*1.96));
    connectivity_conf(1,scan_id) = 1-(mean_connectivity(1)/mean_movie_fc_herit - (sem_connectivity(1)/mean_movie_fc_herit*1.96));
    connectivity_conf(2,scan_id) = 1-(mean_connectivity(1)/mean_movie_fc_herit +(sem_connectivity(1)/mean_movie_fc_herit*1.96));
end

saveas(gcf,[data_dir '/isc_heritability/figs/fig5_areascatter_hyperalignment_fc_profile_linear.fig'])
close all

% FIGURE 5 HEATMAP

figure;
fig_counter = 1;

for scan_id = 1:num_days
    for col = 1:2
        % Determine the grid_out based on column
        if col == 1
            grid_out = vector_to_symmetric_grid(connectivity_herit{1,1}(:,scan_id)-movie_fc_herit(:,scan_id),max(net_ids));
        else
            grid_out = vector_to_symmetric_grid(piecewise_herit{1,1}(:,scan_id)-movie_fc_herit(:,scan_id),max(net_ids));
        end

        % Create subplot
        ax = subplot(2, 2, fig_counter);

        mask           = triu(true(size(grid_out)), 1);
        grid_out(mask) = NaN; % Set upper triangle to NaN
        grid_out       = flipud(grid_out);
        grid_out(num_nets+1,:) = NaN;
        grid_out(:,num_nets+1) = NaN;

        % pcolor plot
        pcolor(grid_out);
        colormap(ax, bluewhitered(256));
        axis square;
        caxis([-.3 .3]);
        set(ax, 'LineWidth', 2, 'TickLength', [0 0], 'FontName', 'Arial');  % Increase line width, remove tick marks, and set font

        % Set X and Y tick labels where necessary
        if col == 1  % First column
            set(ax, 'YTick', 1:max(net_ids), 'YTickLabel', flip(network_names_kong), 'XTick', [], 'XTickLabel', {});
        else
            set(ax, 'YTick', [], 'YTickLabel', {}, 'XTick', [], 'XTickLabel', {});
        end
        if scan_id == 2  % Second row
            set(ax, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
        end

        % Add colorbar for second column
        if col == 2
            c            = colorbar;
            c.Ticks      = [-.3 -.15 0 .15 .3];
            c.TickLabels = {'-.3', '-.15', '0', '.15', '.3'};
            set(c, 'LineWidth', 2, 'FontName', 'Arial'); % Increase colorbar line width and set font
        end

        % Set title
        title(['Day ' num2str(scan_id)]);

        % Increment figure counter
        fig_counter = fig_counter + 1;
    end
end

sgtitle('Movie FC h^2', 'FontName', 'Arial');

saveas(gcf,[data_dir '/isc_heritability/figs/fig5_heatmap_piecewiseconnectivity_fc_profile.fig'])
close all



% FIGURE S3 SCATTER
for scan_id = 1:num_days
    subplot(2,1,scan_id)
    hold on

    % Calculate means and SEMs
    mean_piecewise = squeeze(nanmean(piecewise_fc_magnitude(:,:,scan_id)));
    sem_piecewise  = piecewise_mag_se(:,scan_id);

    mean_connectivity = squeeze(nanmean(connectivity_fc_magnitude(:,:,scan_id)));
    sem_connectivity  = connectivity_mag_se(:,scan_id);

    mean_movie_fc_herit = nanmean(movie_fc_herit_mag(:,scan_id));
    sem_movie_fc_herit  = nanmean(movie_mag_se(:,scan_id));
    mean_piecewise      = [mean_piecewise,mean_movie_fc_herit];
    mean_connectivity   = [mean_connectivity,mean_movie_fc_herit];
    customModel         = fittype('a*x^b + c', 'independent', 'x', 'dependent', 'y');

    options            = fitoptions('Method', 'NonlinearLeastSquares');
    options.Lower      = [-Inf, -Inf, -Inf]; % Adjust lower bounds as necessary
    options.Upper      = [Inf, Inf, Inf]; % Adjust upper bounds as necessary
    options.StartPoint = [-.05, .1, .5]; % Adjust start points as necessary

    % Fit custom model for piecewise data
    [fitresult_piecewise, gof_piecewise] = fit([avg_areas;0], mean_piecewise', customModel, options);
    piecewise_fc_mag_fits{scan_id}       = struct('FitResult', fitresult_piecewise, 'GOF', gof_piecewise);

    % Fit custom model for connectivity data
    [fitresult_connectivity, gof_connectivity] = fit([avg_areas;0], mean_connectivity', customModel, options);
    connectivity_fc_mag_fits{scan_id}          = struct('FitResult', fitresult_connectivity, 'GOF', gof_connectivity);

    % Additional plotting and analysis code remains the same
    % Ensure to replace any instance of plotting or evaluating the fit with the new model

    % For example, calculating new y values for the extended x range using the fit results
    y_extended_piecewise    = feval(fitresult_piecewise, x_extended);
    y_extended_connectivity = feval(fitresult_connectivity, x_extended);

    % Plotting the extended power law lines
    plot(x_extended, y_extended_piecewise, 'LineWidth', 4, 'Color', [254/255 97/255 0/255]);
    plot(x_extended, y_extended_connectivity, 'LineWidth', 4, 'Color', [120/255 94/255 240/255]);

    % boundedline plots
    hl1 = boundedline(avg_areas, mean_piecewise(1:10), sem_piecewise, 'alpha', 'cmap', [254/255 97/255 0/255],'LineWidth',.00001);
    hl2 = boundedline(avg_areas, mean_connectivity(1:10), sem_connectivity, 'alpha', 'cmap', [120/255 94/255 240/255],'LineWidth',.00001);

    % Scatter plots for mean values
    scatter(avg_areas, mean_piecewise(1:10), 100, 'filled', 'MarkerFaceColor', [254/255 97/255 0/255], 'MarkerEdgeColor', 'k');
    scatter(avg_areas, mean_connectivity(1:10), 100, 'filled', 'MarkerFaceColor', [120/255 94/255 240/255], 'MarkerEdgeColor', 'k');
    scatter(0, mean_movie_fc_herit, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    % SEM plot for movie_fc_herit_mag (as a vertical line)
    errorbar(0, mean_movie_fc_herit, sem_movie_fc_herit, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4);

    % Axes settings
    ylim([.2 .5]); xlim([-100 1200]);
    xticks([0 400 800 1200]);
    xlabel('Average Parcel Area (mm^2)'); ylabel('Average Heritability (h^2)');
    yticks([.2 .3 .4 .5]); ax = gca; set(ax, 'TickDir', 'out'); set(ax, 'FontSize', 20); set(ax, 'LineWidth', 3); set(ax, 'Layer', 'bottom');
    grid off; legend({'Response Hyperalignment', 'Connectivity Hyperalignment'}, 'Location', 'best'); pbaspect([2 1 1]); set(gcf, 'position', [100,100,600,1200]);
    piecewise_mag_perc(1,scan_id)    = 1-(mean_piecewise(1)/mean_movie_fc_herit);
    connectivity_mag_perc(1,scan_id) = 1-(mean_connectivity(1)/mean_movie_fc_herit);

    piecewise_mag_conf(1,scan_id)    = 1-(mean_piecewise(1)/mean_movie_fc_herit - (sem_piecewise(1)/mean_movie_fc_herit*1.96));
    piecewise_mag_conf(2,scan_id)    = 1-(mean_piecewise(1)/mean_movie_fc_herit +(sem_piecewise(1)/mean_movie_fc_herit*1.96));
    connectivity_mag_conf(1,scan_id) = 1-(mean_connectivity(1)/mean_movie_fc_herit - (sem_connectivity(1)/mean_movie_fc_herit*1.96));
    connectivity_mag_conf(2,scan_id) = 1-(mean_connectivity(1)/mean_movie_fc_herit +(sem_connectivity(1)/mean_movie_fc_herit*1.96));
end

saveas(gcf,[data_dir '/isc_heritability/figs/fig5_areascatter_hyperalignment_fc_magnitude_linear.fig'])
close all

% FIGURE S3 HEATMAP

figure;
fig_counter = 1;
clear grid_out
for scan_id = 1:num_days
    for col = 1:2
        % Determine the grid_out based on column
        if col == 1
            grid_out = vector_to_symmetric_grid(connectivity_fc_magnitude(:,1,scan_id)-movie_fc_herit_mag(:,scan_id),max(net_ids));
        else
            grid_out = vector_to_symmetric_grid(piecewise_fc_magnitude(:,1,scan_id)-movie_fc_herit_mag(:,scan_id),max(net_ids));
        end

        % Create subplot
        ax = subplot(2, 2, fig_counter);

        % Mask upper triangle
        mask           = triu(true(size(grid_out)), 1);
        grid_out(mask) = NaN; % Set upper triangle to NaN
        grid_out       = flipud(grid_out);
        grid_out(num_nets+1,:) = NaN;
        grid_out(:,num_nets+1) = NaN;

        pcolor(grid_out);

        colormap(ax, bluewhitered(256));
        axis square;
        caxis([-.6 .6]);
        set(ax, 'LineWidth', 2, 'TickLength', [0 0], 'FontName', 'Arial', 'FontSize', 16);  % Increase line width, remove tick marks, set font, and increase font size

        % Set X and Y tick labels where necessary
        if col == 1  % First column
            set(ax, 'YTick', 1:max(net_ids), 'YTickLabel', flip(network_names_kong), 'XTick', [], 'XTickLabel', {});
        else
            set(ax, 'YTick', [], 'YTickLabel', {}, 'XTick', [], 'XTickLabel', {});
        end
        if scan_id == 2  % Second row
            set(ax, 'XTick', 1:max(net_ids), 'XTickLabel', network_names_kong);
        end

        % Add colorbar for second column
        if col == 2
            c            = colorbar;
            c.Ticks      = [-.6 -.3  0 .3 .6];
            c.TickLabels = {'-.6','-.3', '0', '.3','.6'};
            set(c, 'LineWidth', 2, 'FontName', 'Arial', 'FontSize', 16); % Increase colorbar line width, set font, and increase font size
        end

        % Set title
        title(['Day ' num2str(scan_id)], 'FontSize', 16);
        % Increment figure counter
        fig_counter = fig_counter + 1;
    end
end

sgtitle('Movie FC h^2', 'FontName', 'Arial');

saveas(gcf,[data_dir '/isc_heritability/figs/fig5_heatmap_piecewiseconnectivity_fc_magnitude.fig'])
close all


%% Check magnitude difference perm testing
load([data_dir '/isc_heritability/data/solar/fc_perm/results_array_task.mat'])

rest_fc_herit_mag_square_mean  = squeeze(nanmean(rest_fc_herit_mag_square,1));
movie_fc_herit_mag_square_mean = squeeze(nanmean(movie_fc_herit_mag_square,1));

task_diff_fc_herit_mag_square_mean = movie_fc_herit_mag_square_mean - rest_fc_herit_mag_square_mean;
results_array_task_square          = zeros(num_nets,num_nets,1000,num_days);
for scan_id = 1:num_days
    for perm_id = 1:1000
        results_array_task_square(:,:,perm_id,scan_id) = vector_to_symmetric_grid(results_array_task(:,perm_id,scan_id),num_nets);
    end
end

results_array_task_square_mean = squeeze(nanmean(results_array_task_square,1));

results_array_task_pvalue = zeros(num_nets,num_days);
for scan_id = 1:num_days
    for net_id = 1:num_nets
        results_array_task_pvalue(net_id,scan_id) = sum((abs(results_array_task_square_mean(net_id,:,scan_id)) > task_diff_fc_herit_mag_square_mean(net_id,scan_id)))/1000;
    end
end

results_array_task_pvalue_fdr = fdr_bh(results_array_task_pvalue);




