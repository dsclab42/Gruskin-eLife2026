%%%%%%%%
% Author: David Gruskin
% Contact: dcg2153@cumc.columbia.edu
% Project: Heritability of movie-evoked brain activity and connectivity
% Description: This is the main analysis and visualization script for the
% ISC component of the HCP 7T ISC heritability project.
%%%%%%%%

%% Add and set relevant paths
desktop_dir   = '/path/to/desktop';    
downloads_dir = '/path/to/downloads';  
data_dir      = '/path/to/data';       
% ------------------------------------------------------------------------

addpath([data_dir '/isc_heritability/'])
addpath([downloads_dir '/github_repo/bounded_lines'])
addpath([downloads_dir '/github_repo/inpaint_nans'])
addpath(desktop_dir)
addpath([data_dir '/restsync'])
addpath([downloads_dir '/npy-matlab-master/npy-matlab'])
addpath([data_dir '/isc_heritability/matlab_scripts/h2_multi-master/h2_multi/'])
addpath([data_dir '/restsync/inputs/scripts'])
addpath([downloads_dir '/gifti-main-2/'])
addpath([data_dir '/restsync/inputs/scripts/cifti-matlab-master'])
addpath(fullfile(fileparts(mfilename('fullpath')), 'dependencies'))  % helper functions shipped in matlab_code/dependencies
addpath([desktop_dir '/isc_heritability_r'])
addpath([downloads_dir '/palm-alpha119'])

schaefer_400_dscalar      = [desktop_dir '/Schaefer2018_400Parcels_17Networks_order.dscalar.nii'];
schaefer_400_pscalar      = [desktop_dir '/Schaefer2018_400Parcels_17Networks_order.pscalar.nii'];
schaefer_400_pscalar_kong = [desktop_dir '/Schaefer2018_400Parcels_Kong2022_17Networks_order.pscalar.nii'];
schaefer_100_pscalar_kong = [desktop_dir '/Schaefer2018_100Parcels_Kong2022_17Networks_order.pscalar.nii'];
schaefer_1000_pscalar_kong = [desktop_dir '/Schaefer2018_1000Parcels_Kong2022_17Networks_order.pscalar.nii'];

schaefer_400_label = [desktop_dir '/Schaefer2018_400Parcels_17Networks_order.dlabel.nii'];
pscalar_template   = [desktop_dir '/schaefer_400_round2.pscalar.nii'];
medial_mask_path   = [data_dir '/isc_heritability_review/mat_files/supporting_files/medial_mask.csv'];

pheno_table_path            = [data_dir '/restsync/inputs/accessory_files/hcp_7t_pheno.csv'];
pheno_restricted_table_path = [downloads_dir '/inputs/RESTRICTED_<your_hcp_file>.csv'];  % HCP restricted behavioral/kinship CSV; supply your own (filename differs per download)

base_dir = [data_dir '/isc_heritability/data/'];
wb_path  = '/path/to/workbench/bin/wb_command'; 

%% Load some preliminary files
load([data_dir '/isc_heritability_review/mat_files/supporting_files/covari_motion.mat'])
load([data_dir '/isc_heritability_review/mat_files/supporting_files/family_ids.mat'])
load kong_transform.mat

scalar_template = schaefer_400_dscalar;

% align_types to loop over
align_types = {'piecewise';'anatomical'};

% Set up medial mask variables
medial_mask    = csvread(medial_mask_path);
medial_mask_lh = medial_mask(1:32492);
medial_mask_rh = medial_mask(32493:end);

schaefer_dscalar        = ciftiopen(schaefer_400_dscalar);
vox_ids                 = schaefer_dscalar.cdata;
vox_ids(medial_mask==0) = [];
num_trs                 = [1432, 1409];

% Set some basic variables and functions
parc_reses = {'100','200','300','400','500','600','700','800','900','1000'};

hvec        = @(x) x(triu(true(size(x)), 1));
isemptycell = @(x) cellfun(@isempty, x);

num_vox_cortex = 59412;
num_parcs      = 10;
num_days       = 2;
num_fams       = 90;
num_subj       = 184;
num_perm       = 10000;

exclude_rest               = [8,66,124,126,130,135,179,183];
subjs_rest                 = 1:num_subj;
subjs_rest(exclude_rest)   = [];
exclude_movie              = [66,124,126,130,135,183];
subjs_movie                = 1:num_subj;
subjs_movie(exclude_movie) = [];

num_subj_movie = length(subjs_movie);
num_subj_rest  = length(subjs_rest);

pheno_table                = readtable(pheno_table_path);
pheno_restricted           = readtable(pheno_restricted_table_path);
pheno_restricted.Family_ID = string(pheno_restricted.Family_ID);

% Get unique family IDs
unique_famids = unique(pheno_restricted.Family_ID);

% Create a mapping (container) from the unique family IDs to a sequence of integers
famid_mapping = containers.Map(unique_famids, 1:length(unique_famids));

% Apply the mapping to the Family_ID column
% Vectorized approach for efficiency
mapped_famids = arrayfun(@(x) famid_mapping(x), pheno_restricted.Family_ID);

% Update the Family_ID column with the mapped values
pheno_restricted.Family_ID = mapped_famids;
movie_fam_ids              = pheno_restricted.Family_ID(subjs_movie);
rest_fam_ids               = pheno_restricted.Family_ID(subjs_rest);

pheno_restricted_subset          = pheno_restricted(subjs_movie,:);
pheno_restricted_subset.Zygosity = pheno_restricted_subset.ZygosityGT;
pheno_restricted.Zygosity        = pheno_restricted.ZygosityGT;

pheno_restricted_subset.Zygosity(isemptycell(pheno_restricted_subset.Zygosity)) = pheno_restricted_subset.ZygositySR(isemptycell(pheno_restricted_subset.Zygosity));
pheno_restricted.Zygosity(isemptycell(pheno_restricted.Zygosity))               = pheno_restricted.ZygositySR(isemptycell(pheno_restricted.Zygosity));

mztwin_new = cell(num_subj,1);
dztwin_new = cell(num_subj,1);

for subj1 = 1:num_subj
    for subj2 = 1:num_subj
        if subj1~=subj2
            if pheno_restricted.Family_ID(subj1) == pheno_restricted.Family_ID(subj2) && strcmp(pheno_restricted.Zygosity(subj1),'MZ') && strcmp(pheno_restricted.Zygosity(subj2),'MZ')
                mztwin_new{subj1} = pheno_restricted.Family_ID(subj1);
            elseif pheno_restricted.Family_ID(subj1) == pheno_restricted.Family_ID(subj2) && strcmp(pheno_restricted.Zygosity(subj1),'DZ') && strcmp(pheno_restricted.Zygosity(subj2),'DZ')
                dztwin_new{subj1} = pheno_restricted.Family_ID(subj1);
            end
        end
    end
end

flat_array_mz = vertcat(mztwin_new{:});

% Find unique elements and sort them
[unique_elements, ~, idx] = unique(flat_array_mz);

% Create a mapping (integer representation for each unique element)
% Starting index from 1
mapping = 1:length(unique_elements);

% Replace elements in mztwin_new with their corresponding integers
% Initialize the output array with NaNs
MZTWIN = nan(num_subj, 1);

for i = 1:length(mztwin_new)
    if ~isempty(mztwin_new{i})
        MZTWIN(i) = nanmean(mapping(idx(flat_array_mz == mztwin_new{i})));
    end
end

flat_array_dz = vertcat(dztwin_new{:});

% Find unique elements and sort them
[unique_elements, ~, idx] = unique(flat_array_dz);

% Create a mapping (integer representation for each unique element)
% Starting index from 1
mapping = 1:length(unique_elements);

% Replace elements in mztwin_new with their corresponding integers
% Initialize the output array with NaNs
DZTWIN = nan(num_subj, 1);

for i = 1:length(dztwin_new)
    if ~isempty(dztwin_new{i})
        DZTWIN(i) = nanmean(mapping(idx(flat_array_dz == dztwin_new{i})));
    end
end

%% Get dyad ids

kinship = zeros(num_subj,num_subj);
for subj1 = 1:num_subj
    for subj2 = 1:num_subj
        if subj1 == subj2
            kinship(subj1,subj2) = 1;
        elseif MZTWIN(subj1) == MZTWIN(subj2)
            kinship(subj1,subj2) = 1;
        elseif family_ids(subj1,2) == family_ids(subj2,2)
            kinship(subj1,subj2) = 0.5;
        end
    end
end

gender_var = grp2idx(categorical(pheno_table.Gender));
ages       = pheno_restricted.Age_in_Yrs;

% Initialize the arrays and counters
pairwise_mz_id = zeros(1,1);
pairwise_dz_id = zeros(1,1);
pairwise_nz_id = zeros(1,1);

processedPairs = false(num_subj_movie, num_subj_movie); % Logical array to track processed pairs

dz_counter    = 1;
mz_counter    = 1;
nz_counter    = 1;
total_counter = 1;

for subj1 = 1:num_subj_movie
    for subj2 = 1:num_subj_movie
        if subj1 ~= subj2
            % Check if this pair has already been processed
            if ~processedPairs(subj1, subj2) && ~processedPairs(subj2, subj1)
                % Process the pair
                if DZTWIN(subjs_movie(subj1)) == DZTWIN(subjs_movie(subj2)) && ~isnan(DZTWIN(subjs_movie(subj2)))
                    pairwise_dz_id(dz_counter,1) = total_counter;
                    dz_counter                   = dz_counter + 1;
                elseif MZTWIN(subjs_movie(subj1)) == MZTWIN(subjs_movie(subj2)) && ~isnan(MZTWIN(subjs_movie(subj1)))
                    pairwise_mz_id(mz_counter,1) = total_counter;
                    mz_counter                   = mz_counter + 1;
                elseif gender_var(subjs_movie(subj1)) == gender_var(subjs_movie(subj2)) && abs(ages(subjs_movie(subj1))-ages(subjs_movie(subj2)))<1 && family_ids(subjs_movie(subj1)) ~= family_ids(subjs_movie(subj2))
                    pairwise_nz_id(nz_counter,1) = total_counter;
                    nz_counter                   = nz_counter + 1;
                end
                % Mark this pair as processed
                processedPairs(subj1, subj2) = true;
                processedPairs(subj2, subj1) = true;
            end
        end
        % Increment total counter regardless of whether the pair was processed
        total_counter = total_counter + 1;
    end
end

% Initialize the arrays and counters
pairwise_mz_id_rest = zeros(1,1);
pairwise_dz_id_rest = zeros(1,1);
pairwise_nz_id_rest = zeros(1,1);

processedPairs = false(num_subj_rest, num_subj_rest); % Logical array to track processed pairs

dz_counter    = 1;
mz_counter    = 1;
nz_counter    = 1;
total_counter = 1;

for subj1 = 1:num_subj_rest
    for subj2 = 1:num_subj_rest
        if subj1 ~= subj2
            % Check if this pair has already been processed
            if ~processedPairs(subj1, subj2) && ~processedPairs(subj2, subj1)
                % Process the pair
                if DZTWIN(subjs_movie(subj1)) == DZTWIN(subjs_movie(subj2)) && ~isnan(DZTWIN(subjs_movie(subj2)))
                    pairwise_dz_id_rest(dz_counter,1) = total_counter;
                    dz_counter                        = dz_counter + 1;
                elseif MZTWIN(subjs_movie(subj1)) == MZTWIN(subjs_movie(subj2)) && ~isnan(MZTWIN(subjs_movie(subj1)))
                    pairwise_mz_id_rest(mz_counter,1) = total_counter;
                    mz_counter                        = mz_counter + 1;
                elseif gender_var(subjs_movie(subj1)) == gender_var(subjs_movie(subj2)) && abs(ages(subjs_movie(subj1))-ages(subjs_movie(subj2)))<1 && family_ids(subjs_movie(subj1)) ~= family_ids(subjs_movie(subj2))
                    pairwise_nz_id_rest(nz_counter,1) = total_counter;
                    nz_counter                        = nz_counter + 1;
                end
                % Mark this pair as processed
                processedPairs(subj1, subj2) = true;
                processedPairs(subj2, subj1) = true;
            end
        end
        % Increment total counter regardless of whether the pair was processed
        total_counter = total_counter + 1;
    end
end

age_mean = mean(ages(subjs_movie));
age_std  = std(ages(subjs_movie));
age_min  = min(ages(subjs_movie));
age_max  = max(ages(subjs_movie));

hisp_perc    = sum(grp2idx(categorical(pheno_restricted.Ethnicity(subjs_movie)))==1)/num_subj_movie;
white_perc   = sum(grp2idx(categorical(pheno_restricted.Race(subjs_movie)))==4)/num_subj_movie;
black_perc   = sum(grp2idx(categorical(pheno_restricted.Race(subjs_movie)))==2)/num_subj_movie;
asian_perc   = sum(grp2idx(categorical(pheno_restricted.Race(subjs_movie)))==1)/num_subj_movie;
unknown_perc = sum(grp2idx(categorical(pheno_restricted.Race(subjs_movie)))==3)/num_subj_movie;%

nonEmptyCellsCount = sum(~cellfun(@isempty, pheno_restricted.ZygosityGT(subjs_movie)));


%% Run anatomical ISC analyses
if exist(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_herit.mat'),'file') == 2
    load(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_herit.mat'))
    load(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_herit_jack.mat'))
    load(strcat(base_dir,'anatomical/outputs/parc/anatomical_isc_herit_parc.mat'))
    load(strcat(base_dir,'anatomical/outputs/parc/anatomical_isc_herit_parc_perm.mat'))
else
    anatomical_isc_herit      = zeros(num_vox_cortex,2); anatomical_isc_herit_perm = zeros(num_vox_cortex,2); anatomical_isc_herit_jack = zeros(num_vox_cortex,length(unique(movie_fam_ids)),2);
    anatomical_isc_herit_parc = cell(10,1); anatomical_isc_herit_parc_perm = cell(10,1); anatomical_isc_herit_parc_jack = cell(10,1);
    for scan_id = 1:num_days
        [anatomical_isc_herit(:,scan_id), anatomical_isc_herit_perm(:,scan_id),anatomical_isc_herit_jack(:,:,scan_id)]  = runISCHeritability('400', base_dir, kinship, subjs_movie,'anatomical','gray',pairwise_mz_id,pairwise_dz_id,pairwise_nz_id,scan_id,covari_motion,num_perm,[]);
        for parc_id = 4
            curr_parc = parc_reses{parc_id};
            [anatomical_isc_herit_parc{parc_id,1}(:,scan_id), anatomical_isc_herit_parc_perm{parc_id,1}(:,scan_id),anatomical_isc_herit_parc_jack{parc_id,1}(:,:,scan_id)]  = runISCHeritability(curr_parc, base_dir, kinship, subjs_movie,'anatomical','parc',pairwise_mz_id,pairwise_dz_id,pairwise_nz_id,scan_id,covari_motion,num_perm,movie_fam_ids);
        end
    end
    save(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_herit.mat'),'anatomical_isc_herit')
    save(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_herit_jack.mat'),'anatomical_isc_herit_jack')
    save(strcat(base_dir,'anatomical/outputs/parc/anatomical_isc_herit_parc.mat'),'anatomical_isc_herit_parc')
    save(strcat(base_dir,'anatomical/outputs/parc/anatomical_isc_herit_parc_perm.mat'),'anatomical_isc_herit_parc_perm')

end

% Get % of parcels significant on both days and mean/std
% heritability
isc_heritability_parc_pValues      = anatomical_isc_herit_parc_perm{4,1}(:,:);
isc_heritability_parc_pValues_fdr  = zeros(size(isc_heritability_parc_pValues,1),num_days);
for scan_id = 1:num_days
    isc_heritability_parc_pValues_fdr(:,scan_id)  = fdr_bh(isc_heritability_parc_pValues(:,scan_id));
end

isc_heritability_parc_pValues_perc = sum(isc_heritability_parc_pValues_fdr(:,1).*isc_heritability_parc_pValues_fdr(:,2))/400;
isc_heritability_parc_mean         = mean(anatomical_isc_herit_parc{4,1});
isc_heritability_parc_std          = std(anatomical_isc_herit_parc{4,1});

% For paper, get correlation between heritability maps and heritability/ISC
isc_heritability_parc_trt = corr(anatomical_isc_herit_parc{4,1}(:,:),'type','spearman','rows','complete');

if num_perm == 10000
    save(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_herit_perm.mat'),'anatomical_isc_herit_perm')
    save(strcat(base_dir,'anatomical/outputs/parc/anatomical_isc_herit_parc.mat'),'anatomical_isc_herit_parc_perm')
end

[p_val,anatomical_se(:,1),anatomical_se(:,2)] = compareHeritability(nanmean(anatomical_isc_herit(:,1)),nanmean(anatomical_isc_herit_jack(:,:,1),1),nanmean(anatomical_isc_herit(:,2)),nanmean(anatomical_isc_herit_jack(:,:,2),1));

% Correlate heritability and ISC maps
anat_isc_parc   = zeros(400,num_days);
herit_isc_parc  = zeros(400,num_days);
herit_isc_corr  = zeros(1,num_days);
isc_herit_resid = zeros(400,num_days);

for scan_id = 1:num_days
    data                     = ciftiopen(strcat(base_dir,'anatomical/outputs/parc/anatomical_isc_parc_400_scan_',num2str(scan_id),'.pscalar.nii'));
    anat_isc_parc(:,scan_id) = data.cdata;

    data                      = ciftiopen(strcat(base_dir,'anatomical/outputs/parc/anatomical_isc_herit_parc_400_scan_',num2str(scan_id),'.pscalar.nii'));
    herit_isc_parc(:,scan_id) = data.cdata;

    big_ids   = find(anat_isc_parc(:,scan_id)>median(anat_isc_parc(:,scan_id)));
    small_ids = find(anat_isc_parc(:,scan_id)<median(anat_isc_parc(:,scan_id)));

    herit_isc_corr(1,scan_id) = corr(anat_isc_parc(:,scan_id),herit_isc_parc(:,scan_id),'type','spearman');

    mdl                        = fitlm(anat_isc_parc(:,scan_id),herit_isc_parc(:,scan_id));
    isc_herit_resid(:,scan_id) = mdl.Residuals{:,1};
    saveCifti(isc_herit_resid(:,scan_id), strcat([data_dir '/isc_heritability/data/anatomical/outputs/parc/anatomical_isc_herit_raw_residuals_parc_'],curr_parc,'_day',num2str(scan_id),'.pscalar.nii'), wb_path, schaefer_400_pscalar_kong,medial_mask);
end


%% Run piecewise ISC analyses
parpool('Threads')

if exist(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit.mat'),'file') == 2
    load(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit.mat'))
    load(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit_jack.mat'))
else
    piecewise_isc_herit      = zeros(num_vox_cortex,num_parcs,num_days);
    piecewise_isc_herit_perm = zeros(num_vox_cortex,num_parcs,num_days);
    piecewise_isc_herit_jack = zeros(num_vox_cortex,num_fams,num_parcs,num_days);

    for parc_id = 1:length(parc_reses)
        curr_parc = parc_reses{parc_id};
        for scan_id = 1:num_days
            [piecewise_isc_herit(:,parc_id,scan_id), piecewise_isc_herit_perm(:,parc_id,scan_id),piecewise_isc_herit_jack(:,:,parc_id,scan_id)]  = runISCHeritability(curr_parc, base_dir, kinship, subjs_movie,'piecewise','gray',pairwise_mz_id,pairwise_dz_id,pairwise_nz_id,scan_id,covari_motion,num_perm,movie_fam_ids);
            save(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit_jack.mat'),'piecewise_isc_herit_jack')
        end
    end
    save(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit.mat'),'piecewise_isc_herit')
    save(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit_jack.mat'),'piecewise_isc_herit_jack')
end


if num_perm == 10000
    save(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit_perm.mat'),'piecewise_isc_herit_perm')
end

piecewise_isc_se = zeros(num_parcs,num_days);
for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        [~,~,piecewise_isc_se(parc_id,scan_id)] = compareHeritability(0,nanmean(piecewise_isc_herit_jack(:,:,parc_id,scan_id),1),0,nanmean(piecewise_isc_herit_jack(:,:,parc_id,scan_id),1));
    end
end
%% Run connectivity ISC analyses

if exist(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_herit.mat'),'file') == 2
    load(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_herit.mat'))
    load(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_herit_jack.mat'))
else
    connectivity_isc_herit      = zeros(num_vox_cortex,num_parcs,num_days);
    connectivity_isc_herit_perm = zeros(num_vox_cortex,num_parcs,num_days);
    connectivity_isc_herit_jack = zeros(num_vox_cortex,num_fams,num_parcs,num_days);
    for parc_id = 1:num_parcs
        curr_parc = parc_reses{parc_id};
        for scan_id = 1:num_days
            [connectivity_isc_herit(:,parc_id,scan_id), connectivity_isc_herit_perm(:,parc_id,scan_id),connectivity_isc_herit_jack(:,:,parc_id,scan_id)]  = runISCHeritability(curr_parc, base_dir, kinship, subjs_rest,'connectivity','gray',pairwise_mz_id_rest,pairwise_dz_id_rest,pairwise_nz_id_rest,scan_id,covari_motion,num_perm,rest_fam_ids);
            save(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_herit_jack.mat'),'connectivity_isc_herit_jack')
        end
    end
    save(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_herit.mat'),'connectivity_isc_herit')
    save(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_herit_jack.mat'),'connectivity_isc_herit_jack')
end

if num_perm == 10000
    save(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_herit_perm.mat'),'connectivity_isc_herit_perm')
end

connectivity_isc_se = zeros(num_parcs,num_days);
for parc_id = 1:num_parcs
    for scan_id = 1:num_days
        [~,~,connectivity_isc_se(parc_id,scan_id)] = compareHeritability(0,nanmean(connectivity_isc_herit_jack(:,:,parc_id,scan_id),1),0,nanmean(connectivity_isc_herit_jack(:,:,parc_id,scan_id),1));
    end
end

%% Make difference CIFTIs

curr_parc = '400';
for scan_id = 1:num_days
    anatomical_img   = ciftiopen(strcat([data_dir '/isc_heritability/data/anatomical/outputs/gray/anatomical_isc_herit_parc_'],curr_parc,'_day',num2str(scan_id),'.pscalar.nii'));
    piecewise_img    = ciftiopen(strcat([data_dir '/isc_heritability/data/piecewise/outputs/gray/piecewise_isc_herit_parc_'],curr_parc,'_day',num2str(scan_id),'.pscalar.nii'));
    connectivity_img = ciftiopen(strcat([data_dir '/isc_heritability/data/connectivity/outputs/gray/connectivity_isc_herit_parc_'],curr_parc,'_day',num2str(scan_id),'.pscalar.nii'));

    saveCifti(piecewise_img.cdata-anatomical_img.cdata, strcat([data_dir '/isc_heritability/data/piecewise/outputs/parc/piecewise_anatomical_isc_herit_diff_parc_'],curr_parc,'_day',num2str(scan_id),'.pscalar.nii'), wb_path, schaefer_400_pscalar,medial_mask);
    saveCifti(connectivity_img.cdata-anatomical_img.cdata, strcat([data_dir '/isc_heritability/data/connectivity/outputs/parc/connectivity_anatomical_isc_herit_diff_parc_'],curr_parc,'_day',num2str(scan_id),'.pscalar.nii'), wb_path, schaefer_400_pscalar,medial_mask);
end

%% Get average area per parcel

right_areas = gifti([data_dir '/isc_heritability/mat_files/supporting_files/right_areas.func.gii']);
right_areas = right_areas.cdata;

left_areas = gifti([data_dir '/isc_heritability/mat_files/supporting_files/left_areas.func.gii']);
left_areas = left_areas.cdata;

all_areas = cat(1,left_areas,right_areas);
all_areas = all_areas(medial_mask==1);

parc_areas = cell(num_parcs,1);
parc_ids   = zeros(num_vox_cortex,num_parcs);

avg_areas = zeros(num_parcs,1);
for parc_id = 1:length(parc_reses)
    parc_res = parc_reses{parc_id};

    current_parc        = strcat('Schaefer2018_',parc_res,'Parcels_17Networks_order.dscalar.nii');
    parc_ids_temp       = ciftiopen(strcat(desktop_dir, '/', current_parc),wb_path);
    parc_ids(:,parc_id) = parc_ids_temp.cdata(medial_mask==1);
    for parc_val = 1:max(parc_ids(:,parc_id))
        parc_areas{parc_id,1}(parc_val,1) = sum(all_areas(parc_ids(:,parc_id)==parc_val));
    end
    avg_areas(parc_id,1) = mean(parc_areas{parc_id,1});
end

%% Compare hyperalignment heritability difference to area

connectivity_perc = zeros(1,num_days);

hold on
x_min      = 0;
x_max      = max(avg_areas) * 1.5;
x_extended = linspace(x_min, x_max, 10000);

piecewise_fits    = cell(num_days, 1);
connectivity_fits = cell(num_days, 1);

for scan_id = 1:num_days
    subplot(2,1,scan_id)
    hold on

    % Calculate means and SEMs
    mean_movie_herit = nanmean(anatomical_isc_herit(:,scan_id));
    mean_piecewise   = squeeze(nanmean(piecewise_isc_herit(:,:,scan_id)));
    sem_piecewise    = piecewise_isc_se(:,scan_id);

    mean_connectivity = squeeze(nanmean(connectivity_isc_herit(:,:,scan_id)));
    sem_connectivity  = connectivity_isc_se(:,scan_id);

    sem_movie = nanmean(anatomical_se(:,scan_id));

    mean_piecewise    = [mean_piecewise';mean_movie_herit];
    mean_connectivity = [mean_connectivity';mean_movie_herit];

    % Define custom power law model including y-intercept
    customModel = fittype('a*x^b + c', 'independent', 'x', 'dependent', 'y');
    options     = fitoptions('Method', 'NonlinearLeastSquares');

    options.Lower      = [-Inf, -Inf, -Inf]; % Adjust lower bounds as necessary
    options.Upper      = [Inf, Inf, Inf]; % Adjust upper bounds as necessary
    options.StartPoint = [1, 1, 0]; % Adjust start points as necessary

    % Fit custom model for piecewise data
    [fitresult_piecewise, gof_piecewise] = fit([avg_areas;0], mean_piecewise, customModel, options);
    piecewise_fits{scan_id}              = struct('FitResult', fitresult_piecewise, 'GOF', gof_piecewise);

    % Fit custom model for connectivity data
    [fitresult_connectivity, gof_connectivity] = fit([avg_areas;0], mean_connectivity, customModel, options);
    connectivity_fits{scan_id}                 = struct('FitResult', fitresult_connectivity, 'GOF', gof_connectivity);

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

    % Doubling the line thickness in boundedline plots
    set(hl1, 'LineWidth', .001);
    set(hl2, 'LineWidth', .0001);

    % Scatter plots for mean values
    scatter(avg_areas, mean_piecewise(1:10), 150, 'filled', 'MarkerFaceColor', [254/255 97/255 0/255], 'MarkerEdgeColor', 'k','LineWidth',3);
    scatter(avg_areas, mean_connectivity(1:10), 150, 'filled', 'MarkerFaceColor', [120/255 94/255 240/255], 'MarkerEdgeColor', 'k','LineWidth',3);
    scatter(0, mean_movie_herit, 150, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    % SEM plot for movieNetHerit_mag (as a vertical line)
    errorbar(0, mean_movie_herit, sem_movie, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 5);

    % Axes settings
    ylim([.01 .025]); xlim([-100 1200]);

    xlabel('Average Hyperalignment Parcel Area (mm^2)'); ylabel('% of MSM-aligned Heritability');
    yticks([.01 .015 .02 .025]); ax = gca; set(ax, 'TickDir', 'out'); set(ax, 'FontSize', 32); set(ax, 'LineWidth', 5); set(ax, 'Layer', 'bottom');

    xticks([0 400 800 1200]);set(ax, 'TickDir', 'out')
    grid off; legend({'Response Hyperalignment', 'Connectivity Hyperalignment'}, 'Location', 'best'); set(gcf, 'position', [100,100,600,1200]);
    pbaspect([1.5 1 1])
    piecewise_perc(1,scan_id)    = 1-(mean_piecewise(1)/mean_movie_herit);
    connectivity_perc(1,scan_id) = 1-(mean_connectivity(1)/mean_movie_herit);

    piecewise_conf(1,scan_id)    = 1-((mean_piecewise(1) - (sem_piecewise(1)*1.96))/mean_movie_herit);
    piecewise_conf(2,scan_id)    = 1-((mean_piecewise(1) + (sem_piecewise(1)*1.96))/mean_movie_herit);
    connectivity_conf(1,scan_id) = 1-((mean_connectivity(1) - (sem_connectivity(1)*1.96))/mean_movie_herit);
    connectivity_conf(2,scan_id) = 1-((mean_connectivity(1) + (sem_connectivity(1)*1.96))/mean_movie_herit);
end

saveas(gcf,[data_dir '/isc_heritability/figures/fig3_areascatter_power.svg'])
close all

%% Compare hyperalignment heritability difference to area

connectivity_perc = zeros(1,num_days);

hold on
x_min      = 0;
x_max      = max(avg_areas) * 1.5;
x_extended = linspace(x_min, x_max, 10000);

linear_fits_piecewise    = cell(num_days, 1);
linear_fits_connectivity = cell(num_days, 1);

for scan_id = 1:num_days
    subplot(2,1,scan_id)
    hold on

    % Calculate means and SEMs
    mean_movie_herit = nanmean(anatomical_isc_herit(:,scan_id));
    mean_piecewise   = squeeze(nanmean(piecewise_isc_herit(:,:,scan_id)));
    sem_piecewise    = piecewise_isc_se(:,scan_id);

    mean_connectivity = squeeze(nanmean(connectivity_isc_herit(:,:,scan_id)));
    sem_connectivity  = connectivity_isc_se(:,scan_id);

    sem_movie = nanmean(anatomical_se(:,scan_id));

    mean_piecewise    = [mean_piecewise';mean_movie_herit];
    mean_connectivity = [mean_connectivity';mean_movie_herit];

    % Fit a linear model for piecewise data
    [fitresult_piecewise, gof_piecewise] = fit([avg_areas;0], mean_piecewise, 'poly1');
    linear_fits_piecewise{scan_id}       = struct('FitResult', fitresult_piecewise, 'GOF', gof_piecewise);

    % Fit a linear model for connectivity data
    [fitresult_connectivity, gof_connectivity] = fit([avg_areas;0], mean_connectivity, 'poly1');
    linear_fits_connectivity{scan_id}          = struct('FitResult', fitresult_connectivity, 'GOF', gof_connectivity);

    % Plot the extended linear fit lines
    y_extended_piecewise    = feval(fitresult_piecewise, x_extended);
    y_extended_connectivity = feval(fitresult_connectivity, x_extended);

    plot(x_extended, y_extended_piecewise, 'LineWidth', 4, 'Color', [254/255 97/255 0/255]);
    plot(x_extended, y_extended_connectivity, 'LineWidth', 4, 'Color', [120/255 94/255 240/255]);

    % boundedline plots
    hl1 = boundedline(avg_areas, mean_piecewise(1:10), sem_piecewise, 'alpha', 'cmap', [254/255 97/255 0/255],'LineWidth',.00001);

    hl2 = boundedline(avg_areas, mean_connectivity(1:10), sem_connectivity, 'alpha', 'cmap', [120/255 94/255 240/255],'LineWidth',.00001);

    % Doubling the line thickness in boundedline plots
    set(hl1, 'LineWidth', .001);
    set(hl2, 'LineWidth', .0001);

    % Scatter plots for mean values
    scatter(avg_areas, mean_piecewise(1:10), 150, 'filled', 'MarkerFaceColor', [254/255 97/255 0/255], 'MarkerEdgeColor', 'k','LineWidth',3);
    scatter(avg_areas, mean_connectivity(1:10), 150, 'filled', 'MarkerFaceColor', [120/255 94/255 240/255], 'MarkerEdgeColor', 'k','LineWidth',3);
    scatter(0, mean_movie_herit, 150, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    % SEM plot for movieNetHerit_mag (as a vertical line)
    errorbar(0, mean_movie_herit, sem_movie, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 5);

    % Axes settings
    ylim([.01 .025]); xlim([-100 1200]);

    xlabel('Average Hyperalignment Parcel Area (mm^2)'); ylabel('% of MSM-aligned Heritability');
    yticks([.01 .015 .02 .025]); ax = gca; set(ax, 'TickDir', 'out'); set(ax, 'FontSize', 32); set(ax, 'LineWidth', 5); set(ax, 'Layer', 'bottom');

    xticks([0 400 800 1200]);set(ax, 'TickDir', 'out')
    grid off; legend({'Response Hyperalignment', 'Connectivity Hyperalignment'}, 'Location', 'best'); set(gcf, 'position', [100,100,600,1200]);
    pbaspect([1.5 1 1])
    piecewise_perc(1,scan_id)    = 1-(mean_piecewise(1)/mean_movie_herit);
    connectivity_perc(1,scan_id) = 1-(mean_connectivity(1)/mean_movie_herit);

    piecewise_conf(1,scan_id)    = 1-((mean_piecewise(1) - (sem_piecewise(1)*1.96))/mean_movie_herit);
    piecewise_conf(2,scan_id)    = 1-((mean_piecewise(1) + (sem_piecewise(1)*1.96))/mean_movie_herit);
    connectivity_conf(1,scan_id) = 1-((mean_connectivity(1) - (sem_connectivity(1)*1.96))/mean_movie_herit);
    connectivity_conf(2,scan_id) = 1-((mean_connectivity(1) + (sem_connectivity(1)*1.96))/mean_movie_herit);
end

saveas(gcf,[data_dir '/isc_heritability/figures/fig3_areascatter_linear.fig'])
close all

log_fits_piecewise    = cell(num_days, 1);
log_fits_connectivity = cell(num_days, 1);

for scan_id = 1:num_days
    subplot(2,1,scan_id)
    hold on

    % Calculate means and SEMs
    mean_movie_herit = nanmean(anatomical_isc_herit(:,scan_id));
    mean_piecewise   = squeeze(nanmean(piecewise_isc_herit(:,:,scan_id)));
    sem_piecewise    = piecewise_isc_se(:,scan_id);

    mean_connectivity = squeeze(nanmean(connectivity_isc_herit(:,:,scan_id)));
    sem_connectivity  = connectivity_isc_se(:,scan_id);

    sem_movie = nanmean(anatomical_se(:,scan_id));

    mean_piecewise    = [mean_piecewise';mean_movie_herit];
    mean_connectivity = [mean_connectivity';mean_movie_herit];

    % Define a logarithmic model: y = a * log(x) + b
    logModel = fittype('a*log(x) + b', 'independent', 'x', 'dependent', 'y');
    options  = fitoptions('Method', 'NonlinearLeastSquares');
    options.StartPoint = [1, 0];   % Adjust start points as necessary
    options.Lower = [-Inf, -Inf];  % Adjust lower bounds as necessary
    options.Upper = [Inf, Inf];    % Adjust upper bounds as necessary

    % Fit logarithmic model for piecewise data
    [fitresult_piecewise, gof_piecewise] = fit([avg_areas;1], mean_piecewise, logModel, options);
    log_fits_piecewise{scan_id}          = struct('FitResult', fitresult_piecewise, 'GOF', gof_piecewise);

    % Fit logarithmic model for connectivity data
    [fitresult_connectivity, gof_connectivity] = fit([avg_areas;1], mean_connectivity, logModel, options);
    log_fits_connectivity{scan_id}              = struct('FitResult', fitresult_connectivity, 'GOF', gof_connectivity);

    % Plot the extended logarithmic fit lines
    y_extended_piecewise    = feval(fitresult_piecewise, x_extended);
    y_extended_connectivity = feval(fitresult_connectivity, x_extended);

    plot(x_extended, y_extended_piecewise, 'LineWidth', 4, 'Color', [254/255 97/255 0/255]);
    plot(x_extended, y_extended_connectivity, 'LineWidth', 4, 'Color', [120/255 94/255 240/255]);

    % boundedline plots
    hl1 = boundedline(avg_areas, mean_piecewise(1:10), sem_piecewise, 'alpha', 'cmap', [254/255 97/255 0/255],'LineWidth',.00001);
    hl2 = boundedline(avg_areas, mean_connectivity(1:10), sem_connectivity, 'alpha', 'cmap', [120/255 94/255 240/255],'LineWidth',.00001);

    % Doubling the line thickness in boundedline plots
    set(hl1, 'LineWidth', .001);
    set(hl2, 'LineWidth', .0001);

    % Scatter plots for mean values
    scatter(avg_areas, mean_piecewise(1:10), 150, 'filled', 'MarkerFaceColor', [254/255 97/255 0/255], 'MarkerEdgeColor', 'k','LineWidth',3);
    scatter(avg_areas, mean_connectivity(1:10), 150, 'filled', 'MarkerFaceColor', [120/255 94/255 240/255], 'MarkerEdgeColor', 'k','LineWidth',3);
    scatter(1, mean_movie_herit, 150, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    % SEM plot for movieNetHerit_mag (as a vertical line)
    errorbar(1, mean_movie_herit, sem_movie, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 5);

    % Axes settings
    ylim([.01 .025]); xlim([-100 1200]);

    xlabel('Average Hyperalignment Parcel Area (mm^2)'); ylabel('% of MSM-aligned Heritability');
    yticks([.01 .015 .02 .025]); ax = gca; set(ax, 'TickDir', 'out'); set(ax, 'FontSize', 32); set(ax, 'LineWidth', 5); set(ax, 'Layer', 'bottom');

    xticks([0 400 800 1200]); set(ax, 'TickDir', 'out')
    grid off; legend({'Response Hyperalignment', 'Connectivity Hyperalignment'}, 'Location', 'best'); set(gcf, 'position', [100,100,600,1200]);
    pbaspect([1.5 1 1])
    piecewise_perc(1,scan_id)    = 1-(mean_piecewise(1)/mean_movie_herit);
    connectivity_perc(1,scan_id) = 1-(mean_connectivity(1)/mean_movie_herit);

    piecewise_conf(1,scan_id)    = 1-((mean_piecewise(1) - (sem_piecewise(1)*1.96))/mean_movie_herit);
    piecewise_conf(2,scan_id)    = 1-((mean_piecewise(1) + (sem_piecewise(1)*1.96))/mean_movie_herit);
    connectivity_conf(1,scan_id) = 1-((mean_connectivity(1) - (sem_connectivity(1)*1.96))/mean_movie_herit);
    connectivity_conf(2,scan_id) = 1-((mean_connectivity(1) + (sem_connectivity(1)*1.96))/mean_movie_herit);
end

saveas(gcf,[data_dir '/isc_heritability/figures/fig3_areascatter_log.svg'])
saveas(gcf,[data_dir '/isc_heritability/figures/fig3_areascatter_log.fig'])

close all
poly_fits_piecewise    = cell(num_days, 1);
poly_fits_connectivity = cell(num_days, 1);

% Define the polynomial degree
poly_degree = 2;  % For example, degree 2 polynomial (quadratic)

for scan_id = 1:num_days
    subplot(2,1,scan_id)
    hold on

    % Calculate means and SEMs
    mean_movie_herit = nanmean(anatomical_isc_herit(:,scan_id));
    mean_piecewise   = squeeze(nanmean(piecewise_isc_herit(:,:,scan_id)));
    sem_piecewise    = piecewise_isc_se(:,scan_id);

    mean_connectivity = squeeze(nanmean(connectivity_isc_herit(:,:,scan_id)));
    sem_connectivity  = connectivity_isc_se(:,scan_id);

    sem_movie = nanmean(anatomical_se(:,scan_id));

    mean_piecewise    = [mean_piecewise';mean_movie_herit];
    mean_connectivity = [mean_connectivity';mean_movie_herit];

    % Polynomial model (degree 2)
    polynomialModel = fittype('a*x^2 + b*x + c', 'independent', 'x', 'dependent', 'y');
    options = fitoptions(polynomialModel);
    options.StartPoint = [1, 1, 0];  % Set reasonable starting points for a, b, and c

    % Fit polynomial model for piecewise data, including x = 0
    [fitresult_piecewise, gof_piecewise] = fit([avg_areas; 0], mean_piecewise, polynomialModel, options);
    poly_fits_piecewise{scan_id}         = struct('FitResult', fitresult_piecewise, 'GOF', gof_piecewise);

    % Fit polynomial model for connectivity data, including x = 0
    [fitresult_connectivity, gof_connectivity] = fit([avg_areas; 0], mean_connectivity, polynomialModel, options);
    poly_fits_connectivity{scan_id}             = struct('FitResult', fitresult_connectivity, 'GOF', gof_connectivity);

    % Plot the extended polynomial fit lines for x >= 0
    valid_x_extended = x_extended(x_extended >= 0);  % Restrict x to non-negative values
    y_extended_piecewise    = feval(fitresult_piecewise, valid_x_extended);
    y_extended_connectivity = feval(fitresult_connectivity, valid_x_extended);

    plot(valid_x_extended, y_extended_piecewise, 'LineWidth', 4, 'Color', [254/255 97/255 0/255]);
    plot(valid_x_extended, y_extended_connectivity, 'LineWidth', 4, 'Color', [120/255 94/255 240/255]);

    % boundedline plots
    hl1 = boundedline(avg_areas, mean_piecewise(1:10), sem_piecewise, 'alpha', 'cmap', [254/255 97/255 0/255], 'LineWidth', .00001);
    hl2 = boundedline(avg_areas, mean_connectivity(1:10), sem_connectivity, 'alpha', 'cmap', [120/255 94/255 240/255], 'LineWidth', .00001);

    % Doubling the line thickness in boundedline plots
    set(hl1, 'LineWidth', .001);
    set(hl2, 'LineWidth', .0001);

    % Scatter plots for mean values
    scatter(avg_areas, mean_piecewise(1:10), 150, 'filled', 'MarkerFaceColor', [254/255 97/255 0/255], 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    scatter(avg_areas, mean_connectivity(1:10), 150, 'filled', 'MarkerFaceColor', [120/255 94/255 240/255], 'MarkerEdgeColor', 'k', 'LineWidth', 3);
    scatter(0, mean_movie_herit, 150, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    % SEM plot for movieNetHerit_mag (as a vertical line)
    errorbar(0, mean_movie_herit, sem_movie, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 5);

    % Axes settings
    ylim([.01 .025]); xlim([-100 1200]);  % Limit x to non-negative range

    xlabel('Average Hyperalignment Parcel Area (mm^2)'); ylabel('% of MSM-aligned Heritability');
    yticks([.01 .015 .02 .025]); ax = gca; set(ax, 'TickDir', 'out'); set(ax, 'FontSize', 32); set(ax, 'LineWidth', 5); set(ax, 'Layer', 'bottom');

    xticks([0 400 800 1200]); set(ax, 'TickDir', 'out');
    grid off; legend({'Response Hyperalignment', 'Connectivity Hyperalignment'}, 'Location', 'best'); set(gcf, 'position', [100,100,600,1200]);
    pbaspect([1.5 1 1]);

    % Confidence calculations
    piecewise_perc(1,scan_id)    = 1 - (mean_piecewise(1)/mean_movie_herit);
    connectivity_perc(1,scan_id) = 1 - (mean_connectivity(1)/mean_movie_herit);

    piecewise_conf(1,scan_id)    = 1 - ((mean_piecewise(1) - (sem_piecewise(1)*1.96))/mean_movie_herit);
    piecewise_conf(2,scan_id)    = 1 - ((mean_piecewise(1) + (sem_piecewise(1)*1.96))/mean_movie_herit);
    connectivity_conf(1,scan_id) = 1 - ((mean_connectivity(1) - (sem_connectivity(1)*1.96))/mean_movie_herit);
    connectivity_conf(2,scan_id) = 1 - ((mean_connectivity(1) + (sem_connectivity(1)*1.96))/mean_movie_herit);
end

saveas(gcf,[data_dir '/isc_heritability/figures/fig3_areascatter_poly.svg']);
saveas(gcf,[data_dir '/isc_heritability/figures/fig3_areascatter_poly.fig']);

close all;

%% Make ISC MZ/DZ/NZ line plots

% FIGURE 1 LINES
anatomical_parc_dir = [data_dir '/isc_heritability/data/anatomical/outputs/parc/'];

% Initialize anatomical similarity matrices
anatomical_isc_similarity_mz_parc_all = zeros(length(pairwise_mz_id), 400, num_days);
anatomical_isc_similarity_dz_parc_all = zeros(length(pairwise_dz_id), 400, num_days);
anatomical_isc_similarity_nz_parc_all = zeros(length(pairwise_nz_id), 400, num_days);

anatomical_isc_parc_all = zeros(400,num_days);
for scan_id = 1:num_days
    iscDataPath = (strcat(anatomical_parc_dir,'anatomical_isc_scan_',num2str(scan_id),'.mat'));

    isc  = load(iscDataPath);
    name = fieldnames(isc);
    isc  = isc.(name{1});
    % Now separate out for MZ, DZ, and NZ

    for parc = 1:size(isc,3)
        isc_temp                   = squeeze(isc(:,:,parc));
        [m, n]                     = size(isc_temp);
        isc_temp(triu(true(m, n))) = NaN;

        anatomical_isc_similarity_mz_parc_all(:,parc,scan_id) = isc_temp(pairwise_mz_id);
        anatomical_isc_similarity_dz_parc_all(:,parc,scan_id) = isc_temp(pairwise_dz_id);
        anatomical_isc_similarity_nz_parc_all(:,parc,scan_id) = isc_temp(pairwise_nz_id);
    end
    isc(isc==1)= NaN;
    anatomical_isc_parc_all(:,scan_id) = conv_z2r(squeeze(nanmean(nanmean(conv_r2z(isc)))));

end

% Concatenate data
all_xz_isc = cat(2, ...
    conv_z2r(squeeze(nanmean(conv_r2z(anatomical_isc_similarity_mz_parc_all), 1))), ...
    conv_z2r(squeeze(nanmean(conv_r2z(anatomical_isc_similarity_dz_parc_all), 1))), ...
    conv_z2r(squeeze(nanmean(conv_r2z(anatomical_isc_similarity_nz_parc_all), 1))));

% Initialize matrices
mz_std_mat = zeros(size(anatomical_isc_similarity_mz_parc_all, 2), num_days);
dz_std_mat = mz_std_mat;
nz_std_mat = mz_std_mat;
sortIndex1 = mz_std_mat;

% Loop through each scan
for scan_id = 1:num_days
    mz_std_mat(:, scan_id) = conv_z2r(nanstd(conv_r2z(anatomical_isc_similarity_mz_parc_all(:,:,scan_id)))' / sqrt(sum(~isnan(anatomical_isc_similarity_mz_parc_all(:,1,scan_id)))));
    dz_std_mat(:, scan_id) = conv_z2r(nanstd(conv_r2z(anatomical_isc_similarity_dz_parc_all(:,:,scan_id)))' / sqrt(sum(~isnan(anatomical_isc_similarity_dz_parc_all(:,1,scan_id)))));
    nz_std_mat(:, scan_id) = conv_z2r(nanstd(conv_r2z(anatomical_isc_similarity_nz_parc_all(:,:,scan_id)))' / sqrt(sum(~isnan(anatomical_isc_similarity_nz_parc_all(:,1,scan_id)))));

    % Sorting and plotting
    [~, sortIndex1(:, scan_id)] = sort(anatomical_isc_parc_all(:, scan_id), 1);

    sorted_means       = all_xz_isc(sortIndex1(:, scan_id), [scan_id, scan_id + 2, scan_id + 4]);
    sortIndex_curr     = sortIndex1(:, scan_id);
    mz_std_mat_sorted  = mz_std_mat(sortIndex_curr, scan_id);
    dz_std_mat_sorted  = dz_std_mat(sortIndex_curr, scan_id);
    nz_std_mat_sorted  = nz_std_mat(sortIndex_curr, scan_id);
    sorted_std_errors  = cat(2, mz_std_mat_sorted, dz_std_mat_sorted, nz_std_mat_sorted);

    subplot(1, 2, scan_id)

    % Create a new x-axis vector
    x = 1:size(anatomical_isc_parc_all,1);

    % Plot each line with its shaded error region
    color_vec = [[100/255, 143/255, 255/255]; [220/255, 38/255, 127/255]; [255/255, 176/255, 0/255]]; % colors for each line
    for i = 1:size(sorted_means, 2)
        boundedline(x, (sorted_means(:, i)), (sorted_std_errors(:, i)), 'cmap', color_vec(i, :), 'alpha','linewidth',2);
        hold on;
    end
    ylim([0, 1])
    xlim([0, size(anatomical_isc_parc_all,1)])
    xlabel('Parcel Rank')
    ylabel('ISC (r)')
    yticks([0, 0.25, .5, .75, 1])

    % Customize axes
    ax = gca;
    set(ax, 'TickDir', 'out', 'FontSize', 20, 'LineWidth', 3, 'Layer', 'bottom')

    grid off
    legend({'', 'MZ', '', 'DZ', '', 'UR', ''})
    pbaspect([1, 1, 1])
end

% Save the figure
set(gcf, 'position', [100, 100, 1200, 1200])
saveas(gcf, [data_dir '/isc_heritability/figures/fig1_ISClines_avgsort_parc.fig'])


%% ISC group scatter plots (MZ vs DZ, MZ vs UR, DZ vs UR)
%% Group ISC arrays: pull dimensions
P = size(anatomical_isc_similarity_mz_parc_all, 1);   % # parcels
D = size(anatomical_isc_similarity_mz_parc_all, 3);   % # days

%% Ensure DZ/UR arrays are [parcels x subjects x days]
% DZ:
if size(anatomical_isc_similarity_dz_parc_all,1) ~= P && size(anatomical_isc_similarity_dz_parc_all,2) == P
    anatomical_isc_similarity_dz_parc_all = permute(anatomical_isc_similarity_dz_parc_all, [2 1 3]);
end
% UR (nZ):
if size(anatomical_isc_similarity_nz_parc_all,1) ~= P && size(anatomical_isc_similarity_nz_parc_all,2) == P
    anatomical_isc_similarity_nz_parc_all = permute(anatomical_isc_similarity_nz_parc_all, [2 1 3]);
end

%% Parcel-wise mean ISC per group and day
P = size(anatomical_isc_similarity_mz_parc_all,2);
D = size(anatomical_isc_similarity_mz_parc_all,3);  % should be 2 for two days

mz_day = nan(P, D);
dz_day = nan(P, D);
ur_day = nan(P, D);
for d = 1:D
    mz_day(:,d) = conv_z2r( nanmean( conv_r2z(anatomical_isc_similarity_mz_parc_all(:,:,d)), 1 ) );
    dz_day(:,d) = conv_z2r( nanmean( conv_r2z(anatomical_isc_similarity_dz_parc_all(:,:,d)), 1 ) );
    ur_day(:,d) = conv_z2r( nanmean( conv_r2z(anatomical_isc_similarity_nz_parc_all(:,:,d)), 1 ) );
end

%% Scatter grid (rows = days; cols = MZxDZ, MZxUR, DZxUR)
figure;
pairNames = {'MZ','DZ','UR'};
combos     = [1 2; 1 3; 2 3];  % which pairs to plot

for d = 1:D
    for k = 1:3
        idx = (d-1)*3 + k;
        subplot(D, 3, idx);

        % pick the two ISC vectors for this day
        X = eval([lower(pairNames{combos(k,1)}), '_day(:,d)']);
        Y = eval([lower(pairNames{combos(k,2)}), '_day(:,d)']);

        % scatter with black edge
        h = scatter(X, Y, 50, 'filled', 'MarkerEdgeColor', 'k');

        h.MarkerFaceAlpha = 0.8;
        h.LineWidth        = 0.5;
        hold on;

        % unity line
        mn = min([X;Y]); mx = max([X;Y]);
        plot([mn mx], [mn mx], 'k--', 'LineWidth', 1);

        % styling
        axis square;
        xlim([0 .6]); ylim([0 .6]);
        xlabel([pairNames{combos(k,1)}, ' ISC (r)']);
        ylabel([pairNames{combos(k,2)}, ' ISC (r)']);
        title(sprintf('Day %d: %s vs. %s', d, pairNames{combos(k,1)}, pairNames{combos(k,2)}));

        % Pearson r
        R = corr(X, Y, 'Rows', 'complete');
        text(0.05, 0.9, sprintf('r = %.2f', R), 'Units', 'normalized', 'FontSize', 11);

        set(gca, 'FontSize', 12, 'LineWidth', 1, 'TickDir', 'out', 'Layer', 'bottom');
    end
end

set(gcf, 'Position', [200 200 1200 700]);

%% Parcellate Schaefer files and raw data
parc_reses      = {'100','200','300','400','500','600','700','800','900','1000'};
schaefer_dlabel = cell(num_parcs,1);
for parc_res_id = 1:length(parc_reses)
    curr_parc = parc_reses{parc_res_id};

    schaefer_dlabel{parc_res_id} = strcat(desktop_dir, '/Schaefer2018_',curr_parc,'Parcels_17Networks_order.dlabel.nii');
    parc_path                    = strcat(desktop_dir, '/Schaefer2018_',curr_parc,'Parcels_17Networks_order.dscalar.nii');
    parc_command                 = strcat(wb_path,{' '},'-cifti-parcellate',{' '},parc_path,{' '},schaefer_dlabel{parc_res_id},{' '},'COLUMN',{' '} ,strrep(parc_path,'dscalar','pscalar'));
    system(parc_command{1})

    schaefer_dlabel_kong = strcat(desktop_dir, '/Schaefer2018_',curr_parc,'Parcels_Kong2022_17Networks_order.dlabel.nii');
    parc_path            = strcat(desktop_dir, '/Schaefer2018_',curr_parc,'Parcels_Kong2022_17Networks_order.dscalar.nii');
    parc_command         =  strcat(wb_path,{' '},'-cifti-parcellate',{' '},parc_path,{' '},schaefer_dlabel_kong,{' '},'COLUMN',{' '} ,strrep(parc_path,'dscalar','pscalar'));
    system(parc_command{1})
end

iscHeritabilityParcellate([data_dir '/isc_heritability/data/anatomical/inputs/parc/'], 'raw_data_lh_movie_', num_subj_movie,parc_reses,schaefer_400_label,'400');
iscHeritabilityParcellate([data_dir '/isc_heritability/data/anatomical/inputs/parc/'], 'raw_data_lh_rest_', num_subj_movie,parc_reses,schaefer_400_label,'400');

%% Show difference in ISC heritability after hyperalignment (Fig. 4D)

piecewise_diff_out    = zeros(num_vox_cortex,num_days);
connectivity_diff_out = zeros(num_vox_cortex,num_days);

for scan_id = 1:num_days
    for vox = 1:num_vox_cortex
        piecewise_vals    = piecewise_isc_herit(vox, 1, scan_id);
        connectivity_vals = connectivity_isc_herit(vox, 1, scan_id);
        idxs              = parc_ids(vox, :);
        if min(idxs) > 0

            c                 = anatomical_isc_herit(vox, scan_id);
            piecewise_diff    = (piecewise_vals-c);
            connectivity_diff = (connectivity_vals-c);

            % Storing the results
            piecewise_diff_out(vox, scan_id)    = (piecewise_diff);
            connectivity_diff_out(vox, scan_id) = (connectivity_diff);
        end
    end

    % Saving the results
    piecewise_file = fullfile(base_dir, 'piecewise/outputs/gray', ...
        sprintf('piecewise_herit_diff_scan_%d.dscalar.nii', scan_id));
    connectivity_file = fullfile(base_dir, 'connectivity/outputs/gray', ...
        sprintf('connectivity_herit_diff_scan_%d.dscalar.nii', scan_id));

    saveCifti(piecewise_diff_out(:, scan_id)', piecewise_file, wb_path, schaefer_400_dscalar, medial_mask);
    saveCifti(connectivity_diff_out(:, scan_id)', connectivity_file, wb_path, schaefer_400_dscalar, medial_mask);
end

piecewise_diff_trt    = corr(piecewise_diff_out,'type','spearman','rows','complete');
connectivity_diff_trt = corr(connectivity_diff_out,'type','spearman','rows','complete');

piecewise_connectivity_trt(1) = corr(piecewise_diff_out(:,1),connectivity_diff_out(:,1),'type','spearman','rows','complete');
piecewise_connectivity_trt(2) = corr(piecewise_diff_out(:,2),connectivity_diff_out(:,2),'type','spearman','rows','complete');


%% Generate cifti distance matrices

cifti_path        = [desktop_dir '/inflated_brains/left_distances.dconn.nii'];
left_surface      = cifti_read(cifti_path,wb_path);
left_surface_dist = double(left_surface.cdata);
h5create([desktop_dir '/inflated_brains/distances_l.h5'],"/DS1",[32492 32492])
h5write([desktop_dir '/inflated_brains/distances_l.h5'],"/DS1",left_surface_dist)

cifti_path         = [desktop_dir '/inflated_brains/right_distances.dconn.nii'];
right_surface      = cifti_read(cifti_path,wb_path);
right_surface_dist = double(right_surface.cdata);
h5create([desktop_dir '/inflated_brains/distances_r.h5'],"/DS1",[32492 32492])
h5write([desktop_dir '/inflated_brains/distances_r.h5'],"/DS1",right_surface_dist)

%% Generate RSFC matrices for CHA
num_subj_rest = 176;
base_dir      = [data_dir '/isc_heritability/data/anatomical/inputs/gray/'];

for parc_val_id = 1:length(parc_reses)
    parc_val     = parc_reses{parc_val_id};
    current_parc = strcat('Schaefer2018_',parc_val,'Parcels_17Networks_order.dscalar.nii');

    vox_ids = ciftiopen(strcat(downloads_dir, '/ThomasYeoLab CBIG master stable_projects-brain_parcellation_Schaefer2018_LocalGlobal_Parcellations_HCP_fslr32k_cifti/',current_parc));
    vox_ids = vox_ids.cdata(medial_mask==1);

    parc_num = str2double(parc_val);

    for scan_id = 1:num_days
        cha_mat = zeros(num_vox_cortex,parc_num,num_subj_rest);
        for subj = 1:num_subj_rest
            subj_data_lh = readNPY(strcat([data_dir '/isc_heritability/data/anatomical/inputs/gray/raw_data_lh_rest_'],num2str(scan_id),'_',num2str(subj),'.npy'));
            subj_data_rh = readNPY(strcat([data_dir '/isc_heritability/data/anatomical/inputs/gray/raw_data_rh_rest_'],num2str(scan_id),'_',num2str(subj),'.npy'));
            subj_data    = cat(2,subj_data_lh,subj_data_rh);
            parc_data    = zeros(size(subj_data,2),parc_num);
            mean_vec     = zeros(parc_num,size(subj_data,1));
            parfor parc_id = 1:parc_num
                parc_vox_ids        = find(vox_ids == parc_id);
                mean_vec(parc_id,:) = nanmean(subj_data(:,parc_vox_ids)')';
            end
            if parc_num ~=1000
                cha_mat(:,:,subj) = corr(subj_data,mean_vec','rows','complete');
            else
                parfor i = 1:parc_num
                    cha_mat(:,i,subj) = corr(subj_data,mean_vec(i,:)','rows','complete');
                end
            end
        end

        save(strcat(base_dir,'cha_mat_parc_',num2str(parc_num),'_day_',num2str(scan_id),'.mat'),'cha_mat','-v7.3')
    end
end

%% Parcellate piecewise
iscHeritabilityParcellate([data_dir '/isc_heritability/data/piecewise/inputs/parc/'], 'piecewise_parc_', num_subj_movie,parc_reses(10),schaefer_dlabel{4},num_trs);

%% Parcellate connectivity
iscHeritabilityParcellate([data_dir '/isc_heritability/data/connectivity/inputs/parc/'], 'connectivity_parc_', num_subj_rest,parc_reses(10),schaefer_dlabel{4},num_trs);

base_dir = [data_dir '/isc_heritability/data/'];
%% Calculate MNT (movie neural timescale)
if exist(strcat(base_dir,'anatomical/outputs/gray/mnt_gray.mat'),'file') == 2
    load(strcat(base_dir,'anatomical/outputs/gray/mnt_gray.mat'))
else
    mnt_acf = zeros(num_vox_cortex,num_subj_movie,num_days);
    for scan_id = 1:num_days
        for subj = 1:num_subj_movie
            subj_data = readNPY(strcat([data_dir '/isc_heritability/data/anatomical/inputs/gray/anatomical_movie_scan_'], num2str(scan_id), '_subj_', num2str(subj), '.npy'));
            mnt_acf(:,subj,scan_id) = calculate_int(subj_data);
        end
    end
    save(strcat(base_dir,'anatomical/outputs/gray/mnt_gray.mat'),'mnt_acf')
end

if exist(strcat(base_dir,'piecewise/outputs/gray/mnt_gray.mat'),'file') == 2
    load(strcat(base_dir,'piecewise/outputs/gray/mnt_gray.mat'))
else
    mnt_acf = zeros(num_vox_cortex,num_subj_movie,num_days);
    for scan_id = 1:num_days
        for subj = 1:num_subj_movie
            subj_data               = readNPY(strcat([data_dir '/isc_heritability/data/piecewise/inputs/gray/piecewise_parc_100_scan_'], num2str(scan_id), '_subj_', num2str(subj), '.npy'));
            mnt_acf(:,subj,scan_id) = calculate_int(subj_data');
        end
    end
    save(strcat(base_dir,'piecewise/outputs/gray/mnt_gray.mat'),'mnt_acf')
end

%% Are MNT topographies heritable?

for day_id = 1:num_days
    mnt_corr = corr(mnt_acf(:,:,day_id),'rows','complete','type','spearman');
    mnt_topo_herit(day_id) = h2_mat(mnt_corr, kinship(subjs_movie, subjs_movie), [covari_motion(subjs_movie, [1, 2, day_id + 2])], 0,[]);
end

schaefer_parc = ciftiopen(schaefer_400_dscalar);
schaefer_parc = schaefer_parc.cdata(medial_mask==1);
mnt_parc = zeros(400,num_subj_movie,num_days);
for parc_id = 1:400
    mnt_parc(parc_id,:,:) = squeeze(nanmean(mnt_acf(schaefer_parc==parc_id,:,:),1));
end

for scan_id = 1:num_days
    csvwrite(strcat(desktop_dir, '/isc_heritability_r/mnt_movie_scan_',num2str(scan_id),'.csv'),mnt_parc(:,:,scan_id));
end

for day_id = 1:num_days
    mnt_corr = corr(mnt_parc(:,:,day_id),'type','spearman');
    [mnt_topo_herit_parc(day_id), p_perm(day_id), jack_se,h2_jack] = h2_mat(mnt_corr, kinship(subjs_movie, subjs_movie), [covari_motion(subjs_movie, [1, 2, day_id + 2])], 0,[]);
end


for scan_id = 1:num_days
    mnt_mag_herit{1,scan_id} = readtable(strcat([data_dir "/isc_heritability/data/solar/mnt_herit_parcel_scan_"], num2str(scan_id), ".csv"));
end

mnt_mag_herit_mat     = zeros(400,num_days);

for scan_id = 1:num_days
    mnt_mag_herit_mat(:,scan_id)     = mnt_mag_herit{1,scan_id}.Var;
end

%% Is MNT associated with average heritability?

% Get pairwise measures and run ANOVA
mnt_avg         = squeeze(nanmean(mnt_acf,1));
mnt_acf_out_big = zeros(num_subj_movie,num_subj_movie,num_vox_cortex,num_days);
for scan_id = 1:num_days
    mnt_acf_diff                   = pdist2(mnt_avg(:,scan_id),mnt_avg(:,scan_id),'euclidean');
    [m, n]                         = size(mnt_acf_diff);
    mnt_acf_diff(triu(true(m, n))) = NaN;

    mnt_mz(:,scan_id) = mnt_acf_diff(pairwise_mz_id);
    mnt_dz(:,scan_id) = mnt_acf_diff(pairwise_dz_id);
    mnt_nz(:,scan_id) = mnt_acf_diff(pairwise_nz_id);

    mnt_acf_out(:,:,scan_id) = mnt_acf_diff;
    parfor vox = 1:num_vox_cortex
        mnt_acf_out_big(:,:,vox,scan_id) = pdist2(squeeze(mnt_acf(vox,:,scan_id))',squeeze(mnt_acf(vox,:,scan_id))','euclidean');
    end
end

for scan_id = 1:num_days
    outliers_mz{scan_id} = find_outliers(squeeze(mnt_mz(:, scan_id)));
    outliers_dz{scan_id} = find_outliers(squeeze(mnt_dz(:, scan_id)));
    outliers_nz{scan_id} = find_outliers(squeeze(mnt_nz(:, scan_id)));

    mnt_mz(outliers_mz{scan_id}, scan_id) = NaN;
    mnt_dz(outliers_dz{scan_id}, scan_id) = NaN;
    mnt_nz(outliers_nz{scan_id}, scan_id) = NaN;
end

perc_outliers = sum(sum(outliers_nz{1}) + sum(outliers_nz{2}) +sum(outliers_dz{1}) + sum(outliers_dz{2}) + sum(outliers_mz{1}) + sum(outliers_mz{2}))/(2*sum(length(mnt_mz) + length(mnt_dz) + length(mnt_nz)));

% Start MATLAB's parallel pool (use threads)
if isempty(gcp('nocreate'))
    parpool;
end

% Number of regions, days, and permutations
num_perms = 10000;

% Data arrays
dataArrays = {mnt_mz, mnt_dz, mnt_nz};

% Preallocate arrays for t-scores and p-values
Differences = zeros(1, 6); % 6 comparisons (3 pairs x 2 days)
mnt_pvalues = zeros(1, 6);
mnt_pvalues_fdr = zeros(1, 6);
% Parallel computation over regions
for day = 1:num_days
    % Actual data for each group
    dataMZ = dataArrays{1}(:, day);
    dataDZ = dataArrays{2}(:, day);
    dataNZ = dataArrays{3}(:, day);

    % Calculate actual differences
    diffMZDZ = nanmean(dataMZ) - nanmean(dataDZ);
    diffMZNZ = nanmean(dataMZ) - nanmean(dataNZ);
    diffDZNZ = nanmean(dataDZ) - nanmean(dataNZ);

    actualDifferences = [diffMZDZ, diffMZNZ, diffDZNZ];

    % Initialize permutation t-scores
    permDifferences = zeros(num_perms, 3);

    for perm = 1:num_perms
        rng(perm)
        % Shuffle group identities while maintaining proportions
        allData      = [dataMZ; dataDZ; dataNZ];
        shuffledData = allData(randperm(length(allData)));

        numMZ = numel(dataMZ);
        numDZ = numel(dataDZ);

        shuffledMZ = shuffledData(1:numMZ);
        shuffledDZ = shuffledData(numMZ + 1:numMZ + numDZ);
        shuffledNZ = shuffledData(numMZ + numDZ + 1:end);

        % Calculate differences for shuffled data
        permDifferences(perm, 1, day) = nanmean(shuffledMZ) - nanmean(shuffledDZ);
        permDifferences(perm, 2, day) = nanmean(shuffledMZ) - nanmean(shuffledNZ);
        permDifferences(perm, 3, day) = nanmean(shuffledDZ) - nanmean(shuffledNZ);
    end

    % Calculate p-values
    for comp = 1:3
        mnt_pvalues((day - 1) * 3 + comp) = (sum(abs(permDifferences(:, comp, day)) >= abs(actualDifferences(comp)))+1)/(num_perms+1);
    end
    mnt_pvalues_fdr(:,(1:3) + (day-1)*3) = fdr_bh(mnt_pvalues(:,(1:3) + (day-1)*3));

    % Perform FDR correction for p-values corresponding to the current day
    start_idx = (day - 1) * 3 + 1;  % Index for the first component of the current day
    end_idx = day * 3;              % Index for the last component of the current day
    [h(start_idx:end_idx), crit_p(day), ~, adj_p(start_idx:end_idx)] = fdr_bh(mnt_pvalues(start_idx:end_idx));
end

% Calculate means for each group and day
meanMZ = [nanmean(squeeze(mnt_mz(:, 1))), nanmean(squeeze(mnt_mz(:, 2)))];
meanDZ = [nanmean(squeeze(mnt_dz(:, 1))), nanmean(squeeze(mnt_dz(:, 2)))];
meanNZ = [nanmean(squeeze(mnt_nz(:, 1))), nanmean(squeeze(mnt_nz(:, 2)))];

% Calculate standard errors for each group and day
semMZ = [nanstd(squeeze(mnt_mz(:, 1)))/sqrt(sum(~isnan(mnt_mz(:, 1)))), nanstd(squeeze(mnt_mz(:, 2)))/sqrt(sum(~isnan(mnt_mz(:, 2))))];
semDZ = [nanstd(squeeze(mnt_dz(:, 1)))/sqrt(sum(~isnan(mnt_dz(:, 1)))), nanstd(squeeze(mnt_dz(:, 2)))/sqrt(sum(~isnan(mnt_dz(:, 2))))];
semNZ = [nanstd(squeeze(mnt_nz(:, 1)))/sqrt(sum(~isnan(mnt_nz(:, 1)))), nanstd(squeeze(mnt_nz(:, 2)))/sqrt(sum(~isnan(mnt_nz(:, 2))))];

% Create bar plot
figure;
b = bar([meanMZ; meanDZ; meanNZ], 'FaceColor', 'flat','LineWidth',3);
hold on;

% Set colors for different days
b(1).CData = repmat([70/255, 130/255, 180/255], 3, 1); % Color for Day 1
b(2).CData = repmat([107/255, 142/255, 35/255], 3, 1); % Color for Day 2

% Add error bars
ngroups = size([meanMZ; meanDZ; meanNZ], 1);
nbars   = size([meanMZ; meanDZ; meanNZ], 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = b(i).XEndPoints;
    errorbar(x, [meanMZ(i), meanDZ(i), meanNZ(i)], [semMZ(i), semDZ(i), semNZ(i)], 'k', 'linestyle', 'none','LineWidth',3);
end

% Add significance indicators (example for MZ vs NZ)

% Coordinates for line and asterisk
x1 = sum(b(1).XEndPoints(1) + b(2).XEndPoints(1))/2;
x2 = sum(b(1).XEndPoints(3) + b(2).XEndPoints(3))/2;
y  = max([max(meanMZ), max(meanNZ)]) + max([max(semMZ), max(semNZ)]) * 1.5;
plot([x1, x2], [y, y], 'k', 'LineWidth', 3); % Line
plot(mean([x1, x2]), y * 1.05, 'k*', 'LineWidth', 3); % Asterisk

% Set labels and titles
set(gca, 'XTick', 1:3, 'XTickLabel', {'MZ', 'DZ', 'NZ'});
ylabel('BOLD Time Course Similarity');
title('BOLD Time Course Similarity by Group and Day');
legend('Day 1', 'Day 2');
hold off;
% Get the axis handle
ax = gca;
set(ax, 'TickDir', 'out')
set(ax, 'FontSize', 20)

% Set the line width of the x-axis and y-axis
set(ax, 'LineWidth', 3)
set(ax, 'Layer', 'bottom')
hPlots = findobj(ax, 'Type', 'line');
box off
ylim([0 .3])
yticks([0 .1 .2 .3])
pbaspect([2 1 1])
saveas(gcf,[data_dir '/isc_heritability/figures/fig6_barplot.fig'])

%% Control for tbold
% Predefined variables
indices = find(ismember(subjs_movie, subjs_rest));

% Initialization
isc_herit_mnt      = zeros(num_vox_cortex, num_days); % 3D array to store results
isc_herit_mnt_jack = zeros(num_vox_cortex,num_fams, num_days); % 3D array to store results

% Main processing
parpool("Threads",10)

mnt_pval           = zeros(num_vox_cortex,num_days,length(align_types));
num_perm_covari    = 1000;

for align_id = 1:length(align_types)
    align_type = align_types{align_id};
    align_type
    mnt_acf = load(strcat(base_dir,'/',align_type, '/outputs/gray/mnt_gray.mat'),'mnt_acf');
    mnt_acf = mnt_acf.mnt_acf;
    avg_mnt = squeeze(nanmean(mnt_acf,1));
    isc_mnt_corr = zeros(num_vox_cortex,2);
    for scan_id = 1:num_days
        scan_id
        other_scan_id = 3 - scan_id;
        if strcmp(align_type,'anatomical')
            isc = load(sprintf('%s%s/outputs/%s/%s_isc_scan_%s.mat', base_dir, align_type, 'gray', align_type, num2str(scan_id)));
        else
            isc = load(sprintf('%s%s/outputs/%s/%s_isc_parc_100_scan_%s.mat', base_dir, align_type, 'gray', align_type, num2str(scan_id)));
        end
        name = fieldnames(isc);
        isc = isc.(name{1});

        parfor vert = 1:num_vox_cortex
            [isc_herit_mnt(vert, scan_id), ~,~, isc_herit_mnt_jack(vert,:,scan_id,align_id)] = h2_mat(isc(:,:,vert), kinship(subjs_movie, subjs_movie), [covari_motion(subjs_movie, [1, 2, scan_id + 2]), squeeze(mnt_acf(vert, :, other_scan_id))'], num_perm,[]);
        end

        if num_perm_covari == 1000
            perm_order = zeros(num_subj_movie,num_perm_covari);
            for perm_id = 1:num_perm_covari
                rng(perm_id)
                perm_order(:,perm_id) = randperm(num_subj_movie);
            end
            
            parfor vert = 1:num_vox_cortex
                % Temporary variable for this iteration
                temp_null_herit = zeros(num_perm_covari, 1);

                if mod(vert,1000) == 0
                    disp(['vert = ' num2str(vert)])
                    disp(['align_type = ' align_type])
                    disp(['scan_id = ' num2str(scan_id)])
                end

                isc_temp = isc(:,:,vert);
                mnt_temp = squeeze(mnt_acf(vert, :, other_scan_id))';

                for perm_id = 1:num_perm_covari
                    temp_null_herit(perm_id) = h2_mat(isc_temp, kinship(subjs_movie, subjs_movie), [covari_motion(subjs_movie, [1, 2, scan_id + 2]), mnt_temp(perm_order(:,perm_id))], 0,[]);
                end

                null_herit(:,vert,scan_id,align_id) = temp_null_herit;
                mnt_pval(vert,scan_id,align_id)     = (sum(temp_null_herit>=abs(isc_herit_mnt(vert,scan_id)))+1)/(num_perm_covari+1);
            end
            
        end
        isc_herit_mnt(isc_herit_mnt==1) = NaN;

        % Save Cifti
        saveCifti(isc_herit_mnt(:, scan_id)', strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'), wb_path, scalar_template, medial_mask);

        % Make difference ciftis
        if strcmp(align_type,'anatomical')
            no_covari_cifti = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_400_scan_', num2str(scan_id), '.dscalar.nii'));
        else
            no_covari_cifti = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_100_scan_', num2str(scan_id), '.dscalar.nii'));
        end
        mnt_cifti = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'));
        saveCifti(mnt_cifti.cdata(medial_mask == 1)' - no_covari_cifti.cdata(medial_mask == 1)', strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_mnt_diff_scan_', num2str(scan_id), '.dscalar.nii'), wb_path, scalar_template, medial_mask);

        % Now correlate ISC with MNT
        parfor vox = 1:num_vox_cortex
            isc_mnt_corr(vox,scan_id) = corr(hvec(mnt_acf_out_big(:,:,vox,scan_id)),hvec(isc(:,:,vox)),'rows','complete');
        end
        saveCifti(isc_mnt_corr(:, scan_id)', strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_mnt_corr_scan_', num2str(scan_id), '.dscalar.nii'), wb_path, scalar_template, medial_mask);

    end
end


if num_perm_covari == 1000
    save([data_dir '/isc_heritability/null_herit_NT.mat'],'null_herit')
    save([data_dir '/isc_heritability/mnt_pval.mat'],'mnt_pval')
else
    load([data_dir '/isc_heritability/mat_files/supporting_files/mnt_pval.mat'])
    load([data_dir '/isc_heritability/mat_files/supporting_files/null_herit_NT.mat'])

    for align_id = 1:length(align_types)
        align_type = align_types{align_id};

        mnt_acf = load(strcat(base_dir,'/',align_type, '/outputs/gray/mnt_gray.mat'),'mnt_acf');
        mnt_acf = mnt_acf.mnt_acf;
        avg_mnt = squeeze(nanmean(mnt_acf,1));
        isc_mnt_corr = zeros(num_vox_cortex,2);
        for scan_id = 1:num_days
            data = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'));
            if align_id == 2
                nomnt_data = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_scan_', num2str(scan_id), '.dscalar.nii'));
            else
                nomnt_data = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_100_scan_', num2str(scan_id), '.dscalar.nii'));
            end

            isc_herit_m = data.cdata(medial_mask==1);
            isc_herit = nomnt_data.cdata(medial_mask==1);

            null_herit(null_herit==1) = NaN;
            avg_herit_m = nanmean(isc_herit_m);
            avg_herit = nanmean(isc_herit);
            obs_diff_avg = avg_herit_m - avg_herit;
            obs_diff = isc_herit_m - isc_herit;

            parfor vert = 1:num_vox_cortex
                null_diffs = null_herit(:,vert,scan_id,align_id)-isc_herit(vert);
                mnt_pval(vert,scan_id,align_id) = (sum(abs(null_diffs)>=abs(obs_diff(vert,1)))+1)/(num_perm_covari+1);
            end

            avg_null(:,scan_id,align_id) = squeeze(mean(null_herit(:,:,scan_id,align_id),2));
            null_diffs_all=avg_null(:,scan_id,align_id)-avg_herit;
            avg_diff(scan_id,align_id) = (sum(abs(null_diffs_all)>=abs(obs_diff_avg))+1)/(num_perm_covari+1);
        end
    end
end




for scan_id = 1:num_days
    piecewise_pval(:,scan_id)          = fdr_bh(mnt_pval(:,scan_id,1));
end
piecewise_pval_fdr      = piecewise_pval(:,1).*piecewise_pval(:,2);
piecewise_pval_fdr_perc = sum(piecewise_pval_fdr)/length(piecewise_pval_fdr);


for scan_id = 1:num_days
    anatomical_pval(:,scan_id)          = fdr_bh(mnt_pval(:,scan_id,2));
end
anatomical_pval_fdr      = anatomical_pval(:,1).*anatomical_pval(:,2);
anatomical_pval_fdr_perc = sum(anatomical_pval_fdr)/length(anatomical_pval_fdr);

mnt_data_out = zeros(num_vox_cortex,2,2);
for align_id = 1:length(align_types)
    align_type = align_types{align_id};
    for scan_id = 1:num_days
        if align_id == 1
            no_covari_cifti = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_100_scan_', num2str(scan_id), '.dscalar.nii'));
        else
            no_covari_cifti = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_400_scan_', num2str(scan_id), '.dscalar.nii'));
        end
        mnt_cifti       = ciftiopen(strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'));
        nocovari_data   = no_covari_cifti.cdata(medial_mask == 1)';
        mnt_data        = mnt_cifti.cdata(medial_mask == 1)';
        mnt_data(mnt_data==1) = NaN;
        nocovari_data(nocovari_data==1) = NaN;
        mnt_data_out(:,scan_id,align_id) = mnt_data - nocovari_data ;
        mnt_perc(:,scan_id,align_id) = (mnt_data - nocovari_data)./nocovari_data;
        saveCifti(mnt_perc(:,scan_id), strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_mnt_diff_perc_scan_', num2str(scan_id), '.dscalar.nii'), wb_path, scalar_template, medial_mask);
        saveCifti(mnt_data, strcat(base_dir, align_type, '/outputs/gray/', align_type, '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'), wb_path, scalar_template, medial_mask);
    end
end

mnt_med = squeeze(nanmean(mnt_data_out));
% Now make bar plot comparing anatomical/piecewise/int heritability
anatomical_isc_herit_nocovari = zeros(1,num_days);
anatomical_isc_herit_mnt      = zeros(1,num_days);
piecewise_isc_herit_mnt       = zeros(1,num_days);

% First load means
for scan_id = 1:num_days
    data  = ciftiopen(strcat(base_dir, 'anatomical', '/outputs/gray/', 'anatomical', '_isc_herit_gray_400_scan_', num2str(scan_id), '.dscalar.nii'));
    anatomical_isc_herit_nocovari(:,scan_id) = median(data.cdata(medial_mask==1));

    data = ciftiopen(strcat(base_dir, 'anatomical', '/outputs/gray/', 'anatomical', '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'));
    anatomical_isc_herit_mnt(:,scan_id) = median(data.cdata(medial_mask==1));

    data = ciftiopen(strcat(base_dir, 'piecewise', '/outputs/gray/', 'piecewise', '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'));
    piecewise_isc_herit_mnt(:,scan_id) = median(data.cdata(medial_mask==1));
end


% Next get SEs
for scan_id = 1:num_days
    data = ciftiopen(strcat(base_dir, 'anatomical', '/outputs/gray/', 'anatomical', '_isc_herit_gray_400_scan_', num2str(scan_id), '.dscalar.nii'));
    anatomical_isc_herit_nocovari(:,scan_id) = data.cdata(medial_mask==1);

    data = ciftiopen(strcat(base_dir, 'anatomical', '/outputs/gray/', 'anatomical', '_isc_herit_gray_mnt_scan_', num2str(scan_id), '.dscalar.nii'));
    anatomical_isc_herit_mnt(:,scan_id) = data.cdata(medial_mask==1);
end

isc_herit_gray            = zeros(num_vox_cortex,num_days);
isc_herit_gray_perm       = zeros(num_vox_cortex,num_days);
isc_herit_gray_int        = zeros(num_vox_cortex,num_days);
isc_herit_gray_mnt_perm   = zeros(num_vox_cortex,num_days);

isc_herit_gray_mnt_covari_perm   = zeros(num_vox_cortex,num_perm_covari,2);
num_perm_covari = 10000;
for scan_id = 1:num_days
    if scan_id == 1
        other_scan_id = 2;
    else
        other_scan_id = 1;
    end

    isc  = load(sprintf('%s%s/outputs/%s/%s_isc_parc_100_scan_%s.mat', base_dir,'anatomical','gray','anatomical', num2str(scan_id)));
    name = fieldnames(isc);
    isc  = isc.(name{1});

    parfor vert = 1:num_vox_cortex
        [isc_herit_gray_mnt(vert,scan_id), ~]   = h2_multi(isc(:,:,vert), kinship(subjs_movie,subjs_movie), [covari_motion(subjs_movie,[1,2,scan_id+2]), squeeze(mnt_acf(vert,:,other_scan_id))'], 0);
    end
    isc_herit_gray_int(isc_herit_gray_int==1) = NaN;
    isc_herit_gray_mnt(isc_herit_gray_mnt==1) = NaN;
end

% Make difference ciftis

for scan_id = 1:num_days
    no_covari_cifti = ciftiopen(strcat(base_dir,'piecewise/outputs/gray/anatomical_isc_herit_gray_400_scan_',num2str(scan_id),'.dscalar.nii'));
    mnt_cifti       = ciftiopen(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit_gray_mnt_scan_',num2str(scan_id),'.dscalar.nii'));
    saveCifti(mnt_cifti.cdata(medial_mask==1)'-no_covari_cifti.cdata(medial_mask==1)', strcat(base_dir,'piecewise/outputs/gray/isc_herit_gray_mnt_diff_scan_',num2str(scan_id),'.dscalar.nii'), wb_path, scalar_template,medial_mask);

    mnt_otherscan_cifti = ciftiopen(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit_gray_mnt_otherscan_',num2str(scan_id),'.dscalar.nii'));
    saveCifti(mnt_otherscan_cifti.cdata(medial_mask==1)'-no_covari_cifti.cdata(medial_mask==1)', strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_herit_gray_mnt_diff_otherscan_',num2str(scan_id),'.dscalar.nii'), wb_path, scalar_template,medial_mask);
end

mnt_pval_fdr      = fdr_bh(mnt_pval);
mnt_pval_fdr_perc = nansum(mnt_pval_fdr(:,1).*(mnt_pval_fdr(:,2)))/num_vox_cortex;

%% Compare heritability across different anatomical parcellation resolutions

anatomical_herit_parc_multires = zeros(num_parcs,num_days);
anatomical_parc_se             = zeros(num_parcs,num_days);

for parc_id = 1:length(parc_reses)
    for scan_id = 1:num_days
        anatomical_herit_parc_multires(parc_id,scan_id) = nanmean(anatomical_isc_herit_parc{parc_id,1}(:,scan_id));
        [~,anatomical_parc_se(parc_id,scan_id)]         = compareHeritability(0,nanmean(anatomical_isc_herit_parc_jack{parc_id,1}(:,:,scan_id)),0,nanmean(anatomical_isc_herit_parc_jack{parc_id,1}(:,:,scan_id)));
    end
end

% FIGURE 2 SCATTER
hold on
scatter(avg_areas,anatomical_herit_parc_multires(:,1),100,'filled','MarkerFaceColor',[70/255 130/255 180/255],'MarkerEdgeColor','k')
scatter(avg_areas,anatomical_herit_parc_multires(:,2),100,'filled','MarkerFaceColor',[107/255 142/255 35/255],'MarkerEdgeColor','k')

for parc_id = 1:num_parcs
    errorbar(avg_areas(parc_id),anatomical_herit_parc_multires(parc_id,1),anatomical_parc_se(parc_id,1),'Color',[70/255 130/255 180/255],'LineStyle','-','LineWidth',2)
    errorbar(avg_areas(parc_id),anatomical_herit_parc_multires(parc_id,2),anatomical_parc_se(parc_id,2),'Color',[107/255 142/255 35/255],'LineStyle','-','LineWidth',2)
end

subplot(2,1,1)
hold on

% Calculate means and SEMs
mean_isc_parc_day1 = anatomical_herit_parc_multires(:,1);
sem_isc_parc_day1  = anatomical_parc_se(:,1);
mean_isc_parc_day2 = anatomical_herit_parc_multires(:,2);
sem_isc_parc_day2  = anatomical_parc_se(:,2);

% boundedline plots
hl1 = boundedline(avg_areas, mean_isc_parc_day1, sem_isc_parc_day1, 'alpha', 'cmap', [70/255 130/255 180/255]);
hl2 = boundedline(avg_areas, mean_isc_parc_day2, sem_isc_parc_day2, 'alpha', 'cmap', [107/255 142/255 35/255]);

% Doubling the line thickness in boundedline plots
set(hl1, 'LineWidth', 4);
set(hl2, 'LineWidth', 4);

% Scatter plots for mean values
scatter(avg_areas, mean_isc_parc_day1, 100, 'filled', 'MarkerFaceColor', [70/255 130/255 180/255], 'MarkerEdgeColor', 'k');
scatter(avg_areas, mean_isc_parc_day2, 100, 'filled', 'MarkerFaceColor', [107/255 142/255 35/255], 'MarkerEdgeColor', 'k');

% Axes settings
ylim([.05 .08]); xlim([0 1200]);
xlabel('Average Parcel Area (mm^2)'); ylabel('Average Heritability (h^2)');
yticks([.05 .06 .07 .08]); ax = gca; set(ax, 'TickDir', 'out'); set(ax, 'FontSize', 20); set(ax, 'LineWidth', 3); set(ax, 'Layer', 'bottom');
grid off; legend({'Day 1', 'Day 2'}, 'Location', 'best'); pbaspect([2 1 1]); set(gcf, 'position', [100,100,600,1200]);


saveas(gcf,[data_dir '/isc_heritability/figures/fig2_parc_areascatter.fig'])
close all

%% ISC vs neural-timescale correlation: precompute pairwise NT sums (feeds QAP below)

num_perms = 1000;

if exist(strcat(base_dir,'anatomical/outputs/gray/mnt_gray.mat'),'file') == 2
    load(strcat(base_dir,'anatomical/outputs/gray/mnt_gray.mat'))
else
    mnt_acf = zeros(num_vox_cortex,num_subj_movie,2);
    for scan_id = [1,2]
        for subj = 1:num_subj_movie
            subj
            subj_data = readNPY(strcat([data_dir '/isc_heritability/data/anatomical/inputs/gray/anatomical_movie_scan_'], num2str(scan_id), '_subj_', num2str(subj), '.npy'));
            mnt_acf(:,subj,scan_id) = calculate_int(subj_data);
        end
    end
    save(strcat(base_dir,'anatomical/outputs/gray/mnt_gray.mat'),'mnt_acf')
end

mnt_avg = squeeze(nanmean(mnt_acf,1));
mnt_acf_out_big = zeros(num_subj_movie,num_subj_movie,num_vox_cortex,num_days);
for scan_id = 1:num_days
    for vox = 1:num_vox_cortex
        for subj1 = 1:num_subj_movie
            for subj2 = 1:num_subj_movie
                if subj1~=subj2
                    mnt_acf_out_big(subj1,subj2,vox,scan_id) = squeeze(mnt_acf(vox,subj1,scan_id))+squeeze(mnt_acf(vox,subj2,scan_id));
                else
                    mnt_acf_out_big(subj1,subj2,vox,scan_id) = NaN;
                end
            end
        end
    end
end
clear mnt_acf

% Number of permutations
num_perms = 10000;

% Precompute hvec for all voxels and scans
mnt_acf_hvec = cell(num_vox_cortex, 2);
parfor vox = 1:num_vox_cortex
    for scan_id = 1:num_days
        mnt_acf_hvec{vox, scan_id} = hvec(mnt_acf_out_big(:, :, vox, scan_id));
    end
end

%% ISC vs neural-timescale correlation (family-compliant QAP, 1000 permutations)
align_type = 'anatomical';
isc_mnt_corr = zeros(num_vox_cortex, 2);
isc_mnt_pvals = zeros(num_vox_cortex, 2); % To store p-values

% --- SETTINGS & PREPARATION ---
num_perms = 1000;
num_vox_cortex = 59412;
num_active_subjs = length(subjs_movie);
num_pairs = (num_active_subjs * (num_active_subjs - 1)) / 2;

% family_data: num_subjx2 [SubjectID, FamilyID]
[~, loc] = ismember(subjs_movie, family_ids(:, 1));
curr_family_ids = family_ids(loc, 2);

[unique_fams, ~, fam_map] = unique(curr_family_ids);
num_families = length(unique_fams);

% --- 1. PRECOMPUTE HCP-COMPLIANT PERMUTATION INDICES ---
all_perm_indices = zeros(num_active_subjs, num_perms);
for p_id = 1:num_perms
    rng(p_id);
    fam_perm = randperm(num_families);
    curr_p = zeros(num_active_subjs, 1);
    write_ptr = 1;
    for f = 1:num_families
        orig_member_indices = find(fam_map == fam_perm(f));
        within_fam_perm = orig_member_indices(randperm(length(orig_member_indices)));
        num_members = length(within_fam_perm);
        curr_p(write_ptr : write_ptr + num_members - 1) = within_fam_perm;
        write_ptr = write_ptr + num_members;
    end
    all_perm_indices(:, p_id) = curr_p;
end

% --- 2. VOXEL-WISE QAP ANALYSIS ---
for scan_id = 1:num_days
    other_scan_id = 3 - scan_id;
    fprintf('Processing Scan %d (using ISC from Scan %d)...\n', scan_id, other_scan_id);

    isc_data = load(sprintf('%s%s/outputs/%s/%s_isc_scan_%s.mat', ...
        base_dir, 'anatomical', 'gray', 'anatomical', num2str(other_scan_id)));
    name = fieldnames(isc_data);
    isc = isc_data.(name{1});

    % Slice MNT to avoid broadcast overhead in parfor
    current_scan_mnt = mnt_acf_hvec(:, scan_id);

    isc_mnt_corr = zeros(num_vox_cortex, 1);
    isc_mnt_pvals = ones(num_vox_cortex, 1);

    
    parfor vox = 1:num_vox_cortex
        if mod(vox, 10000) == 0, fprintf('Voxel: %d\n', vox); end

        % A. Get and Pre-rank Observed Data
        % We use hvec to get the lower triangle of the 178x178 matrix
        vox_isc_matrix = isc(:, :, vox);
        obs_isc_vec = hvec(vox_isc_matrix);
        obs_mnt_vec = current_scan_mnt{vox};

        % Calculate Observed Spearman Correlation
        isc_mnt_corr(vox) = corr(obs_mnt_vec, obs_isc_vec, 'rows', 'complete', 'type', 'spearman');

        if num_perms > 0
            % B. Optimization: Rank the matrix once
            % We rank the full matrix, then permute, then extract hvec
            ranked_isc_mat = reshape(tiedrank(vox_isc_matrix(:)), num_active_subjs, num_active_subjs);

            % C. Pre-rank and Standardize the Predictor (MNT) for fast Pearson
            % Standardizing allows us to use dot product instead of full corr()
            mnt_rank = tiedrank(obs_mnt_vec);
            mnt_rank_std = (mnt_rank - mean(mnt_rank)) / std(mnt_rank);

            % D. Vectorized Null Generation
            % Instead of a 1000-iteration loop, we build a matrix of permuted ISCs
            perm_isc_mtx = zeros(num_pairs, num_perms);
            for p_id = 1:num_perms
                p_idx = all_perm_indices(:, p_id);
                perm_isc_mtx(:, p_id) = hvec(ranked_isc_mat(p_idx, p_idx));
            end

            % E. Ultra-Fast Matrix-Vector Pearson Correlation
            % Rank and standardize the permuted matrix columns
            perm_isc_mtx = (perm_isc_mtx - mean(perm_isc_mtx, 1)) ./ std(perm_isc_mtx, 0, 1);

            % Resulting null_corrs is 1000x1 vector
            null_corrs = (mnt_rank_std' * perm_isc_mtx) / (num_pairs - 1);

            % F. Compute P-Value (Two-Tailed)
            isc_mnt_pvals(vox) = (sum(abs(null_corrs) >= abs(isc_mnt_corr(vox))) + 1) / (num_perms + 1);
        end
    end
    

    % --- 3. SAVE RESULTS ---
    corr_out = strcat(base_dir, 'anatomical', '/outputs/gray/', 'anatomical', '_isc_mnt_corr_QAP_scan_', num2str(scan_id), '.dscalar.nii');
    pval_out = strcat(base_dir, 'anatomical', '/outputs/gray/', 'anatomical', '_isc_mnt_pvals_QAP_scan_', num2str(scan_id), '.dscalar.nii');

    saveCifti(isc_mnt_corr', corr_out, wb_path, scalar_template, medial_mask);
    saveCifti(isc_mnt_pvals', pval_out, wb_path, scalar_template, medial_mask);
end

for scan_id = 1:num_days
    xx= ciftiopen(strcat(base_dir, 'anatomical', '/outputs/gray/', 'anatomical', '_isc_mnt_pvals_QAP_scan_', num2str(scan_id), '.dscalar.nii'));
    pval_mat(:,scan_id)  = xx.cdata;
    pval_mat_fdr(:,scan_id) = fdr_bh(pval_mat(:,scan_id));
end

for scan_id = 1:num_days
    isc_mnt_pvals_fdr(:,scan_id) = fdr_bh(isc_mnt_pvals(:,scan_id));
end

isc_mnt_pvals_fdr_sig       = isc_mnt_pvals_fdr(:,1).*isc_mnt_pvals_fdr(:,2);
isc_mnt_pvals_fdr_sig_perc  = sum(isc_mnt_pvals_fdr_sig)/num_vox_cortex;
isc_mnt_max                 = max(isc_mnt_corr);
isc_mnt_mean                = mean(isc_mnt_corr);

%% Subsampling (Fig. S3)
% Set the seed for reproducibility
rng(42);  % You can choose any seed value

dataType = 'parc';
taskType = 'movie';
alignType = 'anatomical';
curr_parc = '400';
Subjects = subjs_movie;
num_subjs = num_subj_movie;
num_parcels = 400;

% Define the percentages of families to drop
drop_percentages = [5, 10, 15, 20, 25,30, 40, 50, 75,90];
drop_fractions = drop_percentages / 100;

% Get unique family IDs
unique_families = unique(movie_fam_ids);
num_families = length(unique_families);

% Number of permutations
num_permutations = 100;

% Initialize the results array
isc_herit_subsample = zeros(num_parcels, 2, length(drop_fractions), num_permutations);

% Load subject data (once)
subj_data_full = cell(1,2);
for scan_id = 1:num_days
    % Load data for all subjects
    subj_data = loadData(num_subjs, base_dir, alignType, taskType, dataType, curr_parc, scan_id, num_trs);
    subj_data = permute(subj_data,[2,1,3]);  % [voxels x subjects x time]
    subj_data_full{scan_id} = subj_data;
end

% Calculate ISC for all subjects (once)
isc_full = cell(1,2);
for scan_id = 1:num_days
    isc = zeros(num_subjs, num_subjs, num_parcels);
    subj_data = subj_data_full{scan_id};  % [voxels x subjects x time]
    for vox = 1:num_parcels
        isc(:,:,vox) = corr(squeeze(subj_data(:,vox,:)));
    end
    isc_full{scan_id} = isc;  % Store the ISC matrices
end

for scan_id = 1:num_days
    for grayordinate = 1:num_parcels
        isc_matrix = squeeze(isc_full{scan_id}(:,:,grayordinate));
        [isc_herit_full(grayordinate,scan_id), ~, ~, ~] = h2_mat(squeeze(isc(:,:,grayordinate)), kinship(Subjects,Subjects), covari_motion(Subjects,[1,2,scan_id+2]), 0,[]);
    end
end

for scan_id = 1:num_days
    isc = isc_full{scan_id};  % Precomputed ISC matrices for this scan
    for frac_idx = 1:length(drop_fractions)
        frac = drop_fractions(frac_idx);

        % Calculate the number of families to drop
        num_families_to_drop = round(frac * num_families);

        parfor perm_id = 1:num_permutations
            % Set seed for each permutation
            rng(42 + perm_id);

            % Randomly select families to drop
            families_to_drop = randsample(unique_families, num_families_to_drop);

            % Find indices of subjects to exclude (belonging to dropped families)
            subjects_to_exclude = ismember(movie_fam_ids, families_to_drop);

            % Indices of subjects to include
            included_subjects = ~subjects_to_exclude;

            % Get the indices of included subjects
            idx_included_subjects = find(included_subjects);
            num_subjs_included = length(idx_included_subjects);

            % Get the subject IDs to include
            Subjects_included = Subjects(included_subjects);

            % Extract the relevant ISC subset
            isc_subset = isc(included_subjects, included_subjects, :);  % [subjects x subjects x parcels]

            % Prepare kinship and covariate matrices for included subjects
            kinship_sub = kinship(Subjects_included, Subjects_included);
            covariates_sub = covari_motion(Subjects_included, [1, 2, scan_id+2]);

            % Recalculate heritability
            isc_herit = zeros(num_parcels, 1);

            for grayordinate = 1:num_parcels
                isc_matrix = squeeze(isc_subset(:,:,grayordinate));
                [isc_herit(grayordinate,1), ~, ~, ~] = h2_mat(isc_matrix, kinship_sub, covariates_sub, 0, []);
            end

            % Store the results
            isc_herit_subsample(:, scan_id, frac_idx, perm_id) = isc_herit;
        end
    end
end

% Number of fractions and permutations
num_fractions = length(drop_percentages);
num_permutations = size(isc_herit_subsample, 4);

% Initialize arrays to store MAE and correlation results
MAE_results = zeros(num_days, num_fractions, num_permutations);          % [scan_id x fractions x permutations]
corr_results = zeros(num_days, num_fractions, num_permutations);         % [scan_id x fractions x permutations]

for scan_id = 1:num_days
    % Full-sample heritability for this scan
    isc_herit_full_scan = isc_herit_full(:, scan_id);

    for frac_idx = 1:num_fractions
        for perm_id = 1:num_permutations
            % Subsampled heritability for this permutation
            isc_herit_sub = isc_herit_subsample(:, scan_id, frac_idx, perm_id);

            % Compute absolute error across parcels
            abs_error = abs(isc_herit_sub - isc_herit_full_scan);

            % Compute Mean Absolute Error (MAE)
            MAE = mean(abs_error);
            MAE_results(scan_id, frac_idx, perm_id) = MAE;

            % Compute spatial correlation across parcels
            corr_coef = corr(isc_herit_sub, isc_herit_full_scan(:), 'Type', 'Spearman');
            corr_results(scan_id, frac_idx, perm_id) = corr_coef;
        end
    end
end

% Compute the average MAE and correlation across permutations for each fraction level
mean_MAE = mean(MAE_results, 3);           % [scan_id x fractions]
mean_corr = mean(corr_results, 3);         % [scan_id x fractions]
std_MAE = std(MAE_results, 0, 3);          % Standard deviation across permutations
std_corr = std(corr_results, 0, 3);        % Standard deviation across permutations

% Plotting Mean Absolute Error (MAE) with shading and circles with black edges
figure;
subplot(1,2,1)
hold on;

% Plot shaded area for standard deviation (Scan 1)
fill([drop_percentages, fliplr(drop_percentages)], ...
    [mean_MAE(1, :) + std_MAE(1, :), fliplr(mean_MAE(1, :) - std_MAE(1, :))], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot shaded area for standard deviation (Scan 2)
fill([drop_percentages, fliplr(drop_percentages)], ...
    [mean_MAE(2, :) + std_MAE(2, :), fliplr(mean_MAE(2, :) - std_MAE(2, :))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot mean MAE with circles and black edges
plot(drop_percentages, mean_MAE(1, :), '-o', 'LineWidth', 2, 'Color', 'b', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'DisplayName', 'Scan 1');
plot(drop_percentages, mean_MAE(2, :), '-o', 'LineWidth', 2, 'Color', 'r', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'DisplayName', 'Scan 2');

xlabel('Percentage of Families Dropped');
ylabel('Mean Absolute Error (MAE)');
title('Mean Absolute Error vs. Percentage of Families Dropped');
hold off;
grid off
% Customize axes
ax = gca;
set(ax, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1, 'Layer', 'bottom')
set(ax, 'LineWidth', 1, 'FontName', 'Arial');  % Increase line width, remove k marks, and set font

grid off
legend({'', '', 'Day 1', 'Day 2', ''},'Location','northwest')
pbaspect([1, 1, 1])

% Plotting Average Spatial Correlation with shading and circles with black edges
subplot(1,2,2)
hold on;

% Plot shaded area for standard deviation (Scan 1)
fill([drop_percentages, fliplr(drop_percentages)], ...
    [mean_corr(1, :) + std_corr(1, :), fliplr(mean_corr(1, :) - std_corr(1, :))], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot shaded area for standard deviation (Scan 2)
fill([drop_percentages, fliplr(drop_percentages)], ...
    [mean_corr(2, :) + std_corr(2, :), fliplr(mean_corr(2, :) - std_corr(2, :))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot mean correlation with circles and black edges
plot(drop_percentages, mean_corr(1, :), '-o', 'LineWidth', 2, 'Color', 'b', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'DisplayName', 'Scan 1');
plot(drop_percentages, mean_corr(2, :), '-o', 'LineWidth', 2, 'Color', 'r', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'DisplayName', 'Scan 2');

xlabel('Percentage of Families Dropped');
ylabel('Average Spatial Correlation (\rho)');
title('Spatial Correlation vs. Percentage of Families Excluded', 'FontSize', 20, 'FontName', 'Arial');
grid off;
hold off;
% Customize axes
ax = gca;
set(ax, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1, 'Layer', 'bottom')
set(ax, 'LineWidth', 1, 'FontName', 'Arial');  % Increase line width, remove k marks, and set font

grid off
legend({'', '', 'Day 1', 'Day 2', ''},'Location','northeast')
pbaspect([1, 1, 1])

saveas(gcf,[data_dir '/isc_heritability/figures/figS1_subsample.svg'])
close all

%% Compare hyperaligned and MSM ISC (Figure S6)
mean_piecewise_scaled = mean_piecewise(1:10)./avg_areas;
mean_connectivity_scaled = mean_connectivity(1:10)./avg_areas;

for scan_id = 1:num_days
    data                     = ciftiopen(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_gray_scan_',num2str(scan_id),'.dscalar.nii'));
    anat_isc_gray(:,scan_id) = data.cdata;
    for parc_id = 1:10
        xx = ciftiopen(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_gray_',parc_reses{parc_id},'_scan_',num2str(scan_id),'.dscalar.nii'));
        piecewise_isc_gray(:,parc_id,scan_id)   = xx.cdata;

        xx = ciftiopen(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_gray_',parc_reses{parc_id},'_scan_',num2str(scan_id),'.dscalar.nii'));
        connectivity_isc_gray(:,parc_id,scan_id) =xx.cdata;
        saveCifti(anat_isc_gray(:,scan_id) - piecewise_isc_gray(:,parc_id,scan_id), strcat([data_dir '/isc_heritability/data/anatomical/outputs/gray/anatomical_piecewise_isc_diff_parc_'],parc_reses{parc_id},'_day',num2str(scan_id),'.dscalar.nii'), wb_path, schaefer_400_dscalar,medial_mask);
        saveCifti(anat_isc_gray(:,scan_id) - connectivity_isc_gray(:,parc_id,scan_id), strcat([data_dir '/isc_heritability/data/anatomical/outputs/gray/anatomical_connectivity_isc_diff_parc_'],parc_reses{parc_id},'_day',num2str(scan_id),'.dscalar.nii'), wb_path, schaefer_400_dscalar,medial_mask);

    end
end

for scan_id = 1:num_days
    data                     = ciftiopen(strcat(base_dir,'anatomical/outputs/gray/anatomical_isc_gray_scan_',num2str(scan_id),'.dscalar.nii'));
    anat_isc_gray(:,scan_id) = data.cdata;
    for parc_id = 1:10
        xx = ciftiopen(strcat(base_dir,'piecewise/outputs/gray/piecewise_isc_gray_',parc_reses{parc_id},'_scan_',num2str(scan_id),'.dscalar.nii'));
        piecewise_isc_gray(:,parc_id,scan_id)   = xx.cdata;

        xx = ciftiopen(strcat(base_dir,'connectivity/outputs/gray/connectivity_isc_gray_',parc_reses{parc_id},'_scan_',num2str(scan_id),'.dscalar.nii'));
        connectivity_isc_gray(:,parc_id,scan_id) =xx.cdata;

        piecewise_isc_diff(parc_id,scan_id) = mean(mean(anat_isc_gray(:,scan_id) - piecewise_isc_gray(:,parc_id,scan_id)));
        connectivity_isc_diff(parc_id,scan_id) = mean(mean(anat_isc_gray(:,scan_id) - connectivity_isc_gray(:,parc_id,scan_id)));
    end
end

% --- Data Preparation ---
% Concatenate the MSM baseline with the piecewise and connectivity data
piecewise_data = [nanmedian(anat_isc_gray(:,1)) squeeze(nanmedian(piecewise_isc_gray(:,:,2),1))];
connectivity_data = [nanmedian(anat_isc_gray(:,1)) squeeze(nanmedian(connectivity_isc_gray(:,:,2),1))];

% --- Plotting ---
figure;
hold on; % Allows multiple lines on the same axes

% Plot both lines with 'LineWidth' set to 2 (adjust as needed)
plot(piecewise_data, 'LineWidth', 2, 'DisplayName', 'Piecewise Hyperalignment');
plot(connectivity_data, 'LineWidth', 2, 'DisplayName', 'Connectivity Hyperalignment');

% --- Aesthetics ---
ylim([0.018 0.026]);
ylabel('ISC (r)');
xlabel('Parcellation Resolution');
title('Hyperalignment Comparison');

% Set X-ticks and Labels
xticks(1:11);
xticklabels({'MSM','100','200','300','400','500','600','700','800','900','1000'});

% Remove top and right axes (spines)
set(gca, 'Box', 'off');

% Add Legend
legend('Location', 'northeast');

hold off;

% Define constants
num_scans = 2;
orange_col = [254/255 97/255 0/255];
purple_col = [120/255 94/255 240/255];

figure;
set(gcf, 'position', [100, 100, 800, 1200]);

for scan_id = 1:num_scans
    subplot(2, 1, scan_id)
    hold on

    % --- 1. Data Preparation (Medians) ---

    % MSM (Anatomical Baseline)
    med_msm = double(nanmedian(anat_isc_gray(:, scan_id)));

    % Piecewise Hyperalignment
    med_pw = double(squeeze(nanmedian(piecewise_isc_gray(:,:,scan_id), 1))');

    % Connectivity Hyperalignment
    med_conn = double(squeeze(nanmedian(connectivity_isc_gray(:,:,scan_id), 1))');

    avg_areas_db = double(avg_areas);

    % --- 2. Plotting Connecting Lines ---
    % Since we removed boundedline, we use standard plot calls for the lines
    % Using 'LineWidth', 4 to match the thickness of your original fit lines
    plot(avg_areas_db, med_pw, 'Color', orange_col, 'LineWidth', 4);
    plot(avg_areas_db, med_conn, 'Color', purple_col, 'LineWidth', 4);

    % --- 3. Median Scatter Points ---
    % Piecewise
    hp1 = scatter(avg_areas_db, med_pw, 150, 'filled', ...
        'MarkerFaceColor', orange_col, 'MarkerEdgeColor', 'k', 'LineWidth', 3);

    % Connectivity
    hp2 = scatter(avg_areas_db, med_conn, 150, 'filled', ...
        'MarkerFaceColor', purple_col, 'MarkerEdgeColor', 'k', 'LineWidth', 3);

    % MSM point at x=0
    scatter(0, med_msm, 170, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

    % --- 4. Aesthetics & Formatting ---
    xlim([-100 1200]);
    ylim([0.015 0.03]);
    ylabel('ISC (r)');
    xlabel('Average Parcel Area (mm^2)');

    % Match the heavy axis style from your heritability code
    ax = gca;
    set(ax, 'TickDir', 'out', 'FontSize', 32, 'LineWidth', 5, 'Layer', 'bottom');

    xticks([0 400 800 1200]);
    xticklabels({'MSM', '400', '800', '1200'});

    % Legend only on the top plot
    if scan_id == 1
        legend([hp1, hp2], {'Piecewise', 'Connectivity'}, 'Location', 'best', 'FontSize', 20);
        legend boxoff;
    end

    grid off;
    pbaspect([1.5 1 1]);
end

%% Re-run without GSR (Figure S4)
base_dir = [data_dir '/HCP_7T/'];
num_parcels = 400;

% HCP 7T Movie TR counts (standard lengths)
% Movie 1: 921, Movie 2: 918, Movie 3: 915, Movie 4: 901
tr_counts = [921, 918, 915, 901];

% Get subject list
subjects = dir(fullfile(base_dir, '*'));
subjects = subjects([subjects.isdir]);
subjects = subjects(~strncmp({subjects.name}, '.', 1));
num_subjs = length(subjects);
subjects(2,:) = [];
subjects(186,:) = [];
subjects(185,:) = [];

% --- PRE-ALLOCATION ---
% Creating four 3D matrices: [Parcels x Time x Subjects]
movie1_data = nan(num_parcels, tr_counts(1), num_subjs);
movie2_data = nan(num_parcels, tr_counts(2), num_subjs);
movie3_data = nan(num_parcels, tr_counts(3), num_subjs);
movie4_data = nan(num_parcels, tr_counts(4), num_subjs);

% Store in a cell array for easier indexing during the loop
all_movies = {movie1_data, movie2_data, movie3_data, movie4_data};

% --- LOADING LOOP ---

% Pre-allocate a status tracker: [num_subj subjects x 4 movies]
% 0 = Data OK, 1 = Contains NaNs, 2 = Missing File entirely
subj_nan_status = zeros(num_subj, 4);

tic
for i = 1:num_subj
    subj_id = subjects(i).name;
    results_path = fullfile(base_dir, subj_id, 'MNINonLinear', 'Results');

    % Find all parcellated files for this subject
    ptseries_files = dir(fullfile(results_path, 'tfMRI_MOVIE*', '*hp2000_clean.ptseries.nii'));

    % Track which movies were found for this specific subject
    movies_found = false(1, 4);

    for f = 1:length(ptseries_files)
        fname = ptseries_files(f).name;

        tokens = regexp(fname, 'tfMRI_MOVIE(\d)', 'tokens');
        if isempty(tokens), continue; end
        m_idx = str2double(tokens{1}{1});
        movies_found(m_idx) = true;

        try
            cifti_struct = ciftiopen(fullfile(ptseries_files(f).folder, fname), wb_path);
            current_data = cifti_struct.cdata;

            % --- CHECK FOR NANS ---
            if any(isnan(current_data), 'all')
                subj_nan_status(i, m_idx) = 1;
                fprintf('!! Warning: Subj %s Movie %d contains NaNs\n', subj_id, m_idx);
            end

            % Assign to the correct slice
            all_movies{m_idx}(1:size(current_data,1), 1:size(current_data,2), i) = current_data;

        catch ME
            warning('Failed to load %s: %s', fname, ME.message);
        end
    end

    % Mark missing files
    missing_idx = find(~movies_found);
    if ~isempty(missing_idx)
        subj_nan_status(i, missing_idx) = 2;
        for m = missing_idx
            fprintf('-- Missing: Subj %s Movie %d file not found\n', subj_id, m);
        end
    end
end

% Summary Report
fprintf('\n--- DATA QUALITY SUMMARY ---\n');
fprintf('Subjects with at least one NaN: %d\n', sum(any(subj_nan_status == 1, 2)));
fprintf('Subjects missing at least one run: %d\n', sum(any(subj_nan_status == 2, 2)));

xx = (any(subj_nan_status == 2, 2));
% Unpack cell array back into individual variables
movie1_data = all_movies{1};
movie2_data = all_movies{2};
movie3_data = all_movies{3};
movie4_data = all_movies{4};

% --- 1. PRE-PROCESSING ---
sub_data = cell(4,1);
sub_data{1} = movie1_data(:,:,1:num_subj);
sub_data{2} = movie2_data(:,:,1:num_subj);
sub_data{3} = movie3_data(:,:,1:num_subj);
sub_data{4} = movie4_data(:,:,1:num_subj);

% Define which TRs to cut from each of the 4 movie runs.
% Each movie clip is preceded by 20 seconds of rest. Here, we remove TRs that
% take place in those 20 seconds as well as in the first 20 seconds of each clip
% to account for potential onset transients.
cut_montage = 0;  % 0 keeps the trailing montage clip; reproduces num_trs = [1432, 1409]

% Initialize matrices
rest_trs   = cell(4,1);
movie_cut  = cell(4,1);
movie_mask = cell(4,1);

starts = cell(4,1);
ends = cell(4,1);

% Rest blocks- Run 1
starts{1,1} = [0,264.0833,505.75,713.7917,797.5833,901];
starts{1,1} = floor(starts{1,1} +1);
ends{1,1}   = [19.9583,284.0417,525.7083,733.75,817.5417,920.9583];
ends{1,1}   = floor(ends{1,1} +1);

for i = 1:length(starts{1,1})
    rest_trs{1,1} = [rest_trs{1,1} starts{1,1}(i):ends{1,1}(i)];
end

% Rest blocks- Run 2
starts{2,1} = [0,246.75,525.375,794.625,898];
starts{2,1} = floor(starts{2,1}+1);
ends{2,1}   = [19.9583,266.7083,545.3333,814.5417,917.9583];
ends{2,1}   = floor(ends{2,1}+1);

for i = 1:length(starts{2,1})
    rest_trs{2,1} = [rest_trs{2,1} starts{2,1}(i):ends{2,1}(i)];
end

% Rest blocks- Run 3
starts{3,1} = [0,200.5833,405.125,629.25,791.7917,895];
starts{3,1} = floor(starts{3,1}+1);
ends{3,1}   = [19.9583,220.5417,425.0833,649.2083,811.5417,914.9583];
ends{3,1}   = floor(ends{3,1}+1);

for i = 1:length(starts{3,1})
    rest_trs{3,1} = [rest_trs{3,1} starts{3,1}(i):ends{3,1}(i)];
end

% Rest blocks- Run 4
starts{4,1} = [0,252.3333,502.2083,777.4167,881];
starts{4,1} = floor(starts{4,1}+1);
ends{4,1}   = [19.9583,272.2917,522.1667,797.5417,900.9583];
ends{4,1}   = floor(ends{4,1}+1);

for i = 1:length(starts{4,1})
    rest_trs{4,1} = [rest_trs{4,1} starts{4,1}(i):ends{4,1}(i)];
end

% Movie blocks- Run 1
starts{1,2} = [20,284.0833,525.75,733.7917,817.5833];
ends{1,2}   = [264.0417,505.7083,713.75,797.5417,900.9583];
starts{1,2} = ceil(starts{1,2}+1);
ends{1,2}   = floor(ends{1,2});

if cut_montage == 1
    starts{1,2} = starts{1,2}(:,1:size(starts{1,2},2)-1);
    ends{1,2}   = ends{1,2}(:,1:size(ends{1,2},2)-1);
end

for i = 1:length(starts{1,2})
    movie_cut{1,1}     = [movie_cut{1,1} starts{1,2}(i):starts{1,2}(i)+19];
    movie_mask{1,1}{i} = [starts{1,2}(i)+20 starts{1,1}(i+1)-1];
end

% Movie blocks- Run 2
starts{2,2} = [20,266.75,545.375,814.5833];
ends{2,2}   = [246.7083,525.3333,794.5833,897.9583];
starts{2,2} = ceil(starts{2,2}+1);
ends{2,2}   = floor(ends{2,2});

if cut_montage == 1
    starts{2,2} = starts{2,2}(:,1:size(starts{2,2},2)-1);
    ends{2,2}   = ends{2,2}(:,1:size(ends{2,2},2)-1);
end

for i = 1:length(starts{2,2})
    movie_cut{2,1}     = [movie_cut{2,1} starts{2,2}(i):starts{2,2}(i)+19];
    movie_mask{2,1}{i} = [starts{2,2}(i)+20 starts{2,1}(i+1)-1];
end

% Movie blocks- Run 3
starts{3,2} = [20,220.5833,425.125,649.25,811.5833];
ends{3,2}   = [200.5417,405.0833,629.2083,791.75,894.9583];
starts{3,2} = ceil(starts{3,2}+1);
ends{3,2}   = floor(ends{3,2});

if cut_montage == 1
    starts{3,2} = starts{3,2}(:,1:size(starts{3,2},2)-1);
    ends{3,2}   = ends{3,2}(:,1:size(ends{3,2},2)-1);
end

for i = 1:length(starts{3,2})
    movie_cut{3,1}     = [movie_cut{3,1} starts{3,2}(i):starts{3,2}(i)+19];
    movie_mask{3,1}{i} = [starts{3,2}(i)+20 starts{3,1}(i+1)-1];
end

% Movie blocks- Run 4
starts{4,2} = [20,272.3333,522.2083,797.5833];
ends{4,2}   = [252.2917,502.1667,777.375,880.9583];
starts{4,2} = ceil(starts{4,2}+1);
ends{4,2}   = floor(ends{4,2});

if cut_montage == 1
    starts{4,2} = starts{4,2}(:,1:size(starts{4,2},2)-1);
    ends{4,2}   = ends{4,2}(:,1:size(ends{4,2},2)-1);
end

for i = 1:length(starts{4,2})
    movie_cut{4,1}     = [movie_cut{4,1} starts{4,2}(i):starts{4,2}(i)+19];
    movie_mask{4,1}{i} = [starts{4,2}(i)+20 starts{4,1}(i+1)-1];
end

% Get indices of all frames that need to be removed
all_cut      = cell(4,1);
all_cut{1,1} = unique(cat(2,rest_trs{1,1},movie_cut{1,1}));
all_cut{2,1} = unique(cat(2,rest_trs{2,1},movie_cut{2,1}));
all_cut{3,1} = unique(cat(2,rest_trs{3,1},movie_cut{3,1}));
all_cut{4,1} = unique(cat(2,rest_trs{4,1},movie_cut{4,1}));

for s_id = 1:4
    % FIRST: Remove rest/onset frames
    sub_data{s_id}(:, all_cut{s_id}, :) = [];

    % SECOND: Z-score the remaining movie frames
    for s = 1:num_subj
        % Dimension 2 is Time
        sub_data{s_id}(:,:,s) = zscore(sub_data{s_id}(:,:,s), 0, 2);
    end
end

% Concatenate after cleaning
sub_data_cat{1} = cat(num_days, sub_data{1}, sub_data{2});
sub_data_cat{2} = cat(num_days, sub_data{3}, sub_data{4});

% --- 2. ISC CALCULATION ---
gsr_isc = zeros(num_subj, num_subj, 400, 2);
for d_idx = 1:num_days
    parfor s1 = 1:num_subj
        % Get data for all parcels for subject 1
        data1 = sub_data_cat{d_idx}(:, :, s1);
        for s2 = 1:num_subj
            data2 = sub_data_cat{d_idx}(:, :, s2);
            % Vectorized ISC (faster and avoids the 400x400 matrix)
            % This is the same as diag(corr(data1', data2'))
            gsr_isc(s1, s2, :, d_idx) = mean(data1 .* data2, 2);
        end
    end
end

% --- 3. HERITABILITY ---
for parc_id = 1:400
    for day_id = 1:num_days
        % Use day_id consistently throughout this block
        current_isc = squeeze(gsr_isc(subjs_movie, subjs_movie, parc_id, day_id));
        current_motion = covari_motion(subjs_movie, [1, 2, day_id + 2]);

        isc_gsr_herit(parc_id, day_id) = h2_mat(current_isc, ...
            kinship(subjs_movie, subjs_movie), ...
            current_motion, 0, []);
    end
end

saveCifti(isc_gsr_herit(:,1), strcat([data_dir '/isc_heritability/data/piecewise/outputs/parc/piecewise_anatomical_isc_herit_GSR_parc_'],curr_parc,'_day',num2str(1),'.pscalar.nii'), wb_path, schaefer_400_pscalar_kong,medial_mask);
saveCifti(isc_gsr_herit(:,2), strcat([data_dir '/isc_heritability/data/piecewise/outputs/parc/piecewise_anatomical_isc_herit_GSR_parc_'],curr_parc,'_day',num2str(2),'.pscalar.nii'), wb_path, schaefer_400_pscalar_kong,medial_mask);
