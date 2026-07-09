function [isc_herit, isc_herit_perm,isc_herit_jack] = runISCHeritability(curr_parc, base_dir, kinship, Subjects,alignType,dataType,isc_mz_id,isc_dz_id,isc_nz_id,scan_id,covari_motion,num_perm,fam_ids)

% --- USER CONFIGURATION (edit these for your own system) ---
desktop_dir = '/path/to/desktop';  % originally /Users/davidgruskin/Desktop (Schaefer atlas templates)
data_dir    = '/path/to/data';     % originally /data/data7
% -----------------------------------------------------------

medial_mask = csvread([data_dir '/isc_heritability/medial_mask.csv']);

taskType = 'movie';
if strcmp(dataType,'parc')
    suffix = '.pscalar.nii';
    num_vox_cortex = 400;
    if strcmp(alignType,'anatomical')
        num_vox_cortex = str2double(curr_parc);
    end
else
    suffix = '.dscalar.nii';
    num_vox_cortex = 59412;
end

curr_scalar_template = strcat(desktop_dir, '/Schaefer2018_',curr_parc,'Parcels_17Networks_order',suffix);
if strcmp(curr_parc,'400') && strcmp(alignType,'anatomical')
    curr_scalar_template = strcat(desktop_dir, '/Schaefer2018_',curr_parc,'Parcels_Kong2022_17Networks_order',suffix);
end

wb_path = '/path/to/workbench/bin/wb_command';
numSubj = length(Subjects);
num_tr_list = [1432,1409];

numTRs = num_tr_list;

% Load and process data
iscDataPath = sprintf('%s%s/outputs/%s/%s_isc_parc_%s_scan_%s.mat', base_dir,alignType,dataType,alignType,curr_parc, num2str(scan_id));


if exist(iscDataPath,'file') ~=2

    subj_data = loadData(length(Subjects), base_dir, alignType,taskType,dataType,curr_parc,scan_id,numTRs);

    subj_data = permute(subj_data,[2,1,3]);
    % Calculate ISC
    isc = zeros(numSubj, numSubj, num_vox_cortex);
    for vox = 1:num_vox_cortex
        isc(:,:,vox) = corr(squeeze(subj_data(:,vox,:)));
    end

    % Save ISC data
    save(iscDataPath, 'isc', '-v7.3');
else
    isc = load(iscDataPath);
    name = fieldnames(isc);
    isc = isc.(name{1});
end

% Calculate heritability
isc_herit = zeros(num_vox_cortex, 1);
isc_herit_perm = zeros(num_vox_cortex, 1);
parfor grayordinate = 1:num_vox_cortex
    [isc_herit(grayordinate),isc_herit_perm(grayordinate),~,isc_herit_jack(grayordinate,:)] = h2_mat(squeeze(isc(:,:,grayordinate)), kinship(Subjects,Subjects), covari_motion(Subjects,[1,2,scan_id+2]), num_perm,fam_ids);
end

isc_herit(isc_herit==1) = NaN;
isc(isc==1) = NaN;
isc_mean = conv_z2r(nanmean(squeeze(nanmean(conv_r2z(isc),1))))';

% Save heritability and ISC data


heritDataPath = sprintf('%s%s/outputs/%s/%s_isc_herit_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(isc_herit, heritDataPath, wb_path, curr_scalar_template,medial_mask);

iscDataPath = sprintf('%s%s/outputs/%s/%s_isc_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(isc_mean, iscDataPath, wb_path, curr_scalar_template,medial_mask);


% Now separate out for MZ, DZ, and NZ
for vertex = 1:size(isc,3)
    isc_temp = squeeze(isc(:,:,vertex));
    isc_mz(:,vertex) = isc_temp(isc_mz_id);
    isc_dz(:,vertex) = isc_temp(isc_dz_id);
    isc_nz(:,vertex) = isc_temp(isc_nz_id);
end

% Now create ISC comparison maps
isc_mz(isc_mz==1) = NaN;
isc_dz(isc_dz==1) = NaN;
isc_nz(isc_nz==1) = NaN;

isc_mz_avg = squeeze(nanmean(conv_r2z(isc_mz),1));
isc_dz_avg = squeeze(nanmean(conv_r2z(isc_dz),1));
isc_nz_avg = squeeze(nanmean(conv_r2z(isc_nz),1));

isc_mzdz = conv_z2r(isc_mz_avg - isc_dz_avg);
isc_mznz = conv_z2r(isc_mz_avg - isc_nz_avg);
isc_dznz = conv_z2r(isc_dz_avg - isc_nz_avg);

iscMZDZPath = sprintf('%s%s/outputs/%s/%s_isc_mzdz_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(isc_mzdz', iscMZDZPath, wb_path, curr_scalar_template,medial_mask);

iscMZNZPath = sprintf('%s%s/outputs/%s/%s_isc_mznz_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(isc_mznz', iscMZNZPath, wb_path, curr_scalar_template,medial_mask);

iscDZNZPath = sprintf('%s%s/outputs/%s/%s_isc_dznz_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(isc_dznz', iscDZNZPath, wb_path, curr_scalar_template,medial_mask);

iscMZPath = sprintf('%s%s/outputs/%s/%s_isc_mz_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(conv_z2r(isc_mz_avg'), iscMZPath, wb_path, curr_scalar_template,medial_mask);

iscNZPath = sprintf('%s%s/outputs/%s/%s_isc_nz_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(conv_z2r(isc_nz_avg'), iscNZPath, wb_path, curr_scalar_template,medial_mask);

iscDZPath = sprintf('%s%s/outputs/%s/%s_isc_dz_%s_%s_scan_%s%s', base_dir,alignType,dataType,alignType,dataType,curr_parc, num2str(scan_id),suffix);
saveCifti(conv_z2r(isc_dz_avg'), iscDZPath, wb_path, curr_scalar_template,medial_mask);

end

