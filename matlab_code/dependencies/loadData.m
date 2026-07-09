function subjData = loadData(numSubjects, baseDir, alignType,taskType,dataType,curr_parc,scanIds,numTRs)
% Function to load data for either 'rest' or 'movie' condition

for scanId = scanIds
    if strcmp(taskType,'rest')
        numTR = 1800;
    else
        numTR = numTRs(scanId);
    end

    if strcmp(dataType,'gray')
        subjData = zeros(59412,numTR,numSubjects);
    else
        subjData = zeros(400,numTR,numSubjects);
        if strcmp(alignType,'anatomical')
            subjData = zeros(str2double(curr_parc),numTR,numSubjects);
        end
    end
    parfor subj = 1:numSubjects
        fprintf('Processing subject %d, scan %d\n', subj, scanId);
        % For parcellated data
        if strcmp(dataType,'parc')
            if strcmp(alignType,'anatomical')
                dataPath = sprintf('%s%s/inputs/%s/%s_%s_scan_%d_subj_%d_parc_%s.ptseries.nii', baseDir, alignType,dataType,alignType,taskType,scanId, subj,curr_parc);
                data = squeeze(niftiread(dataPath))';
                subjData(:,:,subj) = data;
            else
                dataPath = sprintf('%s%s/inputs/%s/%s_parc_%s_scan_%d_subj_%d.ptseries.nii', baseDir, alignType,dataType,alignType,curr_parc,scanId, subj);
                data = squeeze(niftiread(dataPath))';
                subjData(:,:,subj) = data;
            end

        else
            if strcmp(alignType,'anatomical')
                dataPath = sprintf('%s%s/inputs/%s/%s_%s_scan_%d_subj_%d.npy', baseDir, alignType,dataType,alignType, taskType, scanId, subj);
            else
                dataPath = sprintf('%s%s/inputs/%s/%s_parc_%s_scan_%d_subj_%d.npy', baseDir, alignType,dataType,alignType,curr_parc,scanId, subj);
            end
            subjData(:,:,subj)  = squeeze(readNPY(dataPath))';
        end

    end
end
end



