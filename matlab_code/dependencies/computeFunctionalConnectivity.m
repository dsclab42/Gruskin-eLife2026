function [fc, fcSquare] = computeFunctionalConnectivity(subjData, numSubjects)

fc = zeros((400*399)/2, numSubjects, 2);
fcSquare = zeros(400, 400, numSubjects, 2);
hvec = @(x) x(triu(true(size(x)), 1));

for scanId = 1:2
    for subj = 1:numSubjects
        fc(:, subj, scanId) = hvec(corr(subjData{scanId}(:,:,subj)', 'rows', 'complete'));
        fcSquare(:,:, subj, scanId) = corr(subjData{scanId}(:,:,subj)', 'rows', 'complete');
    end
end
end
