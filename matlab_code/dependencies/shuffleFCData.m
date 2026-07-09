function [shuffled_rest_fc, shuffled_movie_fc, shuffled_rest_covari, shuffled_movie_covari, swap_indices] = shuffleFCData(rest_fc_square, movie_fc_square, rest_covari, movie_covari,iter)
rng(iter)
% Check if FC data is provided
if exist('rest_fc_square', 'var') && exist('movie_fc_square', 'var')
    [nROIs, nROIs, num_subjects, num_days] = size(rest_fc_square);
    shuffled_rest_fc                       = zeros(size(rest_fc_square));
    shuffled_movie_fc                      = zeros(size(movie_fc_square));
else
    nROIs        = [];
    num_subjects = [];
    num_days     = [];
end

% Check if covariate data is provided
if exist('rest_covari', 'var') && exist('movie_covari', 'var')
    if isempty(num_subjects)
        num_subjects = size(rest_covari, 1);
    end
    shuffled_rest_covari  = zeros(size(rest_covari));
    shuffled_movie_covari = zeros(size(movie_covari));
else
    shuffled_rest_covari  = [];
    shuffled_movie_covari = [];
end

% Generate random swap indices (1: keep original, 2: swap)
swap_indices = randi(2, num_subjects, 1);

% Loop through subjects (if FC data provided)
if ~isempty(nROIs) && ~isempty(num_days)
    for subject_id = 1:num_subjects
        for day_id = 1:num_days
            if swap_indices(subject_id) == 1
                shuffled_rest_fc(:,:,subject_id,day_id)  = rest_fc_square(:,:,subject_id,day_id);
                shuffled_movie_fc(:,:,subject_id,day_id) = movie_fc_square(:,:,subject_id,day_id);
            else
                shuffled_rest_fc(:,:,subject_id,day_id)  = movie_fc_square(:,:,subject_id,day_id);
                shuffled_movie_fc(:,:,subject_id,day_id) = rest_fc_square(:,:,subject_id,day_id);
            end
        end
    end
end

% Loop through subjects (if covariate data provided)
if ~isempty(num_subjects)
    for subject_id = 1:num_subjects
        if swap_indices(subject_id) == 1
            shuffled_rest_covari(subject_id,:)  = rest_covari(subject_id,:);
            shuffled_movie_covari(subject_id,:) = movie_covari(subject_id,:);
        else
            shuffled_rest_covari(subject_id,:)  = movie_covari(subject_id,:);
            shuffled_movie_covari(subject_id,:) = rest_covari(subject_id,:);
        end
    end
end
end
