function int = calculate_int(subj_data)

% Get the size of subj_data
[trs, vertices] = size(subj_data);

% Initialize INT array
int = zeros(vertices, 1);

% Calculate INT for each vertex and each subject

parfor v = 1:vertices
    % Extract the BOLD time series for the kth vertex
    bold_ts = squeeze(subj_data(:, v));

    acf    = autocorr(bold_ts,trs-1);
    neg_id = find(acf<0,1)-1;

    % Find the lag N where ACF first becomes negative
    if isempty(neg_id)
        neg_id = trs - 1; % Use max lag if ACF does not go negative
    end

    % Calculate INT
    int(v, 1) = trapz(acf(1:neg_id));
end
end
