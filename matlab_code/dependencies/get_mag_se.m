function global_se = get_mag_se(directory,task_type, align_type,curr_parc,scan_id)
% Get a list of files matching the pattern
filePattern = strcat('fc_herit_net_',task_type,'_',align_type,'_parc_',curr_parc,'_scan_',num2str(scan_id),'_jackknife_*.csv');
files       = dir(fullfile(directory, filePattern));

% Count the number of files
numFiles = length(files);

% Check if there are any files
if numFiles == 0
    error('No files found matching the specified pattern.');
end

% Initialize a variable to store the number of rows (size of 'Var' column)
numRows = 0;

% Check the size of the 'Var' column in the first file
firstFile = readtable(fullfile(directory, files(1).name));
numRows   = height(firstFile);

% Preallocate an array to store the 'Var' column data from each file
varArray = NaN(numRows, numFiles);  % Using NaN for uninitialized values

% Loop through the files and load the 'Var' column
for k = 1:numFiles
    % Full path to the file
    filePath = fullfile(directory, files(k).name);

    % Load the file as a table
    dataTable = readtable(filePath);

    % Check if the current file has the same number of rows
    if height(dataTable) ~= numRows
        warning('Number of rows in file %s does not match the first file. Skipping this file.', files(k).name);
        continue;
    end

    % Extract the 'Var' column and store it in the array
    varArray(:, k) = dataTable.Var;
end
total_table = strcat('fc_herit_net_',task_type,'_',align_type,'_parc_',curr_parc,'_scan_',num2str(scan_id),'.csv');
filePath    = fullfile(directory, total_table);

varArray(:, all(isnan(varArray), 1)) = [];

% Load the file as a table
dataTable       = readtable(filePath);
[~,global_se,~] = compareHeritability(nanmean(dataTable.Var),nanmean(varArray),nanmean(dataTable.Var),nanmean(varArray));
end
