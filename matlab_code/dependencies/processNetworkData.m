function [netIds, parcelColors] = processNetworkData(filePath, numParcels, networkNames)
    % Function to process network data from a given file

    fileID       = fopen(filePath, 'r');
    netIds       = zeros(numParcels, 1);
    parcelColors = zeros(numParcels, 4);

    for i = 1:numParcels
        % Read the network name line
        networkLine = fgetl(fileID);
        disp(networkLine)
        colorLine = fgetl(fileID);

        % Extract network name
        underscoreIndices = strfind(networkLine, '_');
        networkName       = networkLine(underscoreIndices(2)+1:underscoreIndices(3)-1);

        % Find the index of the network
        networkIndex = find(strcmp(networkNames, networkName));
        netIds(i)    = networkIndex;

        % Extract color values
        colorValues        = sscanf(colorLine, '%d %d %d %d');
        parcelColors(i, :) = colorValues(2:end)';
    end

    fclose(fileID);
end
