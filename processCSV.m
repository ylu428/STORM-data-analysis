function [file_list, folder_name] = processCSV(csvFilePath, imageDir)
    % Open the CSV file for reading
    fileID = fopen(csvFilePath, 'r');

    if fileID == -1
        error('Failed to open file: %s', csvFilePath);
    end

    % Initialize a map to store test IDs and their corresponding particle
    % numbers
    fileMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

    % Read each line from the CSV
    line = fgetl(fileID);
    while ischar(line)
        parts = strsplit(line, ':');
        if numel(parts) < 2
            line = fgetl(fileID);
            continue;
        end
        testID = strtrim(parts{1});
        nums = strtrim(parts{2});
        nums = strsplit(nums, ' ');
        nums = nums(~cellfun('isempty', nums));
        nums = str2double(nums);

        fileMap(testID) = nums;

        line = fgetl(fileID);
    end

    fclose(fileID);

    % Prompt the user to select the folder containing TIFF and CSV files
    folder_name = uigetdir(imageDir, 'Select Folder Containing TIFF and CSV Files');
    if folder_name == 0
        error('No folder selected');
    end

    % Prepare to match files
    tif_files = dir(fullfile(folder_name, '*.tif'));
    csv_files = dir(fullfile(folder_name, '*.csv'));

    file_list = {}; % Initialize an empty cell array to hold file sets

    % Loop through each test ID to find matching TIFF and CSV files
    fileKeys = keys(fileMap);
    for i = 1:length(fileKeys)
        testID = fileKeys{i};
        % Patterns to match file names
        tif_pattern = strcat('*', testID, '*_NP_f1.tif');
        csv_pattern = strcat('*', testID, '*_SMLs.csv');

        % Find files matching the patterns
        matched_tif = dir(fullfile(folder_name, tif_pattern));
        matched_csv = dir(fullfile(folder_name, csv_pattern));

        % If matches are found for both TIFF and CSV
        if ~isempty(matched_tif) && ~isempty(matched_csv)
            tif_full_path = fullfile(folder_name, matched_tif(1).name); % Assuming first match is desired one
            csv_full_path = fullfile(folder_name, matched_csv(1).name);

            % Append this file set to the list
            file_set = {tif_full_path, csv_full_path, fileMap(testID)};
            file_list{end+1} = file_set;
        end
    end

    % Optionally, display selected files and their details
    for i = 1:length(file_list)
        fprintf('TIFF File: %s, CSV File: %s, Test ID: %s, Selected NPs: %s\n', ...
                file_list{i}{1}, file_list{i}{2}, testID, mat2str(file_list{i}{3}));
    end
end