function [folder_name, file_list] = processTextAndFiles()
    % Use UI to get the text file
    [txtFileName, txtPath] = uigetfile('*.txt', 'Select the text file');
    if isequal(txtFileName, 0)
        error('No file selected');
    end
    txtFileFullPath = fullfile(txtPath, txtFileName);

    % Read and parse the text file
    [fileID, errmsg] = fopen(txtFileFullPath, 'r');
    if fileID < 0
        error(errmsg);
    end
    data = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    lines = data{1};

    % Process each line
    file_map = containers.Map;
    for i = 1:length(lines)
        line = lines{i};
        % Splitting line into test name and numbers
        parts = strsplit(line, ':');
        if numel(parts) < 2
            continue; % Skip lines that don't have a ':' character
        end
        test_name = strtrim(parts{1}); % Extract test name
        nums = strtrim(parts{2}); % Extract number string
        nums = strsplit(nums, ' '); % Split by space
        nums = nums(~cellfun('isempty', nums)); % Remove any empty strings

        if length(nums) > 1 % Include line if more than one particle number
            file_map(test_name) = nums;
        end
    end

    % Select folder
    folder_name = uigetdir(txtPath, "Select the folder with _NP_f1.tif files");
    if folder_name == 0
        error('No folder selected');
    end

    % Filter TIFF and CSV files
    tif_files = dir(fullfile(folder_name, '*_NP_f1.tif'));
    csv_files = dir(fullfile(folder_name, '*_SMLs.csv'));

    % Create output list
    file_list = {};
    keys = file_map.keys;
    for i = 1:length(keys)
        key = keys{i};
        tif_name = strcat(key, '_NP_f1.tif');
        csv_name = strcat(key, '_SMLs.csv');

        % Check if both TIFF and CSV files exist
        % if any(strcmp({tif_files.name}, tif_name)) && any(strcmp({csv_files.name}, csv_name))
            file_list{end+1} = {csv_name, tif_name, file_map(key)};
        % end
    end

    % Display the file list in a structured format
    disp('Selected files and PNGs:');
    for i = 1:length(file_list)
        disp(['Test Name: ', file_list{i}{1}]);
        disp(['CSV File: ', file_list{i}{1}]);
        disp(['TIFF File: ', file_list{i}{2}]);
        selectedParticles = strjoin(file_list{i}{3}, ', '); % Remove leading comma
        disp(['Selected Particles: ', selectedParticles]);
        disp('-----------------------------------');
    end
end
