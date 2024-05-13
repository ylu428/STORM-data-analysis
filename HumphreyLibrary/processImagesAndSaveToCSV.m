% This script processes images in a selected directory, allows for manual selection of relevant images, and saves the selected data to a CSV file. 
% It provides a UI for selecting the image directory and displays images with options to exclude them before processing.

function [csvFilePath, imageDir] = processImagesAndSaveToCSV()
    % Use UI to get the directory containing the images
    imageDir = uigetdir('*.*', 'Select the folder holding NP images');
    if imageDir == 0
        % User pressed cancel or closed the dialog
        disp('No folder selected. Exiting...');
        return;
    end
    
    % Display and select images before processing
    allImageFiles = dir(fullfile(imageDir, '*_NP_SML.png')); 
    selectedImages = selectImagesForProcessing(imageDir, {allImageFiles.name});
    
    % Set the CSV file path directly without prompting the user
    csvFilePath = fullfile(imageDir, 'Selected_NP.csv');
    
    % Initialize a container for selections
    selections = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Load and filter selected images
    totalImages = length(selectedImages); % Total number of selected images
    
    % Process each selected image
    try
        for i = 1:totalImages
            filename = selectedImages{i};
            
            % Filter and parse filename
            [isRelevant, testID, particleNumber] = parseAndFilterFilename(filename);
            
            if isRelevant
                % Display image and prompt for user input, include progress in title
                userInput = displayImageAndPrompt(imageDir, filename, i, totalImages);
                
                % If user selects 'Yes', update the selections container and store filename for review
                if userInput == 1
                    if selections.isKey(testID)
                        selections(testID) = [selections(testID), particleNumber];
                    else
                        selections(testID) = [particleNumber];
                    end
                elseif userInput == -1
                    % If 'Done' button is clicked, break the loop
                    break; % Exit the loop early
                end
            end
        end
    catch ME
        if strcmp(ME.identifier, 'ImageProcessing:UserTerminated')
            disp('Image processing terminated early by the user.');
        else
            rethrow(ME); % Handle unexpected errors
        end
    end
    
    % Convert the list of selected and reviewed images back into selections map
    selections = convertSelectionsToMap(selectedImages);
    
    % Update the CSV with the final selections
    updateCSV(csvFilePath, selections);
    return;
end

function selectedFiles = selectImagesForProcessing(imageDir, imageFiles)
    % Creates a figure to display images with an "X" to remove them, excluding certain patterns
    screenSize = get(0, 'ScreenSize');
    figWidth = screenSize(3) * 0.7; % 70% of the screen width
    figHeight = screenSize(4) * 0.7; % 70% of the screen height
    figX = (screenSize(3) - figWidth) / 2; % Center horizontally
    figY = (screenSize(4) - figHeight) / 2; % Center vertically

    fig = figure('Name', 'Select Images for Processing', 'NumberTitle', 'off', 'Toolbar', 'none', 'MenuBar', 'none', ...
                 'Position', [figX, figY, figWidth, figHeight]);

    % Filter imageFiles based on exclusion criteria
    filteredImageFiles = filterImages(imageFiles);
    numImages = length(filteredImageFiles);
    sqrtNum = ceil(sqrt(numImages)); % Calculate grid size for a somewhat square layout
    removeIndices = false(1, numImages);

    for i = 1:numImages
        img = imread(fullfile(imageDir, filteredImageFiles{i}));
        ax = subplot(sqrtNum, sqrtNum, i);
        imshow(img, 'Parent', ax);
        title(ax, sprintf('Image %d', i));

        % Calculate "X" button size relative to the subplot size
        btnSize = min(ax.Position(3:4)) * 0.15;

        % Calculate "X" button position relative to the subplot
        btnPos = [ax.Position(1) + ax.Position(3) - btnSize, ax.Position(2) + ax.Position(4) - btnSize, btnSize, btnSize];

        % Add "X" button or clickable text
        uicontrol('Style', 'pushbutton', 'String', 'X', 'Units', 'normalized', ...
                  'Position', btnPos, 'Callback', {@removeImageCallback, i});
    end

    % Add a "Done" button to finalize the selections
    uicontrol('Style', 'pushbutton', 'String', 'Done', 'Units', 'normalized', ...
              'Position', [0.5 - 0.05, 0.01, 0.1, 0.05], 'Callback', @doneSelectionCallback);

    % Wait for the user to finish with the figure before continuing
    uiwait(fig);
    
    selectedFiles = filteredImageFiles(~removeIndices);
    close(fig);

    function filteredFiles = filterImages(imageFiles)
        % Exclude files ending in 'hist_ori' or matching specific unwanted patterns
        unwantedPatterns = {'hist_ori', 'hist_corr'};
        filteredFiles = {};
        for idx = 1:length(imageFiles)
            skipFile = any(cellfun(@(pattern) endsWith(imageFiles{idx}, pattern), unwantedPatterns));
            if ~skipFile
                filteredFiles{end+1} = imageFiles{idx};
            end
        end
    end

    % Callback function to remove an image
    function removeImageCallback(~, ~, index)
        removeIndices(index) = true; % Mark image for removal
        subplot(sqrtNum, sqrtNum, index); % Select the subplot
        cla reset; % Clear the subplot
        title(sprintf('Image %d Removed', index)); % Indicate removal
    end

    % Callback function for the "Done" button
    function doneSelectionCallback(~, ~)
        uiresume(fig); % Resume execution, indicating selection completion
    end
end
