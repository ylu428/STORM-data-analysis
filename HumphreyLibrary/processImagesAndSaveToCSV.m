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
    
