function [csvFilePath, imageDir] = processImagesAndSaveToCSV()

    % Use UI to get the directory containing the images
    imageDir = uigetdir('*.*', 'Select the folder holding NP images');
    if imageDir == 0
        % User pressed cancel or closed the dialog
        disp('No folder selected. Exiting...');
        return;
    end
    
    % Set the CSV file path directly without prompting the user
    csvFilePath = fullfile(imageDir, 'Selected_NP.csv');
    
    % Initialize a container for selections
    selections = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Initialize a list to store filenames of selected images
    selectedImages = {};
    
    % Load and filter images
    images = dir(fullfile(imageDir, '*_NP_SML.png'));
    totalImages = length(images); % Total number of images
    
    % Process each image
    try
        for i = 1:totalImages
            filename = images(i).name;
            
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
                    selectedImages{end+1} = filename; % Store filename for review
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
    
    % After processing, call reviewSelectedImages to allow user to review and possibly remove selections
    selectedImages = reviewSelectedImages(imageDir, selectedImages);
    % Convert the list of selected and reviewed images back into selections map
    selections = convertSelectionsToMap(selectedImages);
    
    % Update the CSV with the final selections
    updateCSV(csvFilePath, selections);
    return;
end