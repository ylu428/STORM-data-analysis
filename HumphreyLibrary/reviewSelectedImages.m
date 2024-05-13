% This function reviews the selected images with options to remove, skip, or finalize the review process.
% This script creates a user interface displaying all selected images in a grid layout.
% Users can click an "X" button to remove specific images, a "Skip" button to go straight to the review section,
% Selecting the "Done" button to finalize and save the reviewed selection.

function selectedImages = reviewSelectedImages(imageDir, selectedImages)
    % Display all selected images in a single figure with an "X" to remove them,
    % a "Skip" button to ignore the review process, and a "Done" button to finalize the review.
    
    % Get the screen size
    screenSize = get(0, 'ScreenSize');
    
    % Define the figure window size as 70% of the screen size
    figWidth = screenSize(3) * 0.7;
    figHeight = screenSize(4) * 0.7;
    
    % Center the figure window on the screen
    figX = (screenSize(3) - figWidth) / 2;
    figY = (screenSize(4) - figHeight) / 2;

    % Create a figure window for reviewing selected images
    fig = figure('Name', 'Review Selected Images', 'NumberTitle', 'off', 'Toolbar', 'none', 'MenuBar', 'none', ...
                 'Position', [figX, figY, figWidth, figHeight]);

    % Calculate the number of images and determine grid layout
    numImages = length(selectedImages);
    sqrtNum = ceil(sqrt(numImages)); % Grid size for a roughly square layout

    % Initialize an array to keep track of images to be removed
    removeIndices = false(1, numImages);
    
    for i = 1:numImages
        % Load each image from the specified directory
        img = imread(fullfile(imageDir, selectedImages{i}));
        
        % Display the image in a subplot
        ax = subplot(sqrtNum, sqrtNum, i);
        imshow(img, 'Parent', ax);
        title(ax, sprintf('Image %d', i));
        
        % Calculate the "X" button size and position relative to the subplot
        btnSize = min(ax.Position(3:4)) * 0.15;
        btnPos = [ax.Position(1) + ax.Position(3) - btnSize, ax.Position(2) + ax.Position(4) - btnSize, btnSize, btnSize];
        
        % Add the "X" button to remove the image
        btn = uicontrol('Style', 'pushbutton', 'String', 'X', 'Units', 'normalized', ...
                        'Position', btnPos, 'Callback', {@removeImageCallback, i});
    end
    
    % Add a "Skip" button to bypass the review process
    skipBtn = uicontrol('Style', 'pushbutton', 'String', 'Skip Review', 'Units', 'normalized', ...
                        'Position', [0.35, 0.01, 0.1, 0.05], 'Callback', @skipReviewCallback);
    
    % Add a "Done" button to finalize the review process
    doneBtn = uicontrol('Style', 'pushbutton', 'String', 'Done', 'Units', 'normalized', ...
                        'Position', [0.55, 0.01, 0.1, 0.05], 'Callback', @doneReviewCallback);
    
    % Wait for the user to complete the review before proceeding
    uiwait(fig);
    
    % Update the list of selected images, excluding those marked for removal
    selectedImages = selectedImages(~removeIndices);
    
    % Callback function to remove an image
    function removeImageCallback(~, ~, index)
        removeIndices(index) = true; % Mark image for removal
        subplot(sqrtNum, sqrtNum, index); % Select the subplot
        cla reset; % Clear the subplot
        title(sprintf('Image %d Removed', index)); % Indicate removal
        % Check if all images are marked for removal to close the figure
        if all(removeIndices)
            uiresume(fig); % Resume execution if all images are removed
        end
    end
    
    % Callback function for the "Skip" button
    function skipReviewCallback(~, ~)
        uiresume(fig); % Resume execution, skipping the review
        close(fig); % Close the figure window
    end

    % Callback function for the "Done" button
    function doneReviewCallback(~, ~)
        uiresume(fig); % Resume execution, indicating the review is done
        close(fig); % Close the figure window
    end
end
