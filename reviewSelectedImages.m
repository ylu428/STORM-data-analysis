function selectedImages = reviewSelectedImages(imageDir, selectedImages)
    % Display all selected images in a single figure with an "X" to remove them,
    % a "Skip" button to ignore the review process, and a "Done" button to finalize the review.
    
    screenSize = get(0, 'ScreenSize'); % Get the screen size
    figWidth = screenSize(3) * 0.7; % 70% of the screen width
    figHeight = screenSize(4) * 0.7; % 70% of the screen height
    figX = (screenSize(3) - figWidth) / 2; % Center horizontally
    figY = (screenSize(4) - figHeight) / 2; % Center vertically

    % Create a figure window
    fig = figure('Name', 'Review Selected Images', 'NumberTitle', 'off', 'Toolbar', 'none', 'MenuBar', 'none', ...
                 'Position', [figX, figY, figWidth, figHeight]);

    numImages = length(selectedImages);
    sqrtNum = ceil(sqrt(numImages)); % Calculate grid size for a somewhat square layout

    % Keep track of which images to remove
    removeIndices = false(1, numImages);
    
    for i = 1:numImages
        % Load image
        img = imread(fullfile(imageDir, selectedImages{i}));
        
        % Display image in subplot
        ax = subplot(sqrtNum, sqrtNum, i);
        imshow(img, 'Parent', ax);
        title(ax, sprintf('Image %d', i));
        
        % Calculate "X" button size relative to the subplot size
        btnSize = min(ax.Position(3:4)) * 0.15;
        
        % Calculate "X" button position relative to the subplot
        btnPos = [ax.Position(1) + ax.Position(3) - btnSize, ax.Position(2) + ax.Position(4) - btnSize, btnSize, btnSize];
        
        % Add "X" button or clickable text
        btn = uicontrol('Style', 'pushbutton', 'String', 'X', 'Units', 'normalized', ...
                        'Position', btnPos, 'Callback', {@removeImageCallback, i});
    end
    
    % Add a "Skip" button to bypass the review process
    skipBtn = uicontrol('Style', 'pushbutton', 'String', 'Skip Review', 'Units', 'normalized', ...
                        'Position', [0.35, 0.01, 0.1, 0.05], 'Callback', @skipReviewCallback);
    
    % Add a "Done" button to finalize the review process
    doneBtn = uicontrol('Style', 'pushbutton', 'String', 'Done', 'Units', 'normalized', ...
                        'Position', [0.55, 0.01, 0.1, 0.05], 'Callback', @doneReviewCallback);
    
    % Wait for the user to finish with the figure before continuing
    uiwait(fig);
    
    % Remove selected images after review, if not skipped or removed
    selectedImages = selectedImages(~removeIndices);
    
    % Callback function to remove an image
    function removeImageCallback(~, ~, index)
        removeIndices(index) = true; % Mark image for removal
        subplot(sqrtNum, sqrtNum, index); % Select the subplot
        cla reset; % Clear the subplot
        title(sprintf('Image %d Removed', index)); % Indicate removal
        % Check if all images are marked for removal to close the figure
        if all(removeIndices)
            uiresume(fig); % Resume execution if figure is still open
        end
    end
    
    % Callback function for the "Skip" button
    function skipReviewCallback(~, ~)
        uiresume(fig); % Resume execution, closing the figure and skipping the review
        close(fig);
    end

    % Callback function for the "Done" button
    function doneReviewCallback(~, ~)
        uiresume(fig); % Resume execution, indicating review completion
        close(fig);
    end
end
