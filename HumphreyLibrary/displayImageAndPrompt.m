function userInput = displayImageAndPrompt(imageDir, filename, currentImage, totalImages)
    % Define persistent variable to store the user input
    persistent selectionMade;
    selectionMade = false; % Reset for each call
    
    % Initialize user input as NaN (to distinguish between Yes/No and no input)
    userInput = NaN;

    % Get the size of the screen
    screenSize = get(0, 'ScreenSize');
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);
    
    % Calculate the size of the figure to be 2/3 of the screen width and height
    figWidth = screenWidth * 2/3;
    figHeight = screenHeight * 2/3;
    
    % Calculate the position to center the figure on the screen
    figX = (screenWidth - figWidth) / 2;
    figY = (screenHeight - figHeight) / 2;
    
    % Create a figure window with the calculated size and position
    fig = figure('Position', [figX, figY, figWidth, figHeight]);
    
    % Display the image
    img = imread(fullfile(imageDir, filename));
    imshow(img);
    title(sprintf('Do you want to save this image? (%d/%d)', currentImage, totalImages));
    
    % Display the image name in the bottom left corner of the figure
    uicontrol('Style', 'text', 'String', filename,...
        'Position', [10, 10, 300, 20], 'BackgroundColor', 'white', 'FontSize', 13);
    
    % Calculate button size and position based on the figure size
    buttonWidth = 100; % Adjusted button width for fitting
    buttonHeight = 40; % Adjusted button height
    buttonYPosition = 40; % Distance from the bottom of the figure
    buttonSpacing = 10; % Space between buttons
    
    % Position for Yes button
    btnYesXPosition = 10; % Start from the left
    
    % Position for No button
    btnNoXPosition = btnYesXPosition + buttonWidth + buttonSpacing;
    
    % Position for Done button, next to the No button
    btnDoneXPosition = btnNoXPosition + buttonWidth + buttonSpacing;
    
    % Yes Button
    btnYes = uicontrol('Style', 'pushbutton', 'String', 'Yes',...
        'Position', [btnYesXPosition, buttonYPosition, buttonWidth, buttonHeight],...
        'Callback', {@btnCallback, 1});
    
    % No Button
    btnNo = uicontrol('Style', 'pushbutton', 'String', 'No',...
        'Position', [btnNoXPosition, buttonYPosition, buttonWidth, buttonHeight],...
        'Callback', {@btnCallback, 0});
    
    % Done Button
    btnDone = uicontrol('Style', 'pushbutton', 'String', 'Done',...
        'Position', [btnDoneXPosition, buttonYPosition, buttonWidth, buttonHeight],...
        'Callback', {@btnCallback, -1});
    
    % Wait for the user to make a choice or close the figure
    uiwait(fig);
    
    % Check if a selection was made
    if selectionMade
        % Selection was made; proceed normally
    else
        if ishandle(fig)
            close(fig); % Ensure the figure is closed if still open
        end
        error('Image processing was terminated by the user.');
    end

    % Nested function for button callbacks
    function btnCallback(src, event, choice)
        userInput = choice; % Set the user input based on the button pressed
        selectionMade = true; % Mark that a selection has been made
        uiresume(fig); % Resume execution
        if ishandle(fig)
            close(fig); % Close the figure window
        end
    end
end