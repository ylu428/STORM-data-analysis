addpath(genpath('M:\Yi-Han\ONI_Storm\Yi-Han_edited_code\MATLAB_library'))
addpath(genpath('M:\Yi-Han\ONI_Storm\Yi-Han_edited_code\MATLAB_library\HumphreyLibrary'))


FWHMtest= startAutomation();

temp = 1;
for i = 1:length(FWHMtest)
    if FWHMtest(i).pass_FWHMthr_ratio_corr~=1
        % disp(FWHMtest(i).SML_file)
        % disp(strcat('    pass_FWHMthr_ratio_corr =', num2str(FWHMtest(i).pass_FWHMthr_ratio_corr)))
        if length(str2num(FWHMtest(i).pass_FWHMthr_par_num))>=2
            test_result = process1(FWHMtest(i).folder, FWHMtest(i).NP1_fif_file, FWHMtest(i).SML_file, str2num(FWHMtest(i).pass_FWHMthr_par_num));
            FWHMtest2(temp)= test_result;
            temp=temp+1;
            clear test_result
        elseif length(str2num(FWHMtest(i).pass_FWHMthr_par_num))<2
            temp = split(FWHMtest(i).SML_file, '.csv');
            driftcorr_matrix = strcat(FWHMtest(i).folder, '\output\', char(temp(1)), '_driftcorr.mat');
            if exist(driftcorr_matrix,"file") ==2
                delete(driftcorr_matrix)
            end  
        end
    end
end

function FWHMtest = startAutomation()
    
    [csvFilePath, imageDir] = processImagesAndSaveToCSV();
    
    % Check if csvFilePath is not empty (meaning the function didn't early return)
    if isempty(csvFilePath)
        disp('Process was cancelled or an error occurred. Exiting automation...');
        return;
    end
    
    % Call processCSV to read the file and prepare file sets for processing
    [file_list, folder_name] = processCSV(csvFilePath, imageDir);
    
    % Loop through each file set and process
    % for i = 1:length(file_list)
    %     file_set = file_list{i};
    %     tif_file = file_set{1};
    %     csv_file = file_set{2};
    %     select_NP = file_set{3}; % 
        
    for i = 1:length(file_list)
        file_set = file_list{i};
        tif_file = file_set{1};  
        csv_file = file_set{2};
        % Convert cell array of strings to a numeric array
        % select_NP = cellfun(@str2num, file_set{3});
        select_NP = file_set{3};

        % Call the process function for each set
        % process(tif_file, csv_file, select_NP)
        test_result = process1(folder_name, tif_file, csv_file, select_NP);
        FWHMtest(i)= test_result;
        clear test_result
    end
end


% function csvFilePath = processImagesAndSaveToCSV()
% 
%     % Use UI to get the directory containing the images
%     imageDir = uigetdir('', 'Select the folder holding NP images');
%     if imageDir == 0
%         % User pressed cancel or closed the dialog
%         disp('No folder selected. Exiting...');
%         return;
%     end
% 
%     % Set the CSV file path directly without prompting the user
%     csvFilePath = fullfile(imageDir, 'Selected_NP.csv');
% 
%     % Initialize a container for selections
%     selections = containers.Map('KeyType', 'char', 'ValueType', 'any');
% 
%     % Initialize a list to store filenames of selected images
%     selectedImages = {};
% 
%     % Load and filter images
%     images = dir(fullfile(imageDir, '*_NP_SML.png'));
%     totalImages = length(images); % Total number of images
% 
%     % Process each image
%     try
%         for i = 1:totalImages
%             filename = images(i).name;
% 
%             % Filter and parse filename
%             [isRelevant, testID, particleNumber] = parseAndFilterFilename(filename);
% 
%             if isRelevant
%                 % Display image and prompt for user input, include progress in title
%                 userInput = displayImageAndPrompt(imageDir, filename, i, totalImages);
% 
%                 % If user selects 'Yes', update the selections container and store filename for review
%                 if userInput == 1
%                     if selections.isKey(testID)
%                         selections(testID) = [selections(testID), particleNumber];
%                     else
%                         selections(testID) = [particleNumber];
%                     end
%                     selectedImages{end+1} = filename; % Store filename for review
%                 elseif userInput == -1
%                     % If 'Done' button is clicked, break the loop
%                     break; % Exit the loop early
%                 end
%             end
%         end
%     catch ME
%         if strcmp(ME.identifier, 'ImageProcessing:UserTerminated')
%             disp('Image processing terminated early by the user.');
%         else
%             rethrow(ME); % Handle unexpected errors
%         end
%     end
% 
%     % After processing, call reviewSelectedImages to allow user to review and possibly remove selections
%     selectedImages = reviewSelectedImages(imageDir, selectedImages);
%     % Convert the list of selected and reviewed images back into selections map
%     selections = convertSelectionsToMap(selectedImages);
% 
%     % Update the CSV with the final selections
%     updateCSV(csvFilePath, selections);
%     return;
% end
% 
% 
% 
% function [isRelevant, testID, particleNumber] = parseAndFilterFilename(filename)
% 
%     isRelevant = false;  % Default value
%     testID = '';
%     particleNumber = 0;
% 
%     % Check if filename ends with '_NP_SML.png'
%     if endsWith(filename, '_NP_SML.png')
%         % Adjusted regex to capture test ID and particle number
%         % This regex accounts for a prefix, followed by the test ID (with or without a hyphen) and the particle number
%         tokens = regexp(filename, '(.*?)_SMLs?-part-(\d+)_NP_SML\.png', 'tokens');
%         if ~isempty(tokens)
%             token = tokens{1};  % Assuming one match, which is a cell containing the captured groups
%             testID = token{1};  % The test ID is the first capture group
%             particleNumber = str2double(token{2});  % The particle number is the second capture group
%             isRelevant = true;
%         else
%             % Handle the case where the filename does not match the expected pattern
%             isRelevant = false;
%         end
% 
%     end
% end
% 
% function userInput = displayImageAndPrompt(imageDir, filename, currentImage, totalImages)
%     % Define persistent variable to store the user input
%     persistent selectionMade;
%     selectionMade = false; % Reset for each call
% 
%     % Initialize user input as NaN (to distinguish between Yes/No and no input)
%     userInput = NaN;
% 
%     % Get the size of the screen
%     screenSize = get(0, 'ScreenSize');
%     screenWidth = screenSize(3);
%     screenHeight = screenSize(4);
% 
%     % Calculate the size of the figure to be 2/3 of the screen width and height
%     figWidth = screenWidth * 2/3;
%     figHeight = screenHeight * 2/3;
% 
%     % Calculate the position to center the figure on the screen
%     figX = (screenWidth - figWidth) / 2;
%     figY = (screenHeight - figHeight) / 2;
% 
%     % Create a figure window with the calculated size and position
%     fig = figure('Position', [figX, figY, figWidth, figHeight]);
% 
%     % Display the image
%     img = imread(fullfile(imageDir, filename));
%     imshow(img);
%     title(sprintf('Do you want to save this image? (%d/%d)', currentImage, totalImages));
% 
%     % Display the image name in the bottom left corner of the figure
%     uicontrol('Style', 'text', 'String', filename,...
%         'Position', [10, 10, 300, 20], 'BackgroundColor', 'white', 'FontSize', 13);
% 
%     % Calculate button size and position based on the figure size
%     buttonWidth = 100; % Adjusted button width for fitting
%     buttonHeight = 40; % Adjusted button height
%     buttonYPosition = 40; % Distance from the bottom of the figure
%     buttonSpacing = 10; % Space between buttons
% 
%     % Position for Yes button
%     btnYesXPosition = 10; % Start from the left
% 
%     % Position for No button
%     btnNoXPosition = btnYesXPosition + buttonWidth + buttonSpacing;
% 
%     % Position for Done button, next to the No button
%     btnDoneXPosition = btnNoXPosition + buttonWidth + buttonSpacing;
% 
%     % Yes Button
%     btnYes = uicontrol('Style', 'pushbutton', 'String', 'Yes',...
%         'Position', [btnYesXPosition, buttonYPosition, buttonWidth, buttonHeight],...
%         'Callback', {@btnCallback, 1});
% 
%     % No Button
%     btnNo = uicontrol('Style', 'pushbutton', 'String', 'No',...
%         'Position', [btnNoXPosition, buttonYPosition, buttonWidth, buttonHeight],...
%         'Callback', {@btnCallback, 0});
% 
%     % Done Button
%     btnDone = uicontrol('Style', 'pushbutton', 'String', 'Done',...
%         'Position', [btnDoneXPosition, buttonYPosition, buttonWidth, buttonHeight],...
%         'Callback', {@btnCallback, -1});
% 
%     % Wait for the user to make a choice or close the figure
%     uiwait(fig);
% 
%     % Check if a selection was made
%     if selectionMade
%         % Selection was made; proceed normally
%     else
%         if ishandle(fig)
%             close(fig); % Ensure the figure is closed if still open
%         end
%         error('Image processing was terminated by the user.');
%     end
% 
%     % Nested function for button callbacks
%     function btnCallback(src, event, choice)
%         userInput = choice; % Set the user input based on the button pressed
%         selectionMade = true; % Mark that a selection has been made
%         uiresume(fig); % Resume execution
%         if ishandle(fig)
%             close(fig); % Close the figure window
%         end
%     end
% end
% 
% function updateCSV(csvFilePath, selections)
%     fid = fopen(csvFilePath, 'w');
%     keys = selections.keys;
% 
%     for i = 1:length(keys)
%         testID = keys{i};
%         particleNumbers = selections(testID);
%         % Format the line as 'TestID: number number number...'
%         fprintf(fid, '%s: %s\n', testID, num2str(particleNumbers, '%d '));
%     end
% 
%     fclose(fid);
% end
% 
% function file_list = processCSV(csvFilePath)
%     % Open the CSV file for reading
%     fileID = fopen(csvFilePath, 'r');
% 
%     if fileID == -1
%         error('Failed to open file: %s', csvFilePath);
%     end
% 
%     % Initialize a map to store test IDs and their corresponding particle
%     % numbers
%     fileMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
% 
%     % Read each line from the CSV
%     line = fgetl(fileID);
%     while ischar(line)
%         parts = strsplit(line, ':');
%         if numel(parts) < 2
%             line = fgetl(fileID);
%             continue;
%         end
%         testID = strtrim(parts{1});
%         nums = strtrim(parts{2});
%         nums = strsplit(nums, ' ');
%         nums = nums(~cellfun('isempty', nums));
%         nums = str2double(nums);
% 
%         fileMap(testID) = nums;
% 
%         line = fgetl(fileID);
%     end
% 
%     fclose(fileID);
% 
%     % Prompt the user to select the folder containing TIFF and CSV files
%     folder_name = uigetdir('', 'Select Folder Containing TIFF and CSV Files');
%     if folder_name == 0
%         error('No folder selected');
%     end
% 
%     % Prepare to match files
%     tif_files = dir(fullfile(folder_name, '*.tif'));
%     csv_files = dir(fullfile(folder_name, '*.csv'));
% 
%     file_list = {}; % Initialize an empty cell array to hold file sets
% 
%     % Loop through each test ID to find matching TIFF and CSV files
%     fileKeys = keys(fileMap);
%     for i = 1:length(fileKeys)
%         testID = fileKeys{i};
%         % Patterns to match file names
%         tif_pattern = strcat('*', testID, '*_NP_f1.tif');
%         csv_pattern = strcat('*', testID, '*_SMLs.csv');
% 
%         % Find files matching the patterns
%         matched_tif = dir(fullfile(folder_name, tif_pattern));
%         matched_csv = dir(fullfile(folder_name, csv_pattern));
% 
%         % If matches are found for both TIFF and CSV
%         if ~isempty(matched_tif) && ~isempty(matched_csv)
%             tif_full_path = fullfile(folder_name, matched_tif(1).name); % Assuming first match is desired one
%             csv_full_path = fullfile(folder_name, matched_csv(1).name);
% 
%             % Append this file set to the list
%             file_set = {tif_full_path, csv_full_path, fileMap(testID)};
%             file_list{end+1} = file_set;
%         end
%     end
% 
%     % Optionally, display selected files and their details
%     for i = 1:length(file_list)
%         fprintf('TIFF File: %s, CSV File: %s, Test ID: %s, Selected NPs: %s\n', ...
%                 file_list{i}{1}, file_list{i}{2}, testID, mat2str(file_list{i}{3}));
%     end
% end
% 
% 
% 
% function selectedImages = reviewSelectedImages(imageDir, selectedImages)
%     % Display all selected images in a single figure with an "X" to remove them,
%     % a "Skip" button to ignore the review process, and a "Done" button to finalize the review.
% 
%     screenSize = get(0, 'ScreenSize'); % Get the screen size
%     figWidth = screenSize(3) * 0.7; % 70% of the screen width
%     figHeight = screenSize(4) * 0.7; % 70% of the screen height
%     figX = (screenSize(3) - figWidth) / 2; % Center horizontally
%     figY = (screenSize(4) - figHeight) / 2; % Center vertically
% 
%     % Create a figure window
%     fig = figure('Name', 'Review Selected Images', 'NumberTitle', 'off', 'Toolbar', 'none', 'MenuBar', 'none', ...
%                  'Position', [figX, figY, figWidth, figHeight]);
% 
%     numImages = length(selectedImages);
%     sqrtNum = ceil(sqrt(numImages)); % Calculate grid size for a somewhat square layout
% 
%     % Keep track of which images to remove
%     removeIndices = false(1, numImages);
% 
%     for i = 1:numImages
%         % Load image
%         img = imread(fullfile(imageDir, selectedImages{i}));
% 
%         % Display image in subplot
%         ax = subplot(sqrtNum, sqrtNum, i);
%         imshow(img, 'Parent', ax);
%         title(ax, sprintf('Image %d', i));
% 
%         % Calculate "X" button position relative to the subplot
%         btnPos = [ax.Position(1) + ax.Position(3) - 0.05, ax.Position(2) + ax.Position(4) - 0.05, 0.03, 0.03];
% 
%         % Add "X" button or clickable text
%         btn = uicontrol('Style', 'pushbutton', 'String', 'X', 'Units', 'normalized', ...
%                         'Position', btnPos, 'Callback', {@removeImageCallback, i});
%     end
% 
%     % Add a "Skip" button to bypass the review process
%     skipBtn = uicontrol('Style', 'pushbutton', 'String', 'Skip Review', 'Units', 'normalized', ...
%                         'Position', [0.35, 0.01, 0.1, 0.05], 'Callback', @skipReviewCallback);
% 
%     % Add a "Done" button to finalize the review process
%     doneBtn = uicontrol('Style', 'pushbutton', 'String', 'Done', 'Units', 'normalized', ...
%                         'Position', [0.55, 0.01, 0.1, 0.05], 'Callback', @doneReviewCallback);
% 
%     % Wait for the user to finish with the figure before continuing
%     uiwait(fig);
% 
%     % Remove selected images after review, if not skipped or removed
%     selectedImages = selectedImages(~removeIndices);
% 
%     % Callback function to remove an image
%     function removeImageCallback(src, event, index)
%         removeIndices(index) = true; % Mark image for removal
%         subplot(sqrtNum, sqrtNum, index); % Select the subplot
%         cla reset; % Clear the subplot
%         title(sprintf('Image %d Removed', index)); % Indicate removal
%         % Check if all images are marked for removal to close the figure
%         if all(removeIndices)
%             uiresume(fig); % Resume execution if figure is still open
%         end
%     end
% 
%     % Callback function for the "Skip" button
%     function skipReviewCallback(src, event)
%         uiresume(fig); % Resume execution, closing the figure and skipping the review
%         close(fig);
%     end
% 
%     % Callback function for the "Done" button
%     function doneReviewCallback(src, event)
%         uiresume(fig); % Resume execution, indicating review completion
%         close(fig);
%     end
% end
% 
% 
% function selections = convertSelectionsToMap(selectedImages)
%     % Convert the list of selected images back into the selections map
%     % format expected by the updateCSV function.
%     selections = containers.Map('KeyType', 'char', 'ValueType', 'any');
%     for i = 1:length(selectedImages)
%         filename = selectedImages{i};
%         % Assume parseAndFilterFilename is a function you've defined elsewhere
%         [isRelevant, testID, particleNumber] = parseAndFilterFilename(filename);
%         if isRelevant
%             if selections.isKey(testID)
%                 % If testID already exists, append particleNumber to its list
%                 selections(testID) = [selections(testID), particleNumber];
%             else
%                 % Otherwise, create a new entry with this testID and particleNumber
%                 selections(testID) = [particleNumber];
%             end
%         end
%     end
%     return; % Return the map of selections.
% end
% 
% 
% function process(raw_tif_filename, SMLs_filename, select_NP )
%      if ~isnumeric(select_NP)
%         error('select_NP must be a numeric array.');
%      end
%     % Add necessary paths (adjust as needed)
%     addpath(genpath('/Users/humphrey/Desktop/YCC_Matlab copy'));
%     disp(['Raw Tif File: ', raw_tif_filename]);
%     disp(['SML Filename: ', SMLs_filename]);
%     disp('Selected Nano-Particles:');
%     disp(select_NP);
% 
% 
%     % Initialization
% 
%     offset = 0;
%     first_frame_NP = 11;
%     first_frame_NP_after_prebleach = 1011;
%     first_frame_GFP = 1;
%     start_frame = 1010; % GFP 0-9, pre-bleach 10-1009 frames
%     end_frame = 21009; % total 40k frames
%     start_frame = start_frame + offset;
%     end_frame = end_frame + offset;
%     minSMLs_driftcorr = (end_frame - start_frame + 1) * 0.7;
%     nm_per_pixel = 116.999998688698;
%     precision_xy_thr = 20;
%     precision_z_thr = inf;
%     GFP_loc_thr = 200;
%     GFP_SMLs_stdev_thr = 30;
%     GFP_SMLs_stdev_thr_z = 60;
%     FWHMthr = 20;
%     search_radius_nm = 200; % particle size = 200 nm
%     runDBSCAN = false;
%     DBSCAN_eps_1 = 15; % 20 nm cluster size
%     if_2D = true;
%     if_outputimage = true;
%     clus_max_axis_thr = inf;
%     clu1_counts = 0;
%     doparse = true;
% 
%      % Adjusting TopDir and Mypath
%     [TopDir, name, ext] = fileparts(raw_tif_filename); % Extract the directory and filename separately
%     raw_tif_file = [name, ext]; % Reconstruct filename with extension
% 
% 
%     matlaboutbase = '/output/';
%     disp(['TopDir: ', TopDir]);
%     disp(['matlaboutbase: ', matlaboutbase]);
%     mkdir(strcat(TopDir, matlaboutbase, 'image_output'));
%     addpath(genpath(TopDir)); % add the specified folder and subfolders to the path.
% 
%     % Remove .csv extension from SMLs_filename if it exists
%     [~, name, ~] = fileparts(SMLs_filename); % Extract the name without extension
%     SMLs_filename = name;
% 
%     parse_particle_filename = strcat(SMLs_filename, '_NP_Particles');
%     delete_NP=[ ...
%            ];
% select_NP=select_NP
%     % Load particle info (adjust as needed)
%     particleFilePath = strcat(TopDir, matlaboutbase, parse_particle_filename, '.mat');
%     if exist(particleFilePath, 'file')
%         load(particleFilePath);
%     else
%         error('Particle file not found: %s', particleFilePath);
%     end
% 
%     % Load NP low res image
%     Im_NP = imread(strcat(TopDir, '/', raw_tif_file)); % Corrected path
%     % imwrite(Im_NP,strcat(TopDir,raw_tif_filename,'NP_f1.tiff'))
%     % tfilebase_all={'20190917 B-';...
%     %                '20190919 B-'};
% 
% 
% 
%     % output x,y,frame
%     output_ini_f=start_frame;
%     output_fin_f=end_frame;
%     output_total_f=end_frame-start_frame+1;
%     frame=(output_ini_f:1:output_fin_f)';
%     logic_frame=zeros((end_frame-start_frame+1),1);
% 
%     output_pixx=[];
%     output_pixy=[];
%     output_pixx_map=[];
%     output_pixy_map=[];
%     output_xlimspix=[];
%     output_ylimspix=[];
%     output_par_num=[];
% 
%     before_prebleach_pix=[];
% 
%     win_size=20; %x frames average
% 
%     for par=1:size(AllParticles,1)
% 
%     %     if sum(delete_NP==AllParticles(par).partnum)>0
%     %         continue
%     %     end
% 
%         if sum(select_NP==AllParticles(par).partnum)==0
%             continue
%         end
% 
% 
%         toutput_x=nan((end_frame-start_frame+1),1);
%         toutput_y=nan((end_frame-start_frame+1),1);
%         toutput_x_map=nan((end_frame-start_frame+1),1);
%         toutput_y_map=nan((end_frame-start_frame+1),1);
% 
%         tframe=AllParticles(par).SMLLabel.Frame;
% 
%         t_tiff_pix=AllParticles(par).partlocpix;
% 
%         txRaw=AllParticles(par).SMLLabel.XRaw_pix_;
%         tyRaw=AllParticles(par).SMLLabel.YRaw_pix_;
% 
%         tx=AllParticles(par).SMLLabel.X_pix_;
%         ty=AllParticles(par).SMLLabel.Y_pix_;
% 
%         if size(tx,1)<minSMLs_driftcorr
%             continue
%         end
% 
%         for i=1:size(tframe,1)
%             if (tframe(i))<(start_frame) ||... %remove SMLs from prebleaching time
%                   (tframe(i))>end_frame      %remove NP SMLs from GFP illumination time
%                 continue
%             end
%             logic_frame(tframe(i)+1-start_frame)=1;
%             toutput_x(tframe(i)+1-start_frame)=tx(i);
%             toutput_y(tframe(i)+1-start_frame)=ty(i);
%             toutput_x_map(tframe(i)+1-start_frame)=tx(i)-t_tiff_pix(1);
%             toutput_y_map(tframe(i)+1-start_frame)=ty(i)-t_tiff_pix(2);
%         end
% 
%         output_pixx=cat(2,output_pixx,toutput_x);
%         output_pixy=cat(2,output_pixy,toutput_y);
% 
%         output_pixx_map=cat(2,output_pixx_map,nanmean(toutput_x_map));
%         output_pixy_map=cat(2,output_pixy_map,nanmean(toutput_y_map));
% 
%         output_xlimspix=cat(1,output_xlimspix,AllParticles(par).xlimspix);
%         output_ylimspix=cat(1,output_ylimspix,AllParticles(par).ylimspix);
% 
%         output_par_num=cat(1,output_par_num,AllParticles(par).partnum);
% 
%         before_prebleach_pix=cat(1,before_prebleach_pix,t_tiff_pix);
%     end
%     center_pixx=nanmean(output_pixx,1);
%     center_pixy=nanmean(output_pixy,1);
%     local_output_pixx=output_pixx-ones(size(output_pixx,1),1)*center_pixx;
%     local_output_pixy=output_pixy-ones(size(output_pixy,1),1)*center_pixy;
%     mean_pixx=nanmean(local_output_pixx,2);
%     mean_pixy=nanmean(local_output_pixy,2);
% 
%     mean_pixx_map=output_pixx_map;
%     mean_pixy_map=output_pixy_map;
% 
%     plot_merge_drift=figure;
%     meanx_nm=mean_pixx*nm_per_pixel;
%     meany_nm=mean_pixy*nm_per_pixel;
%     pixel_thispart_drift=[meanx_nm, meany_nm];
%     c = linspace(1,10,length(pixel_thispart_drift(:,1)));
%     ax2_drift = axes;
%     plot_clus_drift=scatter(ax2_drift,pixel_thispart_drift(:,1), pixel_thispart_drift(:,2), 1, c, 'filled');
%     xlabel('nm')
%     ylabel('nm')
%     hold on
%     colormap(ax2_drift,'jet')
%     OutDir_partnum_drift = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename);
%     saveas(plot_merge_drift, strcat(OutDir_partnum_drift, '_NP_drift.png'));
%     close(plot_merge_drift);
% 
%     fixed_win_x=zeros(floor(output_total_f/win_size),1);
%     fixed_win_y=zeros(floor(output_total_f/win_size),1);
% 
%     %calculate the drift from the initial xy
%     for frame_win=1:floor(output_total_f/win_size)
%         fixed_win_x(frame_win)=nanmean(mean_pixx((frame_win-1)*win_size+1:(frame_win)*win_size));
%         fixed_win_y(frame_win)=nanmean(mean_pixy((frame_win-1)*win_size+1:(frame_win)*win_size));
%     end
%     ini_x=fixed_win_x(1);
%     ini_y=fixed_win_y(1);
% 
%     %drift correction
%     driftcorr_matrix_x=nan(size(output_pixx,1),1);
%     driftcorr_matrix_y=nan(size(output_pixy,1),1);
%     driftcorr_output_pixx=nan(size(output_pixx));
%     driftcorr_output_pixy=nan(size(output_pixy));
%     for frame_win=1:floor(output_total_f/win_size)
%         tdriftx=fixed_win_x(frame_win)-ini_x;
%         tdrifty=fixed_win_y(frame_win)-ini_y;
% 
%         driftcorr_output_pixx((frame_win-1)*win_size+1:(frame_win)*win_size,:)=...
%             output_pixx((frame_win-1)*win_size+1:(frame_win)*win_size,:)-tdriftx*ones(win_size,size(output_pixx,2));
%         driftcorr_output_pixy((frame_win-1)*win_size+1:(frame_win)*win_size,:)=...
%             output_pixy((frame_win-1)*win_size+1:(frame_win)*win_size,:)-tdrifty*ones(win_size,size(output_pixy,2));
% 
%         driftcorr_matrix_x((frame_win-1)*win_size+1:(frame_win)*win_size,1)= -tdriftx*ones(win_size,1);
%         driftcorr_matrix_y((frame_win-1)*win_size+1:(frame_win)*win_size,1)= -tdrifty*ones(win_size,1);
% 
%     end
%     %remainder from the win_size
%     if rem(output_total_f,win_size)~=0
%       driftcorr_output_pixx((frame_win)*win_size:end,:)=...
%             output_pixx((frame_win)*win_size:end,:)-tdriftx*ones(rem(output_total_f,win_size),size(output_pixx,2));
%       driftcorr_output_pixy((frame_win)*win_size:end,:)=...
%             output_pixy((frame_win)*win_size:end,:)-tdrifty*ones(rem(output_total_f,win_size),size(output_pixy,2));
% 
%       driftcorr_matrix_x((frame_win)*win_size:end,:)= -tdriftx*ones(rem(output_total_f,win_size),1);
%       driftcorr_matrix_y((frame_win)*win_size:end,:)= -tdrifty*ones(rem(output_total_f,win_size),1);
%     end
% 
%     driftcorr_matrix=[driftcorr_matrix_x driftcorr_matrix_y];
%     save(strcat(TopDir, matlaboutbase, SMLs_filename, '_driftcorr.mat'), 'driftcorr_matrix');
% 
%     total_NP_num=0;
%     pass_FWHMthr_num_ori=0;
%     pass_FWHMthr_num_corr=0;
% 
%     Im_low=min(min(Im_NP));
%     Im_high=max(max(Im_NP))*0.4;
% 
%     %%NP loc from Image after prebleach
%     % Parse SMLs with NP presenting; Peak finder to locate centroids in GFP image
%     prebleach_tif_file = strrep(raw_tif_file, '_NP_f1.tif', '_NP_f1001.tif');
%     Im_NP_after_prebleach = imread(strcat(TopDir, '/', prebleach_tif_file));
%     peaksizepix = 5;
%     peakintensitythresh =300;
%     [centr_lowres_NP_after_prebleach] = one_channel_partID_includeSide(Im_NP_after_prebleach, peaksizepix, peakintensitythresh, []);
% 
%     % NP peaks
%     % centr_lowres is
%     % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_after_prebleach;
%     % filter resolvpow>0.8
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_filt_after_prebleach(centr_lowres_NP_filt_after_prebleach(:,7)>0.8,:);
%     % filter perc_backpix>50%
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_filt_after_prebleach(centr_lowres_NP_filt_after_prebleach(:,8)>50,:);
%     % filter SNR>1.5
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_filt_after_prebleach(centr_lowres_NP_filt_after_prebleach(:,6)>5,:);
% 
%     NPpix_after_prebleach = centr_lowres_NP_filt_after_prebleach(:,[1,2]);
% 
%     prebleach_driftcorr_x=[];
%     prebleach_driftcorr_y=[];
%         for par=1:size(output_pixx,2)
% 
% 
% 
%             %% correction of drift in prebleach
%             tsearch_pixx=before_prebleach_pix(par,1)+[-1,1]*search_radius_nm/nm_per_pixel;
%             tsearch_pixy=before_prebleach_pix(par,2)+[-1,1]*search_radius_nm/nm_per_pixel;
% 
%             tpos=NPpix_after_prebleach((NPpix_after_prebleach(:,1) > tsearch_pixx(1)...
%                                       & NPpix_after_prebleach(:,1) < tsearch_pixx(2)...
%                                       & NPpix_after_prebleach(:,2) > tsearch_pixy(1)...
%                                       & NPpix_after_prebleach(:,2) < tsearch_pixy(2)),:);
% 
%             if size(tpos,1)==0
%                 disp('no NP after prebleach step'); continue
%             elseif size(tpos,1)>1
%                 disp('too dense of NP in search field'); continue
%             end
% 
%             prebleach_driftcorr_x=cat(2,prebleach_driftcorr_x,tpos(1)-before_prebleach_pix(par,1));
%             prebleach_driftcorr_y=cat(2,prebleach_driftcorr_y,tpos(2)-before_prebleach_pix(par,2));
%         end
%             prebleach_driftcorr_x
%             prebleach_driftcorr_y
% 
%            %% image output
% 
%         for par=1:size(output_pixx,2)
%             total_NP_num=total_NP_num+1;
% 
%             if if_outputimage
% 
%             %scatter + low res image
%             plot_merge_ori=figure;
%             pixel_thispart_ori=[output_pixx(:,par)-mean_pixx_map(par), output_pixy(:,par)-mean_pixy_map(par)];
%     %          pixel_thispart_ori=[output_pixx(:,par), output_pixy(:,par)];
% 
%             ax1 = axes;
%             plot_SMLs=scatter(ax1,pixel_thispart_ori(:,1), pixel_thispart_ori(:,2), 20, [1,0,0], '.');
%             hold on
%             colormap(ax1,'jet')
%             scale_bar_coord_ori=[output_xlimspix(par,2)-120/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel;...
%                 output_xlimspix(par,2)-20/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel];
%             plot_scalebar=plot(scale_bar_coord_ori(1:2,1),scale_bar_coord_ori(1:2,2),'-c', 'LineWidth', 3);hold on
%     %             pixel_clus = cell(size(clus1,1),1);
%     %             for k = 1: size(clus1,1)
%     %                 pixel_clus = clus1{k,1}.Points*transform_matrix+shift_matrix;
%     %                 plot_clus=plot(alphaShape(pixel_clus(:,1),pixel_clus(:,2)),'FaceAlpha',0.1,'FaceColor',[1,0,0],'EdgeAlpha',0);
%     %                 hold on
%     %             end
% 
%             plot_partnum = imshow(Im_NP, [Im_low, Im_high*0.7]);
%             axisgfp = gca;
%             axisgfp.XLim = output_xlimspix(par,:);
%             axisgfp.YLim = output_ylimspix(par,:);
%             axisgfp.Visible = 'off';
%             hold off
%             alpha(0.5)
% 
%             OutDir_partnum_ori = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_merge_ori, strcat(OutDir_partnum_ori, '_NP_merge_ori.png'));
%             close(plot_merge_ori);
%             end
% 
% 
% 
%             t_pixx=(output_pixx(:,par)-center_pixx(par))*nm_per_pixel;
%             t_pixy=(output_pixy(:,par)-center_pixy(par))*nm_per_pixel;
%             pdx=fitdist(t_pixx,'Normal');
%             pdy=fitdist(t_pixy,'Normal');
% 
%             pdx_fwhm= pdx.sigma*2*sqrt(2*log(2));
%     %             pdx_ci = paramci(pdx);
%     %             pdx_fwhmCI=round(pdx_ci(:,2)'*2*sqrt(2*log(2)),1);
% 
%             pdy_fwhm= pdy.sigma*2*sqrt(2*log(2));
%     %             pdy_ci = paramci(pdy);
%     %             pdy_fwhmCI=round(pdy_ci(:,2)'*2*sqrt(2*log(2)),1);
%             if pdx_fwhm<FWHMthr && pdy_fwhm<FWHMthr
%                pass_FWHMthr_num_ori=pass_FWHMthr_num_ori+1; 
%             end
% 
%             if if_outputimage
%             % histogram of fiducial markers
%             plot_hist_ori = figure;
%             plot_hist_ori.Units = 'inches';
%             plot_hist_ori.Position = [1 1 4 4];
%             plot_hist_ori.PaperPosition = [0 0 4 4];
% 
%             ax3_ori_x = axes('Position', [0.15 0.58 0.8 0.4]);
%             p_histx=histogram(t_pixx,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%     %             xlabel('nm');
%             ylabel('x pdf');
%             str={'FWHM',strcat(num2str(round(pdx_fwhm,1)),'nm')};
%     %                 strcat('CI: ',num2str(pdx_fwhmCI(1)),'-',num2str(pdx_fwhmCI(2)))};
%             peak_ctr_x=max(p_histx.BinCounts)/sum(p_histx.BinCounts)/p_histx.BinWidth;
%             text(pdx.mu,peak_ctr_x-0.01,str,'Color','red','FontSize',12);
%             xgrid = linspace(p_histx.BinEdges(1),p_histx.BinEdges(end),100)';
%             pdfx_Est = pdf(pdx,xgrid);
%             line(xgrid,pdfx_Est)
% 
%             ax3_ori_y = axes('Position', [0.15 0.12 0.8 0.4]);
%             p_histy=histogram(t_pixy,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%             xlabel('nm');
%             ylabel('y pdf');
%             str={'FWHM',strcat(num2str(round(pdy_fwhm,1)),'nm')};
%     %                strcat('CI: ',num2str(pdy_fwhmCI(1)),'-',num2str(pdy_fwhmCI(2)))};
%             peak_ctr_y=max(p_histy.BinCounts)/sum(p_histy.BinCounts)/p_histy.BinWidth;
%             text(pdy.mu,peak_ctr_y-0.01,str,'Color','red','FontSize',12);
%             ygrid = linspace(p_histy.BinEdges(1),p_histy.BinEdges(end),100)';
%             pdfy_Est = pdf(pdy,ygrid);
%             line(ygrid,pdfy_Est)
% 
%             OutDir_partnum_ori = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_hist_ori, strcat(OutDir_partnum_ori, '_NP_hist_ori.png'));
%             close(plot_hist_ori);
%             end
% 
% 
%             if if_outputimage
%             %scatter + low res image
%             plot_merge_corr=figure;
%     %         pixel_thispart_corr=[driftcorr_output_pixx(:,par), driftcorr_output_pixy(:,par)];
%              pixel_thispart_corr=[driftcorr_output_pixx(:,par)-mean_pixx_map(par), driftcorr_output_pixy(:,par)-mean_pixy_map(par)];
% 
%             ax1 = axes;
%             plot_SMLs=scatter(ax1,pixel_thispart_corr(:,1), pixel_thispart_corr(:,2), 20, [1,0,0], '.');
%             hold on
%             colormap(ax1,'jet')
%             scale_bar_coord_corr=[output_xlimspix(par,2)-120/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel;...
%                 output_xlimspix(par,2)-20/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel];
%             plot_scalebar=plot(scale_bar_coord_corr(1:2,1),scale_bar_coord_corr(1:2,2),'-c', 'LineWidth', 3);hold on
%     %             pixel_clus = cell(size(clus1,1),1);
%     %             for k = 1: size(clus1,1)
%     %                 pixel_clus = clus1{k,1}.Points*transform_matrix+shift_matrix;
%     %                 plot_clus=plot(alphaShape(pixel_clus(:,1),pixel_clus(:,2)),'FaceAlpha',0.1,'FaceColor',[1,0,0],'EdgeAlpha',0);
%     %                 hold on
%     %             end
% 
%             plot_partnum = imshow(Im_NP, [Im_low, Im_high*0.7]);
%             axisgfp = gca;
%             axisgfp.XLim = output_xlimspix(par,:);
%             axisgfp.YLim = output_ylimspix(par,:);
%             axisgfp.Visible = 'off';
%             hold off
%             alpha(0.5)
% 
%             OutDir_partnum_corr = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_merge_corr, strcat(OutDir_partnum_corr, '_NP_merge_corr.png'));
%             close(plot_merge_corr);
%             end
% 
% 
% 
% 
%             t_pixx=(driftcorr_output_pixx(:,par)-center_pixx(par))*nm_per_pixel;
%             t_pixy=(driftcorr_output_pixy(:,par)-center_pixy(par))*nm_per_pixel;
%             pdx=fitdist(t_pixx,'Normal');
%             pdy=fitdist(t_pixy,'Normal');
% 
%             pdx_fwhm= pdx.sigma*2*sqrt(2*log(2));
%     %             pdx_ci = paramci(pdx);
%     %             pdx_fwhmCI=round(pdx_ci(:,2)'*2*sqrt(2*log(2)),1);
% 
%             pdy_fwhm= pdy.sigma*2*sqrt(2*log(2));
%     %             pdy_ci = paramci(pdy);
%     %             pdy_fwhmCI=round(pdy_ci(:,2)'*2*sqrt(2*log(2)),1);
%             if pdx_fwhm<FWHMthr && pdy_fwhm<FWHMthr
%                pass_FWHMthr_num_corr=pass_FWHMthr_num_corr+1; 
%             end
% 
%             if if_outputimage
%             % histogram of fiducial markers
%             plot_hist_corr = figure;
%             plot_hist_corr.Units = 'inches';
%             plot_hist_corr.Position = [1 1 4 4];
%             plot_hist_corr.PaperPosition = [0 0 4 4];
% 
%             ax3_corr_x = axes('Position', [0.15 0.58 0.8 0.4]);
%             p_histx=histogram(t_pixx,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%     %             xlabel('nm');
%             ylabel('x pdf');
%             str={'FWHM',strcat(num2str(round(pdx_fwhm,1)),'nm')};
%     %                 strcat('CI: ',num2str(pdx_fwhmCI(1)),'-',num2str(pdx_fwhmCI(2)))};
%             peak_ctr_x=max(p_histx.BinCounts)/sum(p_histx.BinCounts)/p_histx.BinWidth;
%             text(pdx.mu,peak_ctr_x-0.01,str,'Color','red','FontSize',12);
%             xgrid = linspace(p_histx.BinEdges(1),p_histx.BinEdges(end),100)';
%             pdfx_Est = pdf(pdx,xgrid);
%             line(xgrid,pdfx_Est)
% 
%             ax3_corr_y = axes('Position', [0.15 0.12 0.8 0.4]);
%             p_histy=histogram(t_pixy,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%             xlabel('nm');
%             ylabel('y pdf');
%             str={'FWHM',strcat(num2str(round(pdy_fwhm,1)),'nm')};
%     %                strcat('CI: ',num2str(pdy_fwhmCI(1)),'-',num2str(pdy_fwhmCI(2)))};
%             peak_ctr_y=max(p_histy.BinCounts)/sum(p_histy.BinCounts)/p_histy.BinWidth;
%             text(pdy.mu,peak_ctr_y-0.01,str,'Color','red','FontSize',12);
%             ygrid = linspace(p_histy.BinEdges(1),p_histy.BinEdges(end),100)';
%             pdfy_Est = pdf(pdy,ygrid);
%             line(ygrid,pdfy_Est)
% 
%             OutDir_partnum_corr = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_hist_corr, strcat(OutDir_partnum_corr, '_NP_hist_corr.png'));
%             close(plot_hist_corr);
%             end
% 
% 
%         end 
% 
% 
%       pass_FWHMthr_num_ori/total_NP_num
%       pass_FWHMthr_num_corr/total_NP_num
%       pass_FWHMthr_num_corr
% 
% 
%     end
% 
