% updateCSV.m
% This script updates a CSV file with the corresponding images and nano  particiles extracted from the selection process

% Function to update a CSV file
function updateCSV(filename, data)
    % Read the existing CSV data
    existingData = readtable(filename);
    
    % Update the existing data with the new data
    updatedData = [existingData; data];
    
    % Write the updated data back to the CSV file
    writetable(updatedData, filename);
end
