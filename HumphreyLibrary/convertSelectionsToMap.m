% Function to convert a list of selected images back into a selections map.
% This map format is expected by the updateCSV function.

function selections = convertSelectionsToMap(selectedImages)
    % Initialize a container map to store selections with test IDs and particle numbers
    selections = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Loop through each selected image
    for i = 1:length(selectedImages)
        filename = selectedImages{i};
        
        % Assume parseAndFilterFilename is a function you've defined elsewhere
        [isRelevant, testID, particleNumber] = parseAndFilterFilename(filename);
        
        if isRelevant
            if selections.isKey(testID)
                % If testID already exists, append particleNumber to its list
                selections(testID) = [selections(testID), particleNumber];
            else
                % Otherwise, create a new entry with this testID and particleNumber
                selections(testID) = [particleNumber];
            end
        end
    end
    
    % Return the map of selections
    return;
end
