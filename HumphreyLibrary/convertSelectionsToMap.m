function selections = convertSelectionsToMap(selectedImages)
    % Convert the list of selected images back into the selections map
    % format expected by the updateCSV function.
    selections = containers.Map('KeyType', 'char', 'ValueType', 'any');
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
    return; % Return the map of selections.
end