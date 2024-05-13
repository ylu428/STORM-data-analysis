% This function parses and filters filenames based on a specific pattern.
% This script checks if the filename ends with '_NP_SML.png', extracts the test ID and particle number and determines if the filename is relevant based on the pattern.
% '_NP_SML.png' can be changed to whatever type of particle needed for that experiment.


function [isRelevant, testID, particleNumber] = parseAndFilterFilename(filename)
    
    isRelevant = false;  % Default value indicating the filename is not relevant
    testID = '';         % Initialize testID as an empty string
    particleNumber = 0;  % Initialize particleNumber as zero
    
    % Check if the filename ends with '_NP_SML.png'
    if endsWith(filename, '_NP_SML.png')
        % Regex to capture test ID and particle number from the filename
        % The pattern accounts for a prefix, followed by the test ID (with or without a hyphen) and the particle number
        tokens = regexp(filename, '(.*?)_SMLs?-part-(\d+)_NP_SML\.png', 'tokens');
        if ~isempty(tokens)
            token = tokens{1};  % Assuming one match, which is a cell containing the captured groups
            testID = token{1};  % The test ID is the first capture group
            particleNumber = str2double(token{2});  % The particle number is the second capture group
            isRelevant = true;  % Mark the filename as relevant
        else
            % Handle the case where the filename does not match the expected pattern
            isRelevant = false;
        end
    end
end
