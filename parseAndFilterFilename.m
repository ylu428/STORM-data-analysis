function [isRelevant, testID, particleNumber] = parseAndFilterFilename(filename)
    
    isRelevant = false;  % Default value
    testID = '';
    particleNumber = 0;
    
    % Check if filename ends with '_NP_SML.png'
    if endsWith(filename, '_NP_SML.png')
        % Adjusted regex to capture test ID and particle number
        % This regex accounts for a prefix, followed by the test ID (with or without a hyphen) and the particle number
        tokens = regexp(filename, '(.*?)_SMLs?-part-(\d+)_NP_SML\.png', 'tokens');
        if ~isempty(tokens)
            token = tokens{1};  % Assuming one match, which is a cell containing the captured groups
            testID = token{1};  % The test ID is the first capture group
            particleNumber = str2double(token{2});  % The particle number is the second capture group
            isRelevant = true;
        else
            % Handle the case where the filename does not match the expected pattern
            isRelevant = false;
        end

    end
end