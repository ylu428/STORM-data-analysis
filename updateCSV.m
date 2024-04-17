function updateCSV(csvFilePath, selections)
    fid = fopen(csvFilePath, 'w');
    keys = selections.keys;
    
    for i = 1:length(keys)
        testID = keys{i};
        particleNumbers = selections(testID);
        % Format the line as 'TestID: number number number...'
        fprintf(fid, '%s: %s\n', testID, num2str(particleNumbers, '%d '));
    end
    
    fclose(fid);
end