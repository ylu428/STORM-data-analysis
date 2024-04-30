function [output, A]=uniq_Pdist(data)
    A = unique(data(:,1:3), 'rows','stable');
    a = pdist(A,'euclidean')';
    pair = {};
    for i = 1:(size(A,1)-1)
        for j = i+1:size(A,1)
            temp = [i,j];
            pair = cat(1, pair, temp);
            % text = ['dist between points: ', num2str(i), ',', num2str(j)];
            % disp(text);
        end
    
    end
    output = cat(2, a, cell2mat(pair));
end

