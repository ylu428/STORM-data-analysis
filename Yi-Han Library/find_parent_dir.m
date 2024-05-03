function parent_path = find_parent_dir(variable)
% Find the upper-level folder which saves all the folders recorded in
% "variable". "variable" need to be a structure.
    names = fieldnames(variable);
    firstpath = names{1};
    lastpath = names{end};
    if length(names)==1
        parts = strsplit(variable.(lastpath), filesep);
        parent_path = strcat(strjoin(parts(1:length(parts)-2), filesep), '\');
    else
        diff_folder = setdiff(strsplit(variable.(lastpath), filesep), ...
            strsplit(variable.(firstpath), filesep));
        parts = strsplit(variable.(lastpath), filesep);
        parent_path = strcat(strjoin(parts(1:find(strcmp(parts, diff_folder{1}))-1), filesep), '\');
    end

end