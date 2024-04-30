function [tblparse] = findSMLs_mapred_ONI(ds, xlims, ylims, Outputdir)
% extract entries from datastore that lie within x- and ylims
%   Detailed explanation goes here

% parse on 'X Raw pix' and 'Y Raw pix' cols 13 and 14
inbetween = @(data) data{:,13} > xlims(1,1) & data{:,13} < xlims(1,2) &...
                    data{:,14} > ylims(1,1) & data{:,14} < ylims(1,2);

% see
% 'https://www.mathworks.com/help/matlab/import_export/simple-data-subsetting-using-mapreduce.html'
configuredMapper = ...
    @(data, info, intermKVStore) subsettingMapperGeneric(data, info, ...
    intermKVStore, inbetween);

result = mapreduce(ds, configuredMapper, @subsettingReducer,...
    'OutputFolder', strcat(Outputdir,'MapReduceFiles'));

a = readall(result);
tblparse = a.Value{1};

%delete MapReduceFiles\*.mat
end