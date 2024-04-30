function [tblparse] = findSMLs_mapred_ONI_nm_colormapping(ds, xlims, ylims, Outputdir)
% extract entries from datastore that lie within x- and ylims
%   Detailed explanation goes here

% parse on 'X nm' and 'Y nm' cols 3 and 4
% 'X nm' and 'Y nm' are nm loc from 2-color mapping of ONI
inbetween = @(data) data{:,3} > xlims(1,1) & data{:,3} < xlims(1,2) &...
                    data{:,4} > ylims(1,1) & data{:,4} < ylims(1,2);

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