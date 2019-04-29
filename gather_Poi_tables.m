function [poi_table, poi_counts, poi_indices] = gather_Poi_tables(directory, max_fov)
% @param:
%       directory Directory of POI tables
%       max_fov Maximum number of POI tables to load
% @return:
%       poi_table POI table with inputs of all POI tables
%       poi_counts Vector of integers representing number of POI's in each
%                  POI table

% Last Modified: 4/29/2019

%% Extract POI table from directory

poi_counts = zeros(max_fov, 1);
poi_table = zeros(1, 15);
poi_indices = zeros(1, 4);

for iFov = 1:max_fov
    file = fullfile(directory, sprintf('POI_Table_FOV_%d.mat', iFov));
    if ~exist(file, 'file')
        continue;
    end
    load(file, 'POI_Table');
    POI_Table = full(POI_Table);
    [poi_table_i, poi_counts(iFov)] = clean_fov(POI_Table);
    poi_table = [poi_table; poi_table_i];
    fov_number = ones(poi_counts(iFov), 1) * iFov;
    poi_indices = [poi_indices; [fov_number, poi_table_i(1:3:end, 1:3)]];
end
poi_table(1, :) = [];
poi_indices(1, :) = [];


function [poi_table, poi_count] = clean_fov(poi_table)
%
nPoi = size(poi_table, 1) / 3;
dlt = false(3 * nPoi, 1);

for iPoi = 1:nPoi
    iPosition = 3 * (iPoi - 1) + 1;
    [~, idx] = max(poi_table(iPosition:iPosition+2, 4:end));
    matchN = sum(idx == [2, 4, 3, 2, 3, 3, 1, 2, 3, 1, 2, 4]);
    if matchN < 8 || ~within_border(poi_table(iPosition, 2:3))
        dlt(iPosition:iPosition+2) = 1;
    end
end
poi_table(dlt, :) = [];
poi_count = size(poi_table, 1) / 3;

function isWithin = within_border(coordinate)
x = coordinate(1);
y = coordinate(2);
isWithin = x > 5 && x < 1019 && y > 5 && y < 1019;