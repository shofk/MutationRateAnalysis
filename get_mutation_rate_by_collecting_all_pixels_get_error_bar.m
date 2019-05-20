out_dir = 'D:\\MutationRateProject\\_20190520_%s_analysis';
isTaq = input('Taq:1 or BST:0?');


max_fov = 120;
graph_ = false;
save_image = false;


if isTaq
    directory = 'D:\MutationRateProject\Taq_local_max\Point_of_Interest_table';
    out_dir = sprintf(out_dir, 'Taq');
    ratios_low = ones(3, 12) * 0.26;
    ratios_high = ones(3, 12) * 1.50;
    to_gta_cutoff = 0.30;
    to_c_cutoff = 0.4;
    c_low = 0.50;
else % isBST
    directory = 'D:\MutationRateProject\BST1_max_itst\Point_of_Interest_table';
    out_dir = sprintf(out_dir, 'BST');
    ratios_low = ones(3, 12) * 0.20;
    ratios_high = ones(3, 12) * 1.50;
    to_gta_cutoff = 0.25;
    to_c_cutoff = 0.3;
    c_low = 0.25;
end

cutoff_ratio_low = 0.26;
cutoff_ratio_high = 1.5;
c_high = 1.5;

%% Checks if the range of cutoff points are decisive
just_one = true;
len_cutoff_ratio = length(cutoff_ratio_low);
len_c_cutoff = length(c_low);
if len_cutoff_ratio ~= 1 || length(c_low) ~= 1
    just_one = false;
    mutation_gradient = cell(3, 12);
    mutation_gradient_C = cell(1,12);
    for cyc = 1:12
        mutation_gradient_C{cyc} = zeros(len_c_cutoff);
        for ch = 1:3
            mutation_gradient{ch, cyc} = zeros(len_cutoff_ratio);
        end
    end
end

%% Gather initial information
[poi_master, poi_counts, poi_indices] = gather_Poi_tables(directory, max_fov); % Collect the unprocessed intensities
poi = poi_master(:, 4:end);
nPoi = size(poi,1) / 3;

if ~exist(out_dir, 'dir') && save_image
    mkdir(out_dir);
end


cycles = cell(4,1);
cycles{1} = [7, 10];      % G
cycles{2} = [1, 4, 8, 11];% T
cycles{3} = [3, 5, 6, 9]; % A
cycles{4} = [2, 12];      % C

%% Count and collect mutations
image_name = 'mutation_map_FOV%d_cycle%d_mutation_ch%s_to_ch%s.tif';
FOVs = unique(poi_indices(:,1))';
nFov = max(FOVs);
nMutations = cell(nFov, 1);
channel = 1:3;
bases = 'GTAC';
cyc = 1:12;
poi_c = poi_master;
b = {'G', 'T', 'A'};

for ifov = FOVs
    if save_image
        mkdir(fullfile(out_dir, sprintf('FOV%d', ifov)));
    end
    out_directory = fullfile(out_dir, sprintf('FOV%d', ifov));
    poiLocs = find(poi_indices(:,1) == ifov);
    offset = poiLocs(1)-1;
    itst = extract_and_correct_phase2(poi((poiLocs(1)-1)*3+1:(poiLocs(end)-1)*3+3, :), graph_);
    median_itst = get_median_intensities(itst);
    iMut = cell(4, 12);
    for icyc = cyc
        for ich = 1:4
            if sum(ismember(cycles{ich}, icyc))
                chN = ich;
                break;
            end
        end
        citst = reshape(itst(:, icyc), 3, numel(poiLocs))';
        medians = median(citst);
        medians = medians .* (medians < 0.5);
        citst = citst - medians;
        % plot cycle intensities
        if graph_
            figure; boxplot(citst, b);
            ylim([-0.5 2]);
            yticks(-0.5:0.1:2);
            title(sprintf('FOV%d Cycle%d', ifov, icyc));
        end
        if chN == 4 % when C cycle
            for iich = channel
                mut_itst = citst(:, iich);
%                 idx = find(mut_itst > c_low + median(mut_itst) & sum(citst < to_c_cutoff * median_itst, 2) ~= 3 ...
%                     & mut_itst < c_high + median(mut_itst) ...
%                     & sum(citst > c_low + median(mut_itst) & citst < c_high + median(mut_itst), 2) == 1) - 1 + offset;
                idx = find(mut_itst > std(mut_itst) * 1.5 & sum(citst < to_gta_cutoff * median_itst, 2) == 2) - 1 + offset;
                % record the number of mutations
                iMut{iich, icyc} = idx;
                if save_image
                    map = zeros(1024);
                    map(poi_master(idx * 3 + 1, 1)) = 1;
                    imwrite(imdilate(map, strel('disk', 5)),...
                        fullfile(out_directory, sprintf(image_name, ifov, icyc, 'C', bases(iich))));
                end
            end
        else % when G, T, A cycles
            main_itst = citst(:, chN);
            for iich = channel
                if chN ~= iich
                    mut_itst = citst(:, iich);
                    if isTaq
                        idx = find(mut_itst > main_itst * ratios_low(iich, icyc) & sum(citst < to_c_cutoff * median_itst, 2) == 1 ...
                            & mut_itst < main_itst * ratios_high(iich, icyc) & (max(citst, [], 2) == main_itst | max(citst, [], 2) == mut_itst)) - 1 + offset;
                    else
                        idx = find(mut_itst > main_itst * ratios_low(iich, icyc) & sum(citst < to_c_cutoff * median_itst, 2) == 1 ...
                            & mut_itst < main_itst * ratios_high(iich, icyc) & (max(citst, [], 2) == main_itst | max(citst, [], 2) == mut_itst)) - 1 + offset;
                    end
                        % record the number of mutations
                    iMut{iich, icyc} = idx;
                    if save_image
                        map = zeros(1024);
                        map(poi_master(idx * 3 + 1, 1)) = 1;
                        imwrite(imdilate(map, strel('disk', 5)),...
                            fullfile(out_directory, sprintf(image_name, ifov, icyc, bases(chN), bases(iich))));
                    end
                end
            end
            % add mutation to C
            idx = find(sum(citst < to_c_cutoff * median_itst, 2) == 3) - 1 + offset;
            % record the number of mutations
            iMut{4, icyc} = idx;
            if save_image
                map = zeros(1024);
                map(poi_master(idx * 3 + 1, 1)) = 1;
                imwrite(imdilate(map, strel('disk', 5)),...
                    fullfile(out_directory, sprintf(image_name, ifov, icyc, bases(chN), 'C')));
            end
        end
    end
    nMutations{ifov} = iMut;
    poi_c((poiLocs(1)-1)*3+1:(poiLocs(end)-1)*3+3, 4:end) = itst;
    fprintf('FOV%d complete\n', ifov);
end
%%

if just_one
    final_rate = compile_mutation_rate(nMutations, FOVs, nPoi);
    [cycle_mutation_rate, by_channel_mutation_rates] = compile_cycle_mutation_rate(nMutations, FOVs, nPoi, poi_counts);

    compile_results_with_error_bars_8(nMutations, FOVs, out_dir, poi_counts);
    
    xlswrite(fullfile(out_dir, 'mutation_final_rate.xlsx'), final_rate);
    xlswrite(fullfile(out_dir, 'mutation_cycle_rate.xlsx'), cycle_mutation_rate);
end
%%
% for icyc = cyc
%     chN = 4;
%     for ich = channel
%         if ~isempty(intersect(icyc, cycles{ich}))
%             chN = ich;
%         end
%     end
%     if chN ~= 4
%         main_itst = poi(chN:3:end, icyc);
%         for ich = channel
%             if chN ~= ich
%                 mut_itst = poi(ich:3:end, icyc);
%                 idx = find(mut_itst > main_itst * (1/3));
%                 map = zeros(1024);
%                 idx = intersect(idx, find(poi_indices(:, 1)==fov));
%                 map(poi_indices(idx, 2)) = 1;
%                 imwrite(imdilate(map, strel('disk', 10)), fullfile(out_directory, sprintf(image_name, fov, icyc, bases(ich))));
%             end
%         end
%     else
%         for ich = channel
%             mut_itst = poi(ich:3:end, icyc);
%             idx = find(mut_itst > (1/3));
%             map = zeros(1024);
%             idx = intersect(idx, find(poi_indices(:,1) == fov));
%             map(poi_indices(idx, 2)) = 1;
%             imwrite(imdilate(map, strel('disk', 10)), fullfile(out_directory, sprintf(image_name, fov, icyc, 'C')));
%         end
%     end
%     
%     
% end






disp('Complete');

%% Functions
function phasing = create_phasing_matrix(pre, post, ncyc)
    phasing=zeros(ncyc);
    phasing(1,1)=1-pre-post;
    phasing(1,2)=pre;
    for i=2:ncyc
        phasing(i,:) = post*phasing(i-1,:) + ...
            (1-post-pre)*[0 phasing(i-1,1:end-1)] +...
            pre*[0 0 phasing(i-1,1:end-2)]; 
        % assume in each cycle, the same proportion of phasing and pre-phasing molecules occurs.
    end
end

function plot_intensities(poi)
    % get intensities
    g = poi(1:3:end, :);
    t = poi(2:3:end, :);
    a = poi(3:3:end, :);

    % get median
%     gMed = median(g, 1);
%     tMed = median(t, 1);
%     aMed = median(a, 1);

%     figure; hold on;
%     title('Median by cycles');
%     plot(gMed, 'b');
%     plot(tMed, 'g');
%     plot(aMed, 'r');
%     legend('G', 'T', 'A');
%     ylim([0 1]);
%     hold off;

    figure;
    boxplot(g, 'symbol', 'bo', 'outliersize', 1);
    % ylim([0 1]);
    title('G channel intensity flow');

    figure;
    boxplot(t, 'symbol', 'go', 'outliersize', 1);
    % ylim([0 1]);
    title('T channel intensity flow');

    figure;
    boxplot(a, 'symbol', 'ro', 'outliersize', 1);
    % ylim([0 1]);
    title('A channel intensity flow');
end

function [receiver, source] = correct_phasing(receiver, source, ratio)
    itst_passed = zeros(size(source));
    itst_passed(source > 0) = source(source > 0) * ratio;
    source = source + itst_passed;
    receiver(receiver > 0) = receiver(receiver > 0) - itst_passed(receiver > 0);
end

function gathered_Poi_table = gather_poi_tables(directory)
    gathered_Poi_table = [];

    format = 'POI_Table_FOV_%d.mat';

    for iFov = 1:95
        file = fullfile(directory, sprintf(format, iFov));
        if exist(file, 'file')
            load(file, 'POI_Table');
            POI_Table = full(POI_Table);
            % filter out POI's with more than 2 mismatches
            POI_Table = barcode_filter(POI_Table);
            % extract poi locations and add fov number at the beginning
            locations = POI_Table(1:3:end, 1:3);
            nPOI = length(locations);
            addingFOV = ones(nPOI,1) * iFov;
            locations = [addingFOV locations];
            % extract intensities and reform it to a vector
            itst = POI_Table(:, 4:end);
            itst = reshape(itst, 3, nPOI, 12);
            itst = permute(itst,[1, 3, 2]);
            itst = reshape(itst, 1, 36, nPOI);
            itst = permute(itst, [3, 2, 1]);
            itst = squeeze(itst);
            POI_Table = [locations, itst];

            if isempty(gathered_Poi_table)
                gathered_Poi_table = POI_Table;
            else
                gathered_Poi_table = [gathered_Poi_table; POI_Table];
            end
        end
    end
end

function gathered_Poi_locations = gather_poi_indices(directory, max_fov)
    gathered_Poi_locations = [];

    format = 'POI_Table_FOV_%d.mat';

    for iFov = 1:max_fov
        file = fullfile(directory, sprintf(format, iFov));
        if exist(file, 'file')
            load(file, 'POI_Table');
            POI_Table = full(POI_Table);
            % filter out POI's with more than 2 mismatches
            POI_Table = barcode_filter(POI_Table);
            % extract poi locations and add fov number at the beginning
            locations = POI_Table(1:3:end, 1:3);
            nPOI = length(locations);
            addingFOV = ones(nPOI,1) * iFov;
            locations = [addingFOV locations];

            if isempty(gathered_Poi_locations)
                gathered_Poi_locations = locations;
            else
                gathered_Poi_locations = [gathered_Poi_locations; locations];
            end
        end
    end
end

function [poi_c, phaseRatios] = extract_and_correct_phase2(poi, graph_)
    nPoi = size(poi,1) / 3;
    phaseRatios = zeros(3,2);

    cycles = cell(4,1);
    cycles{1} = [7, 10];      % G
    cycles{2} = [1, 4, 8, 11];% T
    cycles{3} = [3, 5, 6, 9]; % A
    cycles{4} = [2, 12];      % C

    % collectingCycles = [7, 4, 3]; % G T A
    collectingCycles = cell(3, 2);
    collectingCycles{1,1} = [6, 9];   % G pre
    collectingCycles{1,2} = [8, 11];  % G post
    collectingCycles{2,1} = [3, 7, 10];
    collectingCycles{2,2} = [5, 9, 12];
    collectingCycles{3,1} = [2, 4, 8];
    collectingCycles{3,2} = [4, 7, 10];

    channelIndex = [1, 2, 3];

    %% acquire phase ratios
    % for i = 1:3
    %     mainCycle = median(poi(channelIndex(i):3:end, collectingCycles(i)));
    %     preCycle = median(poi(channelIndex(i):3:end, collectingCycles(i) - 1));
    %     postCycle = median(poi(channelIndex(i):3:end, collectingCycles(i) + 1));
    %     phaseRatios(i, 1) = preCycle / mainCycle;
    %     phaseRatios(i, 2) = postCycle / mainCycle;
    % end

    x = 1:12;
    for i = 1:3
        preCycle = median(poi(channelIndex(i):3:end, collectingCycles{i, 1}), 1) ./ median(poi(channelIndex(i):3:end, collectingCycles{i, 1} + 1), 1);
        postCycle = median(poi(channelIndex(i):3:end, collectingCycles{i, 2}), 1) ./ median(poi(channelIndex(i):3:end, collectingCycles{i, 2} - 1), 1);
        pPre = polyfit(collectingCycles{i, 1}, preCycle, 1);
        pPost = polyfit(collectingCycles{i, 2}, postCycle, 1);
        phaseRatios(i, 1) = pPre(1) / 2;
        phaseRatios(i, 2) = pPost(1) / 2;
    end
    phaseRatios(phaseRatios < 0) = 0;

    %% create phase correction matrix
    phasing_matricies = cell(3,1);
    for i = 1:3
        phasing_matricies{i} = create_phasing_matrix(phaseRatios(i, 1), phaseRatios(i, 2), 12);
    end

    % Correct phasing
    poi_c = poi';
    poi_c = reshape(poi_c, 12, 3, nPoi);
    poi_c = permute(poi_c, [1, 3, 2]);
    for i = 1:3
        poi_c(:,:,i) = phasing_matricies{i} \ poi_c(:,:,i);
    end
    
    poi_c = permute(poi_c, [1, 3, 2]);
    poi_c = reshape(poi_c, 12, 3 * nPoi);
    poi_c = poi_c';

    %% Plot
    if graph_
        plot_intensities(poi);
        % plot_intensities(poi_c);
    end
    %% Collect medians of signal cycles and normalize
    normalizing_itst = zeros(3, 12);
    cyc = 1:12;
    for ich = 1:3
        normalizing_itst(ich, cycles{ich}) = median(poi_c(channelIndex(ich):3:end, cycles{ich}), 1);
        % set mean as signal median to the rest of the cycles
        normalizing_itst(ich, cyc(~ismember(cyc, cycles{ich}))) = mean(normalizing_itst(ich, cycles{ich}));
    end

    poi_c = poi_c';
    poi_c = reshape(poi_c, 12, 3, nPoi);
    for i = 1:nPoi
        poi_c(:,:,i) = poi_c(:,:,i) ./ normalizing_itst';
    end

    poi_c = reshape(poi_c, 12, 3 * nPoi);
    poi_c = poi_c';

    % plot
    if graph_
        plot_intensities(poi_c);
    end
end

function [poi_c, phaseRatios] = extract_and_correct_phase(poi, graph_)
    nPoi = size(poi,1) / 3;
    % phaseRatios = zeros(3,2);
    phaseRatios = zeros(6, 12);

    cycles = cell(4,1);
    cycles{1} = [7, 10];      % G
    cycles{2} = [1, 4, 8, 11];% T
    cycles{3} = [3, 5, 6, 9]; % A
    cycles{4} = [2, 12];      % C

    % collectingCycles = [7, 4, 3]; % G T A
    collectingCycles = cell(3, 2);
    collectingCycles{1,1} = [6, 9];   % G pre
    collectingCycles{1,2} = [8, 11];  % G post
    collectingCycles{2,1} = [3, 7, 10];
    collectingCycles{2,2} = [5, 9, 12];
    collectingCycles{3,1} = [2, 4, 8];
    collectingCycles{3,2} = [4, 7, 10];

    channelIndex = [1, 2, 3];

    %% acquire phase ratios
    % for i = 1:3
    %     mainCycle = median(poi(channelIndex(i):3:end, collectingCycles(i)));
    %     preCycle = median(poi(channelIndex(i):3:end, collectingCycles(i) - 1));
    %     postCycle = median(poi(channelIndex(i):3:end, collectingCycles(i) + 1));
    %     phaseRatios(i, 1) = preCycle / mainCycle;
    %     phaseRatios(i, 2) = postCycle / mainCycle;
    % end

    x = 1:12;
    for i = 1:3
        preCycle = median(poi(channelIndex(i):3:end, collectingCycles{i, 1}), 1) ./ median(poi(channelIndex(i):3:end, collectingCycles{i, 1} + 1), 1);
        postCycle = median(poi(channelIndex(i):3:end, collectingCycles{i, 2}), 1) ./ median(poi(channelIndex(i):3:end, collectingCycles{i, 2} - 1), 1);
        pPre = polyfit(collectingCycles{i, 1}, preCycle, 1);
        pPost = polyfit(collectingCycles{i, 2}, postCycle, 1);
        phaseRatios(i, :) = polyval(pPre, x);
        phaseRatios(i+3, :) = polyval(pPost, x);
    end
    phaseRatios(phaseRatios < 0) = 0;

    %% create phase correction matrix
    % phasing_matricies = cell(3,1);
    % for i = 1:3
    %     phasing_matricies{i} = create_phasing_matrix(phaseRatios(i, 1), phaseRatios(i, 2), 12);
    % end

    % % Correct phasing
    % poi_c = poi';
    % poi_c = reshape(poi_c, 12, 3, nPoi);
    % poi_c = permute(poi_c, [1, 3, 2]);
    % for i = 1:3
    %     poi_c(:,:,i) = phasing_matricies{i} \ poi_c(:,:,i);
    % end
    % 
    % poi_c = permute(poi_c, [1, 3, 2]);
    % poi_c = reshape(poi_c, 12, 3 * nPoi);
    % poi_c = poi_c';

    % Correct phasing
    poi_c = poi;
    for ich = 1:3
        % get intensities (put the first cycle to the front)
        % (itst:intensities)
        itst = poi_c(ich:3:end, :);

        % correct pre phasing of the first cycle
        [itst(:,1), itst(:,2)] = correct_phasing(itst(:,1), itst(:,2), phaseRatios(ich, 1));

        % correct pre and post phasing of all middles cycles
        for icyc = 2:11
            % correct postphasing
            [itst(:,icyc), itst(:, icyc-1)] = correct_phasing(itst(:, icyc), itst(:, icyc-1), phaseRatios(ich+3, icyc-1));
            % correct prephasing
            [itst(:,icyc), itst(:, icyc+1)] = correct_phasing(itst(:,icyc), itst(:,icyc+1), phaseRatios(ich, icyc+1));
        end

        % correct post phasing of the last cycle
        [itst(:, end), itst(:,end-1)] = correct_phasing(itst(:,end), itst(:,end-1), phaseRatios(ich+3, 12));

        % put back the intensities to the poi table
        poi_c(ich:3:end, :) = itst;
    end
    %% Plot
    if graph_
        plot_intensities(poi);
        % plot_intensities(poi_c);
    end
    %% Collect medians of signal cycles and normalize
    normalizing_itst = zeros(3, 12);
    cyc = 1:12;
    for ich = 1:3
        normalizing_itst(ich, cycles{ich}) = median(poi_c(channelIndex(ich):3:end, cycles{ich}), 1);
        % set mean as signal median to the rest of the cycles
        normalizing_itst(ich, cyc(~ismember(cyc, cycles{ich}))) = mean(normalizing_itst(ich, cycles{ich}));
    end

    poi_c = poi_c';
    poi_c = reshape(poi_c, 12, 3, nPoi);
    for i = 1:nPoi
        poi_c(:,:,i) = poi_c(:,:,i) ./ normalizing_itst';
    end

    poi_c = reshape(poi_c, 12, 3 * nPoi);
    poi_c = poi_c';

    % plot
    if graph_
        plot_intensities(poi_c);
    end
end


%%
% % phaseRatios = zeros(3,2);
% phaseRatios = zeros(6, 12);
% 
% cycles = cell(4,1);
% cycles{1} = [7, 10];      % G
% cycles{2} = [1, 4, 8, 11];% T
% cycles{3} = [3, 5, 6, 9]; % A
% cycles{4} = [2, 12];      % C
% 
% % collectingCycles = [7, 4, 3]; % G T A
% collectingCycles = cell(3, 2);
% collectingCycles{1,1} = [6, 9];   % G pre
% collectingCycles{1,2} = [8, 11];  % G post
% collectingCycles{2,1} = [3, 7, 10];
% collectingCycles{2,2} = [5, 9, 12];
% collectingCycles{3,1} = [2, 4, 8];
% collectingCycles{3,2} = [4, 7, 10];
% 
% channelIndex = [1, 2, 3];
% 
% %% acquire phase ratios
% % for i = 1:3
% %     mainCycle = median(poi(channelIndex(i):3:end, collectingCycles(i)));
% %     preCycle = median(poi(channelIndex(i):3:end, collectingCycles(i) - 1));
% %     postCycle = median(poi(channelIndex(i):3:end, collectingCycles(i) + 1));
% %     phaseRatios(i, 1) = preCycle / mainCycle;
% %     phaseRatios(i, 2) = postCycle / mainCycle;
% % end
% 
% x = 1:12;
% for i = 1:3
%     preCycle = median(poi(channelIndex(i):3:end, collectingCycles{i, 1}), 1);
%     postCycle = median(poi(channelIndex(i):3:end, collectingCycles{i, 2}), 1);
%     pPre = polyfit(collectingCycles{i, 1}, preCycle, 1);
%     pPost = polyfit(collectingCycles{i, 2}, postCycle, 1);
%     phaseRatios(i, :) = polyval(pPre, x);
%     phaseRatios(i+3, :) = polyval(pPost, x);
% end
% phaseRatios(phaseRatios < 0) = 0;
% 
% %% create phase correction matrix
% % phasing_matricies = cell(3,1);
% % for i = 1:3
% %     phasing_matricies{i} = create_phasing_matrix(phaseRatios(i, 1), phaseRatios(i, 2), 12);
% % end
% 
% % % Correct phasing
% % poi_c = poi';
% % poi_c = reshape(poi_c, 12, 3, nPoi);
% % poi_c = permute(poi_c, [1, 3, 2]);
% % for i = 1:3
% %     poi_c(:,:,i) = phasing_matricies{i} \ poi_c(:,:,i);
% % end
% % 
% % poi_c = permute(poi_c, [1, 3, 2]);
% % poi_c = reshape(poi_c, 12, 3 * nPoi);
% % poi_c = poi_c';
% 
% % Correct phasing
% poi_c = poi;
% for ich = 1:3
%     % get intensities (put the first cycle to the front)
%     % (itst:intensities)
%     itst = poi_c(ich:3:end, :);
% 
%     % correct pre phasing of the first cycle
%     itst(:,1) = correct_phasing(itst(:,1), itst(:,2), phaseRatios(ich, 1));
% 
%     % correct pre and post phasing of all middles cycles
%     for icyc = 2:11
%         % correct postphasing
%         itst(:,icyc) = correct_phasing(itst(:, icyc), itst(:, icyc-1), phaseRatios(ich+3, icyc-1));
%         % correct prephasing
%         itst(:,icyc) = correct_phasing(itst(:,icyc), itst(:,icyc+1), phaseRatios(ich, icyc+1));
%     end
% 
%     % correct post phasing of the last cycle
%     itst(:, end) = correct_phasing(itst(:,end), itst(:,end-1), phaseRatios(ich+3, 12));
% 
%     % put back the intensities to the poi table
%     poi_c(ich:3:end, :) = itst;
% end
% %% Plot
% plot_intensities(poi);
% % plot_intensities(poi_c);
% 
% %% Collect medians of signal cycles and normalize
% normalizing_itst = zeros(3, 12);
% cyc = 1:12;
% for ich = 1:3
%     normalizing_itst(ich, cycles{ich}) = median(poi(channelIndex(ich):3:end, cycles{ich}), 1);
%     average = mean(normalizing_itst(ich, cycles{ich}));
%     normalizing_itst(ich, cyc(~ismember(cyc, cycles{ich}))) = average;
% end
% 
% poi_c = poi_c';
% poi_c = reshape(poi_c, 12, 3, nPoi);
% for i = 1:nPoi
%     poi_c(:,:,i) = poi_c(:,:,i) ./ normalizing_itst';
% end
% 
% poi_c = reshape(poi_c, 12, 3 * nPoi);
% poi_c = poi_c';
% 
% % plot
% plot_intensities(poi_c);
% 
% %% Clean up
% clearvars -except poi poi_c poi_master out_directory directory cyc cycles max_fov


function final_rate = compile_mutation_rate(mutation_table, fov_list, nPoi)
    count = zeros(4, 12);
    for iFov = fov_list
        count = count + cellfun(@numel, mutation_table{iFov});
    end
    
    cycles = cell(4,1);
    cycles{1} = [7, 10];      % G
    cycles{2} = [1, 4, 8, 11];% T
    cycles{3} = [3, 5, 6, 9]; % A
    cycles{4} = [2, 12];      % C
    
    final_rate = zeros(4);
    for icyc = 1:4
        cyc = cycles{icyc};
        final_rate(icyc, :) = sum(count(:, cyc), 2) / (nPoi * numel(cyc))';
    end
    
    final_rate = final_rate([3,4,1,2], :);
    final_rate = final_rate(:, [3,4,1,2]);
    
    sum_ = sum(final_rate, 2);
    final_rate(find(eye(4))) = 1 - sum_;
end

function [cycle_mutation_rate, by_channel_mutation_rates] = compile_cycle_mutation_rate(mutation_table, fov_list, nPoi, poi_counts)
    by_channel_mutation_rates = zeros(4, 12, length(poi_counts));
    for iFov = fov_list
        by_channel_mutation_rates(:, :, iFov) = cellfun(@numel, mutation_table{iFov}) / poi_counts(iFov);
    end
    count = zeros(4, 12);
    for iFov = fov_list
        count = count + cellfun(@numel, mutation_table{iFov});
    end
    
    nPoi = numel(fov_list) * nPoi;
    cycle_mutation_rate = zeros(4, 12);
    for icyc = 1:12
        cycle_mutation_rate(:, icyc) = sum(count(:, icyc), 2) / nPoi;
    end
end

function compile_results_with_error_bars_8(counts, fovs, file_dir, poi_counts)
    n_groups = floor(length(fovs) / 8);
    count = cell(n_groups);
    nPoi_list = zeros(n_groups, 1);
    for igroup = 0:n_groups - 1
        ifov = fovs(igroup * 8 + 1 : (igroup + 1) * 8);
        icount = zeros(4, 12);
        for i = ifov
            icount = icount + cellfun(@numel, counts{i});
        end
        count{igroup+1} = icount;
        nPoi_list(igroup+1) = sum(poi_counts(ifov));
    end
    
    rates = cell(n_groups, 1);
    for ig = 1:n_groups
        rates{ig} = count{ig} / nPoi_list(ig);
    end
    compile_results_with_error_bars(rates, 1:n_groups, file_dir);
end

function compile_results_with_error_bars(mutation_rates, Fovs, file_dir)
    cycles = cell(4,1);
    cycles{1} = [7, 10];      % G
    cycles{2} = [1, 4, 8, 11];% T
    cycles{3} = [3, 5, 6, 9]; % A
    cycles{4} = [2, 12];      % C
    folder = fullfile(file_dir, 'errorgraphs');
    mkdir(folder);
    
    
    % A
    mut_to = [];
    for ig = Fovs
        mut_to = [mut_to, mutation_rates{ig}(:, cycles{3})];
    end
    sd_a = std(mut_to, 1) / size(mut_to, 2);
    mn_a = mean(mut_to, 1);
    fig = plot_with_error_bars(mn_a([1,2,4]), sd_a([1,2,4]), {'A to G','A to T','A to C'});
    ylim([0, 2.5e-3]);
    title('Mutation from A Cycle');
    saveas(fig, fullfile(folder, 'Mutation from A Cycle.jpg'));
    % C
    mut_to = [];
    for ig = Fovs
        mut_to = [mut_to, mutation_rates{ig}(:, cycles{4})];
    end
    sd_c = std(mut_to, 1) / size(mut_to, 2);
    mn_c = mean(mut_to, 1);
    fig = plot_with_error_bars(mn_c([1,2,3]), sd_c([1,2,3]), {'C to G','C to T','C to A'});
    ylim([0, 2.5e-3]);
    title('Mutation from C Cycle');
    saveas(fig, fullfile(folder, 'Mutation from C Cycle.jpg'));
    % G
    mut_to = [];
    for ig = Fovs
        mut_to = [mut_to, mutation_rates{ig}(:, cycles{1})];
    end
    sd_g = std(mut_to, 1) / size(mut_to, 2);
    mn_g = mean(mut_to, 1);
    fig = plot_with_error_bars(mn_g([2,3,4]), sd_g([2,3,4]), {'G to T','G to A','G to C'});
    ylim([0, 2.5e-3]);
    title('Mutation from G');
    saveas(fig, fullfile(folder, 'Mutation from G Cycle.jpg'));
    % T
    mut_to = [];
    for ig = Fovs
        mut_to = [mut_to, mutation_rates{ig}(:, cycles{2})];
    end
    sd_t = std(mut_to, 1) / size(mut_to, 2);
    mn_t = mean(mut_to, 1);
    fig = plot_with_error_bars(mn_t([1,3,4]), sd_t([1,3,4]), {'C to G','C to T','C to A'});
    ylim([0, 2.5e-3]);
    title('Mutation from T Cycle');
    saveas(fig, fullfile(folder, 'Mutation from T Cycle.jpg'));
end

function fig = plot_with_error_bars(mean_val, std_val, labels)
    fig = figure('visible', 'off');
    errorbar(1:length(mean_val), mean_val, std_val);
    xlim([0, length(mean_val)+1]);
    xticks(0:length(mean_val)+1);
    xticklabels(cat(2, cat(2, {' '}, labels), {' '}));
    grid on;
end














