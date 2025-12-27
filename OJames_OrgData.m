addpath(genpath('D:\GitHubRepos\Call_ParanodalBridge_2022'));
prop_tbl = readtable('D:\GitHubRepos\Call_ParanodalBridge_2022\OwenUpdatedData\Proportion_cells_with_bridges_updated.xlsx');
prop_tbl = prop_tbl(1:24,:);

nonbridged = prop_tbl.Cells_with_0;
bridged = prop_tbl.Cells_with_1 + prop_tbl.Cells_with_2 + prop_tbl.Cells_with_3;

prop = bridged ./ (nonbridged + bridged);

%% sheaths per cell
data = readtable('D:\GitHubRepos\Call_ParanodalBridge_2022\OwenUpdatedData\iPSC_myelinoid_quantification_sheathlengths_aggregatedpercell.csv');
brg_idx = contains(data.Cell_type,'OLs with bridges');
nonbrg_idx = contains(data.Cell_type,'OLs without bridges');

sheathsPerCell_brg = table2array(data(brg_idx,10));
sheathsPerCell_nonbrg = table2array(data(nonbrg_idx,10));
avg = [mean(sheathsPerCell_nonbrg), mean(sheathsPerCell_brg)];
sem = [calcSEM(sheathsPerCell_nonbrg,1), calcSEM(sheathsPerCell_brg,1)];

[ct,cu] = getFigColors;

figure
plotSpread({sheathsPerCell_nonbrg,sheathsPerCell_brg},'distributionMarker','o','distributionColors',{ct,cu});
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 3])
xticklabels({})
ylim([0 25])
figQuality(gcf,gca,[2.4 2.2])

%% sheath length per cell
sheathLnth_brg = table2array(data(brg_idx,8));
sheathLnth_nonbrg = table2array(data(nonbrg_idx,8));
avg = [mean(sheathLnth_nonbrg), mean(sheathLnth_brg)];
sem = [calcSEM(sheathLnth_nonbrg,1), calcSEM(sheathLnth_brg,1)];

figure
plotSpread({sheathLnth_nonbrg,sheathLnth_brg},'distributionMarker','o','distributionColors',{ct,cu});
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 3])
xticklabels({})
ylim([0 200])
figQuality(gcf,gca,[2.4 2.2])

%% length per cell per sheath type
data = readtable('D:\GitHubRepos\Call_ParanodalBridge_2022\OwenUpdatedData\iPSC_myelinoid_quantification_sheathlengths.csv');

nonbrg_idx = contains(data.Sheath_type,'Regular');
anchr_idx = contains(data.Sheath_type,'Anchored');
brg_idx = contains(data.Sheath_type,'Bridged');

anchr_data = data(anchr_idx,:);
CellIDs = anchr_data.CellID;
uniqCells = unique(anchr_data.CellID,'stable');
mean_anchrLnths = NaN(size(uniqCells));
j = 1;
k = 1;
while j <= length(CellIDs)
    templnths = [];
    while contains(CellIDs{j},uniqCells{k})
        templnths = [templnths; anchr_data.Sheath_length(j)];
%         fprintf([CellIDs{j},'\n',uniqCells{k},'\n'])
        j = j+1;
        if j > length(CellIDs)
%             fprintf('got to here\n')
            break
        end
    end
    mean_anchrLnths(k) = mean(templnths);
    k = k+1;
end
splitData = cellfun(@(x) strsplit(x, '.'), uniqCells, 'UniformOutput', false);
splitAnchrs = vertcat(splitData{:});

nonbrg_data = data(nonbrg_idx,:);
CellIDs = nonbrg_data.CellID;
uniqCells = unique(nonbrg_data.CellID,'stable');
mean_nonbrgLnths = NaN(size(uniqCells));
j = 1;
k = 1;
while j <= length(CellIDs)
    templnths = [];
    while contains(CellIDs{j},uniqCells{k})
        templnths = [templnths; nonbrg_data.Sheath_length(j)];
%         fprintf([CellIDs{j},'\n',uniqCells{k},'\n'])
        j = j+1;
        if j > length(CellIDs)
%             fprintf('got to here\n')
            break
        end
    end
    mean_nonbrgLnths(k) = mean(templnths);
    k = k+1;
end
splitData = cellfun(@(x) strsplit(x, '.'), uniqCells, 'UniformOutput', false);
splitNonbrgs = vertcat(splitData{:});

brg_data = data(brg_idx,:);
CellIDs = brg_data.CellID;
uniqCells = unique(brg_data.CellID,'stable');
mean_brgLnths = NaN(size(uniqCells));
j = 1;
k = 1;
while j <= length(CellIDs)
    templnths = [];
    while contains(CellIDs{j},uniqCells{k})
        templnths = [templnths; brg_data.Sheath_length(j)];
%         fprintf([CellIDs{j},'\n',uniqCells{k},'\n'])
        j = j+1;
        if j > length(CellIDs)
%             fprintf('got to here\n')
            break
        end
    end
    mean_brgLnths(k) = mean(templnths);
    k = k+1;
end
splitData = cellfun(@(x) strsplit(x, '.'), uniqCells, 'UniformOutput', false);
splitBrgs = vertcat(splitData{:});

data2 = readtable('D:\GitHubRepos\Call_ParanodalBridge_2022\OwenUpdatedData\iPSC_myelinoid_quantification_sheathlengths_reorganised_for_Chainlengths.csv');

chain_idx = ~isnan(data2.Chain_length);
data2_chains = data2(chain_idx,:);
uniqCells = unique(data2_chains.CellID,'stable');
CellIDs = data2_chains.CellID;
meanCellChains = NaN(size(uniqCells));
j = 1;
k = 1;
while j <= length(CellIDs)
    templnths = [];
    while contains(CellIDs{j},uniqCells{k})
        templnths = [templnths; data2_chains.Chain_length(j)];
%         fprintf([CellIDs{j},'\n',uniqCells{k},'\n'])
        j = j+1;
        if j > length(CellIDs)
%             fprintf('got to here\n')
            break
        end
    end
    meanCellChains(k) = mean(templnths);
    k = k+1;
end
splitData = cellfun(@(x) strsplit(x, '.'), uniqCells, 'UniformOutput', false);
splitChains = vertcat(splitData{:});

lengths = num2cell([mean_brgLnths;mean_nonbrgLnths;mean_anchrLnths;meanCellChains]);
types = [repmat({'brg'},[length(mean_brgLnths),1]);...
        repmat({'nonbrg'},[length(mean_nonbrgLnths),1]);...
        repmat({'anch'},[length(mean_anchrLnths),1]);...
        repmat({'chain'},[length(meanCellChains),1])];
ids = [splitBrgs;splitNonbrgs;splitAnchrs;splitChains];
concatdata = [lengths types ids];
    
T = cell2table(concatdata, 'VariableNames', {'Length', 'Type', 'Cells', 'Conversion', 'Organoid', 'Cell'});
mdl = fitlme(T,'Length ~ Type + (1|Cells) + (1|Conversion) + (1|Organoid) + (1|Cell)')

sheathLnth_brg = mean_brgLnths;
sheathLnth_nonbrg = mean_nonbrgLnths;
sheathLnth_anchr = mean_anchrLnths;
sheathLnth_chain = meanCellChains;

avg = [mean(sheathLnth_nonbrg), mean(sheathLnth_anchr), mean(sheathLnth_brg), mean(sheathLnth_chain)];
sem = [calcSEM(sheathLnth_nonbrg,1), calcSEM(sheathLnth_anchr,1), calcSEM(sheathLnth_brg,1), calcSEM(sheathLnth_chain,1)];

figure
plotSpread({sheathLnth_nonbrg,sheathLnth_anchr,sheathLnth_brg,sheathLnth_chain},'distributionMarker','o','distributionColors',{ct,[52 75 160]./255,cu,[0.5 0.5 0.5]});
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 5])
xticklabels({})
ylim([0 400])
figQuality(gcf,gca,[2.5 2.2])