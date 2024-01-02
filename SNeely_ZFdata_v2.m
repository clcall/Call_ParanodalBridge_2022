path = 'D:\GitHubRepos\Call_ParanodalBridge_2022\SarahUpdatedData';
cd(path);
files = dir(path);
files = files(3:end,:);
sheathsSNccrev=[];
for i = 1:size(files,1)
    if contains(files(i).name,'.mat')
        continue
    end
    name = files(i).name;
    varname = matlab.lang.makeValidName(name);
    if exist(fullfile(path,[varname '.mat']),'file')
        load([varname '.mat']);
    else
        xmlstruct = parseXML_SingleCell(name);
        save(fullfile(path,[varname '.mat']),'xmlstruct');
    end
    [names,lengths,colors,chains] = parseData(xmlstruct);
    nonbrglnths = lengths(colors==0);
    anchlnths = lengths(colors==1);
    brglnths = lengths(colors==2);
    brgshthlnths = lengths(colors==3);
    sheathsSNccrev = [sheathsSNccrev; ... 
        (length(nonbrglnths)+length(brgshthlnths)) length(brglnths) mean([nonbrglnths;anchlnths;brgshthlnths]) mean(nonbrglnths) mean(anchlnths) mean(brglnths) mean(brgshthlnths) mean(chains)];
        % 1 - total #
        % 2 - # bridges
        % 3 - avg all sheath lnths
        % 4 - avg nonbrg lnth
        % 5 - avg anchoringlnth
        % 6 - avg bridge lnth
        % 7 - avg bridged sheath lnth
        % 8 - avg chain length
end

%% panel d - # sheath initiations (i.e. num nonbrg + num anch)
totSheaths_bridges = [21; sheathsSNccrev(:,1); 14; 19; 27; 21] - [3; sheathsSNccrev(:,2); 1; 2; 4; 2];
totSheaths_noBridges = [20; 11; 14; 7; 9; 19; 12; 9; 20; 8; 13; 18;...
15.5; 14; 17; 13; 19; 27; 21; 11; 10; 10; 9; 12];

avg = [mean(totSheaths_noBridges); mean(totSheaths_bridges)]
sem = [calcSEM(totSheaths_noBridges,1); calcSEM(totSheaths_bridges,1)]
figure
plotSpread({totSheaths_noBridges,totSheaths_bridges},'distributionMarker','o','distributionColors',{[0.2 0.5 1],[1 0.5 0]})
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 3])
ylim([0 25])
% ylim([0 50])
% xticklabels({'OLs without bridges','OLs with bridges'})
xticklabels({})
% ylabel('Number of sheaths per cell')
figQuality(gcf,gca,[1.5,2]);
[p,tbl,stats] = kruskalwallis([totSheaths_bridges;totSheaths_noBridges],[ones(size(totSheaths_bridges));zeros(size(totSheaths_noBridges))])
% [p,~,stats] = ranksum(totSheaths_bridges,totSheaths_noBridges)
% [~,p,tstat] = ttest2(totSheaths_bridges,totSheaths_noBridges)

%% panel e
nonBridgeSheath_avgLnth = [26.24; sheathsSNccrev(:,4); 47.37; 36.25; 38.61; 27.77];

anchSheath_avgLnth = [22.76; sheathsSNccrev(:,5); 38.49; 24.02; 27.95; 22.08];

bridge_avgLnth = [2.94; sheathsSNccrev(:,6); 3.58; 4.52; 5.83; 7.11];

bridgeSheath_avgLnth = [15.16; sheathsSNccrev(:,7); 24.44; 14.05; 21.13; 23.87];

% INCLUDING BRIDGES THEMSELVES
bridgeChain_avgLnth = [40.85; sheathsSNccrev(:,8); 52.46; 51.19; 57.09; 54.85]; 

avg = [mean(nonBridgeSheath_avgLnth); mean(anchSheath_avgLnth); mean(bridge_avgLnth); mean(bridgeSheath_avgLnth); mean(bridgeChain_avgLnth)]
sem = [calcSEM(nonBridgeSheath_avgLnth,1); calcSEM(anchSheath_avgLnth,1); calcSEM(bridge_avgLnth,1); calcSEM(bridgeSheath_avgLnth,1); calcSEM(bridgeChain_avgLnth,1)]
figure
plotSpread({nonBridgeSheath_avgLnth,anchSheath_avgLnth,bridge_avgLnth,bridgeSheath_avgLnth,bridgeChain_avgLnth},'distributionMarker','o','distributionColors',{[0.2 0.5 1],[52 75 160]./255,[1 0 0],[1 0.5 0],[0.4 0.4 0.4]})
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 5])
% xticklabels({'Non-bridged sheaths','Anch sheaths','Bridges','Bridged sheaths','Bridge chain'})
xticklabels({})
% ylabel('Average length ()')
figQuality(gcf,gca,[2.5,2]);

temp = [nonBridgeSheath_avgLnth; anchSheath_avgLnth; bridge_avgLnth; bridgeSheath_avgLnth; bridgeChain_avgLnth];
idx = [repmat({'nonbrg'},size(nonBridgeSheath_avgLnth)); repmat({'anch'},size(nonBridgeSheath_avgLnth)); repmat({'brg'},size(bridge_avgLnth));...
    repmat({'brgshth'},size(bridgeSheath_avgLnth)); repmat({'chain'},size(bridgeChain_avgLnth))];
[p,tbl,stats] = kruskalwallis(temp,idx)

%% ----Local Function----
function [ traceName, traceLength, traceColor, chainLnths] = parseData(xmlstruct)
%get length of paths list
numPaths = size(xmlstruct.paths,2);
traceName = cell(numPaths,1);
traceLength = NaN(numPaths,1);
traceColor = NaN(numPaths,1);
traceID = NaN(numPaths,1);
traceParentID = NaN(numPaths,1);
for i = 1:numPaths
    traceName{i,1} = xmlstruct.paths(i).attribs.name;
    traceLength(i,1) = str2double(xmlstruct.paths(i).attribs.reallength); %not using smoothed because highly sparse nodes in traces
    if contains(xmlstruct.paths(i).attribs.color,'Blue')
        traceColor(i,1) = 1; % meaning anchor
    elseif contains(xmlstruct.paths(i).attribs.color,'Red')
        traceColor(i,1) = 2; % bridge
    elseif contains(xmlstruct.paths(i).attribs.color,'Orange')
        traceColor(i,1) = 3; % bridged sheath
    else        
        traceColor(i,1) = 0;
    end
    traceID(i,1) = str2double(xmlstruct.paths(i).attribs.id);
    if isfield(xmlstruct.paths(i).attribs,'startson')
        traceParentID(i,1) = str2double(xmlstruct.paths(i).attribs.startson);
    end
end
chainLnths = [];
usedIDs = [];
for i = 1:numPaths
    if sum(ismember(usedIDs,i))
        continue
    end
    temp = [];
    usedIDs = [usedIDs; i];
    idx = ismember(traceParentID,i);
    if sum(idx)>0
        temp = traceLength(i);
        children = traceID(idx);
        usedIDs = [usedIDs; children];
        while ~isempty(children)
            idx2 = ismember(traceParentID,children);
            temp = [temp; traceLength(idx2)];
            children = traceID(idx2);
            usedIDs = [usedIDs; children];
        end
        chainLnths = [chainLnths; sum(temp)];
    end
end
end