% NOTE: requires plotbridges_v3.m to be run first to acquire single cell
% data

maindir = 'D:\GitHubRepos\OligodendrocyteAnalysisCode\ParanodalBridges\';
addpath(genpath(maindir));
directory = 'D:\GitHubRepos\OligodendrocyteAnalysisCode\ParanodalBridges\labeledBslnTraces\';
files = dir(directory);
files(1:2) = [];
PV_PNBs = struct;
for k = 1:length(files)
    if contains(files(k).name,'.mat')
        continue
    end
    varname = matlab.lang.makeValidName(files(k).name);
    currXmlstruct = fullfile(files(1).folder,[varname '.mat']);
    if contains(varname,'F_191218w_1')
        branchPts = [7;8];
        branches = [13;25];
    elseif contains(varname,'M_191213w_1')
        branchPts = [4;5;2];
        branches = [3;12;6];
    elseif contains(varname,'F_200129w_4')
        branchPts = [7;6];
        branches = [15;11];
    elseif contains(varname,'F_200129w_5')
        branchPts = [4;4];
        branches = [10;3];
    elseif contains(varname,'M_200129w_1')
        branchPts = [1;1];
        branches = [3;3];
    elseif contains(varname,'M_200129w_2')
        branchPts = [3;6];
        branches = [7;14];
    elseif contains(varname,'M_200129w_3')
        branchPts = [4;7];
        branches = [7;14];
    elseif contains(varname,'M_200129w_6')
        branchPts = [2;1;3];
        branches = [7;3;7];
    elseif contains(varname,'M_200129w_7')
        branchPts = [3;4];
        branches = [6;13];
    end
    tbl = table(branchPts,branches);
    if exist(currXmlstruct,'file')
        load(currXmlstruct)
        PV_PNBs.(varname) = calculatePathsXML_PVpnb(xmlstruct,1,0);
        PV_PNBs.(varname) = [PV_PNBs.(varname) tbl];
    else
        tracesFile = fullfile(directory,files(k).name);                
        fprintf('Parsing traces file %s...\n',num2str(k));
        xmlstruct = parseXML_SingleCell(tracesFile);
        save([directory varname '.mat'],'xmlstruct');
        PV_PNBs.(varname) = calculatePathsXML_PVpnb(xmlstruct,1,0);
        PV_PNBs.(varname) = [PV_PNBs.(varname) tbl];
        fprintf('Done.\n');
    end
%     pause;
end
save([directory 'PV_PNBs.mat'],'PV_PNBs');
%% combine all data
fields = fieldnames(PV_PNBs);
alldata = [];
for f = 1:length(fields)
    alldata = [alldata; PV_PNBs.(fields{f})];
end
%% make figures
data = [alldata.avgBrgSheathLnth,alldata.avgNonBrgSheathLnth];
avgData = mean(data(:,1:2));
sem = calcSEM(data(:,1:2),1);
[~,pvcol] = CTSMcolors;
%% plot average lengths of anchor vs bridged per axon
figure
color = [0.5 0.5 0.5];
x = [1 2];
hold on
for i = 1:size(data,1)
    plot(x,data(i,1:2),'o-','Color',pvcol,'MarkerSize',5,'LineWidth',0.5);
end
xlim([0.5,2.5]);
ylim([0,80]);
xticks([1,2]);
errorbar([1,2],avgData,sem,'ko','MarkerSize',4,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',0)
hold off
set(gca,'XTickLabel',[])
figQuality(gcf,gca,[1.69,2.411]);

[~,p,tstat] = ttest(data(:,1),data(:,2))
%% plot correlation between branch pts and bridge #
figure; 
mdl = fitlm(alldata.branchPts,alldata.numBridges);
h = plot(mdl);
h(1).Marker = 'o';
h(1).MarkerSize = 5;
h(1).Color = 'k';
xlabel('Branch points');
ylabel('Number of bridges');
title('')
xlim([0 10]);
ylim([0 10]);
figQuality(gcf,gca,[3,2.7]);
[rho,pval] = corr(alldata.branchPts,alldata.numBridges)
mdl.Rsquared
%% plot correlation between total myelin and bridge #
figure; 
mdl = fitlm(alldata.totalSheathLength/1000,alldata.numBridges);
h = plot(mdl);
h(1).Marker = 'o';
h(1).MarkerSize = 5;
h(1).Color = 'k';
xlabel('Total myelin length ()');
ylabel('Number of bridges');
title('')
xlim([0 1.5]);
ylim([0 10]);
figQuality(gcf,gca,[3,2.7]);
[rho,pval] = corr(alldata.totalSheathLength,alldata.numBridges)
mdl.Rsquared
%% plot proportion of PV axon bridges vs. generation
%combine bsln and ctrl data
nonpv = [bsln.numbrgs./(bsln.numbrgs + bsln.numNonbrg);...
        ctrl.numbrgs./(ctrl.numbrgs + ctrl.numNonbrg)] .*100;
L23 = l23.numbrgs./(l23.numbrgs + l23.numNonbrg) .* 100;
pv = alldata.pnbSheathRatio .* 100;
plotit(nonpv,L23,pv,'Proportion of sheaths with bridges',0);

[~,p1_23,tstat] = ttest2(nonpv,L23);
[~,pol_pv,tstat] = ttest2(L23,pv);
p1_23*2
pol_pv*2

%% plot non-bridge sheath length data
L1 = [bsln.avgNonbrglength; ctrl.avgNonbrglength];
L23 = l23.avgNonbrglength;
pv = alldata.avgNonBrgSheathLnth;
plotit(L1,L23,pv,'Average non-bridged sheath length',0);
ylim([0 100])

[~,p1_23,tstat] = ttest2(L1,L23);
[~,pol_pv,tstat] = ttest2(L23,pv);
p1_23*2
pol_pv*2
%% plot bridged sheath length data
L1 = [bsln.avgbrglength; ctrl.avgbrglength];
L23 = l23.avgbrglength;
pv = alldata.avgBrgSheathLnth;
plotit(L1,L23,pv,'Average bridged sheath length',0);
ylim([0 120])

[~,p1_23,tstat] = ttest2(L1,L23);
[~,pol_pv,tstat] = ttest2(L23,pv);
p1_23*2
pol_pv*2
%% plot bridge chain length data
L1 = [bsln.avgChainLengthSum; ctrl.avgChainLengthSum];
L23 = l23.avgChainLengthSum;
pv = alldata.avgChainLnth;
plotit(L1,L23,pv,'Average bridge chain length',0);
ylim([0 200])

[~,p1_23,tstat] = ttest2(L1,L23);
[~,pol_pv,tstat] = ttest2(L23,pv);
p1_23*2
pol_pv*2
%% fractions of bridged sheaths that ensheath PV axons (L23 cells)
f = [3 5; 2 5; 0 7; 3 5; 1 1; 4 7; 4 7; 2 3; 6 6; 2 5; 5 7; 8 8; 2 6; 10 11; 6 6];

function [p,tbl,stats] = plotit(nonpv,L23,pv,y_label,do_stats)
[~,pvcol] = CTSMcolors;
avg = [mean(nonpv,'omitnan'),mean(L23),mean(pv)];
sem = [calcSEM(nonpv,1),calcSEM(L23,1),calcSEM(pv,1)];

figure
plotSpread({nonpv,L23,pv},'distributionMarker','o','distributionColors',{'b',[0.5 0.5 0.5],pvcol})
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 4])
xticklabels({'L1 Single OL', 'L2/3 Single OL','Single PV axon'})
ylabel(y_label)
figQuality(gcf,gca,[1.8,2.4]);
% yl = ylim;
% ylim([0 yl(2)])
% yt = yticks;
% ylim([0 (yl(2) + (yt(2) - yt(1)))])
if do_stats
    [p,tbl,stats] = kruskalwallis([nonpv;L23;pv],[repmat({'L1'},size(nonpv));repmat({'L23'},size(L23));repmat({'pv'},size(pv))])
end
end