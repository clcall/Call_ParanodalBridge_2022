ctrlsheaths = [16 17 20 18 27 21 18 21 16 21 11 17 8 8 25 18 8 13 14 12 22 9];
ctrlbridges = [3 2 1 0 5 0 1 1 3 3 4 4 2 0 1 4 1 0 5 0 1 1];
ctrlprop = ctrlbridges./ctrlsheaths;

bontsheaths = [22 14 4 24 12 12 11 10 14 8 9 19 6 8 10 15 14 8 14 15 17 13 13 14 11 8 11];
bontbridges = [0 0 1 4 0 0 0 0 0 0 2 3 0 1 3 0 4 1 0 1 3 0 3 5 3 0 0];
bontprop = bontbridges./bontsheaths;

[ct,cu] = getFigColors;
%% # of sheaths
avg = [mean(ctrlsheaths),mean(bontsheaths)];
sem = [calcSEM(ctrlsheaths),calcSEM(bontsheaths)];
figure
plotSpread({ctrlsheaths,bontsheaths},'distributionMarker','o','distributionColors',{ct,cu});
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 3])
ylim([0 40])
xticklabels({})
figQuality(gcf,gca,[2.1 2.2])

[~,p,ci,stats] = ttest2(ctrlsheaths,bontsheaths);

%% proportion of bridges
avg = [mean(ctrlprop),mean(bontprop)];
sem = [calcSEM(ctrlprop),calcSEM(bontprop)];
figure
plotSpread({ctrlprop,bontprop},'distributionMarker','o','distributionColors',{ct,cu});
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 3])
ylim([0 0.5])
xticklabels({})
figQuality(gcf,gca,[2.1 2.2])

[~,p,ci,stats] = ttest2(ctrlbridges,bontbridges);
