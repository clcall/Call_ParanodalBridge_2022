ctrl = [0 0 1 0 0 1 0 1 0 2 0 1 1 0 0 0 2 2 1 0 1 0 0 0 1 5 0 1 1 0 2 3 1 0 0 0 2 0 0];
bont = [1 0 0 0 2 0 1 3 0 3 3 0 2 0 0 0 1 3 0 0 2 0 1 0 0 1 0 0 0 0 0 0 0];

avg = [mean(ctrl),mean(bont)];
sem = [calcSEM(ctrl),calcSEM(bont)];
[ct,cu] = getFigColors;
figure
plotSpread({ctrl,bont},'distributionMarker','o','distributionColors',{ct,cu});
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 3])
ylim([0 100])
xticklabels({})
figQuality(gcf,gca,[2.1 2.2])