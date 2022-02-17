L1 = [bsln;ctrl];
[~,pvcol] = CTSMcolors;
%% Proportion of bridges
y = {L1.numbrgs./(L1.numbrgs + L1.numNonbrg).*100,l23.numbrgs./(l23.numbrgs + l23.numNonbrg).*100};

y2 = [L1.numbrgs./(L1.numbrgs + L1.numNonbrg).*100; l23.numbrgs./(l23.numbrgs + l23.numNonbrg).*100];

[p,tbl,stats] = kruskalwallis(y2,[repmat({'Layer 1'},size(L1.numbrgs)); repmat({'Layer 2/3'},size(l23.numbrgs))]);

avg = reshape(cellfun(@(y) mean(y,'omitnan'),y),1,2)
sem = reshape(cellfun(@(y) calcSEM(y,1),y),1,2)
sz = cellfun(@length,y);
idx = arrayfun(@(x) ones(x,1), sz, 'UniformOutput', false);
idx2 = [];
for i = 1:length(idx)
    idx2 = [idx2; idx{i} .* i];
end
figure

plotSpread(y2,'distributionIdx',idx2,'distributionMarker','o','distributionColors',{[0.5 0.5 0.5],pvcol})
hold on
errorbar(avg,sem,'ko','MarkerSize',4,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 3])
ylim([0 25])
ax = gca();
xticklabels({'Layer 1' 'Layer 2/3'})
xticklabels({})
figQuality(gcf,gca,[1.8 2.2])
% ytickformat('%.2f')

%% Cell tot length VS bridge number
figure; plot(L1.numbrgs,L1.cellTotLength./1000,'o','MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(l23.numbrgs,l23.cellTotLength./1000,'o','MarkerEdgeColor',pvcol,'LineWidth',1);hold off;
figQuality(gcf,gca,[2.8 2.8])
ylim([0 6])
xlim([0 12])

[R,P,RL,RU] = corrcoef([[L1.numbrgs] [L1.cellTotLength]])

%% Average bridge length VS average sheath length
colors2 = {[0.5 0.5 0.5],[0.5 0.5 0.5],pvcol,pvcol};
xes2 = repmat(1:2,2,1);

y = {L1.avgbrglength, L1.avgNonbrglength,...
     l23.avgbrglength, l23.avgNonbrglength};

plotEmUp(y,xes2,colors2)
xticklabels({'Layer 1' 'Layer 2/3'})
figQuality(gcf,gca,[2.8 2.2])

y2 = [L1.avgbrglength; L1.avgNonbrglength;...
     l23.avgbrglength; l23.avgNonbrglength];
idx = [repmat({'L1_b'},size(l23.cellname)); repmat({'L1_nb'},size(l23.cellname));...
       repmat({'l23_b'},size(L1.cellname)); repmat({'l23_nb'},size(L1.cellname))];
[p,tbl,stats] = kruskalwallis(y2,idx);
multcompare(stats,"CType","dunn-sidak");

%% Average chain length VS average sheath length
y = {L1.avgChainLengthSum, L1.avgNonbrglength,...
    l23.avgChainLengthSum, l23.avgNonbrglength};

plotEmUp(y,xes2,colors2)
xticklabels({'Layer 1' 'Layer 2/3'})
figQuality(gcf,gca,[2.8 2.2])

y2 = [L1.avgChainLengthSum; L1.avgNonbrglength;...
     l23.avgChainLengthSum; l23.avgNonbrglength];
[~,p,tstat] = ttest(L1.avgChainLengthSum,L1.avgNonbrglength);
p*2
[~,p,tstat] = ttest(l23.avgChainLengthSum,l23.avgNonbrglength);
p*2

idx = [repmat({'l1_b'},size(L1.cellname)); repmat({'l1_nb'},size(L1.cellname));...
       repmat({'L23_b'},size(l23.cellname)); repmat({'L23_nb'},size(l23.cellname))];
[p,tbl,stats] = kruskalwallis(y2,idx);
multcompare(stats,"CType","dunn-sidak");