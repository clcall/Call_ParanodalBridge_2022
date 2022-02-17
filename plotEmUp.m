function plotEmUp(y,xes,colors)
avg = reshape(cellfun(@(y) mean(y,'omitnan'),y),[],2)
sem = reshape(cellfun(@(y) calcSEM(y,1),y),[],2)
figure
[b,e] = errorbarbar(xes',avg',sem',{'EdgeColor','none','FaceColor','none'},{'k','LineWidth',1.5,'LineStyle','none','MarkerSize',7,'MarkerFaceColor','k'});
e(1).Marker = 'o';
e(1).CapSize = 0;
e(1).MarkerFaceColor = 'k';
e(1).MarkerSize = 4;

e(2).Marker = 's';
e(2).MarkerFaceColor = 'k';
e(2).MarkerSize = 4.5;
e(2).CapSize = 0;
sz = cellfun(@length,y);
xoff = [get(b(1),'XOffset') get(b(2),'XOffset')];
xoff = repmat(xoff,1,3);
idx = arrayfun(@(x) ones(x,1), sz, 'UniformOutput', false);
idx2 = [];
for i = 1:length(idx)
    idx2 = [idx2; idx{i} .* xoff(i) + ceil(i./2)];
end
hold on
plotSpread(y,'distributionIdx',idx2,'distributionMarkers',repmat({'o','s'},1,2),'distributionColors',colors,'binWidth',1.1);
xticks(xes(1,:))
% row1 = {'< 8wks' '> 8wks' '> 8wks'};
% row2 = {' ' 'Control' 'Cuprizone'};
% labelArray = [row1; row2];
% labelArray = strjust(pad(labelArray),'center'); % 'left'(default)|'right'|'center
% ticklabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels({'Baseline','Control','Cuprizone'})
ax = gca();
% ax.XTickLabel = ticklabels;
hold off
end