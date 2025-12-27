t = readtable('D:\GitHubRepos\Call_ParanodalBridge_2022\PinkyNodeDiameters.csv');
figure;
hold on;
xes = 1:68;
% non-bridged, non-branched
idx1 = t.bridgeStatus==0 & t.branchStatus==0;
plot(xes(idx1),t.diameter(idx1),'ko');
% non-bridged, branched
idx2 = t.bridgeStatus==0 & t.branchStatus==1;
plot(xes(idx2),t.diameter(idx2),'ro');
% bridged, non-branched
idx3 = t.bridgeStatus==1 & t.branchStatus==0;
plot(xes(idx3),t.diameter(idx3),'ksq','MarkerFaceColor','k');
% bridged, branched
idx4 = t.bridgeStatus==1 & t.branchStatus==1;
plot(xes(idx4),t.diameter(idx4),'rsq','MarkerFaceColor','r');
hold off

xticks()
ylim([0 1500])
figQuality(gcf,gca,[4 2])

[h,p,ks2stat] = kstest2(t.diameter(logical(idx1+idx2)), t.diameter(logical(idx3+idx4)))

figure
hold on
histogram(t.diameter(logical(idx1+idx2)));
histogram(t.diameter(logical(idx3+idx4)));
hold off
figQuality(gcf,gca,[4 2])
legend('non-bridged','bridged')
