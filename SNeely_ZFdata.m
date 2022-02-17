%% panel C
totSheaths_bridges = [20
13
27
18
21
21
21
27
14
19
27
21
12];

totSheaths_noBridges = [20
11
14
14
7
12
15
9
22
19
12
28
9
20
8
13
14
18
15.5
15.5
21
17
13
23
11
10
10
9
12];

avg = [mean(totSheaths_noBridges); mean(totSheaths_bridges)]
sem = [calcSEM(totSheaths_noBridges,1); calcSEM(totSheaths_bridges,1)]
figure
plotSpread({totSheaths_noBridges,totSheaths_bridges},'distributionMarker','o','distributionColors',{[0.2 0.5 1],[1 0.5 0]})
hold on
errorbar(avg,sem,'ko','MarkerSize',4,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 3])
ylim([0 25])
% ylim([0 50])
% xticklabels({'OLs without bridges','OLs with bridges'})
xticklabels({})
% ylabel('Number of sheaths per cell')
figQuality(gcf,gca,[1.5,2]);
% [p,tbl,stats] = kruskalwallis([totSheaths_bridges;totSheaths_noBridges],[ones(size(totSheaths_bridges));zeros(size(totSheaths_noBridges))])
% [p,~,stats] = ranksum(totSheaths_bridges,totSheaths_noBridges)
[~,p,tstat] = ttest2(totSheaths_bridges,totSheaths_noBridges)
%% panel D
avgSheathLnth_bridges = [15.9212
26.82738462
19.15692593
28.5345
24.16004762
32.05380952
28.64233333
23.1854074
44.0912857
32.7433684
26.6674074
27.0308095
37.59975];

avgSheathLnth_noBridges = [14.9379
30.81409091
22.42392857
33.41207143
29.15857143
34.19541667
34.67286667
28.756
29.88945455
22.42542105
20.62641667
25.13257143
37.46166667
28.36675
41.433125
36.55338462
41.20992857
34.333944
20.071806
20.071806
18.902524
40.070588
38.876769
26.732435
34.271091
49.7854
44.5586
43.179778
49.279167];

avg = [mean(avgSheathLnth_noBridges); mean(avgSheathLnth_bridges)]
sem = [calcSEM(avgSheathLnth_noBridges,1); calcSEM(avgSheathLnth_bridges,1)]
figure
plotSpread({avgSheathLnth_noBridges,avgSheathLnth_bridges},'distributionMarker','o','distributionColors',{[0.2 0.5 1],[1 0.5 0]})
hold on
errorbar(avg,sem,'ko','MarkerSize',3.8,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 3])
ylim([0 60])
% xticklabels({'OLs without bridges','OLs with bridges'})
xticklabels({})
% ylabel('Average sheath length ()')
figQuality(gcf,gca,[1.5,2]);
% [p,tbl,stats] = kruskalwallis([avgSheathLnth_bridges;avgSheathLnth_noBridges],[ones(size(avgSheathLnth_bridges));zeros(size(avgSheathLnth_noBridges))])
[~,p,tstat] = ttest2(avgSheathLnth_bridges,avgSheathLnth_noBridges)
%% panel E
nonBridgeSheath_avgLnth = [17.1375
27.30181818
19.4794
26.95475
26.2408
34.7955625
29.51731579
24.123455
47.366667
36.2485
28.6051
27.774529
39.5061];

bridgeSheath_avgLnth = [4.9745
24.218
15.126
41.1725
18.95816667
23.2802
20.33
15.8816667
24.439
14.0493333
21.1311429
23.87
28.068];

bridge_avgLnth = [7.722
7.306
1.94
2.986
2.936666667
3.942333333
10.174
6.505
3.58
4.52
5.83425
7.114
2.12];

bridgeChain_avgLnth = [17.671
55.742
32.192
85.331
40.853
64.114
50.834
57.4025
52.458
51.188
57.085
54.854
58.256];

avg = [mean(nonBridgeSheath_avgLnth); mean(bridgeSheath_avgLnth); mean(bridge_avgLnth); mean(bridgeChain_avgLnth)]
sem = [calcSEM(nonBridgeSheath_avgLnth,1); calcSEM(bridgeSheath_avgLnth,1); calcSEM(bridge_avgLnth,1); calcSEM(bridgeChain_avgLnth,1)]
figure
plotSpread({nonBridgeSheath_avgLnth,bridgeSheath_avgLnth,bridge_avgLnth,bridgeChain_avgLnth},'distributionMarker','o','distributionColors',{[0.2 0.5 1],[1 0.5 0],[1 0 0],[0.4 0.4 0.4]})
hold on
errorbar(avg,sem,'ko','MarkerSize',5,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 5])
% xticklabels({'Non-bridged sheaths','Bridged sheaths','Bridge','Bridge chain'})
xticklabels({})
% ylabel('Average length ()')
figQuality(gcf,gca,[2.5,2]);

temp = [nonBridgeSheath_avgLnth; bridgeSheath_avgLnth; bridge_avgLnth; bridgeChain_avgLnth];
idx = [repmat({'nonbrg'},size(nonBridgeSheath_avgLnth)); repmat({'brgshth'},size(bridgeSheath_avgLnth));...
    repmat({'brg'},size(bridge_avgLnth)); repmat({'chain'},size(bridgeChain_avgLnth))];
[p,tbl,stats] = kruskalwallis(temp,idx)
%% panel g 
% NOTE DO NOT INCLUDE F
anchorSheath_avgLnth = [7.553
45.23
25.853
80.629
44.231
40.568
28.714
5.078
37.669
18.962
20.176
11.9
32.4
38.49
24.024
26
6.418
24.024
40.513
25.141
18.199
40.472
11.186
32.971
52.845];

bridgedSheath_avgLnth = [2.396
3.206
4.399
1.716
7.344
11.554
11.946
13.188
6.434
24.946
22.6
2.39
10.388
11.706
7.969
12.939
2.685
21.722
29.601
3.291];

avg = [mean(anchorSheath_avgLnth); mean(bridgedSheath_avgLnth)]
sem = [calcSEM(anchorSheath_avgLnth,1); calcSEM(bridgedSheath_avgLnth,1)]
figure
plotSpread({anchorSheath_avgLnth,bridgedSheath_avgLnth},'distributionMarker','o','distributionColors',{[0.2 0.5 1],[1 0.5 0]})
hold on
errorbar(avg,sem,'ko','MarkerSize',5,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 3])
% xticklabels({'Anchoring sheaths','Bridged sheaths'})
xticklabels({})
% ylabel('Average sheath length ()')
figQuality(gcf,gca,[1.5,2]);
[~,p,tstat] = ttest2(anchorSheath_avgLnth,bridgedSheath_avgLnth)