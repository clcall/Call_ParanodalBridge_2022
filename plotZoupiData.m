function plotZoupiData
a = [7	16; 3	8; 9	10];

b = [3	8;
2	7;
5	15;
1	7;
3	6;
3	6;
1	5;
3	6;
4	7;
1	5;
1	4;
4	7];

c = [10	8;
5	9;
4	5;
3	5;
7	14;
8	23];

d = [2	3;
4	8;
2	5;
3	3;
0	3;
1	3;
2	5;
1	3;
2	2;
4	5;
2	3;
3	9;
2	2;	
1	1;
2	15];

data = [mean(a,1); mean(b,1); mean(c,1); mean(d,1)];
sums = sum(data,2);
data2 = (data(:,1)./sums).*100;
avg = mean(data2);
sem = calcSEM(data2,1);
figure
plotSpread(data2,'distributionMarker','o','distributionColors',[0.5 0.5 0.5])
hold on
errorbar(avg,sem,'ko','MarkerSize',3,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold off
xlim([0 2])
ylim([0 50])
figQuality(gcf,gca,[1.69,2.411]);
xticklabels([])
end