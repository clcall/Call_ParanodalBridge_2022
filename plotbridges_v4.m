cond = fieldnames(ALL);
cond = cond(2:end);
data = struct;
for c = 1:length(cond)
    an = fieldnames(ALL.(cond{c}));
    data.(cond{c}) = struct;
    vars = {'cellname','numbrgs','numNonbrg','avgbrglength','totbrglength'...
            'avgNonbrglength','avganchLength','avgChainLengthSum','cellTotLength',...
            'avgChainNodeRatio','avgNonBrgNodeRatio','xy','z','xy_nb','z_nb'};
    vartypes = ['cellstr', repmat({'double'},1,length(vars)-1)];
    data.(cond{c}).datatable = table('Size',[1, length(vars)],'VariableTypes',vartypes,'VariableNames',vars);
    for a = 1:length(an)
        day = fieldnames(ALL.(cond{c}).(an{a}).raw);
        tempdata14 = ALL.(cond{c}).(an{a}).raw.(day{end});
        if contains(an{a},{'a129' 'a213'}) && contains(cond{c},{'ctrl','L2_3'})
            tempdata00 = tempdata14;
        elseif contains(cond{c},'bsln') % work-around for d00 being empty in bsln struct
            tempdata00 = tempdata14;
        else
            tempdata00 = ALL.(cond{c}).(an{a}).raw.d00;
        end
        brgIdx = tempdata14.Sheath_Info.isbridgecnctd;
        tempdata14.Sheath_Info(isnan(brgIdx),:) = [];
        
        % convert node class to # neighbors
        nodes = tempdata14.Sheath_Info.node_class;
        nodes(nodes==1) = 0;
        nodes(nodes==2) = 1;
        nodes(nodes==3) = 0;
        nodes(nodes==4) = 1;
        nodes(nodes==5) = 2;
        tempdata14.Sheath_Info.node_class = nodes;
        
        brgIdx(isnan(brgIdx)) = [];
        brgIdx = logical(brgIdx);
        brgdata14 = tempdata14.Sheath_Info(brgIdx,:);
        % find original path ID connections at d00
        d00idx = contains(tempdata00.Sheath_Info.intNames, brgdata14.intNames);
        brgdata00 = tempdata00.Sheath_Info(d00idx,:);
        if isempty(brgdata00)
            nonbrgdata = tempdata14.Sheath_Info((~brgIdx & ~contains(tempdata14.Sheath_Info.intNames, anchNames)),:);
            numNonbrg = length(nonbrgdata.intNames);
            cellTotLength = sum(tempdata14.Sheath_Info.intLengths);
            avgNonbrglength = mean(nonbrgdata.intLengths);
            avgNonBrgNodeRatio = mean(nonbrgdata.node_class,'omitnan');
            data.(cond{c}).datatable = [data.(cond{c}).datatable; an(a),0,numNonbrg,NaN,NaN,avgNonbrglength,NaN,NaN,cellTotLength,NaN,avgNonBrgNodeRatio,NaN,NaN,NaN,NaN];
            continue
        end
        flag = [];
        anchNames = cell(size(brgdata00,1),1);
        for i = 1:size(brgdata00,1)
            brgID = brgdata00.intIDs{i};
            temp = brgdata00.intStartsOn{i};
            if ~isempty(temp)
                ec = cellfun('isempty',temp);
                temp(ec) = [];
            end
            brgSO = temp;
            procIdx = find(ismember(cell2mat(tempdata00.Process_Info.procIDs), cell2mat(brgSO)));
            procSO = tempdata00.Process_Info.procStartsOn(procIdx);
            intIDs = cellfun(@cell2mat,tempdata00.Sheath_Info.intIDs,'UniformOutput',false);
            if size(procSO,1) < 1 
                % for whatever reason, bridge path is not connected to an anchoring sheath
                anchdata14.intLengths = NaN;
                anchdata14.node_class = NaN;
            elseif isempty(procSO{1})
                % for whatever reason, bridge path is not connected to an anchoring sheath
                anchdata14.intLengths = NaN;
                anchdata14.node_class = NaN;
            else 
                % change idx and id to sheath names, as ID may have changed between timecourse edits
                anchIdx00 = cell2mat(cellfun(@(x) any(ismember(x,cell2mat(procSO))), intIDs, 'UniformOutput', false));
                anchID00 = tempdata00.Sheath_Info.intIDs(anchIdx00);
                anchName00 = tempdata00.Sheath_Info.intNames(anchIdx00);
                anchdata14 = tempdata14.Sheath_Info(contains(tempdata14.Sheath_Info.intNames, anchName00),:);
                anchNames(i) = anchName00;
            end
            if tempdata00.Sheath_Info.isbridgecnctd(anchIdx00)
                for j = 1:size(brgdata14,1)
                    if contains(brgdata14.intNames{j},anchName00)
                        if j < i % move last bridged sheath in chain to first position e.g. [3* 2* 1*]
                            data.(cond{c}).(an{a}).lengths{j} = [brgdata14.intLengths(i), data.(cond{c}).(an{a}).lengths{j}]; 
                            data.(cond{c}).(an{a}).nodes{j} = [brgdata14.node_class(i), data.(cond{c}).(an{a}).nodes{j}]; 
                        else % store this bridged sheath's data until 'i' iteration meets 'j'
                            flag = [flag; j];
%                             flaganchdata{j} = [tempdata14.Sheath_Info.intLengths(anchIdx00), tempdata14.Sheath_Info.node_class(anchIdx00)];
                            flagbrgdata{j} = [brgdata14.intLengths(i), brgdata14.node_class(i)];
                        end
                    end
                end
            else % anchoring sheath is normal, assign bridged and anchoring lengths
                data.(cond{c}).(an{a}).lengths{i} = [brgdata14.intLengths(i), anchdata14.intLengths(1)];
                data.(cond{c}).(an{a}).anchlengths{i} = anchdata14.intLengths(1);
                data.(cond{c}).(an{a}).nodes{i} = [brgdata14.node_class(i), anchdata14.node_class(1)];
                % add hasmltproc data
                data.(cond{c}).(an{a}).anchMultProcs(i) = anchdata14.hasmltprocs(1);
            end
            % ADD CASE FOR B-N-B
            if any(ismember(flag,i))
                f = flag(ismember(flag,i));
                data.(cond{c}).(an{a}).lengths{i} = [flagbrgdata{j}(1), data.(cond{c}).(an{a}).lengths{i}];
                data.(cond{c}).(an{a}).nodes{i} = [flagbrgdata{j}(2), data.(cond{c}).(an{a}).nodes{i}];
            end
        end
        
        for i = 1:length(data.(cond{c}).(an{a}).lengths)
            data.(cond{c}).(an{a}).lengthSums(i) = sum(data.(cond{c}).(an{a}).lengths{i});
            data.(cond{c}).(an{a}).noderatio(i) = mean(data.(cond{c}).(an{a}).nodes{i},'omitnan');
        end
        
        avglengthSum = mean(data.(cond{c}).(an{a}).lengthSums,'omitnan');
        numbrgs = length(data.(cond{c}).(an{a}).lengths);
        anchNames(cellfun(@isempty,anchNames)) = [];
        nonbrgdata = tempdata14.Sheath_Info((~brgIdx & ~contains(tempdata14.Sheath_Info.intNames, anchNames)),:);
        numNonbrg = length(tempdata14.Sheath_Info.intNames);
        cellTotLength = sum(tempdata14.Sheath_Info.intLengths);
%         cellTotNum = ALL.(cond{c}).(an{a}).num_FinalSheaths(end);
        avganchlength = mean(anchdata14.intLengths);
        avgbrglength = mean(brgdata14.intLengths);
        totbrglength = sum(brgdata14.intLengths,'omitnan');
        avgNonbrglength = mean(nonbrgdata.intLengths);
        avgChainNodeRatio = mean(data.(cond{c}).(an{a}).noderatio,'omitnan');
        avgNonBrgNodeRatio = mean(nonbrgdata.node_class,'omitnan');
        
        doplot = 0;
        [xy,z] = compareTerritories_bridges(tempdata14.xmlstruct,cond{c},1,doplot);
        [xy_nb,z_nb] = compareTerritories_bridges(tempdata14.xmlstruct,cond{c},0,doplot);
        
        data.(cond{c}).datatable = [data.(cond{c}).datatable; an(a),numbrgs,numNonbrg,avgbrglength,totbrglength,avgNonbrglength,avganchlength,avglengthSum,cellTotLength,avgChainNodeRatio,avgNonBrgNodeRatio,xy,z,xy_nb,z_nb];
    end
    data.(cond{c}).datatable(1,:) = [];
end
% DON'T OMIT NANS FOR SUMMARY STATS (within PNB chains)!!!! 

%% Plot standards
ctrl = data.ctrl.datatable;
bsln = data.bsln.datatable;
l23 = data.L2_3.datatable;
[ct,cu] = getFigColors;
colors = {[0.5 0.5 0.5],[0.5 0.5 0.5],ct,ct,cu,cu};
xes = repmat(1:3,2,1);


%% Average lengths for all sheath types
colors2 = {[0.2 0.5 1],[52 75 160]./255,[1 0.5 0],[0.4 0.4 0.4]}; %,[1 0 0]
all = [ctrl;bsln];
 
y = {all.avgNonbrglength,all.avganchLength,all.avgbrglength,all.avgChainLengthSum};
avg = cell2mat(cellfun(@(x) mean(x,'omitnan'),y,'UniformOutput',false));
sem = cell2mat(cellfun(@(x) calcSEM(x,1),y,'UniformOutput',false));
figure
hold on
plotSpread(y,'distributionMarker','o','distributionColors',colors2);
errorbar(avg,sem,'ko','MarkerSize',4,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
xticklabels({})
ylim([0 200])
figQuality(gcf,gca,[2 2.6])

y2 = [all.avgNonbrglength;all.avganchLength;all.avgbrglength;all.avgChainLengthSum];
idx = [repmat({'non'},size(all.avgNonbrglength)); repmat({'anch'},size(all.avganchLength));...
       repmat({'brg'},size(all.avgbrglength)); repmat({'chain'},size(all.avgChainLengthSum))];
[p,tbl,stats] = kruskalwallis(y2,idx);
compresults = multcompare(stats);
%% Average chain length VS average sheath length
y = {bsln.avgChainLengthSum, bsln.avgNonbrglength,...
    ctrl.avgChainLengthSum, ctrl.avgNonbrglength};

plotEmUp(y,xes2,colors2)
xticklabels({})
figQuality(gcf,gca,[2 2.6])
y2 = [bsln.avgChainLengthSum; bsln.avgNonbrglength;...
     ctrl.avgChainLengthSum; ctrl.avgNonbrglength];
[~,p,tstat] = ttest2(bsln.avgChainLengthSum,bsln.avgNonbrglength);
p*2
[~,p,tstat] = ttest2(ctrl.avgChainLengthSum,ctrl.avgNonbrglength);
p*2

%% Proportion of bridges
y = {bsln.numbrgs./(bsln.numbrgs + bsln.numNonbrg),...
    ctrl.numbrgs./(ctrl.numbrgs + ctrl.numNonbrg)};

y2 = [bsln.numbrgs./(bsln.numbrgs + bsln.numNonbrg);...
    ctrl.numbrgs./(ctrl.numbrgs + ctrl.numNonbrg)];
[p,tbl,stats] = kruskalwallis(y2,[repmat({'bsln'},size(bsln.numbrgs));repmat({'ctrl'},size(ctrl.numbrgs))]);

avg = reshape(cellfun(@(y) mean(y,'omitnan'),y),1,2)
sem = reshape(cellfun(@(y) calcSEM(y,1),y),1,2)
sz = cellfun(@length,y);
idx = arrayfun(@(x) ones(x,1), sz, 'UniformOutput', false);
idx2 = [];
for i = 1:length(idx)
    idx2 = [idx2; idx{i} .* i];
end
figure

plotSpread(y2,'distributionIdx',idx2,'distributionMarker','o','distributionColors',{[0.5 0.5 0.5],ct})
hold on
errorbar(avg,sem,'ko','MarkerSize',4,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5);
hold off
xlim([0 3])
% row1 = {'< 8wks' '> 8wks' '> 8wks'};
% row2 = {' ' 'Control' 'Cuprizone'};
% labelArray = [row1; row2];
% labelArray = strjust(pad(labelArray),'center'); % 'left'(default)|'right'|'center
% ticklabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
ax = gca();
% ax.XTickLabel = ticklabels;
xticklabels({})
figQuality(gcf,gca,[1.8 2.6])
ytickformat('%.2f')

%% bsln + ctrl generated VS OLD lost
numbrgs = [bsln.numbrgs; ctrl.numbrgs];
numtot = [bsln.numbrgs + bsln.numNonbrg; + ctrl.numbrgs + ctrl.numNonbrg];
genBrgRatio = numbrgs./numtot;
[out2, lostBrgRatio] = ParseOldMouseTraces;
figure
colors = getPNBFigColors;
plotSpread({genBrgRatio,lostBrgRatio},'distributionMarker',{'o'}, 'distributionColors',{[0.5 0.5 0.5],colors{3}});
figQuality(gcf,gca,[2,2.9])
xticklabels({'Generated','Lost'})
ytickformat('%.1f')
p = kruskalwallis([genBrgRatio;lostBrgRatio],[zeros(size(genBrgRatio));ones(size(lostBrgRatio))])
%% Cell tot length VS bridge number
figure; plot(bsln.numbrgs,bsln.cellTotLength./1000,'o','MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(ctrl.numbrgs,ctrl.cellTotLength./1000,'o','MarkerEdgeColor',ct,'LineWidth',1); hold off;
figQuality(gcf,gca,[2.8 2.8])
ylim([0 6])
xlim([0 12])

[R,P,RL,RU] = corrcoef(bsln.numbrgs,bsln.cellTotLength)
[R,P,RL,RU] = corrcoef(ctrl.numbrgs,ctrl.cellTotLength)
%% Proportion of bridges by length
y = {bsln.totbrglength./(bsln.cellTotLength),...
    ctrl.totbrglength./(ctrl.cellTotLength),...
    cupr.totbrglength./(cupr.cellTotLength)};

avg = reshape(cellfun(@(y) mean(y,'omitnan'),y),1,3);
sem = reshape(cellfun(@(y) calcSEM(y,1),y),1,3);
sz = cellfun(@length,y);
idx = arrayfun(@(x) ones(x,1), sz, 'UniformOutput', false);
idx2 = [];
for i = 1:length(idx)
    idx2 = [idx2; idx{i} .* i];
end
figure
errorbar(avg,sem,'ko','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1.5)
hold on
plotSpread(y,'distributionIdx',idx2,'distributionMarker','o','distributionColors',{[0.5 0.5 0.5],ct,cu})
hold off
xlim([0 4])
row1 = {'< 8wks' '> 8wks' '> 8wks'};
row2 = {' ' 'Control' 'Cuprizone'};
labelArray = [row1; row2];
labelArray = strjust(pad(labelArray),'center'); % 'left'(default)|'right'|'center
ticklabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
ax = gca();
ax.XTickLabel = ticklabels;
figQuality(gcf,gca,[3 2.7])

%% CUPRIZONE COMPARISONS
%% Proportion of bridges
controls = [bsln; ctrl];


y = {ctrl.numbrgs./(ctrl.numbrgs + ctrl.numNonbrg),...
    cupr.numbrgs./(cupr.numbrgs + cupr.numNonbrg)};

y2 = [ctrl.numbrgs./(ctrl.numbrgs + ctrl.numNonbrg);...
    cupr.numbrgs./(cupr.numbrgs + cupr.numNonbrg)];
[p,tbl,stats] = kruskalwallis(y2,[repmat({'control'},size(ctrl.numbrgs));repmat({'remyelinating'},size(cupr.numbrgs))]);

avg = reshape(cellfun(@(y) mean(y,'omitnan'),y),1,2)
sem = reshape(cellfun(@(y) calcSEM(y,1),y),1,2)
sz = cellfun(@length,y);
idx = arrayfun(@(x) ones(x,1), sz, 'UniformOutput', false);
idx2 = [];
for i = 1:length(idx)
    idx2 = [idx2; idx{i} .* i];
end
figure
errorbar(avg,sem,'ko','MarkerSize',5,'MarkerFaceColor','k','CapSize',0,'LineWidth',1.5)
hold on
plotSpread(y,'distributionIdx',idx2,'distributionMarker','o','distributionColors',{[0.5 0.5 0.5],cu})
hold off
xlim([0 3])
% ylim([0 0.25])
ax = gca();
xticklabels({'Control' 'Remyelinating'})
figQuality(gcf,gca,[1.8 2.8])
ax.XTickLabelRotation = 45;
ytickformat('%.2f')

[p,tbl,stats] = kruskalwallis(y2,idx2)

%% Cell tot length VS bridge number
figure; plot(cupr.numbrgs,cupr.cellTotLength./1000,'o','MarkerEdgeColor',cu,'LineWidth',1); hold on;
plot(ctrl.numbrgs,ctrl.cellTotLength./1000,'o','MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',1);hold off;
figQuality(gcf,gca,[2.8 2.8])
ylim([0 6])
xlim([0 12])

[R,P,RL,RU] = corrcoef([[cupr.numbrgs] [cupr.cellTotLength]])

%% Average bridge length VS average sheath length
colors2 = {[0.5 0.5 0.5],[0.5 0.5 0.5],cu,cu};
xes2 = repmat(1:2,2,1);

y = {ctrl.avgbrglength, ctrl.avgNonbrglength,...
     cupr.avgbrglength, cupr.avgNonbrglength};

plotEmUp(y,xes2,colors2)
xticklabels({'Control' 'Remyelinating'})
figQuality(gcf,gca,[2.8 2.2])

y2 = [ctrl.avgbrglength; ctrl.avgNonbrglength;...
     cupr.avgbrglength; cupr.avgNonbrglength];
idx = [repmat({'bsln_b'},size(ctrl.cellname)); repmat({'bsln_nb'},size(ctrl.cellname));...
       repmat({'ctrl_b'},size(cupr.cellname)); repmat({'ctrl_nb'},size(cupr.cellname))];
[p,tbl,stats] = kruskalwallis(y2,idx);
multcompare(stats,"CType","dunn-sidak");
%% Average chain length VS average sheath length
y = {ctrl.avgChainLengthSum, ctrl.avgNonbrglength,...
    cupr.avgChainLengthSum, cupr.avgNonbrglength};

plotEmUp(y,xes2,colors2)
xticklabels({'Control' 'Remyelinating'})
figQuality(gcf,gca,[2.8 2.2])

y2 = [ctrl.avgChainLengthSum; ctrl.avgNonbrglength;...
     cupr.avgChainLengthSum; cupr.avgNonbrglength];
[~,p,tstat] = ttest(ctrl.avgChainLengthSum,ctrl.avgNonbrglength);
p*2
[~,p,tstat] = ttest(cupr.avgChainLengthSum,cupr.avgNonbrglength);
p*2

idx = [repmat({'bsln_b'},size(ctrl.cellname)); repmat({'bsln_nb'},size(ctrl.cellname));...
       repmat({'ctrl_b'},size(cupr.cellname)); repmat({'ctrl_nb'},size(cupr.cellname))];
[p,tbl,stats] = kruskalwallis(y2,idx);
multcompare(stats,"CType","dunn-sidak");

%% Territory comparison
% y = {bsln.xy, bsln.xy_nb, bsln.z, bsln.z_nb,...
%     ctrl.xy, ctrl.xy_nb, ctrl.z, ctrl.z_nb};
% y = {[bsln.xy;ctrl.xy],[bsln.xy_nb;ctrl.xy_nb],[bsln.z;ctrl.z],[bsln.z_nb;ctrl.z_nb]};
% colors2 = {[0.5 0.5 0.5],[0.5 0.5 0.5],ct,ct};
% 
% avg = reshape(cellfun(@(y) mean(y,'omitnan'),y),[],2)
% sem = reshape(cellfun(@(y) calcSEM(y,1),y),[],2)
% sz = cellfun(@length,y);
% idx = arrayfun(@(x) ones(x,1), sz, 'UniformOutput', false);
% idx2 = [];
% for i = 1:length(idx)
%     idx2 = [idx2; idx{i} .* i];
% end
% 
% y2 = [bsln.xy; bsln.xy_nb; bsln.z; bsln.z_nb;...
%     ctrl.xy; ctrl.xy_nb; ctrl.z; ctrl.z_nb];
% 
% y = [[bsln.xy;ctrl.xy],[bsln.xy_nb;ctrl.xy_nb],[bsln.z;ctrl.z],[bsln.z_nb;ctrl.z_nb]];
% figure
% color = [0.5 0.5 0.5];
% x = [1 2; 3 4];
% hold on
% for i = 1:size(y,1)
%     plot(x(1,:),y(i,1:2),'o-','Color',[0.5 0.5 0.5],'MarkerSize',3,'LineWidth',0.5);
% end
% for i = 1:size(y,1)
%     plot(x(2,:),y(i,3:4),'o-','Color',[0.5 0.5 0.5],'MarkerSize',3,'LineWidth',0.5);
% end
% xlim([0.5,4.5]);
% ylim([0,120]);
% xticks([1.5,3.5]);
% errorbar([1,2;3,4],avg',sem','ko','MarkerSize',2,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',0)
% hold off
% set(gca,'XTickLabel',[])
% figQuality(gcf,gca,[2,2.2]);