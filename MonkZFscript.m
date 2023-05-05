cd('D:\GitHubRepos\Call_ParanodalBridge_2022\MonkZFtraces');
% file loop
xmlstruct = parseXML_SingleCell('fish08_cellO_nonbridgedSheath1_230503_no16.traces');
% save struct
%     varname = matlab.lang.makeValidName(contents(k).name);
%     save(fullfile(directory,aClass,[varname '.mat']),'xmlstruct');
%%
useImaginaryEdge = 0;
makeplot = 1;
[d,X,Delta,borders,origin] = calculatePathsXML_ZFbridges(xmlstruct,useImaginaryEdge,makeplot);

%%
smol = min(cell2mat(cellfun(@min,cellfun(@min,Delta,'UniformOutput',false),'UniformOutput',false)));
big = max(cell2mat(cellfun(@max,cellfun(@max,Delta,'UniformOutput',false),'UniformOutput',false)));
rng = smol:big;
Dadj = NaN(size(rng));
Dall = [];
for i = 1:length(Delta)
    if any(Delta{i})
        Dadj(ismember(rng,Delta{i}(:,1))) = movmean(Delta{i}(:,2),50);
        Dall = [Dall;Dadj];
    end
end
% add extra blank rows to be able to get peaks/troughs in first and last rows
Dall = [zeros(1,size(Dall,2));Dall;zeros(1,size(Dall,2))];

%%
% smol = min(cell2mat(cellfun(@min,cellfun(@min,X,'UniformOutput',false),'UniformOutput',false)));
% big = max(cell2mat(cellfun(@max,cellfun(@max,X,'UniformOutput',false),'UniformOutput',false)));
% rng = smol:big;
% drawadj = NaN(size(rng));
% drawall = [];
% for i =1:length(d)
%     if any(d{i})
%         drawadj(ismember(X{i},rng)) = d{i}(:,1);
%         drawall = [drawall; drawadj];
%     end
% end
%%
lmin = islocalmin(Dall,1,"MinProminence",0.5);
lmax = islocalmax(Dall,1,"MinProminence",0.5);
lmin = lmin.*-1;
lmax = double(lmax);
extr = lmin + lmax;
extr(isnan(Dall)) = NaN;
figure
heatmap(extr(2:end-1,:),'Colormap',parula)
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
figQuality(gcf,gca,[11,2],1)

figure; heatmap(Dall(2:end-1,:),'Colormap',parula,'ColorLimits',[-1.5,1.5]); grid off
figQuality(gcf,gca,[11,2],1)
%% calculate extension and keep track of wraps up to then when extended
wraps = zeros(1,size(Dall,2));
for i = 1:size(Dall,1)-1
    if ~any(extr(i,:))
        continue
    end
    
    % continue prev. wraps with extensions
    temp = find(~isnan(Dall(i,:)));
    edges_i = [temp(1),temp(end)];
    temp = find(~isnan(Dall(i-1,:)));
    edges_0i = [temp(1),temp(end)];
    if edges_i(1) < edges_0i(1)
        wraps(edges_i(1):edges_0i(1)) = wraps(edges_i(1):edges_0i(1)) + repmat(mode(wraps(edges_0i(1)+20)),[1,edges_0i(1)-edges_i(1)+1]);
    end
    if edges_i(end) > edges_0i(end)
        wraps(edges_0i(end):edges_i(end)) = wraps(edges_0i(end):edges_i(end)) + repmat(mode(wraps(edges_0i(end)-20)),[1,edges_i(end)-edges_0i(end)+1]);
    end
    
    % get new peaks/troughs
    temp = cummax(cumsum(lmax(1:i-1,:)) + cumsum(lmin(1:i-1,:)));
    peaks_last = find(temp(end,:));
    temp = cummin(cumsum(lmax(1:i-1,:)) + cumsum(lmin(1:i-1,:)));
    troughs_last = find(temp(end,:));
    peaks_curr = find(lmax(i,:));
    troughs_curr = find(lmin(i,:));
    idx_down = peaks_last(ismember(peaks_last,troughs_curr));
    idx_up = troughs_last(ismember(troughs_last,peaks_curr));
    
    wraps(idx_down) = wraps(idx_down) + 0.5;
    wraps(idx_up) = wraps(idx_up) + 0.5;
    
end
figure 
plot(wraps)
hold on
plot(smoothdata(wraps,'movmedian',50))
hold off
figQuality(gcf,gca,[11,2])
box off
%%
data = {};
for i = 1:length(Delta)
if any(Delta{i})
    data{i} = [Delta{i}(1:50:end,1),Delta{i}(1:50:end,2)];
end
end
figure
hold on

