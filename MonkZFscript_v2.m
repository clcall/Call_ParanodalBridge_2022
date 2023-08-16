cd('D:\GitHubRepos\Call_ParanodalBridge_2022\MonkZFtraces');
clear variables
name = 'fish07_cellE_bridgedSheath2_230718_lo.traces';
xmlstruct = parseXML_SingleCell(name);
if contains(name,'lo')
    edgepos = 'lo';
elseif contains(name,'hi')
    edgepos = 'hi';
end
% TO DO:
% file loop
% save struct
%     varname = matlab.lang.makeValidName(contents(k).name);
%     save(fullfile(directory,aClass,[varname '.mat']),'xmlstruct');
%
% convert length back to microns for wrap quant
%
% convert tp to real time
%%
useImaginaryEdge = 0;
makeplot = 1;
[d,X,Delta,borders,origin,firstbrgtp] = calculatePathsXML_ZFbridges(xmlstruct,edgepos,useImaginaryEdge,makeplot);
set(gcf,'color','w');
%%
starts = NaN(size(Delta));
ends = NaN(size(Delta));
for i = 1:length(Delta)
    if ~isempty(Delta{i})
        starts(i) = min(Delta{i}(:,1));
        ends(i) = max(Delta{i}(:,1));
    else
        starts(i) = [];
        ends(i) = [];
    end
end
[rngstart,sI] = min(starts);
[rngend,eI] = max(ends);

globalorigin = origin{sI};
rng = rngstart:rngend;
Dadj = NaN(size(rng));
Dall = [];
for i = 1:length(Delta)
    if any(Delta{i})
        Dadj(ismember(rng,Delta{i}(:,1))) = movmean(Delta{i}(:,2),50);
        Dall = [Dall;Dadj];
    end
end
Dall = [zeros(1,size(Dall,2));Dall;zeros(1,size(Dall,2))]; % adds extra blank rows to be able to get peaks/troughs in first and last rows

lmin = islocalmin(Dall,1,"MinProminence",0.4);
lmax = islocalmax(Dall,1,"MinProminence",0.4);
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

figure; 
heatmap(Dall(2:end-1,:),'Colormap',parula,'ColorLimits',[-1.5,1.5]); grid off
figQuality(gcf,gca,[11,2],1)

%% calculate extension and keep track of wraps up to then when extended
wraps = zeros(1,size(Dall,2));
avgbrgwraps = cell(1,size(Dall,1));
avganchwraps = NaN(1,size(Dall,1));
smallestanch = avganchwraps;
allnonbrg = avganchwraps;
xes = -globalorigin+1 : abs(length(wraps) - globalorigin);

for b = 1:length(borders)
    if any(borders{b}>0)
        rightbord_default(b) = min(borders{b}(borders{b}>0));
        rightbord_default(b) = rightbord_default(b) + globalorigin;
    else
        rightbord_default(b) = NaN;
    end
    if any(borders{b}<0)
        leftbord_default(b) = max(borders{b}(borders{b}<0));
        leftbord_default(b) = leftbord_default(b) + globalorigin;
    else
        leftbord_default(b) = NaN;
    end
end
rightbord_default = min(rightbord_default);
leftbord_default = max(leftbord_default);
% xes = [Delta{end}(1,1)-1; Delta{end}(:,1); Delta{end}(end,1)+1]; % adds 1 extra either side to match Dall
figure
hold on
for i = 1:size(Dall,1)-1
    if ~any(extr(i,:))
        continue
    end
    
    % continue prev. wraps with extensions
    temp = find(~isnan(Dall(i,:))); % curr tp edges
    edges_i = [temp(1),temp(end)];
    temp = find(~isnan(Dall(i-1,:))); % prev tp edges
    edges_0i = [temp(1),temp(end)];
    
    % DOUBLE CHECK THIS SECTION
%     if edges_i(1) < edges_0i(1) % extension to left
%         idx = edges_0i(1)+50; 
%         if idx > edges_0i(end)
%             idx = edges_0i(end);
%         end
%         adj = repmat(mode(wraps(idx)),[1,edges_0i(1)-edges_i(1)+1]);
%         wraps(edges_i(1):edges_0i(1)) = wraps(edges_i(1):edges_0i(1)) + adj;
%     end
%     if edges_i(end) > edges_0i(end) % extension to right
%         idx = edges_0i(end)-50;
%         if idx < 0
%             idx = 1;
%         end
%         adj = repmat(mode(wraps(idx)),[1,edges_i(end)-edges_0i(end)+1]);
%         wraps(edges_0i(end):edges_i(end)) = wraps(edges_0i(end):edges_i(end)) + adj;
%     end
    % deal with retractions
    if edges_i(1) > edges_0i(1) % retraction from left
        wraps(edges_0i(1):edges_i(1)) = 0;
    end
    if edges_i(end) < edges_0i(end) % retraction from right
        wraps(edges_i(end):edges_0i(end)) = 0;
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
    
    % FIND MEAN # WRAPS WITHIN BRIDGE
    nonbrgmask = ones(size(wraps));
    for r = 1:2:length(borders{i})
        startpos = borders{i}(r) + globalorigin;
        endpos = borders{i}(r+1) + globalorigin;
%         startpos = round((border1 + (border2-border1)./2) - 5/0.06); %start pos 1st bridge
%         endpos = round((border1 + (border2-border1)./2) + 5/0.06); %end pos 1st bridge
        %compensate for bridges close to paranodes
%         if startpos < edges_i(1)
%             startpos = edges_i(1);
%         end
%         if endpos > edges_i(end)
%             endpos = edges_i(end);
%         end
        avgbrgwraps{i} = [avgbrgwraps{i}, mean(wraps(startpos:endpos))];
        nonbrgmask(startpos:endpos) = 0;
    end
    
    % FIND MEAN # WRAPS OF ANCHORING SHEATH (from origin to nearest bridge
    % border or other paranode)
    if any(borders{i}>0)
        temp1 = min(borders{i}(borders{i}>0));
        idx1 = globalorigin:temp1+globalorigin;
        rightpart = wraps(idx1);
    else
        idx1 = edges_i(end);
        rightpart = wraps(idx1);
    end
    if any(borders{i}<0)
        temp2 = max(borders{i}(borders{i}<0));
        idx2 = temp2+globalorigin:globalorigin-1;
        leftpart = wraps(idx2);
    else
        idx2 = edges_i(1);
        leftpart = wraps(idx2);
    end
    avganchwraps(i) = mean([leftpart,rightpart]);
    
    %extrapolate back to only include areas never bridged

    if isnan(rightbord_default)
        rightbord = edges_i(end);
    else
        rightbord = rightbord_default;
    end
    if isnan(leftbord_default)
        leftbord = edges_i(1);
    else
        leftbord = leftbord_default;
    end
    
    smallestanch(i) = mean(wraps(leftbord:rightbord),'omitnan');
    
    % get wraps for all bridged sheaths
    idx = ~isnan(Dall(i,:));
    idx(idx1) = 0;
    idx(idx2) = 0;
    idx(leftbord:rightbord) = 0;
    
    allnonbrg(i) = mean(wraps(idx));
    
    % make shaded regions for pnb positions
    if ~isempty(borders{i})
        rectangle('Position',[borders{i}(1), 0, borders{i}(2)-borders{i}(1), 5], 'FaceColor', [0.3 0.3 0.3 0.1], 'EdgeColor', 'none');
        if length(borders{i}) > 2
            rectangle('Position',[borders{i}(3), 0, borders{i}(4)-borders{i}(3), 5], 'FaceColor', [0.3 0.3 0.3 0.1], 'EdgeColor', 'none');
        end
        if length(borders{i}) > 4
            rectangle('Position',[borders{i}(5), 0, borders{i}(6)-borders{i}(5), 5], 'FaceColor', [0.3 0.3 0.3 0.1], 'EdgeColor', 'none');
        end
    end
end

plot(xes,wraps)
plot(xes,smoothdata(wraps,'movmedian',50))

hold off
figQuality(gcf,gca,[11,2])
box off

figure
% plot(avganchwraps,'k-'); 
hold on
plot(smallestanch,'b-');
plot(cell2mat(cellfun(@mean,avgbrgwraps,"UniformOutput",false)),'r-');
plot(allnonbrg,'g-');
hold off
ylabel('Wraps')
xlabel('Timepoint')
legend('anchoring sheath','bridges','all non-anch sheath')
figQuality(gcf,gca,[4,2.5])

