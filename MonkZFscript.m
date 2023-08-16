cd('D:\GitHubRepos\Call_ParanodalBridge_2022\MonkZFtraces');

xmlstruct = parseXML_SingleCell('fish08_cellO_bridgedSheath1_230504_no16.traces');
% TO DO:
% file loop
% save struct
%     varname = matlab.lang.makeValidName(contents(k).name);
%     save(fullfile(directory,aClass,[varname '.mat']),'xmlstruct');
%
% convert length back to microns for wrap quant
%
% convert tp to real time
%
% note locations of bridges per tp
%%
useImaginaryEdge = 0;
makeplot = 1;
[d,X,Delta,borders,origin,firstbrgtp] = calculatePathsXML_ZFbridges(xmlstruct,useImaginaryEdge,makeplot);

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
Dall = [zeros(1,size(Dall,2));Dall;zeros(1,size(Dall,2))]; % adds extra blank rows to be able to get peaks/troughs in first and last rows

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
figQuality(gcf,gca,[11,2])

figure; heatmap(Dall(2:end-1,:),'Colormap',parula,'ColorLimits',[-1.5,1.5]); grid off
figQuality(gcf,gca,[11,2])

%% calculate extension and keep track of wraps up to then when extended
wraps = zeros(1,size(Dall,2));
avgbrgwraps = cell(1,size(Dall,1));
avganchwraps = NaN(1,size(Dall,1));
xes = -origin{end}+1 : abs(length(wraps) - origin{end});
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
    if edges_i(1) < edges_0i(1) % extension to left
        idx = edges_0i(1)+50; 
        if idx > edges_0i(end)
            idx = edges_0i(end);
        end
        adj = repmat(mode(wraps(idx)),[1,edges_0i(1)-edges_i(1)+1]);
        wraps(edges_i(1):edges_0i(1)) = wraps(edges_i(1):edges_0i(1)) + adj;
    end
    if edges_i(end) > edges_0i(end) % extension to right
        idx = edges_0i(end)-50;
        if idx < 0
            idx = 1;
        end
        adj = repmat(mode(wraps(idx)),[1,edges_i(end)-edges_0i(end)+1]);
        wraps(edges_0i(end):edges_i(end)) = wraps(edges_0i(end):edges_i(end)) + adj;
    end
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
    
%     switch i
%         case firstbrgtp-1
%             figure 
%             plot(xes,wraps)
%             title('TP(first bridge) - 1')
%             hold on
%             plot(xes,smoothdata(wraps,'movmedian',50))
%             for r = 1:2:length(borders{i})% make shaded regions for pnb positions IN NEXT TP
%                 rectangle('Position',[borders{i+1}(r), 0, borders{i+1}(r+1)-borders{i+1}(r), 6], 'FaceColor', [0.3 0.3 0.3 0.3], 'EdgeColor', 'none');
%             end
%             hold off
%             figQuality(gcf,gca,[11,2])
%             box off
%             ylim([0,5])
%         case firstbrgtp
%             figure 
%             plot(xes,wraps)
%             title('TP(first bridge)')
%             hold on
%             plot(xes,smoothdata(wraps,'movmedian',50))
%             for r = 1:2:length(borders{i})% make shaded regions for pnb positions
%                 rectangle('Position',[borders{i}(r), 0, borders{i}(r+1)-borders{i}(r), 6], 'FaceColor', [0.3 0.3 0.3 0.3], 'EdgeColor', 'none');
%             end
%             hold off
%             figQuality(gcf,gca,[11,2])
%             box off
%             ylim([0,5])
%         case firstbrgtp+5 % note this is arbitrary, still need to convert tp to real time
%             figure 
%             plot(xes,wraps)
%             title('TP(first bridge) + 5')
%             hold on
%             plot(xes,smoothdata(wraps,'movmedian',50))
%             for r = 1:2:length(borders{i}) % make shaded regions for pnb positions
%                 rectangle('Position',[borders{i}(r), 0, borders{i}(r+1)-borders{i}(r), 6], 'FaceColor', [0.3 0.3 0.3 0.3], 'EdgeColor', 'none');
%             end
%             hold off
%             figQuality(gcf,gca,[11,2])
%             box off
%             ylim([0,5])
%     end
    
    % FIND MEAN # WRAPS WITHIN 10 MICRONS (currently approximate scaling) OF BRIDGE CENTER
    for r = 1:2:length(borders{i})
        border1 = borders{i}(r) + origin{end};
        border2 = borders{i}(r+1) + origin{end};
        startpos = round((border1 + (border2-border1)./2) - 5/0.06); %start pos 1st bridge
        endpos = round((border1 + (border2-border1)./2) + 5/0.06); %end pos 1st bridge
        %compensate for bridges close to paranodes
        if startpos < edges_i(1)
            startpos = edges_i(1);
        end
        if endpos > edges_i(end)
            endpos = edges_i(end);
        end
        avgbrgwraps{i} = [avgbrgwraps{i}, mean(wraps(startpos:endpos))];
    end
    % FIND MEAN # WRAPS OF ANCHORING SHEATH (from origin to nearest bridge
    % border or other paranode)
    % NEED TO EXTRAPOLATE BACKWARDS TO EXCLUDE AREAS THAT BECOME BRIDGED???
    if any(borders{i}>0)
        temp = min(borders{i}(borders{i}>0));
        rightpart = wraps(origin{end}:temp+origin{end});
    else
        rightpart = wraps(origin{end}:edges_i(end));
    end
    if any(borders{i}<0)
        temp = max(borders{i}(borders{i}<0));
        leftpart = wraps(temp+origin{end}:origin{end}-1);
    else
        leftpart = wraps(edges_i(1):origin{end}-1);
    end
    avganchwraps(i) = mean([leftpart,rightpart]);
end
figure 
plot(xes,wraps)
hold on
plot(xes,smoothdata(wraps,'movmedian',50))
% make shaded regions for pnb positions
rectangle('Position',[borders{end}(1), 0, borders{end}(2)-borders{end}(1), 5], 'FaceColor', [0.3 0.3 0.3 0.3], 'EdgeColor', 'none');
if length(borders{end}) > 2
    rectangle('Position',[borders{end}(3), 0, borders{end}(4)-borders{end}(3), 5], 'FaceColor', [0.3 0.3 0.3 0.3], 'EdgeColor', 'none');
end
if length(borders{end}) > 4
    rectangle('Position',[borders{end}(5), 0, borders{end}(6)-borders{end}(5), 5], 'FaceColor', [0.3 0.3 0.3 0.3], 'EdgeColor', 'none');
end
hold off
figQuality(gcf,gca,[11,2])
box off

% WHERE ARE NANS COMING FROM IN AVGWRAPS ???
figure
plot(avganchwraps,'k-'); 
hold on
plot(cell2mat(cellfun(@mean,avgbrgwraps,"UniformOutput",false)),'r-');
hold off
ylabel('Wraps')
xlabel('Timepoint')
legend('anchoring sheath','bridges')
figQuality(gcf,gca,[4,2.5])
%%
data = {};
for i = 1:length(Delta)
if any(Delta{i})
    data{i} = [Delta{i}(1:50:end,1),Delta{i}(1:50:end,2)];
end
end
figure
hold on

