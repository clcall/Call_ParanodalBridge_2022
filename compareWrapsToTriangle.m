function compareWrapsToTriangle(SD)
% NOTE: run MonkZFscript_allsheaths.m FIRST to get SD struct with timelines
path = 'D:\GitHubRepos\Call_ParanodalBridge_2022\MonkZFtraces';
cd(path);
files = dir(path);
files = files(3:end,:);
files = files(~[files.isdir]);
files = files(~contains({files.name},'.mat'));
avgbrgwraps = cell(1,size(files,1));
avganchwraps = cell(1,size(files,1));
smallestanch = cell(1,size(files,1));
allnonbrg = cell(1,size(files,1));
allTL = cell(1,size(files,1));
for f = 1:size(files,1)
    name = files(f).name;
    varname = matlab.lang.makeValidName(name);
    if exist(fullfile(path,[varname '.mat']),'file')
        load([varname '.mat']);
    else
        xmlstruct = parseXML_SingleCell(name);
        save(fullfile(path,[varname '.mat']),'xmlstruct');
    end
    if contains(name,'lo')
        edgepos = 'lo';
    elseif contains(name,'hi')
        edgepos = 'hi';
    end
    
    %run MonkZFscript_allsheaths.m to get SD struct with timelines
    idx = contains({SD.name},name(1:12));
    peak = SD(idx).peakedge;
    TL = SD(idx).adjTL;
    allTL{f} = TL;
    
    % get path data and Delta distance between sheath & lifeact
    makeplot = 0;
    makeplot2 = 0;
    [~,~,Delta,borders,origin,~,framesUsed] = calculatePathsXML_ZFbridges_v2(xmlstruct,edgepos,makeplot,makeplot2,TL,name);
    allTL{f} = allTL{f}(framesUsed);
    
    % smooth Delta and concatenate across timepoints
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
    borders_copy = borders;
    borders = {};
    for i = 1:length(Delta)
        if any(Delta{i})
            Dadj(ismember(rng,Delta{i}(:,1))) = movmean(Delta{i}(:,2),50);
            Dall = [Dall;Dadj];
            borders = [borders, borders_copy(i)];
        end
    end
    % calculate extension
    wraps = abs(Dall);
    allWraps{f} = wraps;
    lastFrameIncl = ~isnan(wraps(end,:));
    allEndWraps{f} = sum(wraps(:,lastFrameIncl),1,'omitnan');
    avgbrgwraps{f} = cell(1,size(Dall,1));
    avganchwraps{f} = NaN(1,size(Dall,1));
    smallestanch{f} = NaN(1,size(Dall,1));
    allnonbrg{f} = NaN(1,size(Dall,1));
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
    
    for i = 2:size(Dall,1)-1
        temp = find(~isnan(Dall(i,:))); % curr tp edges
        edges_i = [temp(1),temp(end)];
        temp = find(~isnan(Dall(i-1,:))); % prev tp edges
        edges_0i = [temp(1),temp(end)];
        
        % FIND MEAN # WRAPS WITHIN BRIDGE
        nonbrgmask = ones(size(wraps(i,:)));
        for r = 1:2:length(borders{i})
            startpos = borders{i}(r) + globalorigin;
            endpos = borders{i}(r+1) + globalorigin;
            avgbrgwraps{f}{i} = [avgbrgwraps{f}{i}, mean(sum(wraps(1:i,startpos:endpos),'omitnan'))];
            nonbrgmask(startpos:endpos) = 0;
        end
        
        % FIND MEAN # WRAPS OF ANCHORING SHEATH (from origin to nearest bridge
        % border or other paranode)
        if any(borders{i}>0)
            temp1 = min(borders{i}(borders{i}>0));
            idx1 = globalorigin:temp1+globalorigin;
            rightpart = sum(wraps(1:i,idx1),'omitnan');
        else
            idx1 = edges_i(end);
            rightpart = sum(wraps(1:i,idx1),'omitnan');
        end
        if any(borders{i}<0)
            temp2 = max(borders{i}(borders{i}<0));
            idx2 = temp2+globalorigin:globalorigin-1;
            leftpart = sum(wraps(1:i,idx2),'omitnan');
        else
            idx2 = edges_i(1);
            leftpart = sum(wraps(1:i,idx2),'omitnan');
        end
        avganchwraps{f}(i) = mean([leftpart,rightpart]);
        
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
        smallestanch{f}(i) = mean(sum(wraps(1:i,leftbord:rightbord),'omitnan'));
        
        % get wraps for all bridged segments
        idx = ~isnan(Dall(i,:));
        idx(idx1) = 0;
        idx(idx2) = 0;
        idx(leftbord:rightbord) = 0;
        allnonbrg{f}(i) = mean(sum(wraps(1:i,idx),'omitnan'));
    end
end

figure
r = 56;
tritest = [flip(1-(1/r):-1/r:0.1), 1, 1-(1/r):-1/r:0.1];
for f = 1:size(files,1)
    [~,idx] = min(abs(allTL{f} - 25)); %25 hrs
    useable = ~isnan(allWraps{f}(idx,:));
    sumwraps = sum(allWraps{f}(1:idx,useable),1,'omitnan');
    sumwraps = rescale(imresize(sumwraps,[1,101]));
    [~,maxidx] = max(sumwraps);
    %     if maxidx > 51
    %         sumwraps = flip(sumwraps);
    %     end
    subplot(4,4,f,'TickDir','out')
    hold on
    [~,p] = kstest2(tritest,sumwraps);
    if p<0.05
        color = 'r';
    else
        color = 'k';
    end
    plot(tritest,'LineWidth',1.5)
    plot(sumwraps,'LineWidth',1.5,'Color',color)
    %     title(['p = ' num2str(p)])
    xlim([1 101])
    xticks([1 51 101])
    xticklabels([-1 0 1])
    yticks([0 0.5 1])
    yticklabels({'0.0' '0.5' '1.0'})
    hold off
    sumwraps_all(f,:) = sumwraps;
end

figure
plot(rescale(mean(sumwraps_all)))
[~,p] = kstest2(tritest,rescale(mean(sumwraps_all)))

%% local function
    function filled_indexed  = averageSheathTLs(cellinput,universalTL,roundedTL)
        TLidx = [];
        idxFL= [];
        indexed = [];
        filled_indexed = [];
        
        [~,TLidx(1,:)] = ismember(universalTL,roundedTL);
        for k = 1:length(universalTL)
            if TLidx(1,k)==0
                indexed(1,k) = NaN;
                continue
            else
                indexed(1,k) = cellinput(TLidx(1,k));
            end
        end
        idxFL(1) = find(~isnan(indexed),1,'first');
        idxFL(2) = find(~isnan(indexed),1,'last');
        filled_indexed = fillmissing(indexed,'linear');
        filled_indexed(1:idxFL(1)) = NaN;
        filled_indexed(idxFL(2):end) = NaN;
        filled_indexed(filled_indexed<0) = 0;
    end
end