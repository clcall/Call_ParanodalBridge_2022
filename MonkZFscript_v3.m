% TO DO:

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
    makeplot = 1;
    [d,X,Delta,borders,origin,firstbrgtp,framesUsed] = calculatePathsXML_ZFbridges_v2(xmlstruct,edgepos,makeplot,TL,name);
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
    for i = 1:length(Delta)
        if any(Delta{i})
            Dadj(ismember(rng,Delta{i}(:,1))) = movmean(Delta{i}(:,2),50);
            Dall = [Dall;Dadj];
        end
    end    
    % calculate extension and keep track of wraps up to then when extended
    wraps = abs(Dall);
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
        
        % get wraps for all bridged sheaths
        idx = ~isnan(Dall(i,:));
        idx(idx1) = 0;
        idx(idx2) = 0;
        idx(leftbord:rightbord) = 0;
        allnonbrg{f}(i) = mean(sum(wraps(1:i,idx),'omitnan'));
    end
end

figure
hold on
for f = 1:size(files,1)
    plot(allTL{f},smallestanch{f},'b-');
    plot(allTL{f},cell2mat(cellfun(@mean,avgbrgwraps{f},"UniformOutput",false)),'r-');
    plot(allTL{f},allnonbrg{f},'g-');
end
hold off
ylabel('Cumulative movement')
xlabel('Time')
legend('anchoring sheath','bridges','all non-anch sheath')
figQuality(gcf,gca,[4,2.5])
%%
alltimeR = NaN(size(files,1),160);
tbl = table;
alltimes = [];
allR = [];
for f = 1:size(files,1)
    R = cell2mat(cellfun(@mean,avgbrgwraps{f},"UniformOutput",false)) ./ smallestanch{f};
    [~,scaledidx] = ismember(round(allTL{f}),-29:130);
    alltimeR(f,scaledidx) = R;
    [a,b] = forceConcat(alltimes,scaledidx',1);
    alltimes = [a b];
    [a,b] = forceConcat(allR,R',1);
    allR = [a b];
    
    tbl = [tbl; table({files(f).name(1:12)},{scaledidx},smallestanch(f),allnonbrg(f),{cell2mat(cellfun(@mean,avgbrgwraps{f},"UniformOutput",false))},{R})];
    
end






