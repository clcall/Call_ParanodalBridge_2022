function SD = MonkZFscript_v4_alignByBridgeIncidence(makeplot,makeplot2)
% NOTE: run MonkZFscript_allsheaths.m FIRST to get SD struct with timelines
SD = MonkZFscript_allsheaths(1);
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
    sheathcolors{f} = SD(idx).cellcolor;
    allTL{f} = TL;
    
    % get path data and Delta distance between sheath & lifeact
    [d,X,Delta,bordersOri,origin,firstbrgtp,framesUsed] = calculatePathsXML_ZFbridges_v2(xmlstruct,edgepos,makeplot,makeplot2,TL,name);

    %% PARSE BRIDGED REGIONS & TIMES
    % get bridged times idx
    brgtime = cell2mat(cellfun(@any,bordersOri,'UniformOutput',false));
    % find consecutive time with bridges
    tmp = find(diff([NaN brgtime NaN] == 1));
    brgtimeidx = reshape(tmp,2,[])'; % reshape in desired form
    brgtimeidx = [brgtimeidx(:,1) brgtimeidx(:,2)-1]; % subtract 1 from second column
    
    %for now, just take longest consecutive timeframe
    [~,row] = max(brgtimeidx(:,2)-brgtimeidx(:,1));
    brgTL{f} = allTL{f}(brgtimeidx(row,1):brgtimeidx(row,2));
    brgtimevect = brgtimeidx(row,1):brgtimeidx(row,2);
    brgDelta = Delta(brgtimevect); % separate Delta for only brgtime
    brgOrigin = origin(brgtimevect); % separate origin for only brgtime
    
    %% smooth brgDelta and concatenate across timepoints
    starts = NaN(size(brgDelta));
    ends = NaN(size(brgDelta));
    for i = 1:length(brgDelta)
        if ~isempty(brgDelta{i})
            starts(i) = min(brgDelta{i}(:,1));
            ends(i) = max(brgDelta{i}(:,1));
        else
            starts(i) = [];
            ends(i) = [];
        end
    end
    [rngstart,sI] = min(starts);
    [rngend,eI] = max(ends);
    globalorigin = brgOrigin{sI};
    rng = rngstart:rngend;
    
    Dadj = NaN(size(rng));
    Dall = [];
    borders_copy = bordersOri(brgtimevect);
    borders = {};
    for i = 1:length(brgDelta)
        if any(brgDelta{i})
            Dadj(ismember(rng,brgDelta{i}(:,1))) = movmean(brgDelta{i}(:,2),50);
            Dall = [Dall;Dadj];
            borders = [borders, borders_copy(i)];
        else
            brgTL{f}(i) = [];
        end
    end
    
    %% calculate extension 
    wraps = abs(Dall);
    allWraps{f} = wraps;
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

%%
elapsedtime = {};
% figure
% hold on
for f = 2:size(files,1) %excluding sheath 1, with only 2 tps
%     anch = smallestanch{f};
%     brg = cell2mat(cellfun(@mean,avgbrgwraps{f},"UniformOutput",false));
%     non = allnonbrg{f};
    elapsedtime{f} = brgTL{f}-brgTL{f}(1);
%     
%     plot(elapsedtime{f},anch-anch(2),'b-');
%     plot(elapsedtime{f},brg-brg(2),'r-');
%     plot(elapsedtime{f},non-non(2),'g-');
end
% hold off
% ylabel('Cumulative movement')
% xlabel('Time')
% legend('anchoring sheath','bridges','all non-anch sheath')
% figQuality(gcf,gca,[4,2.5])

%%
anch_FI = [];
brg_FI = [];
non_FI = [];
universalTL = 0:80;
roundedTL = cellfun(@round,elapsedtime,'UniformOutput',false);

for f = 2:size(files,1) % excluding sheath 1, with only 2 tps
    anch_FI(f-1,:) = averageSheathTLs(smallestanch{f},universalTL,roundedTL{f});
    brg_FI(f-1,:) = averageSheathTLs(cell2mat(cellfun(@mean,avgbrgwraps{f},"UniformOutput",false)),universalTL,roundedTL{f});
    non_FI(f-1,:) = averageSheathTLs(allnonbrg{f},universalTL,roundedTL{f});
    colors(f-1,:) = sheathcolors{f};
end

figure
hold on
avg_anch = mean(anch_FI,1,'omitnan');
avg_brg = mean(brg_FI,1,'omitnan');
avg_non = mean(non_FI,1,'omitnan');
sem_anch = calcSEM(anch_FI,1);
sem_brg = calcSEM(brg_FI,1);
sem_non = calcSEM(non_FI,1);

shadedErrorBar(universalTL(1:30),avg_anch(1:30),sem_anch(1:30),'lineProps',{'Color',[52 75 160]./255});
shadedErrorBar(universalTL(1:30),avg_brg(1:30),sem_brg(1:30),'lineProps',{'r-'});
shadedErrorBar(universalTL(1:30),avg_non(1:30),sem_non(1:30),'lineProps',{'Color',[255 128 0]./255});
hold off
figQuality(gcf,gca,[4 3])

slopes = [];
figure 
hold on
for i = 1:size(anch_FI,1)
    if any(anch_FI(i,1:30))
        useableLI = ~isnan(anch_FI(i,1:30));
        a_pf = polyfit(universalTL(useableLI),anch_FI(i,useableLI),1);
        b_pf = polyfit(universalTL(useableLI),brg_FI(i,useableLI),1);
        n_pf = polyfit(universalTL(useableLI),non_FI(i,useableLI),1);
        slopes = [slopes; a_pf(1) b_pf(1) n_pf(1)];
        plot([a_pf(1); n_pf(1)],'Color',colors(i,:),'LineStyle','-','Marker','o','LineWidth',1)
    end
end
hold off
ylim([-0.05 0.25])
xlim([0 3])
xticks([1 2])
xticklabels({})
figQuality(gcf,gca,[2 3])
        

% LME
Y = [];
T = [];
shthID = [];
cellID = [];
group = [];
for i = 2:length(smallestanch)
    Y = [Y; [smallestanch{i} cell2mat(cellfun(@mean,avgbrgwraps{i},"UniformOutput",false)) allnonbrg{i}]'];
    T = [T; [roundedTL{i} roundedTL{i} roundedTL{i}]'];
    shthID = [shthID; repmat({num2str(i)},length(smallestanch{i})*3,1)];
    cellID = [cellID; repmat({files(i).name(1:12)},length(smallestanch{i})*3,1)];
    group = [group; [repmat({'anch'},size(smallestanch{i})) repmat({'brg'},size(avgbrgwraps{i})) repmat({'non'},size(allnonbrg{i}))]'];
end

TBL = table(Y,T,shthID,cellID,group);
lme_reml = fitlme(TBL,'Y ~ T*group + (1|cellID) + (1|cellID:shthID)','FitMethod','REML');
disp(lme_reml)
end
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