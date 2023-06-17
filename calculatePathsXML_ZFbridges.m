function [d,X,Delta,borders,origin,firstbrgtp] = calculatePathsXML_ZFbridges( xmlstruct, useImaginaryEdge, makeplot)
[traceData_sheaths, traceData_lifeact] = parseData(xmlstruct);
frames_sheaths = cell2mat(traceData_sheaths(:,2));
frames_lifeact = cell2mat(traceData_lifeact(:,2));
d = cell(1,max(frames_sheaths));
X = cell(1,max(frames_sheaths));
Delta = cell(1,max(frames_sheaths));
f1 = figure;
%f2 = figure;
firstbrg = 0;
for i=1:max(frames_sheaths)
    temp = traceData_sheaths(frames_sheaths==i,:);
    if isempty(temp)
        subplot(max(frames_sheaths),1,i)
        ylim([0,1.8])
        xlim([0,1000])
        xticks([])
        xticklabels([])
        axis off
        box off
        continue
    end
    avgxyz = cell2mat(cellfun(@mean,temp(:,6),'UniformOutput',false));
    avgx = avgxyz(:,1);
    [~,I] = sort(avgx,'ascend');
    nfrag = size(temp,1);
    d{i} = NaN(2000,2);
    lnth1 = 1;
    borders = [];
    origin = 0;
    for k = 1:nfrag
        sheathedge = cell2mat(temp(I(k),6));
        if sheathedge(1,1) > sheathedge(end,1)
            sheathedge = flip(sheathedge);
        end
        if cell2mat(temp(I(k),5)) == 1 % is a bridge
            b=1;
            if ~firstbrg
                firstbrgtp = i;
                firstbrg = 1;
            end
            % make straight line always on same edge as sheath traces
            if nargin>1 && useImaginaryEdge
                allpts = cell2mat(temp(:,6));
                startI = knnsearch(allpts,sheathedge(1,:));
                stopI = knnsearch(allpts,sheathedge(end,:));
                start = allpts(startI,:);
                stop = allpts(stopI,:);
                sheathedge = [linspace(start(1),stop(1),size(sheathedge,1))' linspace(start(2),stop(2),size(sheathedge,1))' linspace(start(3),stop(3),size(sheathedge,1))'];
            end
        elseif cell2mat(temp(I(k),5)) == 2 % has sheath origin
            b=2;
        else
            b=0;
        end
        lifeact = traceData_lifeact{frames_lifeact==i,6};
        [~,D] = knnsearch(lifeact,sheathedge);
        lnth0 = lnth1;
        lnth1 = lnth0 + length(D) -1;
        d{i}(lnth0:lnth1,1) = D;
        d{i}(lnth0:lnth1,2) = b;
        if b==1
            borders = [borders,lnth0,lnth1];
        end
        if b==2
            origin = lnth0;
        end
    end
    idx = isnan(d{i}(:,1));
    d{i}(idx,:)=[];
    if origin==0
        X{i} = -lnth1:-1;
    else
        X{i} = -origin:(lnth1-origin-1);
    end
    borders(borders<origin) = -(origin-borders(borders<origin));
    borders(borders>origin) = (borders(borders>origin)-origin);
    
    %% PLOT RECENT CHANGE
    if i>1 && any(X{i-1})
        idx1 = ismember(X{i-1},X{i}); % length of i-1
        idx2 = ismember(X{i},X{i-1}); % length of i
        delta = d{i}(idx2,1)-d{i-1}(idx1,1);
        Delta{i} = [X{i-1}(idx1)',delta];
    end
    %% PLOT ABSOLUTE DISTANCE
    if makeplot
        figure(f1);
        hold on
        subplot(max(frames_sheaths),1,i)
        %         h = plot(X,ones(size(d,1),1),'LineWidth',2);
        h = plot(X{i},d{i}(:,1),'LineWidth',2);
        
        cdo = colormap('parula');
        cd = interp1(linspace(0,1.2,length(cdo)),cdo,d{i}(:,1));
        cd = uint8(cd'*255);
        logdx = all(cd==0);
        if any(logdx)
            cd(1,logdx) = uint8(cdo(end,1)*255);
            cd(2,logdx) = uint8(cdo(end,2)*255);
            cd(3,logdx) = uint8(cdo(end,3)*255);
        end
        cd(4,:) = 255;
        drawnow
        set(h.Edge,'ColorBinding','interpolated','ColorData',cd)
        
        ylim([0,2])
        xlim([-1500,1500])
        xticks([])
        xticklabels([])
        axis off
        box off
        
        if any(borders)
            xline(borders,'Linewidth',1.5);
        end
        xline(0,'Linewidth',1.5,'Color','r');
        hold off
    end
end
end

%% ----Local Function----
function [ parsedData_sheaths, parsedData_lifeact ] = parseData(xmlstruct)
%get length of paths list
numPaths = size(xmlstruct.paths,2);
traceName = cell(numPaths,1);
traceLength = cell(numPaths,1);
traceSWC = cell(numPaths,1);
traceColor = cell(numPaths,1);
traceFrame = cell(numPaths,1);
traceChannel = cell(numPaths,1);
tracePaths = cell(numPaths,1);
for i = 1:numPaths
    traceName{i,1} = xmlstruct.paths(i).attribs.name;
    traceLength{i,1} = xmlstruct.paths(i).attribs.reallength_smoothed;
    traceSWC{i,1} = str2double(xmlstruct.paths(i).attribs.swctype);
    if contains(xmlstruct.paths(i).attribs.color,'Orange')
        traceColor{i,1} = 1; % meaning bridge
    elseif contains(xmlstruct.paths(i).attribs.color,'Red')
        traceColor{i,1} = 2; % meaning sheath frag with first point = process intersection (sheath origin)
    else
        traceColor{i,1} = 0;
    end
    traceFrame{i,1} = str2double(xmlstruct.paths(i).attribs.frame);
    traceChannel{i,1} = xmlstruct.paths(i).attribs.channel;
    tracePaths{i,1} = xmlstruct.paths(i).points.smoothed;
end
%PARSE TRACE INFO
%sheath extraction
index = find(contains(traceChannel,'1')); % Ch1 - binarized myrEGFP, Ch2 - myrEGFP, Ch3 - lifeact
Names = traceName(index);
Lengths = traceLength(index);
SWCs = traceSWC(index);
Colors = traceColor(index);
Paths = tracePaths(index);
Channels = traceChannel(index);
Frames = traceFrame(index);

%create array to send main function name & length info
parsedData_sheaths = [Names, Frames, Lengths, SWCs, Colors, Paths];

index = find(contains(traceChannel,'3')); % Ch1 - binarized myrEGFP, Ch2 - myrEGFP, Ch3 - lifeact
Names = traceName(index);
Lengths = traceLength(index);
SWCs = traceSWC(index);
Colors = traceColor(index);
Paths = tracePaths(index);
Channels = traceChannel(index);
Frames = traceFrame(index);

%create array to send main function name & length info
parsedData_lifeact = [Names, Frames, Lengths, SWCs, Colors, Paths];

end