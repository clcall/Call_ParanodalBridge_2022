function datatable = calculatePathsXML_PVpnb( xmlstruct, smoothed, makeplots )
[parsedAxons,parsedSheaths,parsedPNBs] = parseData(xmlstruct,smoothed,makeplots);
axonNames = parsedAxons.axonNames;
axon_num = max(parsedAxons.axonID); 
totalAxonLength = [];
pnbBranchRatio = [];
pnbSheathRatio = [];
avgNonBrgSheathLnth = [];
avgBrgSheathLnth = [];
totalSheathLength = [];
avgInternodeLength = [];
axonName = [];
numBridges = [];
allSheathLengths = {};
chainLnths = {};
avgChainLnth = [];

for i = 1:axon_num
    axon_idx = parsedAxons.axonID==i;
    currAxData = parsedAxons(axon_idx,:);
    totalAxonLength = [totalAxonLength; sum(currAxData.axonLengths)];
    axonName = [axonName; i];
    
    %EXTRACT BRANCHES
    branchCoords = cell2mat(currAxData.axonStarts);
    if size(branchCoords,1) < (size(currAxData,1) - 1)
        error('At least one branch of axon %s is disconnected.\n',currAxData.axonNames{1}(1:3));
    end
    branchCoords = unique(branchCoords,'rows');
    %GET BRANCH BRIDGES
    pnb_idx = parsedPNBs.pnbAxon==i;
    pnbs = parsedPNBs(pnb_idx,:);
    [~,dists] = knnsearch(branchCoords,cell2mat(pnbs.pnbEndCoords(:,1)),'K',1);
    pnbBranches = dists < 5; 
    numBridges = [numBridges; size(pnbBranches,1)];
    pnbBranchRatio = [pnbBranchRatio; sum(pnbBranches)/size(pnbBranches,1)];
    
    %OVERALL SHEATH DIMENSIONS
    sheath_idx = parsedSheaths.sheathAxon==i;
    internodes = parsedSheaths(sheath_idx,:);
    ids = internodes.sheathID;
    internodes2 = internodes(1,:);
    for j = 1:max(internodes.sheathID)
        temp = internodes(internodes.sheathID==j,:);
        if size(temp,1) == 1
            internodes2(j,:) = temp;
        else
            internodes2(j,1:3) = temp(1,1:3);
            sheathLength = sum(temp.sheathLengths);
            tbl = table(sheathLength);
            internodes2(j,4) = tbl;
            internodes2(j,5:6) = temp(1,5:6);
            tbl = table({temp.sheathEndCoords{1,2} temp.sheathEndCoords{2,2}});
            internodes2(j,7) = tbl;
        end
    end
    
    %PARSE SHEATH CHAIN RELATIONSHIPS
    sheathLengths = internodes2.sheathLengths;
    pnbStarts = cell2mat(pnbs.pnbStarts);
    pnbEnds = cell2mat(pnbs.pnbEndCoords);
    sheathEndsA = cell2mat(internodes2.sheathEndCoords(:,1));
    sheathEndsB = cell2mat(internodes2.sheathEndCoords(:,2));
    flags = [];
    tempChainLnths = [];
    bridgeChainLnth = [];
    for j = 1:size(internodes2,1)
        if ismember(j,flags) %keep track of connected sheaths accounted for so as not to double count
            continue
        end
        bridgeChainLnth = [];
        idx = ismember(pnbStarts, [sheathEndsA(j,:);sheathEndsB(j,:)], 'rows'); %one specific sheath 'j'
        if any(idx)
            idxA = ismember(sheathEndsA,pnbEnds(idx,:),'rows'); %is there a second bridge
            idxB = ismember(sheathEndsB,pnbEnds(idx,:),'rows');
            flags = [flags; j; find(idxA + idxB)]; %keep track of connected sheaths accounted for so as not to double count
            
            bridgeChainLnth = sheathLengths(j) + sum(sheathLengths(find(idxA + idxB)));  %#ok<FNDSB>
            idxAB = [ismember(pnbStarts,sheathEndsA(idxA,:),'rows'), ismember(pnbStarts,sheathEndsB(idxB,:),'rows')]; %check for third bridge
            while any(idxAB)
                idxAB = idxAB(:,1) + idxAB(:,2);
                idxA = ismember(sheathEndsA,pnbEnds(idxAB),'rows');
                idxB = ismember(sheathEndsB,pnbEnds(idxAB),'rows');
                flags = [flags; j; find(idxA + idxB)];
                
                bridgeChainLnth = sheathLengths(j) + sum(sheathLengths(find(idxA + idxB)));  %#ok<FNDSB>
                idxAB = [ismember(pnbStarts,sheathEndsA(idxA,:),'rows'), ismember(pnbStarts,sheathEndsB(idxB,:),'rows')];
            end
        end
        tempChainLnths = [tempChainLnths; bridgeChainLnth];
    end
    chainLnths{i,1} = tempChainLnths;
    avgChainLnth = [avgChainLnth; mean(tempChainLnths)];
    totalSheathLength = [totalSheathLength; sum(internodes2.sheathLengths)];
    avgInternodeLength = [avgInternodeLength; mean(internodes2.sheathLengths)];
    allSheathLengths{i,1} = internodes2.sheathLengths;
    pnbSheathRatio = [pnbSheathRatio; size(pnbBranches,1)/size(internodes2,1)];
    
    bridgedSheaths = internodes2(contains(internodes2.sheathColors,'Orange'),:);
    avgBrgSheathLnth = [avgBrgSheathLnth; mean(bridgedSheaths.sheathLengths)];
    nonBrgSheaths = internodes2(~contains(internodes2.sheathColors,'Orange'),:);
    avgNonBrgSheathLnth = [avgNonBrgSheathLnth; mean(nonBrgSheaths.sheathLengths)];
    clear internodes2
end
percentMyelin = (totalSheathLength ./ totalAxonLength)  * 100;
datatable = table(axonName,totalAxonLength,avgInternodeLength,...
                      totalSheathLength,percentMyelin,allSheathLengths,avgBrgSheathLnth,...
                      avgNonBrgSheathLnth,numBridges,pnbSheathRatio,pnbBranchRatio,chainLnths,avgChainLnth);
end

%----Local Function----
function [parsedAxons, parsedSheaths, parsedPNBs] = parseData(xmlstruct,smoothed,makeplots)
%get length of paths list
numPaths = size(xmlstruct.paths,2);
traceName = cell(numPaths,1);
% traceLength = cell(numPaths,1);
% traceSWC = cell(numPaths,1);
traceColor = cell(numPaths,1);
traceEndCoords = cell(numPaths,2);
for i = 1:numPaths
    traceName{i,1} = xmlstruct.paths(i).attribs.name;
    if smoothed == 1
        lengthtype = 'reallength_smoothed';
    elseif smoothed == 0
        lengthtype = 'reallength';
    end
    if ~ischar(xmlstruct.paths(i).attribs.(lengthtype))
        traceLength(i,1) = xmlstruct.paths(i).attribs.(lengthtype);
    else
        traceLength(i,1) = str2double(xmlstruct.paths(i).attribs.(lengthtype));
    end
    traceSWC(i,1) = str2double(xmlstruct.paths(i).attribs.swctype);
    traceColor{i,1} = xmlstruct.paths(i).attribs.color;
    traceEndCoords{i,1} = [xmlstruct.paths(i).points.xd(1,:), xmlstruct.paths(i).points.yd(1,:), xmlstruct.paths(i).points.zd(1,:)];
    traceEndCoords{i,2} = [xmlstruct.paths(i).points.xd(end,:), xmlstruct.paths(i).points.yd(end,:), xmlstruct.paths(i).points.zd(end,:)];
    if isfield(xmlstruct.paths(i).attribs,'startx')
        traceStart{i,1} = [str2double(xmlstruct.paths(i).attribs.startx),...
                            str2double(xmlstruct.paths(i).attribs.starty),...
                            str2double(xmlstruct.paths(i).attribs.startz)];
    else
        traceStart{i,1} = [];
    end
    if makeplots
        if contains(traceColor{i,1},'Yellow')
            color = [1 0 0];
            mrkr = 16;
            bubr = 3;
            buba = 0.6;
        elseif contains(traceColor{i,1},'Red')
            color = [0 0 0];
            mrkr = 2;
            bubr = 1;
            buba = 1;
        elseif contains(traceColor{i,1},{'Cyan','White','Blue'})
            color = [0 1 1];
            mrkr = 16;
            bubr = 3;
            buba = 0.1;
        elseif contains(traceColor{i,1},'Orange')
            color = [1 0.5 0];
            mrkr = 16;
            bubr = 3;
            buba = 0.1;
        end
        if contains(traceName{i,1},'a01')
            figure(1)
            hold on
            bubbleplot3(xmlstruct.paths(i).points.xd,xmlstruct.paths(i).points.yd,-1.*xmlstruct.paths(i).points.zd,...
                ones(size(xmlstruct.paths(i).points.xd)).*bubr,...
                repmat(color,size(xmlstruct.paths(i).points.xd)),...
                buba);
        elseif contains(traceName{i,1},'a02')
            figure(2)
            hold on
            bubbleplot3(xmlstruct.paths(i).points.xd,xmlstruct.paths(i).points.yd,-1.*xmlstruct.paths(i).points.zd,...
                ones(size(xmlstruct.paths(i).points.xd)).*bubr,...
                repmat(color,size(xmlstruct.paths(i).points.xd)),...
                buba);
        elseif contains(traceName{i,1},'a03')
            figure(3)
            hold on
            bubbleplot3(xmlstruct.paths(i).points.xd,xmlstruct.paths(i).points.yd,-1.*xmlstruct.paths(i).points.zd,...
                ones(size(xmlstruct.paths(i).points.xd)).*bubr,...
                repmat(color,size(xmlstruct.paths(i).points.xd)),...
                buba);
        end
    end
end
hold off
%PARSE TRACE INFO
%axon extraction
axon_index = cellfun(@length,traceName) < 5;
axonNames = traceName(axon_index);
axonLengths = traceLength(axon_index);
axonSWCs = traceSWC(axon_index);
axonColors = traceColor(axon_index);
axonStarts = traceStart(axon_index,:);
%parse out axon number - name format is a##s##[c/f]
for i = 1:length(axonNames)
    s = axonNames{i,1};
    axonID(i,1) = str2double(s(2:3));
end
%sheath extraction
sheath_index = find(contains(traceName,'s')); %finds all indices where 's' for sheath is in the name
sheathNames = traceName(sheath_index);
sheathLengths = traceLength(sheath_index);
sheathSWCs = traceSWC(sheath_index);
sheathColors = traceColor(sheath_index);
sheathEndCoords = traceEndCoords(sheath_index,:);
sheath_indexPos = cell2mat(strfind(traceName,'s'));
for i = 1:length(sheathNames)
    if ~(sheath_indexPos(i) == 4)
        error('Name formatting for %s is invalid.\n',sheathNames{i});
    end   
    s = sheathNames{i};
    sheathAxon(i,1) = str2double(s(2:3));
    sheathID(i,1) = str2double(s(5:6));
end
%paranodal bridge extraction
pnb_index = find(contains(traceColor,'Yellow'));
pnbNames = traceName(pnb_index);
for i = 1:length(pnbNames)
    if ~(sheath_indexPos(i) == 4)
        error('Name formatting for %s is invalid.\n',sheathNames{i});
    end
    s = pnbNames{i};
    pnbAxon(i,1) = str2double(s(2:3));
end
pnbEndCoords = traceEndCoords(pnb_index,end);
pnbStarts = traceStart(pnb_index,:);
%create array to send main function name & length info
parsedAxons = table(axonNames, axonID, axonLengths, axonSWCs, axonColors, axonStarts);
parsedSheaths = table(sheathNames, sheathAxon, sheathID, sheathLengths, sheathSWCs, sheathColors, sheathEndCoords);
parsedPNBs = table(pnbNames, pnbAxon, pnbStarts, pnbEndCoords);
end