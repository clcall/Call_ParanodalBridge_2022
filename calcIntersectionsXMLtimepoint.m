function [IntersectionDistribution, ProcessLengths, sheathIndex, procIndex] = calcIntersectionsXMLtimepoint( xmlstruct , Compiled_Data, smoothed )
if nargin == 1
    Compiled_Data.d00.Sheath_Info.intNames = {};
    smoothed = 1;
end
[traceNames,traceLengths,traceColors,traceSWCs,traceIDs,traceStartsOn] = parseData(xmlstruct, smoothed);
traceLengths = cell2mat(traceLengths);

% makes a sorted list of sheath numbers and finds the largest number
% assigned during tracing (to catch missing sheath numbers)
if ~isempty(Compiled_Data.d00.Sheath_Info.intNames)
    tp1flag = 0;
    name = Compiled_Data.d00.Sheath_Info.intNames{end};
    lastSheathNum = str2double(name(2:end));
    lastSheathName = name;
    sheathIndex = find(~contains(traceNames,'Path'));
else %called for first timepoint
    tp1flag = 1;
    sheathIndex = find(~contains(traceNames,'Path'));
    sheathTraceSize = size(sheathIndex,1);
    sheathNames = traceNames(sheathIndex);
    sheathNumVector = [];
    for i = 1:sheathTraceSize
        namelist = sheathNames{i,1};
        sheathNumVector = [sheathNumVector; str2double(namelist(2:3))];
    end
    lastSheathNum = max(sheathNumVector);
%     lastSheathName = ['s', num2str(lastSheathNum, '%02u')];
end

% PARSE CYTOPLASMIC PROCESSES
procIndex = find(~contains(traceNames,'s'));
procNames = traceNames(procIndex);
procLengths = traceLengths(procIndex);
procIDs = traceIDs(procIndex);
procStartsOn = traceStartsOn(procIndex);

%combines fragments of sheaths together and calculates relative location
%of cytoplasmic intersection
searchName = 's01';
counter = 1;
isFound = ~isempty(find(contains(traceNames,searchName), 1));
while ~isFound & tp1flag
    counter = counter + 1;
    searchName = ['s', num2str(counter, '%02u')];
    isFound = ~isempty(find(contains(traceNames,searchName), 1));
end
if tp1flag
    intNames = {};
else
    intNames = table2cell(Compiled_Data.d00.Sheath_Info(:,1))';
end
while isFound & tp1flag & (counter < lastSheathNum+1)
    intSegs_index = find(contains(traceNames,searchName));
    intNames{end+1} = searchName;
    counter = counter + 1;
    if counter > lastSheathNum
        break
    end
    searchName = ['s', num2str(counter, '%02u')];
    isFound = ~isempty(find(contains(traceNames,searchName), 1));
    while ~isFound
        counter = counter + 1;
        if counter > lastSheathNum
            break
        end
        searchName = ['s', num2str(counter, '%02u')];
        isFound = ~isempty(find(contains(traceNames,searchName), 1));
    end
end

% PREALLOCATE VARIABLES
L = length(intNames);
intersections = NaN(L,1);
intLengths = zeros(L,1);
fLnth = NaN(L,1);
cLnth = NaN(L,1);
intIDs = cell(L,1);
intStartsOn = cell(L,1);
node_class = NaN(L,1);
isbridgecnctd = NaN(L,1);
hasmltprocs = NaN(L,1);
pnode1coords = cell(L,1);
pnode2coords = cell(L,1);

i = 1;
while i < L+1
    while isempty(find(contains(traceNames,intNames{i}), 1))
        i = i + 1;
        if i > length(intNames)
            break
        end
    end
    if i > length(intNames)
        break
    end
    intSegs_index = find(contains(traceNames,intNames{i}));
    intIDs{i} = traceIDs(intSegs_index);
    intStartsOn{i} = traceStartsOn(intSegs_index);
    getpnodecoords; % internal function, below, to find coords of sheath ends, to use in vector plot
    %parses the sheath segments and gets the total sheath length
    if any(contains(traceColors(intSegs_index),{'yellow','Yellow'}))
        if size(intSegs_index,1) == 1
            intLengths(i) = traceLengths(intSegs_index);
        elseif any(contains(traceNames(intSegs_index),'a'))
            intLengths(i) = sum(traceLengths(intSegs_index));
        else
            calc_2SegInts;
        end
        % exclude "clipped" sheaths (but keep length, although
        % underestimated)
        %         elseif any(contains(traceColors(intSegs_index),'black'))
        %             calc_2SegInts;
        %             intersections(end) = NaN;
        %             fLnth(end) = NaN;
        %             cLnth(end) = NaN;
        %             hasmltprocs = [hasmltprocs; 0];
        % determine if sheath is connected by >1  cytoplasmic process
    elseif any(contains(traceColors(intSegs_index),{'white','White'}))
        hasmltprocs(i) = 1;
        if size(intSegs_index,1) == 1
            intLengths(i) = traceLengths(intSegs_index);
        else
            calc_2SegInts;
        end
    elseif size(intSegs_index,1) == 1
        intersections(i) = 0;
        intLengths(i) = traceLengths(intSegs_index);
        hasmltprocs(i) = 0;
    elseif length(intSegs_index) > 1 && ...
            any(contains(traceNames(intSegs_index),'a'))
        intLengths(i) = sum(traceLengths(intSegs_index));
        intersections(i) = 0;
        hasmltprocs(i) = 0;
    else
        calc_2SegInts;
        hasmltprocs(i) = 0;
    end
    % determine if sheath segments are heminodes or contribute to nodes
    %1= at least 1 heminode
    %2= at least 1 node
    %3= isolated sheath
    %4= terminal sheath
    %5= continuous sheath
    node_classes = str2double(traceSWCs(intSegs_index));
    if sum(node_classes) == 3
        node_class(i) = 1;
    elseif sum(node_classes) == 4
        node_class(i) = 2;
    elseif sum(node_classes) == 6 || sum(node_classes == 1)
        node_class(i) = 3;
    elseif sum(node_classes) == 7 || sum(node_classes == 2)
        node_class(i) = 4;
    elseif sum(node_classes) == 8 || sum(node_classes == 5)
        node_class(i) = 5;
    else
        node_class(i) = NaN;
    end
    % determine if sheath is connected by paranodal bridge
    if any(contains(traceColors(intSegs_index),{'orange','#ffc800','Orange'}))
        isbridgecnctd(i) = 1;
    else
        isbridgecnctd(i) = 0;
    end
    i = i + 1;
end
intNames = intNames';
    function calc_2SegInts
        c_index = contains(traceNames(intSegs_index),'c');
        c_seg = sum(traceLengths(intSegs_index(c_index)));
        f_index = contains(traceNames(intSegs_index),'f');
        f_seg = sum(traceLengths(intSegs_index(f_index)));
        segs = [c_seg,f_seg];
        shortestSeg = min(segs);
        longestSeg = max(segs);
        fLnth(i) = f_seg;
        cLnth(i) = c_seg;
        intersections(i) = shortestSeg/(shortestSeg + longestSeg)*100;
        intLengths(i) = shortestSeg + longestSeg;
    end
    function getpnodecoords
        coords = [];
        for p = 1:length(intSegs_index)
            coords = [coords; xmlstruct.paths(intSegs_index(p)).points.smoothed(1,:)];
            coords = [coords; xmlstruct.paths(intSegs_index(p)).points.smoothed(end,:)];
        end
        coords = unique(coords,'rows');
        for p = 1:size(coords,1)
            base = coords(p,:);
            for t = 1:size(coords,1)
                test = coords(t,:);
                dists(p,t) = calcEuclid(base,test);
            end
        end
        [~,I] = max(dists(:));
        [row, col] = ind2sub(size(dists),I);
        pnode1coords(i) = {coords(row,:)};
        pnode2coords(i) = {coords(col,:)};
    end
IntersectionDistribution = table(intNames, intersections, intLengths, fLnth, cLnth, node_class, isbridgecnctd, hasmltprocs, pnode1coords, pnode2coords, intIDs, intStartsOn);
ProcessLengths = table(procNames, procLengths, procIDs, procStartsOn);
end


%----Local Function----
function [ traceName,traceLength,traceColor,traceSWC,traceID,traceStartsOn ] = parseData(xmlstruct, smoothed)

%get length of paths list
numPaths = size(xmlstruct.paths,2);
traceName = cell(numPaths,1);
traceLength = cell(numPaths,1);
traceColor = cell(numPaths,1);
traceSWC = cell(numPaths,1);
traceID = cell(numPaths,1);
traceStartsOn = cell(numPaths,1);
for i = 1:numPaths
    traceName{i} = xmlstruct.paths(i).attribs.name;
    if smoothed == 1
        lengthtype = 'reallength_smoothed';
    elseif smoothed == 0
        lengthtype = 'reallength';
    end
    if ~ischar(xmlstruct.paths(i).attribs.(lengthtype))
        traceLength{i} = xmlstruct.paths(i).attribs.(lengthtype);
    else
        traceLength{i} = str2double(xmlstruct.paths(i).attribs.(lengthtype));
    end
    traceColor{i} = xmlstruct.paths(i).attribs.color;
    traceSWC{i} = xmlstruct.paths(i).attribs.swctype;
    traceID{i} = str2double(xmlstruct.paths(i).attribs.id);
    if isfield(xmlstruct.paths(i).attribs,'startson')
        traceStartsOn{i} = str2double(xmlstruct.paths(i).attribs.startson);
    else
        traceStartsOn{i} = {};
    end
end
end
