clear variables
path = 'D:\GitHubRepos\Call_ParanodalBridge_2022\MonkZFtraces\allsheaths';
cd(path);
fldr = dir(path);
files = {fldr(3:end).name};
n = length(files);
alltimeSheaths = NaN(n,160);
alltimeBrgs = NaN(n,160);

colors = parula;
figure(1);
for f = 1:n
    name = files{f};
    xmlstruct = parseXML_SingleCell(name);
    timeline = gettimeline(name(1:12));

    parsedData_sheaths = parseData(xmlstruct);
    frames = cell2mat(parsedData_sheaths(:,2));
    bridges = cell2mat(parsedData_sheaths(:,3));
    lastframe = max(frames);
    numsheaths = NaN(1,lastframe);
    numbrgs = numsheaths;
    for i = 1:lastframe
        numsheaths(i) = sum(frames==i);
        numbrgs(i) = sum(bridges(frames==i));
        if i>1 && ~(numsheaths(i-1)==0) && numsheaths(i)==0
            numsheaths(i) = numsheaths(i-1);
            numbrgs(i) = numbrgs(i-1);
        end
        if i==lastframe && contains(parsedData_sheaths(end,1),'placeholder')
            numsheaths(i) = numsheaths(i-1);
            numbrgs(i) = numbrgs(i-1);
        end
    end
%     figure
%     plot(timeline,numsheaths,'k-o')
%     hold on
%     plot(timeline,numbrgs,'r-o')
%     hold off
    
    % adjusted timeline to peak of sheaths
    figure(1);
    hold on
    peak = max(numsheaths);
    idx = find(numsheaths==peak);
    peakedge = idx(end);
    adjTL = timeline - timeline(peakedge);
    tempcolor = colors(round(256/n*f),:);
    plot(adjTL,numsheaths,'LineStyle','-','Color',tempcolor)
    plot(adjTL,numbrgs,'LineStyle','--','Color',tempcolor)
    hold off
    [~,scaledidx] = ismember(round(adjTL),-29:130);
    alltimeSheaths(f,scaledidx) = numsheaths;
    alltimeBrgs(f,scaledidx) = numbrgs;
end
        
function [ parsedData_sheaths ] = parseData(xmlstruct)
%get length of paths list
numPaths = size(xmlstruct.paths,2);
traceName = cell(numPaths,1);
traceColor = cell(numPaths,1);
traceFrame = cell(numPaths,1);
for i = 1:numPaths
    traceName{i,1} = xmlstruct.paths(i).attribs.name;
    if contains(xmlstruct.paths(i).attribs.color,'Orange')
        traceColor{i,1} = 1; % meaning bridge
    else
        traceColor{i,1} = 0;
    end
    traceFrame{i,1} = str2double(xmlstruct.paths(i).attribs.frame);
end

%PARSE TRACE INFO
%sheath extraction
parsedData_sheaths = [traceName, traceFrame, traceColor];
end