clear variables
path = 'D:\GitHubRepos\Call_ParanodalBridge_2022\MonkZFtraces\allsheaths';
cd(path);
fldr = dir(path);
files = {fldr(3:end).name};
n = length(files);
alltimeSheaths = NaN(n,160);
alltimeBrgs = NaN(n,160);
colors = parula;
SD = struct; %all sheath data

figure(1);
for f = 1:n
    name = files{f};
    SD(f).name = name(1:12);
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
    
    % adjusted timeline to peak of sheaths
    figure(1);
    hold on
    peak = max(numsheaths);
    idx = find(numsheaths==peak);
    peakedge = idx(end);
    SD(f).peakedge = peakedge;
    adjTL = timeline - timeline(peakedge);
    SD(f).adjTL = adjTL;
    tempcolor = colors(round(256/n*f),:);
    plot(adjTL,numsheaths,'LineStyle','-','Color',tempcolor,'LineWidth',2)
    plot(adjTL,numbrgs,'LineStyle',':','Color',tempcolor,'LineWidth',2)
    hold off
    [~,scaledidx] = ismember(round(adjTL),-29:130);
    alltimeSheaths(f,scaledidx) = numsheaths;
    alltimeBrgs(f,scaledidx) = numbrgs;
end

temp = ~isnan(alltimeSheaths);
idx = cell2mat(arrayfun(@(x) find(temp(x,:),1,'last'),1:size(alltimeSheaths,1),'UniformOutput',false));
filledSheaths = fillmissing(alltimeSheaths(:,30:end),'linear',2);
filledBrgs = fillmissing(alltimeBrgs(:,30:end),'linear',2);
for i=1:length(idx)
    if idx(i)<110
        filledSheaths(i,idx(i)-29:81)=NaN;
        filledBrgs(i,idx(i)-29:81)=NaN;
    end
end

figure(1);
hold on
plot(0:80,mean(filledSheaths(:,1:81),'omitnan'),'LineStyle','-','Color','k','LineWidth',2)
plot(0:80,mean(filledBrgs(:,1:81),'omitnan'),'LineStyle',':','Color','k','LineWidth',2)
figQuality(gcf,gca,[4.6,2.5])
xlim([-40 80])

        
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