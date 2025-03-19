path = 'D:\GitHubRepos\Call_ParanodalBridge_2022\MonkZFtraces\allsheaths_24dpf';
cd(path);
fldr = dir(path);
files = {fldr(3:end).name};
n = length(files);
alltimeSheaths = NaN(n,160);
alltimeBrgs = NaN(n,160);
colors = parula;
SD24 = struct; %all sheath data

for f = 1:n
    name = files{f};
    SD24(f).name = name(1:12);
    xmlstruct = parseXML_SingleCell(name);

    parsedData_sheaths = parseData(xmlstruct);
    SD24(f).totsheaths = length(cell2mat(parsedData_sheaths(:,3)));
    SD24(f).numbridges = sum(cell2mat(parsedData_sheaths(:,3)));
    SD24(f).prop = SD24(f).numbridges./SD24(f).totsheaths;
end


figure;
hold on
plotSpread([SD24.prop]','distributionMarker','o','distributionColors','k');
figQuality(gcf,gca,[2,2.5])
ylim([0 1])
        
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