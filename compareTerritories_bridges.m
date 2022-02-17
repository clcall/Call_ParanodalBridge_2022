function [optRadxy,optRadz] = compareTerritories_bridges(xmlstruct,cond,inclBridge,doplot)
[sheathIndex_noBrg,sheathIdxBrg] = getXMLidx_nobridge(xmlstruct);
if inclBridge
    sheathIndex_all = [sheathIndex_noBrg; sheathIdxBrg];
else
    sheathIndex_all = sheathIndex_noBrg;
end
if doplot
    figure
    subplot(2,1,1)
    plotVolumes([],sheathIndex_noBrg,xmlstruct,'notbridge')
    hold on
    if inclBridge
        q = [];
        for k = 1:length(sheathIdxBrg)
            q = [q; xmlstruct.paths(sheathIdxBrg(k)).points.smoothed];
        end
        c = mean(q,'omitnan');
        q = q-c;
        plot3(q(:,1),q(:,2),q(:,3),'.','Color',[0 85 212]./255)
    end
end
[optRadxy, optRadz, ~] = territoryAlgorithm(sheathIndex_all,xmlstruct);
if doplot
    [x, y, z] = ellipsoid(0,0,0,optRadxy,optRadxy,optRadz,40);
    surf(x, y, z,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none')
    camlight('left')
    lighting gouraud
    xlim([-150 150])
    ylim([-150 150])
    zlim([-50 150])
    view([0 90])
    hold off
    subplot(2,1,2)
    plotVolumes([],sheathIndex_noBrg,xmlstruct,'notbridge')
    hold on
    if inclBridge
        q = [];
        for k = 1:length(sheathIdxBrg)
            q = [q; xmlstruct.paths(sheathIdxBrg(k)).points.smoothed];
        end
        c = mean(q,'omitnan');
        q = q-c;
        plot3(q(:,1),q(:,2),q(:,3),'.','Color',[0 85 212]./255)
    end
    [x, y, z] = ellipsoid(0,0,0,optRadxy,optRadxy,optRadz,40);
    surf(x, y, z,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none')
    camlight('left')
    lighting gouraud
    xlim([-150 150])
    ylim([-150 150])
    zlim([-150 150])
    view([0 0])
    hold off
    figQuality(gcf,gca,([1.5 3]).*1.5)
end
    function [sheathIndex,sheathIdxBrg] = getXMLidx_nobridge(xmlstruct)
        [traceNames,~,traceColors,~,~,~] = parseData(xmlstruct, 1);
        sheathIndex = ~contains(traceNames,'Path');
        orgidx = contains(traceColors,{'orange','#ffc800','Orange'});
        tempidx = sheathIndex & ~orgidx;
        sheathIndex = find(tempidx);
        sheathIdxBrg = find(orgidx);
    end

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
end