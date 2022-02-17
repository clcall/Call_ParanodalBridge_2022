function [TimeSeriesData, optRadxy, optRadz] = Analyze_Time_Series(source_dir, smoothed, useStruct)
%global Compiled_Data
Compiled_Data = struct;
Compiled_Data(:).d00.Sheath_Info.intNames = {};
contents = {source_dir.name};
files = contents(contains(contents,'.traces') | contains(contents,'.xml'));
for i = 1:length(files)
    name = files{i};
    var_name = name(1:3); 
    fprintf('\nCreating xmlstruct for timepoint %s...\n',var_name);
    varname = matlab.lang.makeValidName(name);
    fldr = {source_dir.folder};
    fldr = fldr{1};
    if exist(fullfile(fldr,[varname '.mat']),'file') == 2 && useStruct == 1
        load(fullfile(fldr,[varname '.mat']));
    else
        xmlstruct = parseXML_SingleCell(fullfile(fldr,name));  %changed to singlecell version 20180120 CLC
        save(fullfile(fldr,[varname '.mat']),'xmlstruct');
    end
    [IntersectionDistribution, ProcessLengths, sheathIndex, procIndex] = calcIntersectionsXMLtimepoint(xmlstruct,Compiled_Data,smoothed);
    [Compiled_Data(:).(var_name).Sheath_Info] = deal(IntersectionDistribution);
    [Compiled_Data(:).(var_name).Process_Info] = deal(ProcessLengths);
    [Compiled_Data(:).(var_name).procIndex] = deal(procIndex);
    [Compiled_Data(:).(var_name).sheathIndex] = deal(sheathIndex);
    [Compiled_Data(:).(var_name).xmlstruct] = deal(xmlstruct);
    if contains(var_name,{'d14','d15','d16'})
        figure
        if ispc
            substr = split(fldr,'\');
        else
            substr = split(fldr,'/');
        end
        an = substr{end};
        subplot(2,1,1)
        plotVolumes(an,sheathIndex,xmlstruct,fldr)
        hold on
        [optRadxy, optRadz, coords] = territoryAlgorithm(sheathIndex,xmlstruct);
        Compiled_Data(:).(var_name).sheathCoords = deal(coords);
        [x, y, z] = ellipsoid(0,0,0,optRadxy,optRadxy,optRadz,40);
        surf(x, y, z,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
        camlight('left')
        lighting gouraud
        xlim([-150 150])
        ylim([-150 150])
        zlim([-50 150])
        view([0 90])
        hold off
        subplot(2,1,2)
        plotVolumes(an,sheathIndex,xmlstruct,fldr)
        hold on
        [x, y, z] = ellipsoid(0,0,0,optRadxy,optRadxy,optRadz,40);
        surf(x, y, z,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
        camlight('left')
        lighting gouraud
        xlim([-150 150])
        ylim([-150 150])
        zlim([-150 150])
        view([0 0])
        hold off
        figQuality(gcf,gca,([1.5 3]).*1.5)
    end
    %[Compiled_Data(:).(var_name).DeadProc_Info] = deal(DeadProcessLengths);
end
TimeSeriesData = Compiled_Data;
fprintf('Finished.\n')
end