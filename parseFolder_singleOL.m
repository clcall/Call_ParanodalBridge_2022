% edited 20180531 CLC - added analysis for per-timepoint extensions and
% retractions (avg and num)
function ALL = parseFolder_singleOL(cond_dir, ALL, smoothed, useStruct)
    %   Parse out time series data and assign to struct for each animal, within
    %either cup/ctrl. Then, take that data and analyze length changes and
    %assign within same animal name in the ALL struct.
    sub_folders = {cond_dir.name};
    sub_folders = sub_folders(~contains(sub_folders,'.') & ~contains(sub_folders,'..'));
    dirstr = cond_dir.folder;
    cond_name = dirstr(end-3:end);
    % PREALLOCATE VARIABLES:
    ALL.stats.avgLnthChngPerTP.(cond_name) = []; 
    ALL.stats.extLnthPerTP.(cond_name) = [];
    ALL.stats.extNumPerTP.(cond_name) = [];
    ALL.stats.retLnthPerTP.(cond_name) = [];
    ALL.stats.retNumPerTP.(cond_name) = [];
    ALL.stats.bridges.(cond_name) = [];
    ALL.stats.terrRadii.(cond_name) = [];
    % PARSE FOLDERS
    for k = 1:length(sub_folders)
        cell_name = sub_folders{k}(1:11);
        fprintf('\nParsing folder for cell %s...\n',cell_name);
        if ispc
            dirtemp = [cd,'\',cond_name,'\',sub_folders{k}];
        else
            dirtemp = [cd, '/',cond_name,'/',sub_folders{k}];
        end
        if ~contains(cell_name,{'a609_1_r1c6','a609_1_r1c4','a010_0_r2c5'})
            [ALL.(cond_name).(cell_name).raw, optRadxy, optRadz] = Analyze_Time_Series(dir(dirtemp),smoothed, useStruct);
            ALL.stats.terrRadii.(cond_name) = [ALL.stats.terrRadii.(cond_name); optRadxy, optRadz];
        else
            [ALL.(cond_name).(cell_name).raw] = Analyze_Time_Series(dir(dirtemp),smoothed, useStruct);
        end
        [ALL.(cond_name).(cell_name).sheath_LengthDiffs,...
            ALL.(cond_name).(cell_name).c_LengthDiffs,...
            ALL.(cond_name).(cell_name).f_LengthDiffs,...
            ALL.(cond_name).(cell_name).cfRatioPerTP,...
            ALL.(cond_name).(cell_name).Total_SheathLengths,...
            ALL.(cond_name).(cell_name).proc_LengthDiffs,...
            ALL.(cond_name).(cell_name).sheathLengthPerTP,...
            ALL.(cond_name).(cell_name).nodeClassPerTP] = calcLengthChanges_old(ALL.(cond_name).(cell_name).raw);
        ALL.(cond_name).(cell_name).num_BslineSheaths = size(ALL.(cond_name).(cell_name).raw.d00.Sheath_Info,1);
        tps = fieldnames(ALL.(cond_name).(cell_name).raw);
        ALL.(cond_name).(cell_name).num_FinalSheaths = nnz(ALL.(cond_name).(cell_name).raw.(tps{end}).Sheath_Info.intLengths);
        ALL.(cond_name).(cell_name).mean_BslineSheathLength = ALL.(cond_name).(cell_name).Total_SheathLengths(1) ./ ALL.(cond_name).(cell_name).num_BslineSheaths;
        if isempty(ALL.(cond_name).(cell_name).raw.d00.Sheath_Info.intNames)
            ALL.stats.bridges.(cond_name) = [ALL.stats.bridges.(cond_name);...
                                            NaN, sum(ALL.(cond_name).(cell_name).raw.d00.Sheath_Info.isbridgecnctd,'omitnan')];
        else
            day = fieldnames(ALL.(cond_name).(cell_name).raw);
            ALL.stats.bridges.(cond_name) = [ALL.stats.bridges.(cond_name);...
                                            sum(ALL.(cond_name).(cell_name).raw.d00.Sheath_Info.isbridgecnctd,'omitnan'), sum(ALL.(cond_name).(cell_name).raw.(day{end}).Sheath_Info.isbridgecnctd),'omitnan'];
        end
        % get net length changes 
        ALL.(cond_name).(cell_name).net_sheath_LengthDiffs = sum(ALL.(cond_name).(cell_name).sheath_LengthDiffs,2,'omitnan');
    %% get length changes per timepoint 
        tps = size(ALL.(cond_name).(cell_name).sheath_LengthDiffs,2);
        ALL.stats.avgLnthChngPerTP.(cond_name) = [ALL.stats.avgLnthChngPerTP.(cond_name); mean(ALL.(cond_name).(cell_name).sheath_LengthDiffs,'omitnan'),NaN(1,7-tps)];
        sz = size(ALL.(cond_name).(cell_name).sheath_LengthDiffs);
        % extensions only
        logi = ALL.(cond_name).(cell_name).sheath_LengthDiffs > 0;
        temp = ALL.(cond_name).(cell_name).sheath_LengthDiffs .* logi;
        temp(temp==0) = NaN;
        ALL.(cond_name).(cell_name).extPerTP = [temp, NaN(sz(1),7-tps)];
        avgExts = mean(ALL.(cond_name).(cell_name).extPerTP,'omitnan');
        numExts = sum(~isnan(ALL.(cond_name).(cell_name).extPerTP));
        ALL.stats.extLnthPerTP.(cond_name) = [ALL.stats.extLnthPerTP.(cond_name); avgExts];
        ALL.stats.extNumPerTP.(cond_name) = [ALL.stats.extNumPerTP.(cond_name); numExts];
        % retractions only
        logi = ALL.(cond_name).(cell_name).sheath_LengthDiffs < 0;
        temp = ALL.(cond_name).(cell_name).sheath_LengthDiffs .* logi;
        temp(temp==0) = NaN;
        ALL.(cond_name).(cell_name).retPerTP = [temp, NaN(sz(1),7-tps)];
        avgRets = mean(ALL.(cond_name).(cell_name).retPerTP,'omitnan');
        numRets = sum(~isnan(ALL.(cond_name).(cell_name).retPerTP));
        ALL.stats.retLnthPerTP.(cond_name) = [ALL.stats.retLnthPerTP.(cond_name); avgRets];
        ALL.stats.retNumPerTP.(cond_name) = [ALL.stats.retNumPerTP.(cond_name); numRets];
    end
end