% edited 20180608 CLC - changed total sheath calculation to be specific for
% time point dates
% edited 20180612 CLC - added SheathLengthsPerTP, nodeClassPerTP
% edited 20180618 CLC - added c,f length diffs per TP and c/f ratio per TP

function [sheath_LengthDiffs, c_LengthDiffs, f_LengthDiffs, cfRatioPerTP,...
    Total_SheathLengths, proc_LengthDiffs,...
    sheathLengthPerTP, nodeClassPerTP]  = calcLengthChanges_old(TimeSeriesData)
%find trace names
days = fieldnames(TimeSeriesData)';
sheathNames = TimeSeriesData.(days{1}).Sheath_Info.intNames;
procNames = TimeSeriesData.(days{1}).Process_Info.procNames;

sheath_LengthDiffs = [];
c_LengthDiffs = [];
f_LengthDiffs = [];
cfRatioPerTP = [];
proc_LengthDiffs = [];
Total_SheathLengths = [];
sheathLengthPerTP = [];
nodeClassPerTP = [];

%compare timepoints 2-1, 3-2, 4-3, etc
i=0;
j=1;
skipflag = 0;
while j<8
    si1 = ['d' num2str(i*2,'%02d')];
    si2 = ['d' num2str(i*2+1,'%02d')];
    si3 = ['d' num2str(i*2+2,'%02d')];
    si4 = ['d' num2str(i*2+3,'%02d')];
    if j+1>length(days)
        sheath_LengthDiffs(1:length(sheathNames),j) = NaN;
        i=i+1;
        j=j+1;
    elseif length(days)==2
        nextTimepointNames = TimeSeriesData.(days{end}).Sheath_Info.intNames;
        nextTimepointLengths = TimeSeriesData.(days{end}).Sheath_Info.intLengths;
        thisTimepointNames = TimeSeriesData.d00.Sheath_Info.intNames;
        thisTimepointLengths = TimeSeriesData.d00.Sheath_Info.intLengths;
        for k = 1:length(sheathNames)
            if any(contains(nextTimepointNames,sheathNames{k}))
                idxN = ~cellfun(@isempty,strfind(nextTimepointNames,sheathNames{k}));
                idxT = ~cellfun(@isempty,strfind(thisTimepointNames,sheathNames{k}));
                sheath_LengthDiffs(k,j+skipflag) = nextTimepointLengths(idxN) - thisTimepointLengths(idxT);
            else
                idxT = ~cellfun(@isempty,strfind(thisTimepointNames,sheathNames{k}));
                sheath_LengthDiffs(k,j+skipflag) = - thisTimepointLengths(idxT);
            end
        end
        break
    elseif (contains(days{j},si1) | contains(days{j},si2)) & (contains(days{j+1},si3) | contains(days{j+1},si4))
        nextTimepointNames = TimeSeriesData.(days{j+1}).Sheath_Info.intNames;
        nextTimepointLengths = TimeSeriesData.(days{j+1}).Sheath_Info.intLengths;
        thisTimepointNames = TimeSeriesData.(days{j}).Sheath_Info.intNames;
        thisTimepointLengths = TimeSeriesData.(days{j}).Sheath_Info.intLengths;
        for k = 1:length(sheathNames)
            if any(contains(nextTimepointNames,sheathNames{k}))
                idxN = ~cellfun(@isempty,strfind(nextTimepointNames,sheathNames{k}));
                idxT = ~cellfun(@isempty,strfind(thisTimepointNames,sheathNames{k}));
                sheath_LengthDiffs(k,j+skipflag) = nextTimepointLengths(idxN) - thisTimepointLengths(idxT);
            else
                sheath_LengthDiffs(k,j+skipflag) = - thisTimepointLengths(idxT);
            end
        end
        i=i+1;
        j=j+1;
    elseif (contains(days{j},si1) | contains(days{j},si2)) & ~(contains(days{j+1},si3) | contains(days{j+1},si4))
        sheath_LengthDiffs = [sheath_LengthDiffs, NaN(length(sheathNames),2)];
        i=i+2;
        j=j+1;
        skipflag = skipflag + 1;
    end
end
% calculate total cell cohort length per time point, fixed for different timelines
i=0;
j=1;
while j<9
    si1 = ['d' num2str(i*2,'%02d')];
    si2 = ['d' num2str(i*2+1,'%02d')];
    si3 = ['d' num2str(i*2+2,'%02d')];
    si4 = ['d' num2str(i*2+3,'%02d')];
    if contains(si1,'d16')
        j=9;
        continue
    elseif length(days)==2
        Total_SheathLengths = [sum(TimeSeriesData.d00.Sheath_Info.intLengths,'omitnan'), sum(TimeSeriesData.(days{end}).Sheath_Info.intLengths,'omitnan')];
        sheathLengthPerTP = [TimeSeriesData.d00.Sheath_Info.intLengths, TimeSeriesData.(days{end}).Sheath_Info.intLengths];
        nodeClassPerTP = [TimeSeriesData.d00.Sheath_Info.node_class, TimeSeriesData.(days{end}).Sheath_Info.node_class];
        break
    elseif j>length(days)
        Total_SheathLengths = [Total_SheathLengths, NaN];
        sheathLengthPerTP = [sheathLengthPerTP, NaN(size(sheathLengthPerTP,1),1)];
        nodeClassPerTP = [nodeClassPerTP, NaN(size(nodeClassPerTP,1),1)];
        i=i+1;
        j=j+1;
    elseif contains(days{j},si1) | contains(days{j},si2)
        Total_SheathLengths = [Total_SheathLengths, sum(TimeSeriesData.(days{j}).Sheath_Info.intLengths,'omitnan')];
        sheathLengthPerTP = [sheathLengthPerTP, TimeSeriesData.(days{j}).Sheath_Info.intLengths];
        nodeClassPerTP = [nodeClassPerTP, TimeSeriesData.(days{j}).Sheath_Info.node_class];
        i=i+1;
        j=j+1;
    elseif contains(days{j},si3) | contains(days{j},si4)
        Total_SheathLengths = [Total_SheathLengths, NaN, sum(TimeSeriesData.(days{j}).Sheath_Info.intLengths,'omitnan')];
        sheathLengthPerTP = [sheathLengthPerTP, NaN(size(sheathLengthPerTP,1),1), TimeSeriesData.(days{j}).Sheath_Info.intLengths];
        nodeClassPerTP = [nodeClassPerTP, NaN(size(nodeClassPerTP,1),1), TimeSeriesData.(days{j}).Sheath_Info.node_class];
        i=i+2;
        j=j+1;
    end
end
i=0;
j=1;
skipflag = 0;
while j<7
    si1 = ['d' num2str(i*2,'%02d')];
    si2 = ['d' num2str(i*2+1,'%02d')];
    si3 = ['d' num2str(i*2+2,'%02d')];
    si4 = ['d' num2str(i*2+3,'%02d')];
    if length(days)==2 % for bsln to rec5wk comparison
            nextTimepointNames = TimeSeriesData.(days{end}).Process_Info.procNames;
            nextTimepointLengths = TimeSeriesData.(days{end}).Process_Info.procLengths;
            thisTimepointNames = TimeSeriesData.d00.Process_Info.procNames;
            thisTimepointLengths = TimeSeriesData.d00.Process_Info.procLengths;
            for k = 1:length(procNames)
                if any(contains(nextTimepointNames,procNames{k}))
                    idxN = ~cellfun(@isempty,strfind(nextTimepointNames,procNames{k}));
                    idxT = ~cellfun(@isempty,strfind(thisTimepointNames,procNames{k}));
                    proc_LengthDiffs(k) = nextTimepointLengths(idxN) - thisTimepointLengths(idxT);
                    if proc_LengthDiffs(k) > 0
                        proc_LengthDiffs(k) = 0;
                    end
                elseif any(contains(thisTimepointNames,procNames{k})) 
                    idxT = ~cellfun(@isempty,strfind(thisTimepointNames,procNames{k}));
                    proc_LengthDiffs(k) = - thisTimepointLengths(idxT);
                else
                    proc_LengthDiffs(k) = NaN;
                end
            end
            break
    elseif j+1>length(days)
        proc_LengthDiffs(1:length(procNames),j) = NaN;
        i=i+1;
        j=j+1;
    elseif (contains(days{j},si1) | contains(days{j},si2)) & (contains(days{j+1},si3) | contains(days{j+1},si4))
            nextTimepointNames = TimeSeriesData.(days{j+1}).Process_Info.procNames;
            nextTimepointLengths = TimeSeriesData.(days{j+1}).Process_Info.procLengths;
            thisTimepointNames = TimeSeriesData.(days{j}).Process_Info.procNames;
            thisTimepointLengths = TimeSeriesData.(days{j}).Process_Info.procLengths;
            for k = 1:length(procNames)
                if j>=2
                    lastTimepointNames = TimeSeriesData.(days{j-1}).Process_Info.procNames;
                    if contains(lastTimepointNames,procNames{k}) & ~contains(thisTimepointNames,procNames{k}) & contains(nextTimepointNames,procNames{k})
                        sprintf('In %s process %s is missing\n',days{j},procNames{k})
                        break
                    end
                end
                if any(contains(nextTimepointNames,procNames{k}))
                    idxN = ~cellfun(@isempty,strfind(nextTimepointNames,procNames{k}));
                    idxT = ~cellfun(@isempty,strfind(thisTimepointNames,procNames{k}));
                    proc_LengthDiffs(k,j+skipflag) = nextTimepointLengths(idxN) - thisTimepointLengths(idxT);
                    if proc_LengthDiffs(k,j+skipflag) > 0
                        proc_LengthDiffs(k,j+skipflag) = 0;
                    end
                elseif any(contains(thisTimepointNames,procNames{k})) 
                    idxT = ~cellfun(@isempty,strfind(thisTimepointNames,procNames{k}));
                    proc_LengthDiffs(k,j+skipflag) = - thisTimepointLengths(idxT);
                else
                    proc_LengthDiffs(k,j+skipflag) = NaN;
                end
            end
            i=i+1;
            j=j+1;
    elseif (contains(days{j},si1) | contains(days{j},si2)) & ~(contains(days{j+1},si3) | contains(days{j+1},si4))
        proc_LengthDiffs = [proc_LengthDiffs, NaN(size(proc_LengthDiffs,1),2)];
        i=i+2;
        j=j+1;
        skipflag = skipflag + 1;
    end
end
% calculate length changes between time points for f and c portions of sheaths
i=0;
j=1;
skipflag = 0;
while j<8
    si1 = ['d' num2str(i*2,'%02d')];
    si2 = ['d' num2str(i*2+1,'%02d')];
    si3 = ['d' num2str(i*2+2,'%02d')];
    si4 = ['d' num2str(i*2+3,'%02d')];
    if j+1>length(days)
        c_LengthDiffs(1:length(sheathNames),j) = NaN;
        f_LengthDiffs(1:length(sheathNames),j) = NaN;
        i=i+1;
        j=j+1;
    elseif length(days)==2
        nextTimepointNames = TimeSeriesData.(days{end}).Sheath_Info.intNames;
        thisTimepointNames = TimeSeriesData.d00.Sheath_Info.intNames;
        CnextTimepointLengths = TimeSeriesData.(days{end}).Sheath_Info.cLnth;
        CthisTimepointLengths = TimeSeriesData.d00.Sheath_Info.cLnth;
        FnextTimepointLengths = TimeSeriesData.(days{end}).Sheath_Info.fLnth;
        FthisTimepointLengths = TimeSeriesData.d00.Sheath_Info.fLnth;
        for k = 1:length(sheathNames)
            if any(contains(nextTimepointNames,sheathNames{k}))
                idxN = ~cellfun(@isempty,strfind(nextTimepointNames,sheathNames{k}));
                idxT = ~cellfun(@isempty,strfind(thisTimepointNames,sheathNames{k}));
                c_LengthDiffs(k,j+skipflag) = CnextTimepointLengths(idxN) - CthisTimepointLengths(idxT);
                f_LengthDiffs(k,j+skipflag) = FnextTimepointLengths(idxN) - FthisTimepointLengths(idxT);
            else
                c_LengthDiffs(k,j+skipflag) = - CthisTimepointLengths(idxT);
                f_LengthDiffs(k,j+skipflag) = - FthisTimepointLengths(idxT);
            end
        end
        break
    elseif (contains(days{j},si1) | contains(days{j},si2)) & (contains(days{j+1},si3) | contains(days{j+1},si4))
        nextTimepointNames = TimeSeriesData.(days{j+1}).Sheath_Info.intNames;
        thisTimepointNames = TimeSeriesData.(days{j}).Sheath_Info.intNames;
        CnextTimepointLengths = TimeSeriesData.(days{j+1}).Sheath_Info.cLnth;
        CthisTimepointLengths = TimeSeriesData.(days{j}).Sheath_Info.cLnth;
        FnextTimepointLengths = TimeSeriesData.(days{j+1}).Sheath_Info.fLnth;
        FthisTimepointLengths = TimeSeriesData.(days{j}).Sheath_Info.fLnth;
        for k = 1:length(sheathNames)
            if any(contains(nextTimepointNames,sheathNames{k}))
                idxN = ~cellfun(@isempty,strfind(nextTimepointNames,sheathNames{k}));
                idxT = ~cellfun(@isempty,strfind(thisTimepointNames,sheathNames{k}));
                c_LengthDiffs(k,j+skipflag) = CnextTimepointLengths(idxN) - CthisTimepointLengths(idxT);
                f_LengthDiffs(k,j+skipflag) = FnextTimepointLengths(idxN) - FthisTimepointLengths(idxT);
            else
                c_LengthDiffs(k,j+skipflag) = - CthisTimepointLengths(idxT);
                f_LengthDiffs(k,j+skipflag) = - FthisTimepointLengths(idxT);
            end
        end
        i=i+1;
        j=j+1;
    elseif (contains(days{j},si1) | contains(days{j},si2)) & ~(contains(days{j+1},si3) | contains(days{j+1},si4))
        c_LengthDiffs = [c_LengthDiffs, NaN(length(sheathNames),2)];
        f_LengthDiffs = [f_LengthDiffs, NaN(length(sheathNames),2)];
        i=i+2;
        j=j+1;
        skipflag = skipflag + 1;
    end
end
cfRatioPerTP = c_LengthDiffs./f_LengthDiffs;
end