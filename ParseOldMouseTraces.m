function [out2, brgRatio] = ParseOldMouseTraces
maindir = 'D:\GitHubRepos\Call_ParanodalBridge_2022\OldMOBPTraces';
addpath(genpath(maindir))
main = dir(maindir);
files = {main.name};
files = files(3:end);
pnbdata = cell(17,1);
pnbdataS = struct;
names = [];
for f = 1:length(files)
    idx = strfind(files{f},'old');
    an = ['old' files{f}(idx+3)];
    idx = strfind(files{f},'Reg');
    reg = ['reg' files{f}(idx+3)];
    idx = strfind(files{f},'Quad');
    quad = ['quad' str2double(files{f}(idx+3))];
    names = [names; an reg quad];
    xmlstruct = parseXML_SingleCell(files{f});
    pnbdata{f} = arrayfun(@(x) convertCharsToStrings(x.attribs.name),xmlstruct.paths);
    if isfield(pnbdataS,an)
        if isfield(pnbdataS.(an),reg)
            pnbdataS.(an).(reg) = [pnbdataS.(an).(reg) arrayfun(@(x) convertCharsToStrings(x.attribs.name),xmlstruct.paths)];
        else
            pnbdataS.(an).(reg) = arrayfun(@(x) convertCharsToStrings(x.attribs.name),xmlstruct.paths);
        end    
    else
        pnbdataS.(an).(reg) = arrayfun(@(x) convertCharsToStrings(x.attribs.name),xmlstruct.paths);
    end
end
out = [];
for i = 1:14
    out = [out [sum(contains(pnbdata{i},'both')); sum(contains(pnbdata{i},'pnb')); ...
           sum(contains(pnbdata{i},'non')) - sum(contains(pnbdata{i},'short'))]];
end

out(out==-1) = 0;
out2 = out';
out2 = [out2 sum(out2,2)];
out2 = sortrows(out2,4,'descend');
figure; 
b = bar(out2(:,1:3),'stacked');
colors = getPNBFigColors;
b(1).FaceColor = colors{3};
b(2).FaceColor = colors{2};
b(3).FaceColor = colors{1};
figQuality(gcf,gca,[4 3])

figure
plotSpread(out2(:,1:2),'distributionMarker',{'o'}, 'distributionColors',{colors{2} colors{3}})
ylim([0 7])
xlim([.5 2.5])
figQuality(gcf,gca,[2 3])
xticklabels([])

WSRp = signrank(out2(:,1),out2(:,2))

brgRatio = sum(out2(:,1:2),2) ./ out2(:,4);
%%
% an = fieldnames(pnbdata);
% out = cell(1, length(an));
% for a = 1:length(an)
%     reg = fieldnames(pnbdata.(an{a}));
%     for r = 1:length(reg)
%         out{a} = [out{a} [sum(contains(pnbdata.(an{a}).(reg{r}),'pnb')); sum(contains(pnbdata.(an{a}).(reg{r}),'both'));...
%                   sum(contains(pnbdata.(an{a}).(reg{r}),'non')) - sum(contains(pnbdata.(an{a}).(reg{r}),'short'))]];
%     end
% end
%%
% out2 = arrayfun(@(x) riddance(x{:}), out, 'UniformOutput',false);
% 
% function y = riddance(x)
%     y = [x(1:2,:);
%         x(3,:) - x(4,:)];
% end
