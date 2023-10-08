addpath(genpath('D:\GitHubRepos\Call_ParanodalBridge_2022'))
data = readcell('D:\ParanodalBridges\loopReconForFigs\skels\Node_03_Loops_uncollected.csv');
fpath = 'D:\ParanodalBridges\loopReconForFigs\Node_03\';
dataorig = data;
headers = cellfun(@ismissing,data(:,7),'UniformOutput',false);
headers = cellfun(@sum,headers);
data(logical(headers),:) = [];
IDs = unique(data(:,7));
separated = cell(1,length(IDs));
for i = 1:length(IDs)
    idx = contains(data(:,7),IDs(i));
    seg1 = data(idx,1);
    seg2 = data(idx,2);
    segcoords = cellfun(@(x) strsplit(x,["(",", ",")"]),seg1,'UniformOutput',false);
    lastcoord = cellfun(@(x) strsplit(x,["(",", ",")"]),seg2(end),'UniformOutput',false);
    segcoords2 = NaN(length(segcoords),3);
    for j = 1:length(segcoords)
        segcoords2(j,:) = [str2double(segcoords{j}(2)) str2double(segcoords{j}(3)) str2double(segcoords{j}(4))];
    end
    segcoords2 = [segcoords2; [str2double(lastcoord{1}(2)) str2double(lastcoord{1}(3)) str2double(lastcoord{1}(4))]];
    separated{i} = segcoords2;
end 

% get maxes mins for xyz, reset smallest coord to origin
max_x = max(cell2mat(cellfun(@(x) max(x(:,1)), separated,'UniformOutput',false))) + 100;
max_y = max(cell2mat(cellfun(@(x) max(x(:,2)), separated,'UniformOutput',false))) + 100;
max_z = max(cell2mat(cellfun(@(x) max(x(:,3)), separated,'UniformOutput',false))) + 10;

min_x = min(cell2mat(cellfun(@(x) min(x(:,1)), separated,'UniformOutput',false))) - 100;
min_y = min(cell2mat(cellfun(@(x) min(x(:,2)), separated,'UniformOutput',false))) - 100;
min_z = min(cell2mat(cellfun(@(x) min(x(:,3)), separated,'UniformOutput',false))) - 10;

bbox = [max_x,max_y,max_z; min_x,min_y,min_z];

sep_ori = cell(size(separated));
for i = 1:length(separated)
    sep_ori{i} = separated{i} - [min_x,min_y,min_z];
%     
%     loopsSm = smoothdata(sep_ori{i},"gaussian",3);
%     
%     img = zeros([max_x-min_x, max_y-min_y, max_z-min_z]);
%     for j = 1:size(sep_ori{i},1)
%         img(sep_ori{i}(j,1),sep_ori{i}(j,2),sep_ori{i}(j,3)) = 1;
%     end
%     nhood = [0 1 0; 1 1 1; 0 1 0];
%     img = imdilate(img,nhood);
%     img = imdilate(img,nhood);
%     img = imdilate(img,nhood);
%     img = imdilate(img,nhood);
%     slices = size(img,3);
%     loopTraceNum = num2str(j);
%     for s = 1:slices
%         num = num2str(s,'%03.f');
%         imgname = fullfile(fpath,['trace',loopTraceNum,'_s',num,'.tif']);
%         imwrite(img(:,:,s),imgname);
%     end
end

axonmask = load_3D_gray(fullfile(fpath,['MASK_med30_erode7x','.tif']));
axonmask = double(axonmask);
ind = find(axonmask);
[i1, i2, i3] = ind2sub(size(axonmask), ind);


%%

A = alphaShape(i1.*3.58, i2.*3.58, i3.*40, 100, 'HoleThreshold',10000);
% [f,v] = isosurface(axonmask,3);
% p = patch('Vertices',v,'Faces',f);
% isonormals(axonmask,p);
% p.FaceColor = 'red';
% p.EdgeColor = 'blue';

%%
figure
hold on
lines = struct;
for i = 1:length(sep_ori)
    loops = sep_ori{i};
    loopsSm = smoothdata(loops,"gaussian",3);
    loopsSm = [loops(1,:); loopsSm; loops(end,:)];
    lines(i).index = i;
    lines(i).coords = [loopsSm(:,2).*3.58,loopsSm(:,1).*3.58,loopsSm(:,3).*(40)];
    switch i
        case 1
            lines(i).color = ones(length(loopsSm),1).*0.4;
        case 2
            lines(i).color = ones(length(loopsSm),1).*0.1;
        case 3
            lines(i).color = ones(length(loopsSm),1).*0.8;
    end
%     plot3(loopsSm(:,2),loopsSm(:,1),loopsSm(:,3).*(40/3.58),'-','LineWidth',8,'LineJoin','round');
end
if i<3
    lines(3).index = 3;
    lines(3).coords = repmat(zeros(length(loopsSm),1),1,3);
    lines(3).color = ones(length(loopsSm),1).*0.8;
end
h = RenderLines2Tubes(lines,40,40);
plot(A,'EdgeColor','none','FaceAlpha',0.4,'FaceColor',[0.2,0.2,0.2])

axis tight
grid on
hold off
%% PNB_04
xlim([0 2000])
ylim([1000 5000])
zlim([0 5000])

xticks('auto')
yticks('auto')
zticks('auto')
xticklabels('auto')
yticklabels('auto')
zticklabels('auto')

lgt = camlight('headlight');
lighting gouraud
xticklabels([])
yticklabels([])
zticklabels([])

lgt.Position = [lgt.Position(1) lgt.Position(2).*2.5 lgt.Position(3)];
%% PNB_15
set(gca, 'YDir','reverse') 
set(gca, 'ZDir','reverse')

xlim([1500 6500])
ylim([0 4000])
zlim([0 2000])

xticks(1500:1000:6500)
yticks(0:2000:4000)
zticks(0:1000:2000)

lgt = camlight('headlight');
lighting gouraud
xticklabels([])
yticklabels([])
zticklabels([])

%% Node_00
xlim([1000 5000])
ylim([1000 5000])
zlim([0 2500])

xticks(1000:1000:5000)
yticks(1000:1000:5000)

lgt = camlight('headlight');
lighting gouraud
xticklabels([])
yticklabels([])
zticklabels([])
%% Node_03
set(gca, 'ZDir','reverse')

xlim([500 4500])
ylim([500 4500])
zlim([-500 2500])

zticks(-500:1000:2500)
xticks(500:1000:4500)
yticks(500:1000:4500)

lgt = camlight('headlight');
lighting gouraud
xticklabels([])
yticklabels([])
zticklabels([])

%% theoretical schematic
figure
hold on


t = linspace(1.59,10*pi,1000);
r1 = 1 + (t./20);
r2 = ones(1,1000);
r = [r2(51:end) r1(1:50)];
rs = smoothdata(r,"gaussian",10);
x = rs.*cos(t);
y = rs.*sin(t);
z = t/10;
% plot3(x,z-0.2,y);
% plot3(x,flip(z+3.2),y);
% plot3([1.05 1.05],[-1 7],[0 0])


linesT = struct;
linesT(1).index = 1;
linesT(1).color = ones(1000,1).*0;
linesT(1).coords = [x',(z-0.2)',y'];

linesT(2).index = 2;
linesT(2).color = ones(1000,1).*-0.8;
linesT(2).coords = [x',flip(z+3.2)',y'];

linesT(3).index = 3;
linesT(3).color = ones(1000,1).*0.8;
linesT(3).coords = [ones(1000,1).*1.2, linspace(-1,7,1000)', zeros(1000,1)];

hT = RenderLines2Tubes(linesT,0.07,12);


[x,y,z] = cylinder(0.85,1000);
h = surf(x,(z.*8),y);
h.EdgeColor = 'none';
h.FaceColor = [0.1 0.1 0.1];
h.FaceAlpha = 0.4;
h2 = surf(x,(z.*-2),y);
h2.EdgeColor = 'none';
h2.FaceColor = [0.1 0.1 0.1];
h2.FaceAlpha = 0.4;
% h3 = surf((z.*-2),x+3,y);
% h3.EdgeColor = 'none';
% h3.FaceColor = [0.1 0.1 0.1];
% h3.FaceAlpha = 0.4;
hold off

axis tight

lgt = camlight('headlight');
lgt.Position = [4 6 -6];
lighting gouraud
xticklabels([])
yticklabels([])
zticklabels([])