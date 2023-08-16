data = readcell('D:\ParanodalBridges\loopReconForFigs\skels\PNB_16_Loops_uncollected.csv');
fpath = 'D:\ParanodalBridges\loopReconForFigs\PNB_16\';
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

A = alphaShape(i1, i2, i3.*(40/3.58), 20, 'HoleThreshold',10000);
% [f,v] = isosurface(axonmask,3);
% p = patch('Vertices',v,'Faces',f);
% isonormals(axonmask,p);
% p.FaceColor = 'red';
% p.EdgeColor = 'blue';

figure
hold on
lines = struct;
for i = 1:length(sep_ori)
    loops = sep_ori{i};
    loopsSm = smoothdata(loops,"gaussian",3);
    loopsSm = [loops(1,:); loopsSm; loops(end,:)];
    lines(i).index = i;
    lines(i).coords = [loopsSm(:,2),loopsSm(:,1),loopsSm(:,3).*(40/3.58)];
    switch i
        case 1
            lines(i).color = ones(length(loopsSm),1).*0.4;
        case 2
            lines(i).color = ones(length(loopsSm),1).*0.8;
        case 3
            lines(i).color = ones(length(loopsSm),1).*0.1;
    end
%     plot3(loopsSm(:,2),loopsSm(:,1),loopsSm(:,3).*(40/3.58),'-','LineWidth',8,'LineJoin','round');
end
h = RenderLines2Tubes(lines,10,10);
plot(A,'EdgeColor','none','FaceAlpha',0.4,'FaceColor',[0.2,0.2,0.2])

axis tight
grid on
%%
xlim([-2000 -400])
ylim([500 1700])
zlim([50 650])

xticks('auto')
yticks('auto')
zticks('auto')
xticklabels('auto')
yticklabels('auto')
zticklabels('auto')
%%
camlight
lighting gouraud
xticklabels([])
yticklabels([])
zticklabels([])
