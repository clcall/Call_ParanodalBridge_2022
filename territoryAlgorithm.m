% get soma center, subtract center from all proc/sheath coords
% calc euclidian distance from origin (all points are now vectors), find furthest point from soma center
% 
% (x/a)2 + (y/b)2 + (z/c)2 = 1
% if < 1, the point is inside the ellipsoid
%
function [optRadxy, optRadz, coords] = territoryAlgorithm(index,xmlstruct)
coords=[]; 
for i = 1:length(index)
    coords = [coords; xmlstruct.paths(index(i)).points.smoothed];
end
c = mean(coords);
coords = coords - c;
ubx = max(coords(:,1));
uby = max(coords(:,2)); 
ubxy = max([ubx uby]); %since there is no polarity, x and y are equivalent
ubz = max(abs(coords(:,3)));
lb = 2;

xytest = lb:1:ubxy;
ztest = lb:1:ubz;
x = coords(:,1);
y = coords(:,2);
z = coords(:,3);
R = NaN(length(xytest), length(ztest));
for i = 1:length(xytest)
    a = xytest(i);
    b = xytest(i);
    for j = 1:length(ztest)
        c = ztest(j);
        in = (x./a).^2 + (y./b).^2 + (z./c).^2 <= 1;
        propin = sum(in) ./ size(coords,1);
        v = (4/3)*pi*a*b*c;
        if propin >= .8
            R(i,j) = propin./v;
        else
            R(i,j) = NaN;
        end
    end
end
[~,I] = max(R(:));
[optxyIdx, optzIdx] = ind2sub(size(R),I);
optRadxy = xytest(optxyIdx);
optRadz = ztest(optzIdx);
end
