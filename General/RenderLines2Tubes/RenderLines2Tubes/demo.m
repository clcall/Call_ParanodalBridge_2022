% RenderLines2Tubes renders the input curve line cluster into a set of tubes.
% This program is built upon 'tubeplot' developed by Janus H. Wesenberg in 2016, 
% but with the newly organized data structure, it allows to render multiple curves at a single run, 
% which improves the effiency of the original project (tubeplot) significantly.
% Copyright (c) 2020, Junpeng Wang
% All rights reserved. Please read the "license.txt" for license terms.	

% clear all
clc

Lines = ReadData('D:\GitHubRepos\OligodendrocyteAnalysisCode\General\RenderLines2Tubes\RenderLines2Tubes\data\dataset_2.dat');
hd = RenderLines2Tubes(Lines, 0.005, 16);
%set(hd, 'FaceColor', [0 1 0]); 
axis equal; axis off;
lighting gouraud;
[az el] = view;
camlight(az, el);

	
function val = ReadData(fileName)
	fid = fopen(fileName, 'r');
	numLines = fscanf(fid, '%d', 1);
	val = LinesStruct();
	val = repmat(val, numLines, 1);
	for ii=1:1:numLines
		ilth = fscanf(fid, '%d', 1);
		tmp = fscanf(fid, '%f %f %f %f', [4 ilth])';
		val(ii).coords = tmp(:,1:3);
		val(ii).color = tmp(:,4);
	end
	fclose(fid);
end

function val = LinesStruct()
	val = struct(...
		'index',	0,	...
		'coords',	[],	...
		'color',	[]	...
	);
end