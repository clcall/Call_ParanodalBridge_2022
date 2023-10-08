function usrHandle = RenderLines2Tubes(lines, r, n, ct)	
	% RenderLines2Tubes renders the input curve line cluster into a set of tubes.
	% This program is built upon 'tubeplot' developed by Janus H. Wesenberg in 2016, 
	% but with the newly organized data structure, it allows to render multiple curves at a single run, 
	% which improves the effiency of the original project (tubeplot) significantly.
	% Copyright (c) 2020, Junpeng Wang
	% All rights reserved. Please read the "license.txt" for license terms.	
	%
	% Arguments:
	% lines: lines data including the Cartesian coordinates and their attributes, the latter is used 
	%		 to color code the lines
	% r      the radius of the tube
	% n      number of points to use on circumference. Defaults to 8
	% ct     threshold for collapsing points. Defaults to r/2	
	if nargin<3 || isempty(n), n=8;
		if nargin<2, error('Give at least lines and radius');
        end
    end
	if nargin<4 || isempty(ct)
		ct=0.5*r;
	end	
	numLines = length(lines);
	gridXYZ = zeros(3,n+1,1);
	gridC = zeros(n+1,1);
	for ii=1:1:numLines		
		curve = lines(ii).coords';
		npoints = size(curve,2);
		%deltavecs: average for internal points. first strecth for endpoitns.		
		dv = curve(:,[2:end,end])-curve(:,[1,1:end-1]);		
		%make nvec not parallel to dv(:,1)
		nvec=zeros(3,1); [buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;
		xyz=repmat([0],[3,n+1,npoints+2]);
		%precalculate cos and sing factors:
		cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
		sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
		%Main loop: propagate the normal (nvec) along the tube
		xyz = zeros(3,n+1,npoints+2);
		for k=1:npoints
			convec=cross(nvec,dv(:,k));
			convec=convec./norm(convec);
			nvec=cross(dv(:,k),convec);
			nvec=nvec./norm(nvec);
			%update xyz:
			xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1]) + cfact.*repmat(r*nvec,[1,n+1]) + sfact.*repmat(r*convec,[1,n+1]);
		end;
		%finally, cap the ends:
		xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
		xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
		gridXYZ(:,:,end+1:end+npoints+2) = xyz;	
		color = lines(ii).color';	
		c = [color(1) color color(end)];
		c = repmat(c, n+1, 1);
		gridC(:,end+1:end+npoints+2) = c;
	end
	gridX = squeeze(gridXYZ(1,:,:));
	gridY = squeeze(gridXYZ(2,:,:));
	gridZ = squeeze(gridXYZ(3,:,:));
	usrHandle = surf(gridX,gridY,gridZ,gridC); shading interp;
end