directories = newArray("D:\\ParanodalBridges\\OldMice\\MOBPM_190417_old3_Reg2\\");
for ( k=0; k<directories.length; k++ ) { 
	dir = directories[k];
	savename = File.getName(dir);
	print(savename);
	Quad1dir = dir + "\Quad1\\";
	File.makeDirectory(Quad1dir); 
	Quad2dir = dir + "\Quad2\\";
	File.makeDirectory(Quad2dir);
	Quad3dir = dir + "\Quad3\\";
	File.makeDirectory(Quad3dir);
	Quad4dir = dir + "\Quad4\\";
	File.makeDirectory(Quad4dir);
	list = getFileList(dir);
	for ( i=0; i<list.length; i++ ) { 
		list = getFileList(dir);
		for ( i=0; i<list.length; i++ ) { 
			print(dir + list[i]);
			open( dir + list[i] ); 
		} 
		getDimensions(w,h,channels,slices,frames);
		slices = toString(slices);
		run("Concatenate...", "all_open title=[Concatenated Stacks] open");
		tp = toString(list.length - 4);
		//calculate stack size (number of frames)
		run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices=" + slices + " frames=" + tp + " display=Color");
		title = getTitle;
		b = w/2;
		makeRectangle(b,0,b,b);
		run("Duplicate...", "title=Quad1 duplicate");
		saveAs("Tiff", Quad1dir + savename + "_Quad1.tif");
		run("Close");
		
		selectWindow(title);
		makeRectangle(0,0,b,b);
		run("Duplicate...", "title=Quad2 duplicate");
		saveAs("Tiff", Quad2dir + savename + "_Quad2.tif");
		run("Close");
		
		selectWindow(title);
		makeRectangle(0,b,b,b);
		run("Duplicate...", "title=Quad3 duplicate");
		saveAs("Tiff", Quad3dir + savename + "_Quad3.tif");
		run("Close");
		
		selectWindow(title);
		makeRectangle(b,b,b,b);
		run("Duplicate...", "title=Quad4 duplicate");
		saveAs("Tiff", Quad4dir + savename + "_Quad4.tif");
		run("Close");
		//REGISTER QUAD 1
		run("Close All");
		open(Quad1dir + savename + "_Quad1.tif");
		title = getTitle();
		run("Split Channels");
		imageCalculator("Subtract create stack", "C1-"+title, "C2-"+title);
		run("Correct 3D drift", "channel=1 multi_time_scale only=0 lowest=10 highest=30");
		run("Median 3D...","X radius=1 Y radius=1 Z radius=1");		
		filepath = Quad1dir + savename + "_Quad1_3D1MedFilt" + ".tif";
		saveAs("Tiff", filepath);
		Split(Quad1dir,filepath);
		
		//REGISTER QUAD 2
		run("Close All");
		open(Quad2dir + savename + "_Quad2.tif");
		title = getTitle();
		run("Split Channels");
		imageCalculator("Subtract create stack", "C1-"+title, "C2-"+title);
		run("Correct 3D drift", "channel=1 multi_time_scale only=0 lowest=10 highest=30");
		run("Median 3D...","X radius=1 Y radius=1 Z radius=1");
		filepath = Quad2dir + savename + "_Quad2_3D1MedFilt" + ".tif";
		saveAs("Tiff", filepath);
		Split(Quad2dir,filepath);

		//REGISTER QUAD 3
		run("Close All");
		open(Quad3dir + savename + "_Quad3.tif");
		title = getTitle();
		run("Split Channels");
		imageCalculator("Subtract create stack", "C1-"+title, "C2-"+title);
		run("Correct 3D drift", "channel=1 multi_time_scale only=0 lowest=10 highest=30");
		run("Median 3D...","X radius=1 Y radius=1 Z radius=1");
		filepath = Quad3dir + savename + "_Quad3_3D1MedFilt" + ".tif";
		saveAs("Tiff", filepath);
		Split(Quad3dir,filepath);

		//REGISTER QUAD 4
		run("Close All");
		open(Quad4dir + savename + "_Quad4.tif");
		title = getTitle();
		run("Split Channels");
		imageCalculator("Subtract create stack", "C1-"+title, "C2-"+title);
		run("Correct 3D drift", "channel=1 multi_time_scale only=0 lowest=10 highest=30");
		run("Median 3D...","X radius=1 Y radius=1 Z radius=1");
		filepath = Quad4dir + savename + "_Quad4_3D1MedFilt" + ".tif";
		saveAs("Tiff", filepath);
		Split(Quad4dir,filepath);

		run("Close All");
	}
	
function Split(dir,filepath) {
			splitDir=dir + "\TimeSeriesSplit\\";
			File.makeDirectory(splitDir); 
			open(filepath);
			imgName=getTitle();
			run("Stack Slicer", "split_timepoints stack_order=XYCZT");
			close(imgName);
			ids=newArray(nImages); 
			for (i=0;i<nImages;i++) { 
				selectImage(i+1);
				title = getTitle;
				print(title);
			    ids[i]=getImageID; 
			    saveAs("Tiff", splitDir + title + ".tif"); 
			} 
			run("Close All");
}