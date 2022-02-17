//directories = newArray(getDirectory("Choose a Directory"));
basedir = "\\\\MPICoreData.win.ad.jhu.edu\\Data\\Bergles lab\\Cody\\OldMOBP\\MOBPM_190417_old4_Reg1\\";
directories = newArray(basedir + "Quad2\\",basedir + "Quad3\\",basedir + "Quad4\\");
for ( k=0; k<directories.length; k++ ) { 
	dir = directories[k];
	print(dir);
	savename = File.getName(dir);
	print(savename);
	list = getFileList(dir);
	for ( i=0; i<list.length; i++ ) { 
		open( dir + list[i] );
	} 
	getDimensions(w,h,channels,slices,frames);
	title = getTitle();
	a = (w/2)-(w/4);
	b = w/2;
	makeRectangle(a, a, b, b);
	run("Correct 3D drift", "channel=1 only=50 lowest=20 highest=100");
	selectWindow(title);
	run("Close");
	
	//save the registered stack and get its name
	saveAs("Tiff", dir + savename + ".tif");
	
	//dir = getDirectory("Choose a Directory");
	//savename = File.nameWithoutExtension;
	
	filepath = dir + savename + ".tif";
	Split(dir,filepath);
	function Split(dir,filepath) {
		splitDir=dir + "\TimeSeriesSplit\\";
		File.makeDirectory(splitDir); 
		//open(filepath);
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
	print("DONE");
}