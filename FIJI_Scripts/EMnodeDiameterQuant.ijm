run("Set Measurements...", "perimeter redirect=None decimal=3");
path = "D:\\OneDrive - Oregon Health & Science University\\Pictures\\Screenshots\\";
dirList = getFileList(path);
for (i = 0; i < dirList.length; i++) {
	open(dirList[i]);
	title = getTitle();
	print(title);
	run("Set Scale...", "distance=115 known=300 unit=nm global");
	setTool("freehand");
	waitForUser;
	run("Measure");
	close();
}