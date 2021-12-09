//Simple macro for separating 3 channels with maximum projection
//requirements: 3 folders C1 C2 and Merged must be made inside the main file folder
//Files to be processed must be inside the main file folder
//First prompt will ask for the main file folder
//Second prompt will ask for the file to be analyzed

dir=getDirectory("Choose a Directory")
dirA=dir+"Cropped"+File.separator;


Images=getFileList(dir);
//Array.sort(Images);
for(i=0; i<Images.length; i++) {
	if (endsWith(Images[i], "tif")) {
		open(Images[i]);
		name=getTitle();
		setTool("polygon");
		waitForUser("center the rectangle then press ok");
	    roiManager("Add");
		run("Measure");
		run("Close All");
	}
}

