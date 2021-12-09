//for projecting the images to sum of all slices prior to analysis

dir=getDirectory("Choose a Directory")
dirA=dir+"Sum"+File.separator;


Images=getFileList(dir);
Array.sort(Images);
for(i=0; i<Images.length; i++) {
	if (endsWith(Images[i], "tif")) {
		open(Images[i]);
		name=getTitle();
		run("Z Project...", "projection=[Sum Slices] all");
		saveAs("Tiff", dirA + "C1-SUM_Result of C1" +name);
		run("Close All");
	}
}
