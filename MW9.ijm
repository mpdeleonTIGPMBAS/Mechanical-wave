Main=getDirectory("Folder to process"); // main file folder
dir=getDirectory("Choose a Directory"); // choose label(resize)
list = getFileList(dir); //unsorted!!!! REMEMBER!!!
Array.sort(list);
for (i = 0; i < list.length ; i++) {
	if (endsWith(list[i], "tif")) {
		print(list[i]);	
		open(dir+list[i]);
		name=getTitle();
		print(name);
		selectWindow(name);
		run("Bin...", "x=5 y=5 z=1 bin=Max");
		setThreshold(1, 99999);
		run("Convert to Mask", "method=Default background=Dark calculate black");	
		run("Create Selection");
		newImage("Untitled", "8-bit black", 120, 1024, 1);
		selectWindow("Untitled");
		run("Restore Selection");
		setForegroundColor(255, 255, 255);
		run("Draw", "slice");
		dir1= Main + "Seg-outline_origsize" + File.separator;
		File.makeDirectory(dir1);
		saveAs("tiff", dir1 + i+1);
 		run("Close All");
	}} // end for