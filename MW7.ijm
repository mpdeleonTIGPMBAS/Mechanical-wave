//This macro will make a reslice of a whole tailfin image for the analysis of thickness
//Marco 20May06

dir=getDirectory("Choose a Directory");
name=File.getName(dir);
list = getFileList(dir);
Array.sort(list);
for (m = 0; m < list.length; m++) {
	if(endsWith(list[m], "tif") &&  startsWith(list[m], "F") ){ 
	print(list[m]);
	open(list[m]);
	selectWindow(list[m]);
	run("Show Info...");
	saveAs("Text", dir + name);
	run("Duplicate...", "title=duplicate duplicate");
	selectWindow("duplicate");
	run("Z Project...", "projection=[Max Intensity]");
	selectWindow("MAX_duplicate");
	setTool("polyline");
	waitForUser("Please select center points for all areas of interest. Click OK when done");
	saveAs("TIFF", dir + "MAX_Fused");
	run("Create Mask");
	selectWindow("Mask");
	saveAs("Tiff", dir + "Mask");
	run("Create Selection");
	selectWindow("duplicate");
	run("Restore Selection");
	run("Crop");
	selectWindow("duplicate");
	saveAs("tiff", dir + name);
	setTool("line");
	//run("XY Coordinates...");
	saveAs("XY Coordinates", dir + name + "_linecoordinates");
	waitForUser("Please select center points for all areas of interest. Click OK when done");
	run("Reslice [/]...", "output=20.000 slice_count=1 avoid");
	run("Scale...", "x=1 y=20 width=1024 height=1000 interpolation=None average create");
	run("Rotate 90 Degrees Right");
	saveAs("tiff", dir + name + "_Reslice");
	run("Close All");
	}}
