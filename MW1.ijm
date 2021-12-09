// clear log
print("\\Clear")
// close all images first
run("Close All");


//select main folder
dir=getDirectory("Choose a Directory"); 
print("Dir:"+dir); 

list = getFileList(dir);
for (i = 0; i < list.length ; i++) {
	if (endsWith(list[i], "tif")) {
		print(list[i]);	

		//open raw image
		open(dir+list[i]);
		
		name=getTitle();
		print(name);
		
		selectWindow(name);

		//create output folder for the current raw img
		dirOut= dir + "" + name;		
		print("before trimming: " + dirOut);  
		dirOut=substring(dirOut, 0, lengthOf(dirOut)-4); // removes the 4 characters ".tif" 
		print("after trimming: " + dirOut);  

		run("Show Info...");
		saveAs("Text",dirOut+"Info.txt");
		run("Close");
		
		selectWindow(name);
 		run("Z Project...", "projection=[Max Intensity]");
		//run("Properties...", "channels=1 slices=25 frames=1 unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1 frame=[120.00 sec]");
		run("Set Scale...", "distance=1 known=1 pixel=1 unit=pixel");
		run("Select None");
		selectWindow("MAX_Raw.tif");
		setTool("multipoint");
		waitForUser("Please select center points for all areas of interest. Click OK when done");
		saveAs("XY Coordinates",dirOut+"points.txt");
		getSelectionCoordinates(xCoordinates,yCoordinates);
		saveAs("Tiff",dirOut+"_Max");
		close();
		
		// in Aligned.tif image, create rect and save as ray#
		for(w=0;w<lengthOf(xCoordinates);w++) {
			selectWindow("Raw.tif");
			//x=xCoordinates[i]/pixelWidth;
			makeRectangle(xCoordinates[w]-60,0,120,1024);
			run("Duplicate...", "duplicate");
			//run("Tiff...");
			makeRectangle(23, 924, 77, 91);
			saveAs("Tiff",dirOut+"Raw_"+w+"");
			//run("Close")
			run("Estimate Drift", "time=1 max=0 reference=[first frame (default, better for fixed)] apply choose=["+dirOut+".njt]");
			saveAs("Tiff",dir+"Img_nanocor_"+w+"");
			

			name=getTitle();
			print(name);
			dirA= dir +"Img_nanocor_" +w; 
			File.makeDirectory(dirA);  
			print(dirA); 
			//save as Image Sequence
			name1= substring(name, 0, lengthOf(name)-4);
			file = dirA + File.separator +name1+ "_";

			
			run("Image Sequence... ", "format=TIFF name=" +name1+ "_0 digits=2 save=["+ file +".tif]");
			print(file);
			selectWindow(name);
			close;
		}
		Tscale=20;
			File.append("TimeScale="+Tscale,dirA+"calibration.txt");
			run("Close All");
	}
}



	