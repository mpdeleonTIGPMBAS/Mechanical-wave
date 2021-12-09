print("\\Clear");
run("Close All");
setBatchMode(true);

mainDir = getDirectory("Choose a Directory")

files = getFileList(mainDir)
Array.sort(files)

for (m = 0; m < files.length; m++) {
	//if(endsWith(files[m], "analyzed/")) { // get all analyzed folder
	if(endsWith(files[m], "/") &&  startsWith(files[m], "M") ) { // get all analyzed folder		
		print(files[m]);		

		// access subfolders of analyzed subfolder
		dir1 = mainDir + files[m];
		print(dir1);
		file2 = getFileList(dir1);
		Array.sort(file2);
		
		for (p = 0; p < file2.length; p++) {
			if( endsWith( file2[p], "/")  &&  startsWith(file2[p], "Img_nano") ) { // get all folders
				print("nano "+file2[p]);

				// this is now the sub-subfolder Img_nano directory
				subsubDir = dir1 + file2[p];
				name=dir1 + substring(file2[p], 4, lengthOf(file2[p])-1);
				dir2=  name +"_whitewave"+File.separator;				
				File.makeDirectory(dir2);
				print(dir2);
				
				//The resultant images are saved in the "Difference_Images" folder (=dir3) that is created in dir2.
				
				dir3 = dir2 + "Difference_Images" + File.separator;
				File.makeDirectory(dir3);
				print(dir3);
				
				list = getFileList(subsubDir); //accessed all Image Sequence
				Array.sort(list);
				for(i=0; i<list.length-1; i++) {
					if (endsWith(list[i], "tif")) {
						//print(list[i]);
						path1 = subsubDir+list[i];
						open(path1);				
						rename("A");
						path2 = subsubDir+list[i+1];
						open(path2);
						rename("B");
						imageCalculator("Difference create", "A","B");
						//run("Image Calculator...", "image1=A operation=Difference image2=B create");
						v = "Diff_" +  substring(path1, lengthOf(path1)-7, lengthOf(path1)-4);
						saveAs("tiff", dir3 + v);
						run("Close All");
					}
					
				}
				run("Image Sequence...", "open=["+dir3+"Diff_000.tif] sort use");
				run("Reslice [/]...", "output=1 start=Left avoid");
				//run("Z Project...", "projection=[Max Intensity]");
				run("Z Project...", "projection=[Sum Slices]");
				run("8-bit");
				run("Scale...", "x=1 y=20 width=1024 height=1000 interpolation=None average create");
				run("Flip Horizontally");
				run("Rotate 90 Degrees Left");
				saveAs("tiff", name);
				run("Show Info...");
				saveAs("Text", name);
				run("Close All");


				
				
			} // end of if p
		}//end for p
				
	} // end of if m
	
}//end for m

print("\\Update:Done!")



