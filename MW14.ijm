setTool("polygon");

waitForUser("crop the interray");
setBackgroundColor(0, 0, 0);
run("Clear Outside");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Close-");
saveAs("Tiff", "V:/Marco/Wound Healing/Whole tailfin for width measurement/From tailfin thickness/M382A-mask.tif"); //
selectImage("M382A-mask.tif"); //
run("Reslice [/]...", "output=1.000 start=Top avoid");
run("Analyze Particles...", "display add stack");
saveAs("Results", "V:/Marco/Wound Healing/Whole tailfin for width measurement/From tailfin thickness/M382A-widthprofile.csv"); //
run("Close All");

