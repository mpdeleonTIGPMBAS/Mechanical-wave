name=getTitle();
run("Brightness/Contrast...");
waitForUser("Please select center points for all areas of interest. Click OK when done");
run("Apply LUT", "stack");
saveAs("Tiff","V:/Marco/Wound Healing/Recoil speed/Reanalysis/blue channel/Cut/" + name + "intenadj");
makeRectangle(453, 297, 150, 250);
waitForUser("Please select center points for all areas of interest. Click OK when done");
roiManager("Add");
run("Duplicate...", "duplicate");
run("Reslice [/]...", "output=1.0 start=Left avoid");
run("Z Project...", "projection=[Min Intensity]");
run("Scale...", "x=10 y=10 width=1500 height=110 interpolation=Bilinear average create");
name=getTitle();
saveAs("Tiff","V:/Marco/Wound Healing/Recoil speed/Reanalysis/blue channel/Cut/" + name);
getStatistics(area, median);
factor = 1.5
thres = median * factor
setThreshold(thres, 255);
run("Convert to Mask");
saveAs("Tiff","V:/Marco/Wound Healing/Recoil speed/Reanalysis/blue channel/Cut/" + name + "binary"+"-"+factor);
run("Close All");

