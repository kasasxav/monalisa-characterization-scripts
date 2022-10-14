 //get the profile, and do a fit, try to save the data in RESULT and in a folder
 //do the same for all 8 slices in an ABCD series

 //here write the save name
 svnm="line-profiles"
 savename="D:\\OneDrive - KTH\\2021_full\\2022-01-07\\fullchip_analysis\\"+svnm;
  //name of the orignal image
   parent=getImageID;
  run("Clear Results");
//use the ROI manager to select the lines
rois=roiManager("count")

  
for (k=0; k<rois; k++){
selectImage(parent);
roiManager("select",k);

run("Plot Profile");
Plot.getValues(x, y);
close();
 for (i=0; i<x.length; i++){
//      setResult("x"+k, i, x[i]);
      setResult("y"+k, i, y[i]);
  updateResults;
 }
}
   saveAs("Results", savename+".csv");