macro "Measure distance to line [q]" {
  requires("1.49k");
  path = getDirectory("image");  // used for outputting the results
  title = getTitle();                          // used for outputting the results
  title = replace(title, ".tif", "");  title = replace(title, ".lsm", "");  title = replace(title, ".czi", "");   title = replace(title, ".jpg", "");  
  title = replace(title, ".TIF", "");  title = replace(title, ".LSM", "");  title = replace(title, ".CZI", ""); title = replace(title, ".JPG", "");  
  run("Remove Overlay");  // cleans up any previous overlays -- may be deleted or commented out
  getPixelSize(unit, pixelWidth, pixelHeight);
  run("Set Scale...", "distance=1 known=1 unit=pixel");
  roiManager("select", 0);
  if (selectionType() < 5) closedShape = true; else closedShape = false;
  run("Interpolate", "interval=2");  //  half as fast and arguably more precise, change to 1
  getSelectionCoordinates(x, y);
  Roi.setStrokeColor("#0080ff");
  run("Add Selection...");
  run("Set Measurements...", "  centroid redirect=None decimal=1");
  output = File.open(path + title + "_results.txt");
  print("selection#  \t  units \t distance \t units");
  print(output, "selection#  \t  units \t distance \t units");
  for (selection=1;  selection<roiManager("count");  selection++){
      roiManager("select", selection);
      run("Measure");
      xc = getResult("X", nResults()-1);          // get centroid
      yc = getResult("Y", nResults()-1);
      min = 9999999999999;
      for (line=0; line<x.length; line++) {
        dx = xc - x[line];
        dy = yc - y[line];
        d = sqrt( (dx*dx) + (dy*dy));
        if (d < min) { min = d;  lx = x[line];  ly = y[line]; }
      }  // for each point on the line
      // make negative if inside the big shape.
     // select the big shape

     color1 = "#00ff";   
     if (closedShape) {
       roiManager("select", 0);
       inside = Roi.contains(xc, yc);  //Returns "1" if the point x,y is inside the current selection or "0" if it is not. 
       if (inside > 0) {
         min = min * (-1);
         color1 = "#ff00";
       } // if > 0
     }  // if closedShape
     print((selection+1) + " \t " + min + " \t pixels\t "+ (min*pixelWidth) + " \t " + unit);
     print(output, (selection+1) + " \t " + min + " \t pixels\t "+ (min*pixelWidth) + " \t " + unit);
     makeLine(xc, yc, lx, ly);
     Roi.setStrokeColor(color1 + "ff");
     run("Add Selection...");

  }  // for each selection

  run("Select None");
  run("Set Scale...", "distance=1 known="+pixelWidth+"  unit="+unit);  // restore scale to image; assumes square pixels
  run("Clear Results");
  File.close(output);
  selectWindow("Log");
}