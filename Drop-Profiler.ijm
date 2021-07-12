////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//FUNCTIONS USED IN SCRIPT WHICH STARTS @ LINE ~275:

////////////////////////////////////////////////////////////////////////
function setCustomThreshold() { 
// Fitting Gaussian to background peak (left-most peak in LUT histogram),
// then sets lower threshold to a user-defined number of standard
// deviations to the right of the background peak mean.

	print(" \n## START setCustomThreshold() ## ");
	
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	threshold_upper = max;
	pixel_depth = max-min;
	print("pixel_depth = "+ pixel_depth);
	nBins = pixel_depth/pow(2, -floor(-pixel_depth*0.0004));
	print("nBins = "+ nBins);
	
	image_selected = getImageID();
	
	getHistogram(values, counts, nBins);
	plot_values = values;
	plot_counts = counts;
	
	Array.getStatistics(counts, min, max, mean, stdDev);
	print("mean_counts = "+mean);
	
	// Find maxima...
	maxima_list = Array.findMaxima(counts, mean*0.1);
	maxima_list_sorted = Array.sort(maxima_list);
	dummy_array = newArray(nBins-1, nBins-1, nBins-1, nBins-1, nBins-1, nBins-1, nBins-1, nBins-1, nBins-1, nBins-1, nBins-1);
	 maxLoc = Array.concat(maxima_list,dummy_array);
	
	print("maxima locations (bins): "); Array.print(maxLoc);
	print("maxLoc[0] = " + maxLoc[0]);
	print("values[maxLoc[0]] = " + values[maxLoc[0]]);
	print("counts[maxLoc[0]] = " + counts[maxLoc[0]]);
	
	//Exclude "premature" maxima (outliers) and define bkgrnd peak location (bin#)...
	if (maxLoc[0] == 0) {
		peak_x_location = maxLoc[1];
	} else {
		peak_x_location = maxLoc[0];
	}
	
	for (i = 0; i < 10; i++) {
		if ((peak_x_location == maxLoc[i]) == true
				&&	(maxLoc[i+1]-maxLoc[i] <= nBins*0.1) 
				&& (counts[maxLoc[i+1]] >= counts[maxLoc[i]]*0.9) 
				== true) {
			peak_x_location = maxLoc[i+1];
	};
	};
	
	 print("bkgnd maximum (peak_x_location): " + peak_x_location);
	 peak_x = values[peak_x_location]; 
	 print("bkgnd maximum, pixel value (peak_x) = " + peak_x);
	 maxCount = counts[peak_x_location];
	peak_x_bin = plot_values[peak_x_location];
	
	// The maximum cutoff for the fitting of the gaussian to the 
	// background peak in the data.
	
	//max_cutoff = peak_x_location * 1.2; //(includes left half of distribution only) 
	max_cutoff = peak_x_location * 3;
	
	// make trimmed dataset...
	new_x_values = Array.slice(values, nBins*0.00, max_cutoff);
	new_y_values = Array.slice(counts, nBins*0.00, max_cutoff);
	
	// Fit gaussian to the trimmed dataset...
	
	//using built-in functions:
	SD_initialguess = peak_x*0.3; print("SD_initialguess = "+ SD_initialguess);
	initialGuesses = newArray(maxCount, peak_x, SD_initialguess);
	Fit.doFit("Gaussian (no offset)", new_x_values, new_y_values, initialGuesses); //function is "y = a*exp(-(x-b)*(x-b)/(2*c*c)))"
	
	//DISPLAY PLOT (for trouble-shooting):
	//Fit.plot();
	
	// Extract fitted parameters from the gaussian fit...
	fit_max = Fit.p(0); //'a' in "Gaussian (no offset)"
	fit_mean = Fit.p(1);
	fit_sd = abs(Fit.p(2));
	
	print("SD_fit = " + fit_sd);
	
	// **Set the new threshold for the image ("user value" is arbitrary)...
	threshold_lower = fit_mean + fit_sd * user_value; 
	
	selectImage(image_selected);
	setThreshold(threshold_lower, threshold_upper);
	
	////DISPLAY PLOT (for trouble-shooting):
	//// below just creates data from the fitted gaussian for the creation of
	//// the plot.
	//
	//fitted_counts = Array.copy(plot_counts);;
	//
	//for (i = 0; i < plot_values.length; i++) {
	//	fitted_counts[i] = Fit.f(values[i]);
	//}
	//
	//// Creation of the plot for visual represenation of the data.
	//
	// Plot.create("Total Image Histogram", "Pixel Values", "Counts");
	// Plot.setColor("red");
	// Plot.setLineWidth(5);
	// Plot.add("dot", plot_values, plot_counts);
	// Plot.setLineWidth(2);
	// Plot.setColor("cyan");
	// Plot.add("line", plot_values, fitted_counts);
	// Plot.setColor("red");
	//
	// Plot.drawLine(threshold_lower, maxCount/2, threshold_lower, 0);
	// Plot.addText("Threshold Value", threshold_lower, maxCount/2);
	// 
	// Plot.setColor("black");
	// Plot.show;

	print(" *Lower Threshold = "+threshold_lower);
	print(" *Upper Threshold (max) = "+threshold_upper+" \n## END setCustomThreshold() ## \n ");
};
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
function findPlateSurface() { 
// Finding slice that best corresponding to plate surface (bottom of droplets):
// **Requires "setCustomThreshold()" function**
	
	print(" \n## START findPlateSurface() ## ");
	
	//First find slice with max average intensity (to estimate base of droplet):
	Zstack = getImageID();
	getVoxelSize(width, height, depth, unit);
	Voxel_depth = depth;
	Voxel_unit = unit;
	run("Plot Z-axis Profile");
	Plot.getValues(z_micron, Imean);
	run("Close");
	Array.getStatistics(Imean, min, max, mean, stdDev);
	maxLoc = Array.findMaxima(Imean, max/3);
	
	//find x (slice in micron) where y (Imean) = max
	z_base_intensity = Imean[maxLoc[0]];
	z_base_micron = z_micron[maxLoc[0]];
	z_base_slice = maxLoc[0]+1;
	
	print("highest intensity slice = "+z_base_slice);
	print("mean intensity = "+z_base_intensity);
	print("Z coord ("+Voxel_unit+") = "+z_base_micron);
	
	//check for slices immediately below this that have larger drop section area: 
	print("\nSearch for plate surface (largest total droplet area)...");
	selectImage(Zstack);
	setSlice(z_base_slice);
	setCustomThreshold();
 //*** custom function defined above^^
	getThreshold(lower, upper);
 //use these values in for loop below... 
	
	run("Create Selection");
 
	getStatistics(area, mean, min, max, std, histogram);
	drop_area = area;
	
	print(" Slice \\ Drop_Area");
	print(" "+z_base_slice+"    \\ "+drop_area);
	
	for (i = z_base_slice-1; i >= 1; i--) {
	    setSlice(i);
	    setThreshold(lower, upper);
	    run("Create Selection");
	    getStatistics(area, mean, min, max, std, histogram);
	    drop_area_temp = area;
	    	if ((drop_area_temp > drop_area) == true) {
	    		drop_area = drop_area_temp;
	    		z_base_slice = getSliceNumber();
   			    print(" "+z_base_slice+"    \\ "+drop_area);
	    	} else {
	    		print(" "+z_base_slice-1+"    \\ "+drop_area_temp);
	    		i = 0;
	    	};
	};

	setSlice(z_base_slice);
	run("Select None");
	print("plate surface (slice) = "+z_base_slice);
	z_base_micron = Voxel_depth*(z_base_slice-1);
	print("Z coord ("+Voxel_unit+") = "+z_base_micron);
	print("## END findPlateSurface() ##\n ");
};
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
function findContactAngle() { 
	// This will fit a circle to the surface of a  droplet and calculate the 
	// width of the droplet and droplet contact angle (i.e. the inside
	// angle at the base of the droplet between the tangent to the droplet surface
	// and the surface of the plate/slide).
	// Writes to Results Table.
	// *REQUIRES A PRE-THRESHOLDED IMAGE, WITH THE 
	// THRESHOLD DEFINING THE CURVED SURFACE.
	// *REQUIRES DROPLET TO BE DOWNWARD FACING,	WITH TOP OF IMAGE
	// DEFINING SURFACE OF THE SLIDE/PLATE.
	
	run("Create Selection");
	
	getSelectionCoordinates(xpoints, ypoints);
	
	makeSelection("point", xpoints, ypoints);
	
	run("Fit Circle");
	
	//find radius...
	getSelectionBounds(x, y, w, h);
	radius = w*0.5;
	print("Radius: "+radius);
	
	//find centre...
	x_centre = x+radius;
	y_centre = y+radius;
	print("centre: ("+x_centre+", "+y_centre+")");
	
	print("find droplet edges (find x where y=0)...");
	print("Circle equation is:");
	print("  (x - x_centre)^2 + (y - y_centre)^2 = radius^2");
	print("when y=0...");
	print("  (x - x_centre)^2 + y_centre^2 = radius^2");
	print("rearrange to standard form...");
	print("   ax^2 +       bx       +           c            = 0");
	print("  (x^2) + (-2*x_centre*x) + (x_centre^2 + y_centre^2 - radius^2) = 0");
	
	a=1; 
	print("a = "+a);
	b=(-2*x_centre); 
	print("b = "+b);
	c=(x_centre*x_centre + y_centre*y_centre - (radius*radius));
	print("c = "+c);
	
	print("solve for x using Quadratic formula...");
	print("  x = (-b +/- sqrt(b^2 - 4ac))/2a");
	edge_left = (-b - sqrt(b*b - 4*a*c))/2*a;
	edge_right = (-b + sqrt(b*b - 4*a*c))/2*a;
	
	print("edge_left: "+edge_left);
	print("edge_right: "+edge_right);
	
	droplet_width = edge_right-edge_left;
	print("droplet_width = "+droplet_width);
	
	
	// Use isosceles triangle formed by radii and chord (droplet_width) to calculate
	// angle (lambda) between chord and radius (trigonometry)...
	lambda = atan((2 * sqrt(radius*radius - ((droplet_width*droplet_width)/4)) ) / droplet_width);
	
	// ... multiply by 180/pi to get degrees...
	lambda = lambda * 180/PI;
	
	// ... the angle(s) between the chord and tangent(s) (i.e the Droplet Contact 
	// Angle) is equal to 90 minus lambda...
	ContactAngle = 90 - lambda;
	
	print("Contact Angle = "+ContactAngle);
	
		row = nResults;
		setResult("Width_px", row, droplet_width);
		setResult("Contact_Angle", row, ContactAngle);
		updateResults();
};
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

////////////////////////-- START --/////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//Z-STACK MUST BE BOTTOM->TOP (i.e. slice 1 is below plate surface)

//open Zstack and get details
if (nImages==0) {
	waitForUser("Please open Z-stack of droplet images, then hit OK.");
} else {
	waitForUser("Ensure that Z-stack of droplet images is the *only* image window open, then hit OK. \n"+
	"(If other images are open, close them now)");
};

//set "user value" for thresholding
Dialog.create("User input");
Dialog.addNumber("Droplet thresholding:", 15);
Dialog.addMessage("^ Number of positive standard deviations from mean of background peak of \n blank-subtracted image ('user value' - Wang et al, 'A Molecular Grammar...', Cell, 2018).\n Increasing this will increase the lower intensity threshold for defining pixels as condensed phase. \n (use ~10-15 for confocal images)");
Dialog.addCheckbox("Run median filter?", true);
Dialog.addNumber("Median filter radius:", 2, 1, 5, "pixels");
Dialog.addCheckbox("Select droplets manually?", false);
Dialog.addMessage("^ If unchecked, droplets will be selected for analysis automatically using the parameters below...");
Dialog.addMessage("--v--(auto droplet picking)--v--");
Dialog.addNumber("Droplet XY size, min:", 0, 1, 5, "microns^2");
Dialog.addNumber("Droplet XY size, max:", 1000, 1, 5, "microns^2");
Dialog.addSlider("Droplet *circularity*", 0, 1, 0.99);
Dialog.addMessage("(* Specifies how circular (in XY plane) a droplet must be to be included (1 = perfect circle).");
Dialog.addCheckbox("Interpolate data?", false);
Dialog.addMessage("--^--(auto droplet picking)--^--");

Dialog.show();

user_value = Dialog.getNumber();
filter_status = Dialog.getCheckbox();
filter_radius = Dialog.getNumber();
Manual_drop_select_status = Dialog.getCheckbox();
auto_min = Dialog.getNumber();
auto_max = Dialog.getNumber();
auto_circularity = Dialog.getNumber();
auto_interpolate_status = Dialog.getCheckbox();

//run median filter:
if (filter_status==true) {
	run("Median...", "radius="+filter_radius+" stack");
}

//add scalebar:
run("Scale Bar...", "width=10 height=4 font=14 color=White background=None location=[Lower Right] bold overlay");

//sample details...
selectImage(1);
run("In [+]");
run("Scale to Fit");
Zstackfilename = getTitle();
//XYZ = getImageID();
XYZ = getTitle();
dir = getDirectory("image");
print("Directory = " + dir);
print("Z-stack = " + Zstackfilename);

//create output directory:
output_dir = dir+Zstackfilename+"_analysis";
File.makeDirectory(output_dir);
print("Output Directory = " + output_dir);


	getVoxelSize(Vx_width, Vx_height, Vx_depth, Vx_unit);
		Vx_width=Vx_width;
		Vx_height=Vx_height;
		Vx_depth=Vx_depth;
		Vx_unit=Vx_unit;

	// write voxel_dimensions.csv to file:
	run("Clear Results");
	setResult("Vx_width", 0, Vx_width);
	setResult("Vx_height", 0, Vx_height);
	setResult("Vx_depth", 0, Vx_depth);
	setResult("Vx_unit", 0, Vx_unit);
	updateResults();
	selectWindow("Results");
	saveAs("results", output_dir+"\\voxel_dimensions.csv");
	run("Clear Results");
	run("Close");
	


print("\n##################################################################");
print("### Set master threshold using slice that corresponds \n### best to average plate surface... ");

selectImage(XYZ);
findPlateSurface();
plate_surface_slice_average = getSliceNumber();

setCustomThreshold();

//Save master threshold values for reusing when defining droplet surface 
//(after extracting side-on profiles)  
getThreshold(lower, upper);
master_threshold_lower = lower;
//master_threshold_upper = lower*1.8;
master_threshold_upper = lower*1.2;

print("Master threshold values (lower/upper):");
print(master_threshold_lower);
print(master_threshold_upper);

print("##################################################################\n");


Dialog.create("Select Orthogonal Views");
Dialog.addMessage("Which side-on views would you like to include in analysis?");
Dialog.addCheckbox("XZ", true);
Dialog.addCheckbox("YZ", true);
Dialog.show();

XZrunstatus = Dialog.getCheckbox();
YZrunstatus = Dialog.getCheckbox();


////////////////////////////////////////////////////////////////////////////////////////////
///// MANUAL DROPLET PICKING ///////
////////////////////////////////////////////////////////////////////////////////////////////
if (Manual_drop_select_status == true) {
		
	for (j = 1; j < 1000; j++) {
		
		if (j > 1 && nImages > 0) {
			selectImage(XYZ);
			close();
			print("good to 4");
			selectImage("XZ"+(j-1)+"_dropsurface.tif");
			print("good to 5");
			close();
			print("good to 6");
			selectImage("YZ"+(j-1)+"_dropsurface.tif");
			close();
			print("sdf");
		}
		
		print("\n##################################################################");
		print("### DROPLET #"+j+"... ");
		
		//re-open Zstack upon iteration
		if (nImages==0) {
			open(dir+Zstackfilename);
			run("In [+]");
			run("Scale to Fit");
			XYZ = getImageID();
			setSlice(plate_surface_slice_average);
			setThreshold(master_threshold_lower, master_threshold_upper);
		};
		
		//set slice to average plate surface:
		selectImage(XYZ);
		setSlice(plate_surface_slice_average);
		//run("Plot Z-axis Profile");
		//Plot.getValues(z_micron, Imean);
		//run("Close");
		//Array.getStatistics(Imean, min, max, mean, stdDev);
		//maxLoc = Array.findMaxima(Imean, max/2);
		//max_intensity_slice = maxLoc[0]+1;
		//setSlice(max_intensity_slice);
		
		////////////////
		// USER INPUT:
		waitForUser("draw selection box around a single droplet");
		
		getSelectionBounds(x, y, width, height);
		print("Droplet#"+j+" location and size (user-defined selection box):");
		print(" location (top left coords) = ("+x+", "+y+") pxls");
		print(" size (W x H) = "+width+" x "+height+" pxls"); 
		
		run("Crop");
		run("In [+]");
		run("In [+]");
		run("In [+]");
		run("In [+]");
		run("Scale to Fit");
		
		//re-find slice that best represents plate surface 
		//(in case plate is not perfectly flat):
		findPlateSurface();
		
		//remove all slices below this...
		//plate_surface_slice_drop = getSliceNumber() - 1; //(EXCLUSIVE)
		plate_surface_slice_drop = getSliceNumber();; //(INCLUSIVE)
		print("plate surface slice (for Droplet#"+j+") = "+plate_surface_slice_drop);
		
		selectImage(XYZ);
		run("Slice Remover", "first=1 last="+plate_surface_slice_drop+" increment=1");
		
		run("In [+]");
		run("In [+]");
		
		run("In [+]");
		
		run("Orthogonal Views");
		
		Stack.getOrthoViewsIDs(XY, YZ, XZ);
		selectImage(XYZ);
		Stack.setSlice(1);
		resetMinAndMax;
		
		// USER INPUT:
		waitForUser("click on centre of droplet");
		
		
	print("Threshold for defining droplet surface: \n "+
			"  Lower = "+master_threshold_lower+" \n "+
			"  Upper = "+master_threshold_upper);
			
		//////////////////////////////////////////////
		
		if (XZrunstatus == true) {
		
			selectImage(XZ);
			saveAs("Tiff", output_dir+"\\XZ"+j+".tif");
			open(output_dir+"\\XZ"+j+".tif");
			for (i = 0; i < 6; i++) {
				run("In [+]");
			};
			run("Scale to Fit");
	
			setThreshold(master_threshold_lower, master_threshold_upper);
		
		
			run("Create Selection");
			run("Save XY Coordinates...", "save=["+output_dir+"\\"+j+"_XZ_SurfacePx.csv]");
		
			saveAs("Tiff", output_dir+"\\"+"XZ"+j+"_dropsurface.tif");
		
			print(" \n**Output file containing XZ coordinates \n  of Droplet#"+j+" Surface can be found here: \n "+
			""+output_dir+"\\SurfacePx_XZ"+j+".csv");
	
			findContactAngle();
			
			setResult("slice_direction", nResults-1, "XZ");
			updateResults();
		};
		
		//////////////////////////////////////////////
		
		if (YZrunstatus == true) {
			
			selectImage(YZ);
			run("Rotate 90 Degrees Right");
			saveAs("Tiff", output_dir+"\\YZ"+j+".tif");
			open(output_dir+"\\YZ"+j+".tif");
			for (i = 0; i < 6; i++) {
				run("In [+]");
			};
			run("Scale to Fit");
	
			setThreshold(master_threshold_lower, master_threshold_upper);
		
			run("Create Selection");
			run("Save XY Coordinates...", "save=["+output_dir+"\\"+j+"_YZ_SurfacePx.csv]");
		
			saveAs("Tiff", output_dir+"\\"+"YZ"+j+"_dropsurface.tif");
		
			print(" \n**Output file containing YZ coordinates \n  of Droplet#"+j+" Surface can be found here: \n "+
			""+output_dir+"\\SurfacePx_YZ"+j+".csv");
			
			findContactAngle();
			
			setResult("slice_direction", nResults-1, "YZ");
			updateResults();
		};
		
		//////////////////////////////////////////////
		
		waitForUser("Inspect images");
		
		Dialog.create("Keep analysing?");
		Dialog.addRadioButtonGroup("Would you like to analyse another droplet from this Z-stack?", newArray("Yes","No"), 2, 1, "Yes");
		Dialog.addMessage("*If 'Yes', please WAIT for Z-stack to re-open automatically before clicking anywhere!");
		Dialog.show();
		
		continue_status = Dialog.getRadioButton();
		print(continue_status);
		if (continue_status == "No") {
			print("\n##################################################################");
			print("SUMMARY:");
			print("Number of droplets analysed = "+j);
			print("Slice direction(s) (side-on droplet profile(s)):");
			if (XZrunstatus==true) {print("  XZ");};
			if (YZrunstatus==true) {print("  YZ");};
			print("Output Directory: "+output_dir);

			selectWindow("Log");
			saveAs("text", output_dir+"\\log.txt");
			close("Log");
			selectWindow("Results");
			saveAs("results", output_dir+"\\results.csv");
			exit
		};

	};

	
} else {
////////////////////////////////////////////////////////////////////////////////////////////
///// AUTOMATED DROPLET PICKING ///////
////////////////////////////////////////////////////////////////////////////////////////////

	roiManager("reset");
	
	run("Analyze Particles...", "size="+auto_min+"-"+auto_max+" circularity="+auto_circularity+"-1.00 exclude add slice");
	
	n = roiManager('count');
	print("no. of drops to be analysed = "+n);

	if (n == 0) {
		Dialog.create("No Droplets!");
		Dialog.addMessage("No droplets meet specified criteria");
		Dialog.show();
	} else {
		Dialog.create("Automated droplet analysis");
		Dialog.addMessage(""+n+" droplets will be analysed.");
		Dialog.addRadioButtonGroup("Continue?", newArray("Yes","No"), 1, 2, "Yes");
		Dialog.show();

		Auto_continue_status = Dialog.getRadioButton();
		if (Auto_continue_status == "No") {
			exit
			}
	};


setBatchMode(true);
		
	for (k = 0; k < n; k++) {
		
		selectImage(XYZ);
	    roiManager("select", k);
	    Roi.getBounds(x, y, width, height);
	    makeRectangle(x-2, y-2, width+4, height+4);
	    run("Duplicate...", "duplicate");
			
		//re-find slice that best represents plate surface 
		//(in case plate is not perfectly flat):
		findPlateSurface();
		
		//remove all slices below this...
		plate_surface_slice_drop = getSliceNumber() - 1; //(EXCLUSIVE)
//		plate_surface_slice_drop = getSliceNumber();; //(INCLUSIVE)
		print("plate surface slice (XY) for Droplet#"+(k+1)+" = "+plate_surface_slice_drop);
			
		run("Slice Remover", "first=1 last="+plate_surface_slice_drop+" increment=1");
		setSlice(1);

		saveAs("tiff", output_dir+"\\"+(k+1)+"_XY_drop.tif");
//		close();
	}



	for (k = 0; k < n; k++) {

		open(output_dir+"\\"+(k+1)+"_XY_drop.tif");
		XYZ_crop = getImageID();

		if (XZrunstatus == true) {

			//Take side-on slices (XZ) of droplet (1 slice per pixel): 
		    run("Reslice [/]...", "output="+Vx_width+" start=Top avoid");
			
			//Find slice with max average intensity (to estimate middle of droplet):
			XYZ_crop_XZreslice = getImageID();
			run("Plot Z-axis Profile");
			Plot.getValues(z_micron, Imean);
			Array.getStatistics(Imean, min, max, mean, stdDev);
			maxLoc = Array.findMaxima(Imean, max/20, 1);
			
			//find x (XZ slice) where y (Imean) = max
			mid_slice_number = maxLoc[0]+1;
			
			selectImage(XYZ_crop_XZreslice);
			print("\n Droplet#"+(k+1));
			print("number of XZ slices = "+nSlices);		
			print("highest intensity slice = "+mid_slice_number);
			setSlice(mid_slice_number);
			
			run("Duplicate...", "use");
			XZdrop = getImageID();
			selectImage(XYZ_crop_XZreslice);
			close();
			selectImage(XZdrop);
			
			//Set Threshold and select to define surface of droplet:
			selectImage(XZdrop);
			setThreshold(master_threshold_lower, master_threshold_upper);
			run("Create Selection");
			run("Save XY Coordinates...", "save=["+output_dir+"\\"+(k+1)+"_XZ_SurfacePx.csv]");
			
			run("RGB Color");
			setForegroundColor(255, 0, 0);
			fill();
			run("Select None");
			saveAs("tiff", output_dir+"\\"+(k+1)+"_XZ_drop.tif");

			if (auto_interpolate_status == true) {
					
				//Convert pixel coordinates to micron coordinates ("post-interpolation"):
				open(output_dir+"\\"+(k+1)+"_XZ_SurfacePx.csv");
				IJ.renameResults((k+1)+"_XZ_SurfacePx.csv","Results");
				for (i = 0; i < nResults(); i++) {
				    X = getResult("X", i);
				    setResult("X_micron", i, X*Vx_width);
				    Y = getResult("Y", i);
				    setResult("Y_micron", i, Y*Vx_depth);
				}
				updateResults();
				saveAs("Results", output_dir+"\\"+(k+1)+"_XZ_SurfacePx.csv");
				run("Close");
			}
		}


		if (YZrunstatus == true) {
			
			selectImage(XYZ_crop);
			//Check correct window selected...
			if (nSlices == 1) {
				selectImage(XYZ_crop);
			}
			
			//Take side-on slices (YZ) of droplet (1 slice per pixel): 
		    run("Reslice [/]...", "output="+Vx_width+" start=Left avoid");
			
			//Find slice with max average intensity (to estimate middle of droplet):
			XYZ_crop_YZreslice = getImageID();
			run("Plot Z-axis Profile");
			Plot.getValues(z_micron, Imean);
			Array.getStatistics(Imean, min, max, mean, stdDev);
			maxLoc = Array.findMaxima(Imean, max/1000, 1);

			//find x (YZ slice) where y (Imean) = max
			mid_slice_number = maxLoc[0]+1;
			
			selectImage(XYZ_crop_YZreslice);
			print("number of YZ slices = "+nSlices);		
			print("highest intensity slice = "+mid_slice_number);
			setSlice(mid_slice_number);

			run("Duplicate...", "use");
			YZdrop = getImageID();
			selectImage(XYZ_crop_YZreslice);
			close();
			selectImage(YZdrop);
			
			//Set Threshold and select to define surface of droplet:
			setThreshold(master_threshold_lower, master_threshold_upper);
			run("Create Selection");
			run("Save XY Coordinates...", "save=["+output_dir+"\\"+(k+1)+"_YZ_SurfacePx.csv]");

			run("RGB Color");
			setForegroundColor(255, 0, 0);
			fill();
			run("Select None");
			saveAs("tiff", output_dir+"\\"+(k+1)+"_YZ_drop.tif");


			if (auto_interpolate_status == true) {
	
				//Convert pixel coordinates to micron coordinates ("post-interpolation"):
				open(output_dir+"\\"+(k+1)+"_YZ_SurfacePx.csv");
				IJ.renameResults((k+1)+"_YZ_SurfacePx.csv","Results");
				for (i = 0; i < nResults(); i++) {
				    X = getResult("X", i);
				    setResult("X_micron", i, X*Vx_width);
				    Y = getResult("Y", i);
				    setResult("Y_micron", i, Y*Vx_depth);
				}
				updateResults();
				saveAs("Results", output_dir+"\\"+(k+1)+"_YZ_SurfacePx.csv");
				run("Close");
			}
		}
	}


/////////////////////////////////////////////////////////

print("\n##################################################################");
print("SUMMARY:");
print("Voxel:\n  width = "+Vx_width+"\n  height = "+Vx_height+"\n  depth = "+Vx_depth+"\n  units = "+Vx_unit);
print("Number of droplets analysed = "+n);
print("Slice direction(s) (side-on droplet profile(s)):");
	if (XZrunstatus==true) {print("  XZ");};
	if (YZrunstatus==true) {print("  YZ");};
print("Output Directory: "+output_dir);

// Save Log to output directory:
setBatchMode(false);
selectWindow("Log");
saveAs("text", output_dir+"\\log");

	Dialog.create("Done!");
	Dialog.addMessage(""+n+" droplets have been analysed.");
	Dialog.addMessage("Results can be found here:");
	Dialog.addMessage(output_dir);
	Dialog.show();

