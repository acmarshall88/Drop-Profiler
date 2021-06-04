////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//FUNCTIONS:

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
	Array.getStatistics(Imean, min, max, mean, stdDev);
	maxLoc = Array.findMaxima(Imean, max/2);
	
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
	setCustomThreshold(); //*** custom function defined above^^
	getThreshold(lower, upper); //use these values in for loop below... 
	
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
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//open Zstack and get details
if (nImages==0) {
	waitForUser("Please open Z-stack of droplet images, then hit OK.");
} else {
	waitForUser("Ensure that Z-stack of droplet images is the *only* image window open, then hit OK. \n"+
	"(If other images are open, close them now)");
};

//set "user value" for thresholding
Dialog.create("User input");
Dialog.addNumber("(5) Droplet thresholding:  ", 15);
Dialog.addMessage("^ Number of positive standard deviations from mean of background peak of \n blank-subtracted image ('user value' - Wang et al, 'A Molecular Grammar...', Cell, 2018).\n Increasing this will increase the lower intensity threshold for defining pixels as condensed phase. \n (use ~10 for confocal images)");
Dialog.show();

user_value = Dialog.getNumber();

//sample details...
selectImage(1);
run("In [+]");
run("Scale to Fit");
Zstackfilename = getTitle();
XYZ = getImageID();
dir = getDirectory("image");
print("Directory = " + dir);
print("Z-stack = " + Zstackfilename);

//create output directory (need to put before loop!!)
output_dir = dir+Zstackfilename+"_analysis";
File.makeDirectory(output_dir);
print("Output Directory = " + output_dir);


print("\n##################################################################");
print("### Set master threshold using slice that corresponds \n### best to average plate surface... ");

findPlateSurface();
plate_surface_slice_average = getSliceNumber();

setCustomThreshold();

//Save master threshold values for reusing when defining droplet surface 
//(after extracting side-on profiles)  
getThreshold(lower, upper);
master_threshold_lower = lower;
master_threshold_upper = lower*1.8;

print("Master threshold values (lower/upper):");
print(master_threshold_lower);
print(master_threshold_upper);
print("##################################################################\n");


Dialog.create("Select Orthogonal Views");
Dialog.addMessage("Which orthogonal views would you like to include in analysis?");
Dialog.addCheckbox("XZ", true);
Dialog.addCheckbox("YZ", true);
Dialog.show();

XZrunstatus = Dialog.getCheckbox();
YZrunstatus = Dialog.getCheckbox();


/////////////////////////////////////////
///////...PUT LOOP START HERE...////////////
/////////////////////////////////////////


for (j = 1; j < 1000; j++) {
	
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
	//setSlice(1);
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
		run("Save XY Coordinates...", "save=["+output_dir+"\\SurfacePx_XZ"+j+".csv]");
	
		saveAs("Tiff", output_dir+"\\"+"XZ"+j+"_dropsurface.tif");
	
		print(" \n**Output file containing XZ coordinates \n  of Droplet#"+j+" Surface can be found here: \n "+
		""+output_dir+"\\SurfacePx_XZ"+j+".csv");
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
		run("Save XY Coordinates...", "save=["+output_dir+"\\SurfacePx_YZ"+j+".csv]");
	
		saveAs("Tiff", output_dir+"\\"+"YZ"+j+"_dropsurface.tif");
	
		print(" \n**Output file containing YZ coordinates \n  of Droplet#"+j+" Surface can be found here: \n "+
		""+output_dir+"\\SurfacePx_YZ"+j+".csv");
	};
	
	waitForUser("Inspect images");
	
	Dialog.create("Keep analysing?");
	Dialog.addRadioButtonGroup("Would you like to analyse another droplet from this Z-stack?", newArray("Yes","No"), 2, 1, "Yes");
	Dialog.show();
	
	continue_status = Dialog.getRadioButton();
	
	if (continue_status == "Yes") {
		close("*");
	} else {
		selectWindow("Log");
		saveAs("text", output_dir+"\\log.txt");
		close("Log");
			exit
	};
};
