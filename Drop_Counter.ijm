////////////////////////////////////////////////////////////////////////////////////////////
//// DROP COUNTER //////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

// INPUT PARAMETERS (**UPDATE THESE**):
	//sample protein concentration:
	protein_uM = 13; 
	
	//pathway to blank images (NB: these images must have file names: "##uM_Gblur30.tif", where ## is an integer):
	blank_directory = "//uniwa.uwa.edu.au/userhome/staff7/00101127/My Documents/LLPS results/20201112_gfp-sfpq(1-265)/BLANKS_1hrPlateII_nospin_B1-8_05peg/"+round(protein_uM)+"uM_Gblur30.tif";
	
	//'tolerance' for finding maxima in background peak of raw sample image (see 'Array.findMaxima()', LINE 102):
		//(increase this value if too many premature maxima are found in histogram) 
	tolerance = 0.1;
	
	//droplet threshold value (number of bankground peak standard deviations) ("user value", Wang et al 2018):
	user_value = 3;
		//(Increasing user_value will increase intensity threshold for defining pixels as condensed phase... shouldn't have to change)

////////////////////////////////////////////////////////////////////////////////////////////
// This macro is for calculating the relative amount 
//	of condensed protein in liquid-liquid phase separated 
//	samples imaged using fluorescence microscopy.

// MACRO PROCESS OVERVIEW:
// 1. Extract histogram/LUT (x,y = values,counts) of sample image.
// 2. Find the mean background pixel value (x) by fitting gaussian to background peak in histogram.
// 3. Use gaussian blur of image of control sample (same [protein], but no LLPS - e.g. in high [salt]) as 'blank' to subtract from sample image.
// 4. Normalise blank by multipyling it by: 0.99*sample_background_mean/blank_mean. ('mean' is mean pixel value; 99% factor is to leave some background signal to fit gaussian)
// 5. Subtract normalised blank from sample image.
// 6. Extract histogram (x,y = values,counts) of background-subtracted image.
// 7. Fit gaussian to background peak of background-subtracted image to set a reproducible/consistent threshold for all sample images.
// 8. Use threshold to define droplet mask.
// 9. Relative condensed protein is calculated using two methods:
	//	Method#1) total integrated intensity inside droplets divided by total integrated intensity outside droplets (Idroplet/Imedia)
		//	*assumes no saturated pixels... also assumes background subtraction won't affect(?) 
	//	Method#2) multiplying average droplet volume (caluculated from average droplet area) by total number of droplets 
		//	*assumes droplets are approximately hemispherical

////////////////////////////////////////////////////////////////////////////////////////////
// BEFORE RUNNING MACRO YOU NEED TO MAKE 'BLANK' IMAGES FIRST AND SPECIFY LOCATION ABOVE (LINE 10)
// *take images of samples at ~ respective [protein] but in conditions prohibitive of LLPS (at least with no big drops)
	// 1. open image of desired [protein] sample in series (should have no large droplets)
	// 2. Process -> Filters -> Gaussian Blur... 
	// 3. Set 'Sigma (Radius)' = 30.0
	// 4. Save As... TIFF (INCLUDE [protein] in filename)
	// 5. repeat for samples to cover [protein] range
	// 6. Enter file location above (LINE 10)


////////////////////////////////////////////////////////////////////////////////////////////
// TO RUN MACRO, FIRST OPEN SAMPLE IMAGE, THEN HIT 'RUN' (ctrl+R)...
////////////////////////////////////////////////////////////////////////////////////////////
//selectImage(1);

// starts log:
print(" ");
print("###########################################################################");
print("################################# NEW SAMPLE ##############################");
print("###########################################################################");
print("Sample [protein]: "+protein_uM+" uM");
print(" ");

// Assigns 'sample' to open (selected) image:
sample = getTitle();


////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 1. Extract histogram/LUT of raw sample image... ###                        ");
////////////////////////////////////////////////////////////////////////////////////////////

// Defines number of bins for histogram:
getRawStatistics(nPixels, mean, min, max, std, histogram);
pixel_depth = max-min;
print("pixel_depth = "+ pixel_depth);
nBins = pixel_depth/pow(2, -floor(-pixel_depth*0.0001)); //(nBins should be ~500-1000 and a factor of pixel_depth)
print("nBins = "+ nBins);

// Gets values from the LUT histogram. 
// This returns two arrays "values" and "counts".
// The histogram is split into into nBins number of bins. 

getHistogram(values, counts, nBins);
plot_values = values;
plot_counts = counts;

Plot.create(sample+" histogram", "values", "counts", plot_values, plot_counts);

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 2. Find the mean background pixel value (x) by fitting gaussian to ###");
print("###     background peak (left-most peak) in histogram of raw image... ###");
////////////////////////////////////////////////////////////////////////////////////////////

// "Array.getStatistics()" gets statistics from the 
// arrays that were extracted using getHistogram.
// Only really useful for the counts array.

// First, this removes first and last values in array 
// to remove these potential 'maxima'...

values = Array.slice(values, 1, nBins-1);
counts = Array.slice(counts, 1, nBins-1);
 
Array.getStatistics(counts, min, max, mean, stdDev);
print("mean_counts = "+mean);
print("max_counts = "+max);
print("min_counts = "+min);

// Need to locate max value of background peak to define region to fit gaussian... 
// Finds maxima:
maxima_list = Array.findMaxima(counts, max*tolerance);
maxima_list_sorted = Array.sort(maxima_list);
dummy_array = newArray(nBins-4, nBins-4, nBins-4, nBins-4, nBins-4, nBins-4, nBins-4, 
	nBins-4, nBins-4, nBins-4, nBins-4, nBins-4, nBins-4, nBins-4, nBins-4, nBins-4);
 maxLoc = Array.concat(maxima_list,dummy_array);
		//(^adding dummy array is so that for loop (below) will work even if 
		// no premature maxima are defined)
		
print("maxima locations (bins): ") Array.print(maxLoc);

//Excludes "premature" maxima (outliers) and define peak location (bin#)...

if (maxLoc[0] == 0) {
	peak_x_location = maxLoc[1];
} else {
	peak_x_location = maxLoc[0];
}

for (i = 0; i < 15; i++) {
	if ((peak_x_location == maxLoc[i]) == true
			&&	(maxLoc[i+1]-maxLoc[i] <= nBins*0.05) 
			&& (counts[maxLoc[i+1]] >= counts[maxLoc[i]]*0.9) 
			== true) {
		peak_x_location = maxLoc[i+1];
};
};

 print("bkgnd maximum bin (peak_x_location): " + peak_x_location);
 peak_x = values[peak_x_location]; 
 print("bkgnd maximum, pixel value (peak_x) = " + peak_x);
 maxCount = counts[peak_x_location];
 print("bkgnd maximum, counts (maxCount) = " + maxCount);
peak_x_bin = plot_values[peak_x_location];

// Defines maximum cutoff for the fitting of the gaussian 
// to the left-most peak in the histogram:

// max_cutoff = peak_x_location * 1.1; 
	//(^includes left half of distribution only)
max_cutoff = peak_x_location * 1.5; 
	//(^excludes right shoulder of distribution) 

// makes trimmed dataset...
bg_x_values = Array.slice(values, 0, max_cutoff);
bg_y_values = Array.slice(counts, 0, max_cutoff);

// Fits a given function ("Gaussian (no offset)" in 
// this case) to the trimmed dataset using "Fit.doFit()".
// Fit.plot() then prints the fitted plot into
// a new window...

//using built-in function, "Gaussian (no offset)"... function is "y = a*exp(-(x-b)*(x-b)/(2*c*c)))"...
SD_initialguess = peak_x*0.03; 
print("SD_initialguess = "+ SD_initialguess);
initialGuesses = newArray(maxCount, peak_x, SD_initialguess);
Fit.doFit("Gaussian (no offset)", bg_x_values, bg_y_values, initialGuesses);
Fit.plot();

rename(sample+" histogram background peak fit raw");

// Extracts fitted parameters from the gaussian fit using "Fit.p()"...
fit_max = Fit.p(0); //'a' in "Gaussian (no offset)"
sample_background_mean = Fit.p(1); //'b' in "Gaussian (no offset)"
fit_sd = abs(Fit.p(2)); //'c' in "Gaussian (no offset)"

print("Raw Image Background Fit (Gaussian):");
print("sample_background_max = "+fit_max);
print("*sample_background_mean* = "+sample_background_mean);
print("sample_background_SD = "+fit_sd);

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 3. Use gaussian blur of image of control sample (same [protein], but no ###");
print("###     LLPS - e.g. in high [salt]) as 'blank' to subtract from sample image... ###");
////////////////////////////////////////////////////////////////////////////////////////////

// Opens 'blank' image representative of background to subtract: 
open(blank_directory);

blank = getTitle();
getRawStatistics(nPixels, mean, min, max, std, histogram);
blank_mean = mean;
print("blank_mean = "+blank_mean);

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 4. Normalise blank to sample_background_mean...  ###");
////////////////////////////////////////////////////////////////////////////////////////////

normalisation_factor = 0.99*(sample_background_mean/blank_mean); 
	//(^99% factor is so that there remains a background normal distribution to fit)
print("normalisation_factor = " + normalisation_factor);

run("Multiply...", "value="+ normalisation_factor);

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 5. Subtract normalised blank from sample image...  ###");
////////////////////////////////////////////////////////////////////////////////////////////

imageCalculator("Subtract create",sample,blank);

// removes noise:
run("Median...", "radius=2");

// Assigns 'bgsubtracted_sample_ID' to new image:
bgsubtracted_sample_ID = getImageID();

rename(""+sample+" after subtracting blank");
print("outputs new image called: "+getInfo("window.title"));

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 6. Extract histogram (x,y = values,counts) of background-subtracted image...###");
////////////////////////////////////////////////////////////////////////////////////////////

// Defines number of bins for histogram:
getRawStatistics(nPixels, mean, min, max, std, histogram);
pixel_depth = max-min;
print("pixel_depth (bgsubtracted) = "+ pixel_depth);
nBins = pixel_depth/pow(2, -floor(-pixel_depth*0.0001)); 
	//(^nBins should be ~500-1000 and a factor of pixel depth)
print("nBins (bgsubtracted) = "+ nBins);

// Gets values from the LUT histogram. 
// This returns two arrays "values" and "counts".
// The histogram is split into into nBins number of bins. 

getHistogram(values, counts, nBins);

plot_values = Array.slice(values, 0, nBins);
plot_counts = Array.slice(counts, 0, nBins);

Plot.create(bgsubtracted_sample_ID+" histogram", "values", "counts", plot_values, plot_counts);

// Gets the statistics from the 
// arrays that were extracted using getHistogram.
// Only really useful for the counts array.
Array.getStatistics(counts, min, max, mean, stdDev);
//print("bg_sub_mean_counts = "+mean);

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 7. Fit Gaussian to background peak of background-subtracted image to  ###");
print("###     set a reproducible/consistent threshold for all sample images...  ###");
////////////////////////////////////////////////////////////////////////////////////////////

// Makes trimmed dataset using fitted SD from original gaussian
// (pre-subtracted image - line 166) to estimate appropriate max cutoff 
// for fitting Gaussian to BG peak in histogram of BG-subtracted image
//... (*exclude ZEROth point*):
new_bg_x_values = Array.slice(values, 1, fit_sd); //'values' array from line 
new_bg_y_values = Array.slice(counts, 1, fit_sd);

// Fits a given function ("customGaussian" in this case) to the 
// trimmed dataset using "Fit.doFit()".
// Fit.plot() then prints the fitted plot into a new window...

	//Using previously fitted max and SD parameters from above (re-fit mean only):
		initialGuesses = newArray(2);
			//(^sets initial guesses for offset 'a' and mean 'b' as 0 and 0)
		//customGaussian = "y = "+ fit_max +"*exp(-(x-a)*(x-a)/(2*"+ fit_sd +"*"+ fit_sd +"))";
			//(^no y offset; i.e. y min = 0)
		customGaussian = "y = a + ("+ fit_max +"-a)*exp(-(x-b)*(x-b)/(2*"+ fit_sd +"*"+ fit_sd +"))"; 
			//(^WITH OFFSET; i.e. y min = 'a')
		Fit.doFit(customGaussian, new_bg_x_values, new_bg_y_values, initialGuesses);
		
	// Extract fitted parameters...
	fit_offset = Fit.p(0);	//('a' above)
	refit_mean = Fit.p(1);	//('b' above)

// re-make trimmed dataset... (*exclude ZEROth point*)
	new_bg_x_values = Array.slice(values, 1, refit_mean+fit_sd*1.5);
	new_bg_y_values = Array.slice(counts, 1, refit_mean+fit_sd*1.5);

	//Do second round of fitting to tweak max and SD (mean is fixed):
		initialGuesses = newArray(fit_offset, fit_max, fit_sd);
		customGaussian3 = "y = a + (b-a)*exp(-(x-("+ refit_mean +"))*(x-("+ refit_mean +"))/(2*c*c))";
			//(^WITH OFFSET; i.e. y min = 'a')
		Fit.doFit(customGaussian3, new_bg_x_values, new_bg_y_values, initialGuesses);

	// Extract fitted parameters...
	refit_offset = Fit.p(0); 	//('a' above)
	refit_max = Fit.p(1); 		//('b' above)
	refit_sd = Fit.p(2); 		//('c' above)

Fit.plot();

rename(sample+" histogram background peak fit after blank subtraction");

print("Y-OFFSET_fit = " + refit_offset);
print("MAX_fit = " + refit_max);
print("MEAN_fit = " + refit_mean);
print("SD_fit = " + refit_sd);

// Sets the new threshold for the image ("user value" (arbitrary) is defined in line 17)
// using mean and SD parameters of Gaussian fit above...
new_thresh = refit_mean + abs(refit_sd)*user_value; 
print("threshold = " + new_thresh);

// Creates data from the fitted gaussian for creation of plot...
fitted_counts = Array.copy(plot_counts);

for (i = 0; i < plot_values.length; i++) {
	fitted_counts[i] = Fit.f(values[i]);
}

// Creation of the plot for visual represenation of the data...
 Plot.create("Total Image Histogram", "Pixel Values", "Counts");
 Plot.setColor("red");
 Plot.setLineWidth(5);
 Plot.add("dot", plot_values, plot_counts);
 Plot.setLineWidth(2);
 Plot.setColor("cyan");
 Plot.add("line", plot_values, fitted_counts);
 Plot.setColor("red");
Array.getStatistics(new_bg_y_values, min, max, mean, stdDev);
maxCount = max;
 Plot.drawLine(new_thresh, maxCount/2, new_thresh, 0);
 Plot.addText("Threshold Value", new_thresh, maxCount/2);
 Plot.setColor("black");
 Plot.show;

rename(sample+" thresholded histogram");

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 8. Use threshold to define droplet mask...  ###");
////////////////////////////////////////////////////////////////////////////////////////////

selectImage(bgsubtracted_sample_ID);

getStatistics(area, mean, min, max, std, histogram);
setThreshold(new_thresh, max);

// Creates selection from the threshold for the droplets...
run("Create Selection");

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 9. Calculate amount of condensed protein...  ###");
////////////////////////////////////////////////////////////////////////////////////////////

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 9a. Method #1 - compare integrated intensities (relative)...  ###");
////////////////////////////////////////////////////////////////////////////////////////////

// Gets statistics for the condensed phase (inside droplets)...
getStatistics(area, mean, min, max, std);
	Idroplet=mean*area; // integrated intensity inside droplets
	drop_area=area;

// Changes selection to the background (everything else)...
run("Make Inverse");

// Get statistics for the dilute phase (outside droplets)...
getStatistics(area, mean, min, max, std);
	Imedia=mean*area; // integrated intensity outside droplets
	back_area=area;

// proportion of condensed protein BY AREA =
proportion_cond_area=drop_area/(back_area+drop_area);

// First, so that new results will be appended to current results table... 
//selectWindow("Method#1 (compare integrated intensities)"); IJ.renameResults("Results"); 
///////////////////////////////////////////
//***COMMENT OUT ^THIS LINE^ (ABOVE) FIRST TIME***
///////////////////////////////////////////

row=nResults;

setResult("prot_conc_uM", row, protein_uM);
setResult("condensed_area", row, drop_area);
setResult("dilute_area", row, back_area);
setResult("proportion_condensed_area", row, proportion_cond_area);
setResult("Idroplet", row, Idroplet);
setResult("Imedia", row, Imedia);
setResult("Idroplet/Imedia", row, Idroplet/Imedia);
 
 IJ.renameResults("Method#1 (compare integrated intensities)");

print(" ");
////////////////////////////////////////////////////////////////////////////////////////////
print("###########################################################################");
print("### 9b. Method #2 - calculate total condensed volume (estimate of absolute)... ###");
////////////////////////////////////////////////////////////////////////////////////////////

// Clears selection... (image remains thresholded)...
run("Select None");

// "run("Analyze Particles...")" counts droplets and calculates average area...

run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 display summarize");

//run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Outlines display summarize");
// IJ.renameResults("Particle list for "+protein_uM+"uM sample (Method#2)");
	//(^use this to display list of particles in results window)

// Calculates theoretical condensed volume at well bottom
// using average droplet area ("Average Size") to calculate
// average droplet volume (ASSUMING A HEMISPHERE), and then
// multiplying by number of droplets ("Count")...
selectWindow("Summary"); IJ.renameResults("Results"); 

// Gets "Count" and "Average Size" values from last row of Results table...
Count = getResult("Count", nResults-1);
Average_Size = getResult("Average Size", nResults-1);

// Rearrange formula for volume-of-a-sphere to calculate 
// vol (in cubic microns) from area ("Average Size")... 
Av_calculated_sphere_Vol = (4/(3*sqrt(PI)))*pow(Average_Size, 3/2);
// Halve for hemisphere assumption...
Av_calculated_hemisphere_Vol = 0.5*Av_calculated_sphere_Vol;
// Multiply by number of droplets...
Total_calculated_Vol_um3 = Av_calculated_hemisphere_Vol*Count;

// Convert to nanoliters...

Total_calculated_Vol_nL = Total_calculated_Vol_um3/pow(10,6);
	//(^ 1 nanoliter = 1,000,000 cubic microns)

// Interior dimensions of well bottom are ~ 3.3 x 3.3 mm.
print("Well bottom dimensions ~ 3.3 x 3.3 mm");
well_area_um2 = pow(3300,2);
print("Well bottom area = "+well_area_um2+" microns^2");

// Calculates area covered by image...
selectImage(bgsubtracted_sample_ID);
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);
image_width = width*pixelWidth;
image_height = height*pixelHeight;
print("Image dimensions = "+image_width+" x "+image_height+" "+unit);
image_area_um2 = image_width*image_height;
print("Image area = "+image_area_um2+" "+unit+"^2");

// Calculates factor to multiply calculated Volume by
// to extrapolate to total condensed volume in sample...
extrapolation_factor = well_area_um2/image_area_um2;
print("extrapolation_factor = "+extrapolation_factor);

Total_calculated_Vol_in_sample_nL = Total_calculated_Vol_nL*extrapolation_factor;

// Append to Results Table...
setResult("Condensed Volume (um^3)", nResults-1, Total_calculated_Vol_um3);
setResult("Cond Vol in image (nL)", nResults-1, Total_calculated_Vol_nL);
setResult("Total Cond Vol in sample (nL)", nResults-1, Total_calculated_Vol_in_sample_nL);

IJ.renameResults("Method#2 (calculate total condensed volume)");

print(" "); 