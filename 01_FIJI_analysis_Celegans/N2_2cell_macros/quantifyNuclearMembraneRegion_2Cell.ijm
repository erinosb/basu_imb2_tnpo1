//////////////////////////////////////
///    SELECT YOUR PARAMETERS      ///
//////////////////////////////////////

// DATE:     Started May 27, 2026; modified July 2, 2026
// AUTHOR:   Erin Osborne Nishimura
// SCRIPT:   quantifyNuclearMembraneRegion.ijm
// PURPOSE:  This script takes as input a folder of .dv microscopy imgage files. Files have 3 or 4 channels and multiple z-projections. 
//           This script then does the following:
//              - Creates a z-projection of a fixed number of slices (user input required to set a start slice for each image)
//              - Rotates the image (user input required)
//              - Cropped to a given rectangle size (user input required to place the rectangle)
// 	            - Color split 
//              - Brightness adjusted by user specified values or optimal values determined from the first image analyzed
//              - Color merge based on user parameters
//              - Scale bar added
//              - Saves .tif, .png, and inverted.png as output
//              - Writes a logfile for each run

// Please change this part to your input and output directories and desired contact sheet file name

// DIRECTORIES: Enter the input and output directories. 

  // Please replace <insertInputDirHere> with your own input directory. Remove the < and > symbols, too. 
  // This directory must exist and must contain a series of image files that match the image file type in Line 127.
inputDir="/Users/ambikabasu/Nishimura_lab/Macros/QuantificationNuclearMembrane /2Cell/01_input";

  // Please replcce <insertOutputDirHere> with your own output directory. Remove the < and > symbols, too.
  // This directory does not need to exist. It will be created.
outputDir="/Users/ambikabasu/Nishimura_lab/Macros/QuantificationNuclearMembrane /2Cell/03_output";

// Z-PROJECT: Specify how many z-slices will be selected for z-projection. Default is 35:
zSliceDepth = 7;

// CHANNELS: How many channels do you have? Default is 4, but 3 can also work. Change this value to 3 if necessary:
channelNo=3;

// SAMPLES: Specify what is being measured in each channel
channel1_sample = "imb-2";
channel2_sample = "set-3";
channel3_sample = "";
channel4_sample = "";

// MAX/MIN Brightness - OPTIONAL: Set the channel max and min to adjust brightness if you like. Leave at 0 if you want it automatically determined. 
// Automatically determined brightness settings will determine min and max settings that yield a contrast of 0.35 in the first image
// All downstream images will use that same min and max channel setting. 
ch1min = 0;
ch1max = 0;
ch2min = 0;
ch2max = 0;
ch3min = 0;
ch3max = 0;
ch4min = 0;
ch4max = 0;

//////////////////////////////////////////////
///     PARAMETERS SELECTION COMPLETE      ///
//////////////////////////////////////////////



//////////////////////////
///    START CODE      ///
//////////////////////////

// Initialize Logfile
print("\\Clear");

// Get todays date:
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
month = month + 1;
todaysDate = d2s(year, 0) + "-" + d2s(month, 0) + "-" + d2s(dayOfMonth, 0);
timeNow  = d2s(hour, 0) + d2s(minute, 0);

// Create the output directory
File.makeDirectory(outputDir);

// Initialize logfile:
// Choose output directory
logPath = outputDir + "/" + todaysDate + "_" + timeNow + "_" + "logfile.txt";
dataPath = outputDir + "/" + todaysDate + "_" + timeNow + "_" + "datafile.txt";

// set the working directory
File.setDefaultDir(inputDir);

// get the list of files in this directory
fileList = getFileList(inputDir);

// Record time
logfile = File.open(logPath);
print("Initiaiting macro: quantifyNuclearMembraneRegion.ijm");
File.append("Initiaiting macro: quantifyNuclearMembraneRegion.ijm", logPath);
print("\nDate: " + todaysDate);
print("Time: " + timeNow);
File.append("\nDate: " + todaysDate, logPath);
File.append("Time: " + timeNow, logPath);

// Record the input and output directory. 
print("\nInput Directory: " + inputDir);
print("Output Directory: " + outputDir + "\n");
File.append("\nInput Directory: " + inputDir, logPath);
File.append("Output Directory: " + outputDir + "\n", logPath);


// restrict list to R3D_D3D.dv files
dv_array = newArray(0);

// Save .dv files in a dv_array:
for (i=0; i<fileList.length; i++) {
    if(endsWith(fileList[i], "R3D.dv"))
    {
    	dv_array = Array.concat(dv_array, fileList[i]);
    }
}

// Print the list of dv files to process:
print("\nWill process the following DV files:");
File.append("\nWill process the following DV files:", logPath);
for (i=0; i<dv_array.length; i++) {
    print("\t" + dv_array[i]);
    File.append("\t" + dv_array[i], logPath);
}

// cycle through the images in the directory
for (i=0; i<dv_array.length; i++) {
	
	// Get file name 
    full_dv = inputDir + "/" + dv_array[i];
    
    // Report 
    print("\nNow processing file:\t" + dv_array[i]);
    File.append("\nNow processing file:\t" + dv_array[i], logPath);
    
    // Open the .dv file
    run("Bio-Formats Importer", "open=full_dv");
    
    // set the zoom
    run("Set... ", "zoom=75 x=512 y=512");
    
    // Enhance the contrast
    run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(2);
    run("Enhance Contrast", "saturated=0.35");
    Stack.setChannel(3);
    run("Enhance Contrast", "saturated=0.35");
    if (channelNo == 4) {
    	Stack.setChannel(4);
    	run("Enhance Contrast", "saturated=0.35");
    }
    
    // Z project
    print("\t--Starting Z-projection");
    File.append("\t--Starting Z-projection", logPath);
    
    // Get Z slices
    getDimensions(width, height, channels, slices, frames);
    print("\t--There are " + slices + " slices");
    File.append("\t--There are " + slices + " slices", logPath);
    
    // Wait for user to determine a good Z-range:
    waitForUser("Action Required", "By default, 7 slices will be selected for Z-projection. \nPlease look through the image and select the middle slice. \nClick OK when ready.");
    
    // Create a dialog box for user input
    Dialog.create("Specify Z Projection Range");
    Dialog.addNumber("Start slice:", 1);
    Dialog.addString("Projection Method (e.g., Max, Average):", "Sum Slices");
    Dialog.show();

    // Get the user input values
    midSlice = Dialog.getNumber();
    startSlice = midSlice - 3;
    stopSlice = midSlice + 3;
    method = Dialog.getString();
    SliceDepth = stopSlice - startSlice;
    
    // Determine whether to take 35 slices or the maximum number of slices
    print("\t--Splice Depth is: " + zSliceDepth);
    File.append("\t--Splice Depth is: " + zSliceDepth, logPath);
    
    // Run the Z Project function with user-defined parameters
    run("Z Project...", "start=" + startSlice + " stop=" + stopSlice + " projection=[" + method + "]");

    // set the zoom
    run("Set... ", "zoom=75 x=512 y=512");
    
    // Report Z Projection
    print("\t--Z-projection complete from " + startSlice + " to " + stopSlice + " using " + method + " method.");
	File.append("\t--Z-projection complete from " + startSlice + " to " + stopSlice + " using " + method + " method.", logPath);
	
     if (i == 0){
    	print("\t--Evaluating Min and Max Thresholds");
    	File.append("\t--Evaluating Min and Max Thresholds", logPath);
		
		// If ch1min or ch1max are set to 0, get values from auto brightness
		if (ch1min == 0){
    		Stack.setChannel(1);
    		run("Enhance Contrast", "saturated=0.35");
    		getMinAndMax(min, max);
    		ch1min = min;
    		print("\t--Setting ch1min:" + ch1min);
    		File.append("\t--Setting ch1min:" + ch1min, logPath);
    	}
    	if (ch1max == 0){
    		Stack.setChannel(1);
    		run("Enhance Contrast", "saturated=0.35");
    		getMinAndMax(min, max);
    		ch1max = max;
    		print("\t--Setting ch1max:" + ch1max);
    		File.append("\t--Setting ch1max:" + ch1max, logPath);
    	}
		
		// If ch2min or ch2max are set to 0, get values from auto brightness
    	if (ch2min == 0){
    		Stack.setChannel(2);
    		run("Enhance Contrast", "saturated=0.35");
    		getMinAndMax(min, max);
    		ch2min = min;
    		print("\t--Setting ch2min:" + ch2min);
    		File.append("\t--Setting ch2min:" + ch2min, logPath);
    	}
    	if (ch2max == 0){
    		Stack.setChannel(2);
    		run("Enhance Contrast", "saturated=0.35");
    		getMinAndMax(min, max);
    		ch2max = max;
    		print("\t--Setting ch2max:" + ch2max);
    		File.append("\t--Setting ch2max:" + ch2max, logPath);
    	}
		
		// If ch3min or ch3max are set to 0, get values from auto brightness
    	if (ch3min == 0){
    		Stack.setChannel(3);
    		run("Enhance Contrast", "saturated=0.35");
    		getMinAndMax(min, max);
    		ch3min = min;
    		print("\t--Setting ch3min:" + ch3min);
    		File.append("\t--Setting ch3min:" + ch3min, logPath);
    	}
    	if (ch3max == 0){
    		Stack.setChannel(3);
    		run("Enhance Contrast", "saturated=0.35");
    		getMinAndMax(min, max);
    		ch3max = max;
    		print("\t--Setting ch3max:" + ch3max);
    		File.append("\t--Setting ch3max:" + ch3max, logPath);
    	}
		
		// Check if there is a channel4. If there is, check if the ch1min and ch1max are set to 0. 
		if (channelNo == 4){
    	    if (ch4min == 0){
    			Stack.setChannel(4);
    			run("Enhance Contrast", "saturated=0.35");
    			getMinAndMax(min, max);
    			ch4min = min;
    			print("\t--Setting ch4min:" + ch4min);
    			File.append("\t--Setting ch4min:" + ch4min, logPath);
    		}
    		if (ch4max == 0){
    			Stack.setChannel(4);
    			run("Enhance Contrast", "saturated=0.35");
    			getMinAndMax(min, max);
    			ch4max = max;
    			print("\t--Setting ch4max:" + ch4max);
    			File.append("\t--Setting ch4max:" + ch4max, logPath);
    		}
		}
	}
    

   // Set brightness in channel1
    Stack.setChannel(1);
    setMinAndMax(ch1min, ch1max);
    print("\t--Setting Channel1 min: " + ch1min);
	File.append("\t--Setting Channel1 min: " + ch1min, logPath);
    print("\t--Setting Channel1 max: " + ch1max);
	File.append("\t--Setting Channel1 max: " + ch1max, logPath);

    // Set brightness in channel2
    Stack.setChannel(2);
    setMinAndMax(ch2min, ch2max);
    print("\t--Setting Channel2 min: " + ch2min);
	File.append("\t--Setting Channel2 min: " + ch2min, logPath);
    print("\t--Setting Channel2 max: " + ch2max);
	File.append("\t--Setting Channel2 max: " + ch2max, logPath);

    // Set brightness in channel3
    Stack.setChannel(3);
    setMinAndMax(ch3min, ch3max);
    print("\t--Setting Channel3 min: " + ch3min);
	File.append("\t--Setting Channel3 min: " + ch3min, logPath);
    print("\t--Setting Channel3 max: " + ch3max);
	File.append("\t--Setting Channel3 max: " + ch3max, logPath);
    
    // If there is a channel4, set its brightness
    if (channelNo == 4);
    {	
    	Stack.setChannel(4);
    	setMinAndMax(ch4min, ch4max);
    	print("\t--Setting Channel4 min: " + ch4min);
    	File.append("\t--Setting Channel4 min: " + ch4min, logPath);
    	print("\t--Setting Channel4 max: " + ch4max);
    	File.append("\t--Setting Channel4 max: " + ch4max, logPath);
	}
    


    // Determine nuclear axis from user-clicked middle-of-nucleus and plasma membrane
    // (image pixel data is never rotated - only the rectangle SELECTION is rotated,
    // so raw intensity values are preserved exactly, unaffected by interpolation)
    print("\t--Starting nuclear membrane axis determination");  
    File.append("\t--Starting nuclear membrane axis determination", logPath);
    
    // This whole block (point-clicking through rectangle placement) repeats until
    // the user confirms the rectangle looks right. Only this block redoes - Z-projection
    // and brightness settings above are NOT reset between redo attempts.
    apConfirmed = false;
    redoAttempt = 0;
    do {
    	redoAttempt = redoAttempt + 1;
    	if (redoAttempt > 1) {
    		print("\t--Redo attempt #" + redoAttempt + " for nuclear midpoint/rectangle");
    		File.append("\t--Redo attempt #" + redoAttempt + " for nuclear midpoint/rectangle", logPath);
    	}
    	
    	setTool("multipoint");
    	waitForUser("Action Required", "Using the Point tool, click ONCE on the middle of the nucleus, then ONCE on the plasma membrane in the direction of interest (2 points total, in that order). Click OK when both points are placed.");
    	
    	getSelectionCoordinates(apXpoints, apYpoints);
    	if (apXpoints.length != 2) {
    		exit("Expected exactly 2 points (nuclear midpoint, plasma membrane) but got " + apXpoints.length + ". Please re-run and place exactly 2 points.");
    	}
    	
    	anteriorX = apXpoints[0];
    	anteriorY = apYpoints[0];
    	posteriorX = apXpoints[1];
    	posteriorY = apYpoints[1];
    	
    	// Angle of the AP axis, for logging only (not used to build the rectangle below,
    	// which is built directly from the two clicked points via makeRotatedRectangle)
    	myangle = atan2(posteriorY - anteriorY, posteriorX - anteriorX) * 180 / PI;
    	
    	// Report angle and clicked points
    	print("\t--Nuclear midpoint: " + anteriorX + ", " + anteriorY);
    	File.append("\t--Nuclear midpoint: " + anteriorX + ", " + anteriorY, logPath);
    	print("\t--PM point: " + posteriorX + ", " + posteriorY);
    	File.append("\t--PM point: " + posteriorX + ", " + posteriorY, logPath);
    	print("\t--Axis angle: " + myangle + " degrees (image pixel data not rotated)");
    	File.append("\t--Axis angle: " + myangle + " degrees (image pixel data not rotated)", logPath);
    	
    	// Build a rotated rectangle SELECTION (not a rotated image) with the same
    	// dimensions as the original fixed rectangle (333 x 69), centered on the AP
    	// axis defined by the two clicked points. makeRotatedRectangle(x1,y1,x2,y2,width)
    	// takes a centerline from (x1,y1) to (x2,y2) as the rectangle's LONG axis, and
    	// 'width' as the perpendicular (short) dimension - so the centerline length
    	// must equal the original long dimension (333) and width must equal the
    	// original short dimension (69), regardless of how far apart the two clicked
    	// points actually are.
    	print("\t--Get Region of Interest (ROI) Selection"); 
    	File.append("\t--Get Region of Interest (ROI) Selection", logPath);
    	rectLength = 117; // long dimension, same as original rectangle width
    	rectThickness = 46; // short dimension, same as original rectangle height
    	
    	// Unit vector along the AP axis (from anterior to posterior)
    	apDeltaX = posteriorX - anteriorX;
    	apDeltaY = posteriorY - anteriorY;
    	apDist = sqrt(apDeltaX*apDeltaX + apDeltaY*apDeltaY);
    	apUnitX = apDeltaX / apDist;
    	apUnitY = apDeltaY / apDist;
    	
    	// Centerline of length rectLength, centered at the midpoint of the two clicked points,
    	// oriented along the AP axis
    	midX = (anteriorX + posteriorX) / 2;
    	midY = (anteriorY + posteriorY) / 2;
    	lineX1 = midX - (rectLength/2) * apUnitX;
    	lineY1 = midY - (rectLength/2) * apUnitY;
    	lineX2 = midX + (rectLength/2) * apUnitX;
    	lineY2 = midY + (rectLength/2) * apUnitY;
    	
    	makeRotatedRectangle(lineX1, lineY1, lineX2, lineY2, rectThickness);
    	waitForUser("Action Required", "Move the rotated rectangle to capture the nuclear membrane in the middle. Click OK to measure.");
    	
    	// Report selection bounds
    	// NOTE: getSelectionBounds returns the axis-aligned BOUNDING BOX of the rotated
    	// rectangle, not its true length/thickness. The rectangle's actual dimensions remain
    	// rectLength x rectThickness (333 x 69) regardless of rotation; bounding-box x/y/width/height
    	// below are logged for reference only and should not be interpreted as the rectangle's
    	// own dimensions when myangle is not a multiple of 90.
    	getSelectionBounds(x, y, width, height);
    	print("\t--Selection of ROI complete. Rectangle dimensions: " + rectLength + " x " + rectThickness + ", AP angle " + myangle + " degrees. Bounding box: " + x + ", " + y + ", " + width + ", " + height); 
    	File.append("\t--Selection of ROI complete. Rectangle dimensions: " + rectLength + " x " + rectThickness + ", AP angle " + myangle + " degrees. Bounding box: " + x + ", " + y + ", " + width + ", " + height, logPath);
    	
    	// Confirm before proceeding. "Cancel" exits the whole macro (getBoolean default
    	// behavior); "No" redoes the point-clicking and rectangle placement above;
    	// "Yes" continues on to saving/measurement below.
    	apConfirmed = getBoolean("Does the rectangle look correctly placed and oriented along the nuclear membrane?", "Yes, continue", "No, redo points/rectangle");
    	if (!apConfirmed) {
    		print("\t--User requested redo of points/rectangle");
    		File.append("\t--User requested redo of points/rectangle", logPath);
    	}
    } while (!apConfirmed);
    

	// Save a .tif of the full image with the rotated ROI box drawn on it as an overlay.
	// run("Add Selection...") burns the rectangle into the image's overlay (does not
	// change pixel values). saveAs() renames the active window to match the saved
	// filename, so the original title is restored immediately after via rename(),
	// since the SUM_/MAX_/AVERAGE_ title is needed further down for Plot Profile.
	print("\t--Saving ROI overlay image (.tif)");
	File.append("\t--Saving ROI overlay image (.tif)", logPath);
	
	originalTitle = getTitle();
	run("Add Selection...");
	baseName = File.getNameWithoutExtension(dv_array[i]);
	tifOutPath = outputDir + "/" + baseName + "_ROI.tif";
	saveAs("Tiff", tifOutPath);
	rename(originalTitle);
	print("\t--Saved: " + tifOutPath);
	File.append("\t--Saved: " + tifOutPath, logPath);


	// Get values from Ch01
    print("\t--Select Channel 01"); 
    File.append("\t--Select Channel 01", logPath);
	Stack.setChannel(1);
	roiManager("Add");
	run("Plot Profile");
	
	xarray = newArray(0);
	yarray = newArray(0);
	Plot.getValues(xarray, yarray);
	
	buffer = "";

	// Now you can iterate over the arrays and save the data in the datafile
	//File.close(logfile);
	//datafile = File.open(dataPath);
	
	// Save distances if this is the first loop
	if (i == 0){
		buffer = "file\tchannel\tsample\t";
		for (j=0; j<xarray.length; j++) {
			buffer = buffer + xarray[j] + "\t";
		}
		File.append(buffer, dataPath);
	}
	
	// Save ch1 values
    print("\t--Reporting Channel 01 values: " + channel1_sample); 
    File.append("\t--Reporting Channel 01 values: " + channel1_sample, logPath);
	buffer = "";
	buffer = dv_array[i] + "\tch1\t" + channel1_sample + "\t";
	for (j=0; j<xarray.length; j++) {
    	buffer = buffer + yarray[j] + "\t";
	}
	File.append(buffer, dataPath);

	// Close plot
	selectWindow(getTitle());
	close();
	
    print("\t--Select Channel 02"); 
    File.append("\t--Select Channel 02", logPath);
    
	if (method == "Sum Slices"){
		selectWindow("SUM_" + dv_array[i]);
	} else if ( method == "Max"){
		selectWindow("MAX_" + dv_array[i]);		
	} else if ( method == "Average"){
		selectWindow("AVERAGE_" + dv_array[i]);		
	}

	Stack.setChannel(2);
	makeRectangle(x, y, width, height);

	roiManager("Add");
	run("Plot Profile");

	ch2xarray = newArray(0);
	ch2yarray = newArray(0);
	Plot.getValues(ch2xarray, ch2yarray);
	
	buffer = "";

	// Save ch2 values
    print("\t--Reporting Channel 02 values: " + channel2_sample); 
    File.append("\t--Reporting Channel 02 values: " + channel2_sample, logPath);
	buffer = "";
	buffer = dv_array[i] + "\tch2\t" + channel2_sample + "\t";
	for (j=0; j<ch2yarray.length; j++) {
    	buffer = buffer + ch2yarray[j] + "\t";
	}
	File.append(buffer, dataPath);


    
   
    // close file
    close(dv_array[i]);

	// Close plot
	selectWindow(getTitle());
	close();
	
    // close files 
    //selectWindow(dv_array[i]);
    //close();
    
}

exit("end of script");

//////////////////////////
///    END CODE        ///
//////////////////////////
