
// SELECT YOUR PARAMETERS. Please change this part to your input and output directories and desired contact sheet file name


inputDir="/Users/ambikabasu/Nishimura_lab/Macros/260527_N2_imb2_set3";
// note - this directory must exist already and contain .dv files 
outputDir="/Users/ambikabasu/Nishimura_lab/Macros/260527_N2_imb2_set3/output";
// You can specify an output directory, and it will create a new directory or add to an existing one

// Defaults: Please modify these as desired:
// * Default is to look for R3D.dv files. To modify, edit line #32
// * Default z-projection range is 1 - 20. To modify, edit line 70
// * Default brightness max and min are specified between lines 78 & 90. Please run the program once, select your favorite settings, and modify those lines appropriately
// * Default way to merge channels into specific colors is in line 91. Please modify as desired
// * Default exports are in lines 102 - 110

//////////////////////////////////////////////////////////

// Create the output directory
File.makeDirectory(outputDir);

// set the working directory
File.setDefaultDir(inputDir);

// print the input directory
print(inputDir);

// get the list of files in this directory
fileList = getFileList(inputDir);

// restrict list to R3D_D3D.dv files
dv_array = newArray(0);

// Save .dv files in a dv_array:
for (i=0; i<fileList.length; i++) {
    if(endsWith(fileList[i], "R3D.dv"))
    {
    	dv_array = Array.concat(dv_array, fileList[i]);
    }
}

//print the files;
print("\nWill process the following DV files:");
for (i=0; i<dv_array.length; i++) {
    print("\t" + dv_array[i]);
}


// cycle through the images in the directory
for (i=0; i<dv_array.length; i++) {
    full_dv = inputDir + "/" + dv_array[i];
    
    // Report 
    print("\nNow processing file:\t" + dv_array[i]);
    
    open(dv_array[i]);
    
    //select the first image
    selectWindow(dv_array[i]);
    
    // set the zoom
    run("Set... ", "zoom=75 x=512 y=512");
    
    // adjust the contrast
	run("Enhance Contrast", "saturated=0.35");
	run("Enhance Contrast", "saturated=0.35");
	run("Enhance Contrast", "saturated=0.35");
	
	// z project
	run("Z Project...", "start=1 stop=20 projection=[Max Intensity]");
	
	// split channels
    run("Split Channels");
    
    // get channel names
    c1title = "C1-MAX_" + dv_array[i];
    c2title = "C2-MAX_" + dv_array[i];
    c3title = "C3-MAX_" + dv_array[i];
    //c4title = "C4-MAX_" + dv_array[i];
    
    // set the brightness - c1
    selectWindow(c1title);
    setMinAndMax(385, 5524); 
    
    // set the brightness - c2
    selectWindow(c2title);
    setMinAndMax(385, 5524);   
    
    // set the brightness - c3
    selectWindow(c1title);
    setMinAndMax(2583, 9935);       

    // set the brightness - c4
    selectWindow(c4title);
    setMinAndMax(2583, 9935);  
    
    run("Merge Channels...", "c2=[" + c2title + "] c3=[" + c4title + "] c6=[" + c1title + "] create");
    
    // set the zoom
    run("Set... ", "zoom=75 x=512 y=512");
    
    // save image
    rootname = File.getNameWithoutExtension(dv_array[i]);
    print(rootname);
    
    
    // Save the .tif output
    tiffoutput = outputDir + "/" + rootname + ".tif";
    print(tiffoutput);
    saveAs("Tiff", tiffoutput);
    
    // Save the .jpg output
    jpgoutput = outputDir + "/" + rootname + ".jpg";
    print(jpgoutput);
    saveAs("jpeg", jpgoutput);    
   
    // close file
    close(dv_array[i]);
    
    //saveAs("Tiff", "MAX_image_7_R3D.tif");
	//saveAs("Jpeg", "/Volumes/NESS/260205_Ambika_Stress/251101_wLPW019_O2stress/stress_test/MAX_image_7_R3D.jpg");
    
}

// run the script makeContactSheet.sh to create the total contact sheet - this doesn't work yet
// Usage: bash makeContactSheet.sh <contactsheetname.jpg>
//File.setDefaultDir(outputDir);
//exec("/Volumes/NESS/260206_macroBuilding/makeContactSheet.sh", contactsheetname);
