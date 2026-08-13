
// SELECT YOUR PARAMETERS. Please change this part to your input and output directories and desired contact sheet file name


inputDir="/Users/ambikabasu/Nishimura_lab/Macros/250806_set3_imb2_rep1";
// note - this directory must exist already and contain .dv files 
outputDir="/Users/ambikabasu/Nishimura_lab/Macros/250806_set3_imb2_rep1/output";
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
	run("Z Project...", "start=10 stop=30 projection=[Max Intensity]");
	
	// split channels
    run("Split Channels");
    
    // get channel names
    c1title = "C1-MAX_" + dv_array[i];
    c2title = "C2-MAX_" + dv_array[i];
    c3title = "C3-MAX_" + dv_array[i];
    //c4title = "C4-MAX_" + dv_array[i];
    
//    // set the brightness - c1 red
//    selectWindow(c1title);
//    setMinAndMax(711, 1158); 
    
    // set the brightness - c2 magenta
    selectWindow(c1title);
    setMinAndMax(600, 2500);   
    
    // set the brightness - c3 green
    //selectWindow(c2title);
    //setMinAndMax(450, 3700);       

    // set the brightness - c4 blue
    selectWindow(c3title);
    setMinAndMax(200, 412);  
    
    //run("Merge Channels...", "c6=[" + c1title + "] c2=[" + c2title + "] c3=[" + c3title + "] create");
    run("Merge Channels...", "c6=[" + c1title + "] c3=[" + c3title + "] create");
    
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
    
}

// Generate <foldername>_README.txt in the output folder, listing all images produced above
makeReadme(outputDir);

// run the script makeContactSheet.sh to create the total contact sheet - this doesn't work yet
// Usage: bash makeContactSheet.sh <contactsheetname.jpg>
//File.setDefaultDir(outputDir);
//exec("/Volumes/NESS/260206_macroBuilding/makeContactSheet.sh", contactsheetname);

// =====================================================================
// README scaffolding
// Writes <date>_<foldername>_README.txt into the given folder.
//   * Date  : the date the macro is run (yymmdd), via getDateAndTime().
//   * Folder PATH : the output folder path.
//   * Annotation  : one '<image_name>: ' line per image in the folder,
//             de-duplicated by basename so a paired .tif/.jpg saved from
//             the same source appears only once.
// The whole file is built as one string and written with a single
// File.saveString() call, so a half-written/empty file can never result.
// To keep file extensions on annotation lines: change `stripExtension(name)`
// to `name` in makeReadme(). To change date format: edit `currentDate()`.
// =====================================================================

  // True if 'name' ends with a recognized image extension.
  // The extension list is local (self-contained) on purpose - no globals.
  function isImage(name) {
      exts = newArray(".tif", ".tiff", ".jpg", ".jpeg", ".png");
      lower = toLowerCase(name);
      for (e = 0; e < exts.length; e++) {
          if (endsWith(lower, exts[e])) return true;
      }
      return false;
  }

  function stripExtension(name) {
      idx = lastIndexOf(name, ".");
      if (idx < 0) return name;
      return substring(name, 0, idx);
  }

  // Current date (the day the macro is run) as a yymmdd string, e.g. "260519".
  // Uses the built-in IJ.pad for zero-padding, and concatenates string
  // variables (not chained function calls) so the macro interpreter always
  // treats the result as a string. To change format, edit the return line:
  //   yy-mm-dd:    return "" + yy + "-" + mm + "-" + dd;
  //   20yy-mm-dd:  return "20" + yy + "-" + mm + "-" + dd;
  function currentDate() {
      getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
      // getDateAndTime returns month as 0-11, so add 1.
      yy = IJ.pad(year % 100, 2);
      mm = IJ.pad(month + 1, 2);
      dd = IJ.pad(dayOfMonth, 2);
      return "" + yy + mm + dd;
  }

  function makeReadme(folderPath) {
      // Ensure trailing separator for path operations
      if (!endsWith(folderPath, File.separator))
          folderPath = folderPath + File.separator;

      trimmed = substring(folderPath, 0, lengthOf(folderPath) - 1);
      folderName = File.getName(trimmed);

      // Date: the day the macro is run (yymmdd)
      date = currentDate();

      // README filename: <date>_<foldername>_README.txt, but avoid doubling
      // the date if the output folder is already prefixed with today's date.
      if (!startsWith(folderName, date))
          readmeName = date + "_" + folderName + "_README.txt";
      else
          readmeName = folderName + "_README.txt";

      readmePath = folderPath + readmeName;
      if (File.exists(readmePath)) {
          print("README already exists, skipping: " + readmePath);
          return 0;
      }

      // Gather image files, sorted
      files = getFileList(folderPath);
      Array.sort(files);

      // Build the entire README as one string.
      content = "Date: " + date + "\n";
      content = content + "\n";
      content = content + "Name: Ambika Basu \n"; // Change this to user name 
      content = content + "\n";
      content = content + "Folder PATH: " + trimmed + "\n";
      content = content + "\n";
      content = content + "Project:  Ambika imb-2 vs set-3 quant \n"; // Change this to project title
      content = content + "\n";
      content = content + "Purpose: Make contact sheet so we can quantify signal across nuclear membrane  \n"; // Change purpose here
      content = content + "\n";
      content = content + "Notes: \n"; // can edit this section to autopopulate notes in readme
      content = content + "\n";
      content = content + "C1 = imb-2 <- cy5 # Not used \n";
      content = content + "C2 = set-3 <- mCherry-Magenta \n";
      content = content + "C3 =  DAPI <- DAPI-blue \n";
      //content = content + "C4 = DNA <- DAPI-Blue \n";
      content = content + "\n";
      content = content + "Annotation: \n";
      content = content + "\n";

      // One '<image_name>: ' line per image, de-duplicated by basename so
      // paired .tif/.jpg outputs only generate one annotation entry.
      seen = newArray(0);
      imgCount = 0;
      for (k = 0; k < files.length; k++) {
          name = files[k];
          if (endsWith(name, "/")) continue;       // subdirectory
          if (!isImage(name)) continue;
          base = stripExtension(name);

          dup = false;
          for (j = 0; j < seen.length; j++) {
              if (seen[j] == base) { dup = true; break; }
          }
          if (dup) continue;
          seen = Array.concat(seen, base);

          content = content + base + ": \n";
          imgCount++;
      }

      // Single atomic write - no file handle, no buffer to lose.
      File.saveString(content, readmePath);
      print("Wrote README (" + imgCount + " images listed): " + readmePath);
      return 1;
  }
