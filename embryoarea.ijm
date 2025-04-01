function embryoarea(input, filename) {
    // Open the image
    open(input + "/" + filename);
    wait(500); // Wait for the file to open

    run("RGB Color");
    
    // Set black background option
    setOption("BlackBackground", true);
    
    // Convert to Mask
    run("Convert to Mask");
    
    // Fill Holes
    run("Fill Holes");
    
    // Perform watershed segmentation
    run("Watershed", "stack");
    
    // Measure the regions of interest (ROIs)
    run("Set Measurements...", "area mean standard fit shape redirect=None decimal=3");
    run("Analyze Particles...", "size=10000-Infinity circularity=0.00-1.00 show=Masks display clear add");
    
    // Add the filename to the Results table for each measurement
    var currentResults = nResults;
    var startRow = currentResults - nResults;

    for (var j = startRow; j < currentResults; j++) {
        setResult("Image", j, filename); // Add the filename to the new column for each result
    }

    // Save the individual results for this image
    saveAs("Results", input + "/Results_" + filename.replace(".tif", ".csv"));
    
    // Clear the Results table to prepare for the next file
    run("Clear Results");

    // Close all open image windows
    while (nImages > 0) {
        selectImage(nImages);
        close();
    }
}

// Define the input directory
input = "/Volumes/eels/Crab project /Embryo_data/Brood_3/Low/cage 8/2. 1-31-25/C";

// Clear the Results table before starting
run("Clear Results");

// Loop through the files G1.tif to G12.tif
for (var i = 1; i <= 12; i++) {
    filename = "A" + i + ".tif";
    if (File.exists(input + "/" + filename)) {
        print("Processing " + filename);
        embryoarea(input, filename);
    } else {
        print("File not found: " + filename);
    }
}

// Combine all individual CSV files into one final CSV file
//combinedResultsFile = input + "/Combined_Results.csv";
//run("Clear Results");

// Read and combine each individual CSV file
//for (var i = 1; i <= 12; i++) {
   // filename = "Results_G" + i + ".csv";
    //if (File.exists(input + "/" + filename)) {
        //open(input + "/" + filename);
        //wait(500); // Wait for the file to open
        
        // Append the results to the combined file
        //if (i == 1) {
            // For the first file, include the header
         //   run("Table... ", "save=[" + combinedResultsFile + "]");
        //} else {
            // For subsequent files, append without the header
        //    run("Table... ", "append=[" + combinedResultsFile + "]");
        //}
        
        // Close the opened CSV
      //  close();
   // } else {
    //    print("File not found: " + filename);
  //  }
//}

//print("All individual CSV files have been combined into " + combinedResultsFile);


