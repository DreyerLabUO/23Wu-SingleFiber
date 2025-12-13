//-------------------------------------------------------------
//  23_Wu Analysis Macro
//  - Processes .lif files with Bio-Formats
//  - Saves background-subtracted TIF stacks
//  - Generates SD projections and skeletons
//  - Output organized as: [SubjectID]/[Timepoint]/[Leg]/[TIFs|STDIP|Skel]
//-------------------------------------------------------------

inputRoot  = getDirectory("Choose folder with 23_Wu .lif files");
if (inputRoot == "") exit("No input folder chosen.");

outputRoot = getDirectory("Choose root folder for 23_Wu results");
if (outputRoot == "") exit("No output folder chosen.");

setBatchMode(true);
run("Bio-Formats Macro Extensions");

lifFiles = getFileList(inputRoot);

for (f = 0; f < lifFiles.length; f++) {
    if (!endsWith(lifFiles[f], ".lif")) continue;

    lifPath = inputRoot + lifFiles[f];
    print("Processing " + lifFiles[f]);

    // Parse 23_Wu filename: 23_Wu_[SubjectID]_[Timepoint]_[Leg].lif
    baseName = replace(lifFiles[f], ".lif", "");
    parts = split(baseName, "_");
    if (parts.length < 5) {
        print("   Skipped (filename not in expected format): " + lifFiles[f]);
        continue;
    }

    subjectRaw = parts[2];
    tpToken    = parts[3];
    legToken   = parts[4];

    // Normalize tokens
    subject = subjectRaw;
    tpLower = toLowerCase(tpToken);
    if (tpLower == "bl") timePoint = "BL";
    else if (tpLower == "d14" || tpLower == "14d") timePoint = "D14";
    else if (tpLower == "w8") timePoint = "W8";
    else timePoint = tpToken;  // fallback to raw

    leg = legToken;

    // Build output directories
    subjDir = outputRoot + subject;
    tpDir   = subjDir + File.separator + timePoint;
    legDir  = tpDir   + File.separator + leg;
    tifDir  = legDir  + File.separator + "TIFs";
    stdDir  = legDir  + File.separator + "STDIP";
    skelDir = legDir  + File.separator + "Skel";

    File.makeDirectory(subjDir);
    File.makeDirectory(tpDir);
    File.makeDirectory(legDir);
    File.makeDirectory(tifDir);
    File.makeDirectory(stdDir);
    File.makeDirectory(skelDir);

    // Open the LIF with Bio-Formats
    Ext.setId(lifPath);
    Ext.getSeriesCount(nSeries);

    for (s = 0; s < nSeries; s++) {
        run("Bio-Formats Importer", 
            "open=[" + lifPath + "] " +
            "autoscale " +
            "color_mode=Default " +
            "view=Hyperstack " +
            "stack_order=XYCZT " +
            "series_" + s);

        title = getTitle();
        if (indexOf(title, "Merged") == -1) {
            close();
            continue;
        }

        // --- Background subtraction ---
        run("Subtract Background...", "rolling=30 stack");

        // Save background-subtracted stack as TIF
        mergedTitle = replace(title, " ", "_");
        saveAs("Tiff", tifDir + File.separator + mergedTitle + ".tif");

        // --- Standard Deviation Projection ---
        run("Z Project...", "projection=[Standard Deviation]");
        run("Median...", "radius=1");
        setAutoThreshold("Otsu dark no-reset");
        run("Threshold...");
        setOption("BlackBackground", true);
        run("Convert to Mask");

        projTitle = "STDIP_" + replace(title," ","_");
        rename(projTitle);
        saveAs("Tiff", stdDir + File.separator + projTitle + ".tif");

        // --- Skeleton workflow ---
        run("Duplicate...", "title=SkelImage");
        selectWindow("SkelImage");
        run("Gaussian Blur...", "sigma=55 scaled");
        setAutoThreshold("Otsu dark");
        setOption("BlackBackground", true);
        run("Convert to Mask");
        run("8-bit");
        run("Make Binary");
        run("Open");

        run("Analyze Particles...", "size=70000-Infinity clear include add");
        roiCnt = roiManager("count");

        if (roiCnt > 0) {
            newImage("TempMask", "8-bit black", getWidth(), getHeight(), 1);
            selectWindow("TempMask");
            roiManager("Select", 0);
            setForegroundColor(255, 255, 255);
            roiManager("Fill");

            run("Convert to Mask");
            run("8-bit");

            skelTitle = "Skel_" + replace(title," ","_");
            rename(skelTitle);
            run("Skeletonize");
            saveAs("Tiff", skelDir + File.separator + skelTitle + ".tif");
            close(skelTitle);
        } else {
            print("   No large ROI in " + title + " (skipped)");
        }
        roiManager("Reset");
        roiManager("deselect");

        // Cleanup
        close("SkelImage");
        close(projTitle);
        close(title);
        run("Close All");
    }

    Ext.close();
    call("java.lang.System.gc");
    print("   Finished " + lifFiles[f]);
}

setBatchMode(false);
print("All 23_Wu projects processed. Results in: " + outputRoot);