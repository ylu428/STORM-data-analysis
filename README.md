# STORM-data-analysis


### Step 1: 
- Read raw STORM data from NanoImager (Unprocessed tif image files and SMLs.csv).
- Utilize the user interface to select and read **tif files** for all fields within a single sample.
- A folder containing "SMLs.csv" and single-frame images for each field will be created/saved.
- For each field of view (FOV), three single-frame images will be extracted from the raw image. They are:
  - first_frame_GFP=3: Taken from frame #3, used for locating GFPs and NPs in subsequent steps.
  - first_frame_NP=13: Taken from frame #13, used for locating NPs in subsequent steps.
  - start_frame_NP=1013: Taken from frame #1013, used for NP-dependent drift correction.
- Cut the entire created folder and paste it into a personal folder.
- DO NOT need the functions in our libraries.
- Original filename: dSTORM_ONI_lowres_imageoutput_Yihan.m

### Step 2: 
- Read raw image data from DeltaVision. Skip this step if no DV data was taken.
- Utilize the user interface to select and read **dv files** for all fields within a single sample.
- Automatically find all the particles by analyzing the intensities of GFP and AF647.
- Plot the cumulative distribution of GFP and AF647 intensities.
- Save images with particles found and cumulative distribution in the folder "output"
- Save particle information, including coordinates and intensities for each channel, into an Excel file.
- Need libraries in "YCC_Matlab"
- Original filename: A_GFP_Vpr_IF_2dyes_YiHan.m

### Step 3:
- Read files created in step 1 ("_NP.tif" and "_SMLs.csv").
- Utilize the user interface to select and read all **NP_f1.tif** files for a single sample.
- Set the **start frame** and **end frame** correctly based on the light program you used in imaging acquisition step.
- Mapreduce process takes a long time. It could be 1-10 hours or more depends on number of fields and how many NPs per field.
- Two folders will be created: "MapReduceFiles" and "output".
  - MapReduceFiles: Results of mapreduce are saved here.
  - output:
    1. "*_SMLs_NP_Particles.mat": MATLAB data file. Info of all nanoparticles.
    2. "*_SMLsNP_part-detection.png": Detected NPs.
    3. "*_SMLsNP_part-filtering.png": Cumulative distribution function of particles.  Upper left: corrected sum of signal; upper right: signal-to-noise ratio; lower left: resolving power; lower right: percentage of background pixel.
    4. **image_output**: a folder where NP localization distribution plots are saved.
- Need libraries in "YCC_Matlab".
- Original filename: one_dSTORM_ONI_parse_NP.m

### Step 4:
- Drift correction processing. Automatically read all of the **NP images** saved in the "**image_output**" folder created in Step 3.
- Set the **start frame** and **end frame** correctly based on the light program you used in imaging acquisition step.
- Select the **"image_output"** folder.
- ![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/e5bb4ef1-6708-4693-8d16-258fdb8ace83)

- Need libraries in "YCC_Matlab" and "HumphreyLibrary"

