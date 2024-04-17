# STORM-data-analysis


### Step 1: 
- Read raw STORM data from NanoImager (Unprocessed tif image files and SMLs.csv).
- Utilize the user interface to select and read TIF files for all fields within a single sample.
- A folder containing "SMLs.csv" and single-frame images for each field will be created/saved.
- For each field of view (FOV), three single-frame images will be extracted from the raw image. They are:
  - first_frame_GFP=3: Taken from frame #3, used for locating GFPs and NPs in subsequent steps.
  - first_frame_NP=13: Taken from frame #13, used for locating NPs in subsequent steps.
  - start_frame_NP=1013: Taken from frame #1013, used for NP-dependent drift correction.
- DO NOT need the functions in our library.
- Original filename: dSTORM_ONI_lowres_imageoutput_Yihan.m

### Step 2: 
- Read raw image data from DeltaVision.
- Utilize the user interface to select and read dv files for all fields within a single sample.
- Automatically find all the particles by analyzing intensities of GFP and AF647.
- Plot the cumulative distribution of GFP and AF647 intensities.
- Save images with particles found and cumulative distribution in the folder "output"
- Save particle information, including coordinates and intensities for each channel, into an Excel file.
- Original filename: A_GFP_Vpr_IF_2dyes_YiHan.
