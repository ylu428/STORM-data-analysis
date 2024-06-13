# STORM-data-analysis


### Step 1: Extract the necessary information for data analysis from the raw data.
- The program will read raw STORM data from NanoImager (Unprocessed tif image files and SMLs.csv).
- Utilize the user interface to select and read **tif files** for all fields within a single sample.
- A folder containing "SMLs.csv" and single-frame images for each field will be created/saved.
- For each field of view (FOV), three single-frame images will be extracted from the raw image. They are:
  - first_frame_GFP=3: Taken from frame #3, used for locating GFPs and NPs in subsequent steps.
  - first_frame_NP=13: Taken from frame #13, used for locating NPs in subsequent steps.
  - start_frame_NP=1013: Taken from frame #1013, used for NP-dependent drift correction.
- Cut the entire created folder and paste it into a personal folder.
- DO NOT need the functions in our libraries.
- Original filename: **dSTORM_ONI_lowres_imageoutput_Yihan.m**

### Step 2: Examine the labeling efficiency from the DeltaVision images.
- The program will read raw image data from DeltaVision. Skip this step if no DV data was taken.
- Utilize the user interface to select and read **dv files** for all fields within a single sample.
- Automatically find all the particles by analyzing the intensities of GFP and AF647.
- Plot the cumulative distribution of GFP and AF647 intensities.
- Save images with particles found and cumulative distribution in the folder "output"
- Save particle information, including coordinates and intensities for each channel, into an Excel file.
- Need libraries in "*YCC_Matlab*"
- Original filename: **A_GFP_Vpr_IF_2dyes_YiHan.m**

### Step 3: Locate the nanoparticles in each field and identify the single molecule localizations originating from the nanoparticles.
- The program will read files created in step 1 ("_NP.tif" and "_SMLs.csv").
- Utilize the user interface to select and read all **NP_f1.tif** files for a single sample.
- Set the **start frame** and **end frame** correctly based on the light program you used in imaging acquisition step.
- Mapreduce process takes a long time. It could be 1-10 hours or more depends on number of fields and how many NPs per field.
- Two folders will be created: "MapReduceFiles" and "output".
  - MapReduceFiles: Results of mapreduce are saved here.
  - output:
    1. "*_SMLs_NP_Particles.mat": MATLAB data file. Info of all nanoparticles.
    2. "*_SMLsNP_part-detection.png": Detected NPs.
![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/ed82d10c-31f1-4db7-b3dc-6b6fa3dc73fe)
    3. "*_SMLsNP_part-filtering.png": Cumulative distribution function of particles.  
![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/a76db8e5-d70a-434a-acf9-255d0e4b2d99)
      - Upper left: corrected sum of signal;
      - upper right: signal-to-noise ratio;
      - lower left: resolving power;
      - lower right: percentage of background pixel.
    4. **image_output**: a folder where NP localization distribution plots are saved.
- Need libraries in "*YCC_Matlab*".
- Original filename: **one_dSTORM_ONI_parse_NP.m**

### Step 4: Drift correction processing.
- The program can automatically read all of the **NP images** saved in the "**image_output**" folder created in Step 3.
- Set the **start frame** and **end frame** correctly based on the light program you used in imaging acquisition step.
- Select the **"image_output"** folder.
- All of the NP scatter plots saved in the **"image_output"** folder will appear sequentially for user selection. (Scale bar: 100 nm)
![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/e5bb4ef1-6708-4693-8d16-258fdb8ace83)
  - Select "Yes" to choose this NP.
  - Select "No" to skip this NP.
  - Select "Done" to finish the NP selection. Any remaining NPs will be skipped.
- After completing the NP selection or terminating by clicking 'Done,' a review window will appear:
![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/47908e07-81fb-4b4f-89ac-8a7651eab915)
  - Click "X" to remove a specific NP from the selection.
  - Select "Skip Review" or "Done" to complete selection.
- Another open file dialog will appear. Select the folder where the "\*_NP_f1.tif" and "\*_SMLs.csv" files are saved.
- Drift correction is starting. The program will automatically process all fields that contain more than two valid NPs.
  - "\*_SMLs_NP_drift.png", "\*_NP_hist_corr.png", "\*NP_hist_ori.png", "\*_NP_merge_corr.png", "\*_NP_merge_ori.png" will be created and saved for each valid NP.
  - After drift correction complete, "Selected_NP.csv" will be saved in the "image_output" folder.
  - "\*_SMLs_driftcorr.mat" files will be saved in the "output" folder for fields that has been drift-corrected.
- Need libraries in "*YCC_Matlab*", "*\Yi-Han_edited_code\MATLAB_library*" and "*\Yi-Han_edited_code\MATLAB_library\HumphreyLibrary*".
- Original filename: **two_dSTORM_ONI_bead_driftcorr_Step4_v2.m**.

### Step 5: Locate the virion in each field and identify the single molecule localizations originating from the immunofluorescence.
- The program will read files created in step 1 ("_GFP.tif" and "_SMLs.csv").
- Utilize the user interface to select and read all **GFP_f1.tif** files for a single sample.
- Set the **start frame** and **end frame** correctly based on the light program you used in imaging acquisition step.
- Mapreduce process takes a long time. It could be 1-10 hours or more depends on number of fields, number of channels and how many virions per field.
- Results of mapreduce are saved in the "MapReduceFiles" folder created in step 3.
- The following MATLAB data will be saved in the "output" folder created in step 3.
  - MapReduceFiles: Results of mapreduce are saved here.
  -  output:
     1. "*_SMLs_GFP_Particles.mat": MATLAB data file. Info of all AF647 SMLs. (1-color STORM)
        "*_SMLs_GFP_Particles_0.mat" and "*_SMLs_GFP_Particles_1.mat": MATLAB data file for CF/AF568 and AF647, respectively. (2-color STORM)
     2. "*_SMLsGFP_part-detection.png": Detected NPs.
     3. "*_SMLsGFP_part-filtering.png"
     4. "*_SMLs_GFP_Particles_GFPSMLs.mat" MATLAB data file for all GFP signals.
- Need libraries in "*YCC_Matlab*" and "*\Yi-Han_edited_code\MATLAB_library*".
- Original filename: **three_dSTORM_ONI_parse_GFP_nm_2color_SMLs.m**. (this script works for both 1-color and 2-color STORM) or **three_dSTORM_ONI_parse_GFP_nm_SMLs.m** (only works for 1-color STORM)

### Step 6: DBSCAN for SMLs
- The program will read files created in step 1 ("_GFP.tif" and "_SMLs.csv") and step 5 (["*_SMLs_GFP_Particles.mat"] or ["*_SMLs_GFP_Particles_0.mat" + "*_SMLs_GFP_Particles_1.mat"])
- Utilize the user interface to select and read all **GFP_f1.tif** files for a single sample.
- Support doing DBSCAN for 1-color or 2-color STORM. Set variable "two_color" to "true" for 2-color, "false" for 1-color.
- Support multi-sample analysis by adjusting the number of files in "Sample", "stand_SML_n" and "Curr_SML_n"
![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/486d877c-6670-40d8-af0e-67c50daff03a)
- Adjust the number in "stand_SML_n" and "Curr_SML_n" to perform sample size adjustment.
- Adjust or add numbers to "minSMLs_1" to do DBSCAN for clusters with different minimum size.
- After finishing the analysis, a folder "analysis" will be created to save DBSCAN results, clustering images, pairwise info and other information.
- Clustering imaging: (Scale bar: 100 nm)
![SVA_58ABC_58C_SMLs-part-122 DB_15nm_20SMLs_GFP_2colorSMLs](https://github.com/ylu428/STORM-data-analysis/assets/41119470/70922ed3-21c3-44af-b37d-83daf5135aa6)
- The following MATLAB files will be saved in "\analysis\output\clusinfo" folder.
![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/4de368c0-1db6-4d91-b80c-6110dd480eed)
- An excel file with SMLs and cluster's information will be saved in the "process" folder.
- Need libraries in "*YCC_Matlab*" and "*\Yi-Han_edited_code\MATLAB_library*".
- Original filename: **dSTORM_ONI_DB_removeGFPclosetoNP_20240606.m**

### Step 7: Pairwise analysis
- The program will read files created in step 6 ("*_DB15nm_20SMLs_pdist_cf568_cf568.mat" and "*_DB15nm_20SMLs_pdist_af647_af647.mat").
- Run the script, select "*_DB15nm_20SMLs_pdist_cf568_cf568.mat" first then "*_DB15nm_20SMLs_pdist_af647_af647.mat" for 2-color STORM.
- After the figure presented, modify the legend/title/axis labels by double click on them. Zoom in or out using mouse scroll wheel.
![image](https://github.com/ylu428/STORM-data-analysis/assets/41119470/13a6fd06-0bee-4178-918b-6da8fe4ce19c)
- Save the plot after finishing label modification.
- - Need libraries in "*YCC_Matlab*"
- Original filename: **YiHan_pw_pdf_plot_ONI_EnvEnv_pc_5ABCD_2G12_AF647.m**
