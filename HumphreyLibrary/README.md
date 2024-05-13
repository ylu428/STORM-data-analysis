# README

## Overview
This repository contains a set of MATLAB scripts for processing and reviewing images of nanoparticles. The scripts allow users to select images from a directory, review and filter the images, and save the selected images' details to a CSV file. The main functions included are:

1. `processImagesAndSaveToCSV`: Handles the initial selection of images and saves the results to a CSV file.
2. `reviewSelectedImages`: Provides a GUI for users to review the selected images.
3. `processCSV`: Reads a CSV file and matches its data with corresponding TIF and CSV files in a directory.
4. `updateCSV`: Updates a CSV file with the corresponding images and nano  particiles extracted from the selection process
5. `parseAndFilterFilename`: Parses filenames to extract test IDs and particle numbers.
6. `displayImageAndPrompt`: Displays an image and prompts the user for input.
7. `convertSelectionsToMap`: Converts a list of selected images into a selections map.

## Getting Started
To use these scripts, ensure you have MATLAB installed on your system. Follow the steps below to process and review your NP images.

### Step 1: Process Images and Save to CSV
The `processImagesAndSaveToCSV.m` script is the main scirpt starter and it used to process images in a selected directory and save the results to a CSV file. This function is what is called in Step4_V1.m

**Purpose**: This script allows users to select a folder, process images within it, and save the details of selected images to a CSV file.

**Usage**:
1. Follow the prompts to select a folder containing NP images.
2. The script will process each image, allow user input to save or skip images, and save the selected images' details to a CSV file.

### Review Selected Images
After the initial selection, use the `reviewSelectedImages.m` script to review the selected images.

**Purpose**: This script provides a user interface displaying all selected images in a grid layout. Users can remove specific images, skip the review, or finalize and save the reviewed selection.

**Usage**:
- The `processImagesAndSaveToCSV.m` script will call this function after the initial image selection.
- Users can interact with the GUI to review their selections.

###  Process the CSV File
The `processCSV.m` script reasd the CSV file and match its data with corresponding TIF and CSV files in a selected directory.

**Purpose**: This script reads test IDs and particle numbers from the CSV file, matches them with files in the directory, and outputs a list of matched file sets.

**Usage**:
1. Open MATLAB and run the `processCSV.m` script.
2. Follow the prompts to select the CSV file and the folder containing the corresponding TIF and CSV files.

## Helper Functions

### 1. `parseAndFilterFilename`
**Purpose**: This function parses filenames to extract test IDs and particle numbers and determines if the filename is relevant based on a specific pattern.

**Usage**: This function is called by other scripts to filter and extract information from filenames.

### 2. `displayImageAndPrompt`
**Purpose**: This function displays an image and prompts the user for input. It shows an image and provides Yes, No, and Done buttons for user interaction.

**Usage**: This function is called by `processImagesAndSaveToCSV.m` to display images and collect user input.

### 3. `convertSelectionsToMap`
**Purpose**: This function converts a list of selected images into a selections map format expected by the `updateCSV` function.

**Usage**: This function is called by `processImagesAndSaveToCSV.m` after the review process to prepare the final selections for saving to the CSV file.

## Example Workflow

1. **Run `processImagesAndSaveToCSV`** in Step4_v2:
    - This script prompts the user to select a folder containing NP images.
    - It processes each image, allows user input to save or skip images, and saves the selected images' details to a CSV file.

2. **Review Selected Images**:
    - After the initial selection, the `reviewSelectedImages` function provides a GUI to review the images.
    - Users can remove specific images, skip the review, or finalize the review.

3. **Process the CSV File**:
    - The `processCSV` function reads the CSV file and matches its data with corresponding TIF and CSV files in a selected directory.

