# In situ seq
This repository contains scripts for in situ sequencing processing pipeline.
Mats Nilsson Lab, Stockholm University
Xiaoyan Qian, 2018

Download the repository, add lib to MATLAB path. 
Except MATLAB, no additional Mathworks product is required.
Tested on R2016b and R2018a.

InSituSequencing.m is the top-level script.

Example images from the original Nature Methods paper can be found here https://figshare.com/articles/Slide_A_stitched_images/7246952

## This folder contains:
- **InSituSequencing.m**: script for performing either **TILING** previously to Blob detection or **GENE DECODING** and visualizing the output 
- **lib**: subset _of functions_ commonly used in the main repository scripts
- **InputExample/NatMeth**: example of blobs.csv files and taglist file. 
- **LICENSE**: license for use. Not needed for running anything
- **.gitignore**: default git-associated file
