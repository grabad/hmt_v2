# Histone-Modification-Toolbox
Contains all processing code from the paper "A Toolbox to Quantify the Spatial Distribution of Post-translational Histone Modifications from 3D Spectrscopic SMLM Images". This code is intended for use on multicolor images of nuclear proteins such as histone modifications or RNA polymerase
1.	Process all 3D multicolor SMLM datasets and save the localization lists as CSV files with the x-locations in the 3rd column, the y-locations in the 4th column, and the z-location in the 12th column.
2.	Split the results into the two color channels (1 and 2) with the channel 1 filename appended as “_C1” and the channel 2 filename appended as “_C2”.
3.	Find the largest file and place it and its “_C1” and “_C2” files in the “Initialization” folder.
4.	Open the MATLAB file “denstest_1.m” and In lines 8-9, change the “_C1” and “_C2” filenames to match the CSV files you placed in the initialization folder. Run the code.
5.	Open the file “simulation_optimization¬_2.m” and change line one to match the filename of original file you placed in the Initialization folder (i.e. not the “_C1” or “_C2” appended files). Run the code
6.	Open and run the code “sim_densities_3.m”
7.	Open and run “split_script_4.m”
8.	Open and run “param_opt2_5.m”
9.	Transfer the output files “methyl_epsilon_surround2” and “acetyl_epsilon_surround2” to the “Analysis” folder
10.	Place all 2 color images and “_C1” and “_C2” appended files for a single experimental condition in the analysis folder.
11.	Open “analyze_STORM.m”
12.	 In line 3 and line 6 change the list of names (saved to the variable “simnames”) to the names of the csv files.
13.	Run the code. View the resulting figures and ensure that the nuclear binarization looks correct.
14.	If the binarized nucleus includes too much signal, lower the threshold in line 155 (set to 0.17). If it includes too little, raise it.
15.	If no objects appear, lower the threshhold set in line 158 (set to 500000). If too many objects appear, raise it.
16.	Open the file “compile_data.m” and change each “load” and “save” filename to match the file names of the localization csv files. Adjust the number of variables in lines 128-170 accordingly.

