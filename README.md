# popgen.analysis.pipeline
 
/aDNA_analysis contains all the scripts used for analysing ancient DNA. There are separate workflows for logistic regression and MARS. The workflows will fit the models, generate the partial dependence plots and compute the FIRM scores for the models. 

/data contains the dataframe produced by computing the summary statistics for the data in chapter 4. 

/Data_cleaning contains 2 scripts to clean the modern data (ch 5) and aDNA data (ch 6) respectively. 

/DNA_analysis contains the workflow for analysing modern data (ch 5).

/Exploratory_DA contains scripts for doing parallel coordinates plots, PCA and PLS. 

/Plots_raw_stats has scripts for parallel coordinates plots and PCA.

/simulation contains the scripts to simulate the modern and ancient data for ch5 and 6 respectively. 

/vignette contains short tutorials for the project.

/HPC contains job scripts to compute our dataframes on the cluster at 
Adelaide University. 