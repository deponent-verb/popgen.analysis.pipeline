#!/bin/bash
#SBATCH -p batch           	                                # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (no MPI, so we only use a single node)
#SBATCH -n 7              	                                # number of cores
#SBATCH --time=4:00:00    	                                # walltime allocation, which has the format (D-HH:MM:SS), here set to 1 hour
#SBATCH --mem=5GB         	                                # memory required per node (here set to 4 GB)

# Notification configuration 
#SBATCH --mail-type=END					    	# Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL   					# Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=a1708050@adelaide.edu.au  	# Email to which notifications will be sent

# Execute the program
# (The example here is a sequential bash script; use a suitable program for your case.)
R CMD BATCH cpop_aDNA.R
cat cpop_aDNA.Rout