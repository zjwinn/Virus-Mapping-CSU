#!/bin/bash
#SBATCH --job-name=virus-gwas           # Name of the job
#SBATCH --partition=batch               # Run on the 'batch' partition
#SBATCH --nodes=1                       # Request one whole node
#SBATCH --ntasks=1                      # Request one single-process task
#SBATCH --cpus-per-task=40              # Use 40 threads/cores for this single task
#SBATCH --mem=250G                      # Request 250 Gigabytes of memory
#SBATCH --time=7-00:00:00               # Set the time limit to 7 days
#SBATCH --output=%x_%j.out              # %x=job_name, %j=job_id
#SBATCH --error=%x_%j.err               # Same thing except for errors
#SBATCH --mail-type=START,FAIL,END      # Says when the job starts, fails, or ends
#SBATCH --mail-user=zjwinn@colostate.edu

# Run R script
Rscript Analysis.R
