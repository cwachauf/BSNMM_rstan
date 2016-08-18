args=commandArgs(trailingOnly=TRUE)
source("SNMM_functions.R")
SNMM_Stan_Start_Job_By_ID("snmm_job_specifications.csv",as.numeric(args[1]))