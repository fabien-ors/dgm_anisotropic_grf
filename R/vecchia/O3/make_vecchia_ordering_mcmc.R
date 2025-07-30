############################################################################
# RUN MCMC WITH VECCHIA MODEL
############################################################################

library(parallel)

#----------------------------------------------------------------

## Load function to create Vecchia model
source("setupVecchiaTemp.R")

#----------------------------------------------------------------

# Define MCMC parameters
tau2obs=0.05**2
path_files="./test_folder/"

timed=FALSE

#--------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(sum(is.na(args[1:3]))>0){
  stop("Three numerical argumuments are expected for the first three arguements of the R script")
}
id_sim=as.numeric(args[1]) #  Id of the full simulation
id_mask=as.numeric(args[2]) # Id of the mask
kVecchia=as.numeric(args[3])
subsamp=NA
datafilename=args[4]

#--------------------------------------------------------------
## File paths

## Id file and export folder
id_file=paste0("run_",datafilename,"_",id_sim,id_mask,"_v",kVecchia)
export_folder=paste0(path_files,id_file,"/")
if(!dir.exists(export_folder)){
  dir.create(export_folder)
}
id_file=paste0(id_file,"_")

#--------------------------------------------------------------

cat("MCMC data: Data Id=",id_sim,"/ Mask Id=",id_mask,"/ Vecchia Approx=",kVecchia,"/ Subsamp=",subsamp,"\n")

#--------------------------------------------------------------

sgvOrdering=createVecchiaModel(datafilename,path_files,id_file,id_sim,id_mask,kVecchia,tau2obs,subsamp,onlyOrd=TRUE)

#----------------------------------------------------------------
