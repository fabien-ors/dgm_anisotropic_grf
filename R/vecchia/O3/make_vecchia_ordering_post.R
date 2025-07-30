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

## Ordering
orderingFile=paste0(path_files,id_file,"obs_ordering.rds")
if(file.exists(orderingFile)){
  obsOrdering=readRDS(orderingFile)
}else{
  stop("Make mcmc ordering first")
}

res=createVecchiaModel(datafilename,path_files,id_file,id_sim,id_mask,kVecchia,tau2obs,subsamp,
                       sgvOrdering=obsOrdering,onlyOrd = FALSE,constInit = TRUE)
Rmodel=res$model
Yobs=res$data
obsCoord=res$obsCoord
nodeMat=res$allCoord
sel_mask=res$indObs
sel_mask=res$indObs
z_sim=res$image
param_sim=res$param_sim
mask=res$mask
knot_blocks=res$knot_blocks
xknots=res$xknots
ordObs=res$ord
mObs=res$mObs
sdObs=res$sdObs
rm(res)

#----------------------------------------------------------------

## Prediction locations
sel_pred=setdiff(1:nrow(nodeMat),sel_mask)
Xpred=nodeMat[sel_pred,]

#----------------------------------------------------------------

orderingFile=paste0(path_files,id_file,"pred_ordering.rds")
if(file.exists(orderingFile)){
  sgvOrdering=readRDS(orderingFile)
}else{
  message("\nOrdering the prediction locations and determining neighbors/conditioning sets for SGV (this may take a minute).\n")
  sgvOrdering <- sgvSetup(coords = obsCoord[ordObs,], coords_pred = Xpred, k = kVecchia, 
                          seed=sample(1e5,1),order_coords = FALSE)
  saveRDS(sgvOrdering,orderingFile)
}



