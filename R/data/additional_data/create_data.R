library(terra)
library(rhdf5)
library(minigst)
library(fields)
library(bestNormalize) ## For z-score transform
library(rnaturalearth) ## To plot
library(sf) ## To plot

setwd("~/Documents/Work/dev/dgm_anisotropic_grf/R/data/rev/")

## World map
world <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world)

setwd("~/Documents/Work/dev/dgm_anisotropic_grf/R/data/rev/")

## Size
Ngd=256

rast_to_mat<-function(r){
  t(matrix(values(r)[,1],nrow = nrow(r),byrow = T))[1:Ngd,Ngd:1]
}

# Load NetCDF as a raster (works well for LST, NDVI, etc.)
r <- terra::rotate(rast("./raw/era5_O3_cloud.nc")) ## u10, v10, Total Column Ozone, Cloud
rb <- terra::rotate(rast("./raw/era5_slhf.nc")) ## Surface Latent Heat Flux
## Dates: Jan 1st, 2025 at 12:00

## Id of variables
idu10=1
idv10=2
idO3=3 # In 
idCloud=4
idHF=1

# Resolution
df=res(r)*Ngd 

## Selected zones
# tab_x=rbind(c(-10,-50)-df/2,
#             c(80,-50)-df/2)

## Selected zones
tab_x=rbind(c(-10,-50)-df/2,
            c(80,-25)-df/2,
            c(-140,-50)-df/2)

## Plot
# png(filename = "map_O3.png",width=8,height = 5,res=600,units = "in")
plot((r[[idO3]]))
lines(world_vect, col = "white", lwd = 1)
rect(tab_x[,1],tab_x[,2],tab_x[,1]+df[1],tab_x[,2]+df[2],lwd=2)
text(tab_x[,1]+df[1]/2,tab_x[,2]+df[2]/2,-1+1:nrow(tab_x))
# dev.off()


# png(filename = "map_SLHF.png",width=8,height = 5,res=600,units = "in")
plot((rb[[idHF]]))
lines(world_vect, col = "white", lwd = 1)
rect(tab_x[,1],tab_x[,2],tab_x[,1]+df[1],tab_x[,2]+df[2],lwd=2)
text(tab_x[,1]+df[1]/2,tab_x[,2]+df[2]/2,-1+1:nrow(tab_x))
# dev.off()

#---------------------------------------


## Name of the export file
filename_radix="data_O3"
list_transfo=list()
filename=paste0(filename_radix,".h5")
if (file.exists(filename)) {
  file.remove(filename)
}
h5createFile(filename)
h5createDataset(filename, "data", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "data_norm", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask_sub", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
## Fill data
for(i in 1:nrow(tab_x)){
  xmin=tab_x[i,1]
  xmax=xmin+df[1]
  ymin=tab_x[i,2]
  ymax=ymin+df[2]
  
  window(r)=ext(c(xmin,xmax,ymin,ymax))
  
  # Plot the first layer
  #plot(r[[idO3]],main=paste0(xmin,",",ymax,",",xmax,",",ymin))
  dat=rast_to_mat(r[[idO3]])
  print(dim(dat))
  image.plot(dat)
  h5write(dat, file=filename, name="data", index=list(1:Ngd,1:Ngd,i))
  
  
  dat_cloud=rast_to_mat(r[[idCloud]])<0.5+0
  image.plot(dat_cloud)
  h5write(dat_cloud, file=filename, name="mask", index=list(1:Ngd,1:Ngd,i))
  
  
  set.seed(42+i)
  sel_mask=sample(which(dat_cloud==1),10**3)
  dat_cloud_bis=dat_cloud*0
  dat_cloud_bis[sel_mask]=1
  image.plot(dat_cloud_bis)
  h5write(dat_cloud_bis, file=filename, name="mask_sub", index=list(1:Ngd,1:Ngd,i))
  
  
  ## Transform data
  bn <- orderNorm(as.vector(dat[which(dat_cloud==1)]))
  y_trans <- matrix(predict(bn,dat),Ngd,Ngd)
  h5write(y_trans, file=filename, name="data_norm", index=list(1:Ngd,1:Ngd,i))
  list_transfo[[i]]=bn
}
saveRDS(list_transfo,paste0(filename_radix,"_normTransform.rds"))

#---------------------------

## Name of the export file
filename_radix="data_u10"
list_transfo=list()
filename=paste0(filename_radix,".h5")
if (file.exists(filename)) {
  file.remove(filename)
}
h5createFile(filename)
h5createDataset(filename, "data", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "data_norm", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask_sub", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)

## Fill data
for(i in 1:nrow(tab_x)){
  xmin=tab_x[i,1]
  xmax=xmin+df[1]
  ymin=tab_x[i,2]
  ymax=ymin+df[2]
  
  window(r)=ext(c(xmin,xmax,ymin,ymax))
  
  # Plot the first layer
  plot(r[[idu10]],main=paste0(xmin,",",ymax,",",xmax,",",ymin))
  dat=rast_to_mat(r[[idu10]])
  print(dim(dat))
  image.plot(dat)
  h5write(dat, file=filename, name="data", index=list(1:Ngd,1:Ngd,i))
  
  
  dat_cloud=rast_to_mat(r[[idCloud]])<0.5+0
  image.plot(dat_cloud)
  h5write(dat_cloud, file=filename, name="mask", index=list(1:Ngd,1:Ngd,i))
  
  
  set.seed(42+i)
  sel_mask=sample(which(dat_cloud==1),10**3)
  dat_cloud_bis=dat_cloud*0
  dat_cloud_bis[sel_mask]=1
  image.plot(dat_cloud_bis)
  h5write(dat_cloud_bis, file=filename, name="mask_sub", index=list(1:Ngd,1:Ngd,i))
  
  
  ## Transform data
  bn <- orderNorm(as.vector(dat[which(dat_cloud==1)]))
  y_trans <- matrix(predict(bn,dat),Ngd,Ngd)
  h5write(y_trans, file=filename, name="data_norm", index=list(1:Ngd,1:Ngd,i))
  list_transfo[[i]]=bn
}
saveRDS(list_transfo,paste0(filename_radix,"_normTransform.rds"))




#---------------------------

## Name of the export file
filename_radix="data_v10"
list_transfo=list()
filename=paste0(filename_radix,".h5")
if (file.exists(filename)) {
  file.remove(filename)
}
h5createFile(filename)
h5createDataset(filename, "data", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "data_norm", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask_sub", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
## Fill data
for(i in 1:nrow(tab_x)){
  xmin=tab_x[i,1]
  xmax=xmin+df[1]
  ymin=tab_x[i,2]
  ymax=ymin+df[2]
  
  window(r)=ext(c(xmin,xmax,ymin,ymax))
  
  # Plot the first layer
  plot(r[[idv10]],main=paste0(xmin,",",ymax,",",xmax,",",ymin))
  dat=rast_to_mat(r[[idv10]])
  print(dim(dat))
  image.plot(dat)
  h5write(dat, file=filename, name="data", index=list(1:Ngd,1:Ngd,i))
  
  
  dat_cloud=rast_to_mat(r[[idCloud]])<0.5+0
  image.plot(dat_cloud)
  h5write(dat_cloud, file=filename, name="mask", index=list(1:Ngd,1:Ngd,i))
  
  
  set.seed(42+i)
  sel_mask=sample(which(dat_cloud==1),10**3)
  dat_cloud_bis=dat_cloud*0
  dat_cloud_bis[sel_mask]=1
  image.plot(dat_cloud_bis)
  h5write(dat_cloud_bis, file=filename, name="mask_sub", index=list(1:Ngd,1:Ngd,i))
  
  
  ## Transform data
  bn <- orderNorm(as.vector(dat[which(dat_cloud==1)]))
  y_trans <- matrix(predict(bn,dat),Ngd,Ngd)
  h5write(y_trans, file=filename, name="data_norm", index=list(1:Ngd,1:Ngd,i))
  list_transfo[[i]]=bn
}
saveRDS(list_transfo,paste0(filename_radix,"_normTransform.rds"))




#---------------------------

## Name of the export file
filename_radix="data_SLHF"
list_transfo=list()
filename=paste0(filename_radix,".h5")
if (file.exists(filename)) {
  file.remove(filename)
}
h5createFile(filename)
h5createDataset(filename, "data", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "data_norm", rev(c(nrow(tab_x),Ngd,Ngd)),
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
h5createDataset(filename, "mask_sub", rev(c(nrow(tab_x),Ngd,Ngd)), 
                chunk=rev(c(1,Ngd,Ngd)),native=TRUE)
## Fill data
for(i in 1:nrow(tab_x)){
  xmin=tab_x[i,1]
  xmax=xmin+df[1]
  ymin=tab_x[i,2]
  ymax=ymin+df[2]
  
  window(rb)=ext(c(xmin,xmax,ymin,ymax))
  
  # Plot the first layer
  plot(rb[[idHF]],main=paste0(xmin,",",ymax,",",xmax,",",ymin))
  dat=rast_to_mat(rb[[idHF]])
  print(dim(dat))
  image.plot(dat)
  h5write(dat, file=filename, name="data", index=list(1:Ngd,1:Ngd,i))
  
  
  window(r)=ext(c(xmin,xmax,ymin,ymax))
  dat_cloud=rast_to_mat(r[[idCloud]])<0.5+0
  image.plot(dat_cloud)
  h5write(dat_cloud, file=filename, name="mask", index=list(1:Ngd,1:Ngd,i))
  
  set.seed(42+i)
  sel_mask=sample(which(dat_cloud==1),10**3)
  dat_cloud_bis=dat_cloud*0
  dat_cloud_bis[sel_mask]=1
  image.plot(dat_cloud_bis)
  h5write(dat_cloud_bis, file=filename, name="mask_sub", index=list(1:Ngd,1:Ngd,i))
  
  
  ## Transform data
  bn <- orderNorm(as.vector(dat[which(dat_cloud==1)]))
  y_trans <- matrix(predict(bn,dat),Ngd,Ngd)
  h5write(y_trans, file=filename, name="data_norm", index=list(1:Ngd,1:Ngd,i))
  list_transfo[[i]]=bn
}
saveRDS(list_transfo,paste0(filename_radix,"_normTransform.rds"))

