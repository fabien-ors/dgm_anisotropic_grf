###########################################
# Prediction comparison
###########################################


library(NSGP)
library(Matrix)
library(rhdf5)
library(scoringRules)
library(fields)

# Package for viridis color palette
library(viridis)
viridis_colors <- viridis(128)

plot_sim=function(x,zlim=NULL,...){
  if(!is.null(zlim)){
    xplt=x
    xplt[x<zlim[1]]=zlim[1]
    xplt[x>zlim[2]]=zlim[2]
  }else{
    xplt=x
  }
  image.plot(matrix(xplt,nrow = Ngd)[,Ngd:1],zlim=zlim,...)
}

path_files="rev/"

method_name=c("O3","u10","v10","SLHF")

export=F

#--------------------------------------------------------------

Ngd=256
Xgd=expand.grid(seq(from=0,to=1,length.out=256),seq(from=0,to=1,length.out=256))

#------------------------------------

n_cases=3
crps_tab=data.frame(matrix(NA,nrow = 3,ncol=2*length(method_name)))
names(crps_tab)=c(paste0(method_name,"_mean"),paste0(method_name,"_sd"))
es_tab=vs_tab=vsw_tab=crps_tab

#-------------------------------------

nsamples=100
n_cval=32

for(id_gen in 1:length(method_name)){
  
  dat=h5dump(paste0(path_files,"data/data_",method_name[id_gen],".h5"))
  
  ttl=paste0("~/Documents/Work/doc/diff_post_gauss_paper/images/app_crps_",method_name[id_gen],".png")
  if(export) png(ttl,width = 10,height = 4,res=600,units = "in",bg = "transparent")
  par(mfrow=c(1,3))
  
  for(id_case in 1:3){
    
    z_complete=c(dat$data_norm[, ,id_case])
    sel_obs=which(c(dat$mask[, ,id_case])==1)
    sel_pred=which(c(dat$mask[, ,id_case])==0)
    
    set.seed(42)
    shulfled_indices=sample(1:length(sel_pred),length(sel_pred))
    cval_indices=list()
    invDist_weigths=list()
    npairs=list()
    for (id_val in 0:(n_cval-1)) {
      ipp=floor(seq(from=0,to=length(sel_pred),length.out=n_cval+1))
      cval_indices[[id_val+1]]=shulfled_indices[(ipp[id_val+1]+1):(ipp[id_val+2])]
      invDist=1/as.matrix(dist(Xgd[cval_indices[[id_val+1]],]))
      diag(invDist)=min(invDist)*1e-6
      invDist_weigths[[id_val+1]]=invDist/max(invDist)
      npairs[[id_val+1]]=length(cval_indices[[id_val+1]])*(length(cval_indices[[id_val+1]])-1)/2
    }
    
    z_mask=rep(NA,length(z_complete)); z_mask[sel_obs]=1
    z_obs=z_complete[sel_obs]
    z_unobs=z_complete[sel_pred]
    
    cat("Number of observations:",length(z_obs),"\n")
    

    cat("Method",method_name[id_gen],"/ Case",id_case,"\n")
    fnn=paste0(path_files,"MGDM/data_",method_name[id_gen],"_",id_case-1,".h5")
    print(h5ls(fnn))
    dat_pred=h5dump(fnn)$samples
    
    y_complete=(dat_pred[, ,1,])
    y_complete=array(y_complete,dim=c(Ngd*Ngd,dim(y_complete)[length(dim(y_complete))]))
    
    id_na=apply(y_complete,2,function(x){any(is.na(x))})
    cat("Number of samples without NaNs: ",sum(!id_na),"\n")
    sel_samples=which(!id_na)
    y_complete=y_complete[,sel_samples]
    cat("---------------\n")
    
    y_obs=y_complete[sel_obs,]
    y_unobs=y_complete[sel_pred,]
    
    sco=do.call(c,lapply(cval_indices, function(l){mean(crps_sample(z_unobs[l],y_unobs[l,]))}))
    crps_tab[id_case,id_gen]=mean(sco)
    crps_tab[id_case,length(method_name)+id_gen]=sd(sco)
    
    sco=do.call(c,lapply(cval_indices, function(l){(es_sample(z_unobs[l],y_unobs[l,]))}))
    es_tab[id_case,id_gen]=mean(sco)
    es_tab[id_case,length(method_name)+id_gen]=sd(sco)
    
    # sco=do.call(c,lapply(cval_indices, function(l){(vs_sample(z_unobs[l],y_unobs[l,],p=2))}))
    sco=sapply(1:length(cval_indices), function(l){(vs_sample(z_unobs[cval_indices[[l]]],y_unobs[cval_indices[[l]],],p=2)/(npairs[[l]]))})
    vs_tab[id_case,id_gen]=mean(sco)
    vs_tab[id_case,length(method_name)+id_gen]=sd(sco)
    
    sco=sapply(1:length(cval_indices), function(l){(vs_sample(z_unobs[cval_indices[[l]]],y_unobs[cval_indices[[l]],],w_vs = invDist_weigths[[l]],p=2)/(npairs[[l]]))})
    vsw_tab[id_case,id_gen]=mean(sco)
    vsw_tab[id_case,length(method_name)+id_gen]=sd(sco)
    
    ## Create map
    all_crps=crps_sample(z_unobs,y_unobs)
    crps_map=rep(NA, Ngd*Ngd)
    crps_map[sel_pred]=all_crps
    plot_sim(crps_map,col=viridis_colors,main=method_name[id_gen],zlim=c(0,3.8))
    
  }
  par(mfrow=c(1,1))
  if(export) dev.off()
}


## Create CSV
ndig=3
col.names=c("","img_id",method_name)
n_methods=length(method_name)

crsp_val=matrix(NA,nrow = n_cases,ncol=n_methods+2)
crsp_val[,2]=0:(n_cases-1)
crsp_val[,1]=0
for(i in 1:n_cases){
  crsp_val[i,(1:n_methods)+2]=paste0(round(crps_tab[i,1:n_methods],ndig)," (",round(crps_tab[i,n_methods+1:n_methods],ndig),") |")
}
colnames(crsp_val)=col.names
if(export) write.csv(crsp_val,"~/Documents/Work/doc/diff_post_gauss_paper/data/crps_O3_post_temp.csv",row.names = F,quote = F)


es_val=matrix(NA,nrow = n_cases,ncol=n_methods+2)
es_val[,2]=0:(n_cases-1)
es_val[,1]=0
for(i in 1:n_cases){
  es_val[i,(1:n_methods)+2]=paste0(round(es_tab[i,1:n_methods],ndig)," (",round(es_tab[i,n_methods+1:n_methods],ndig),") |")
}
colnames(es_val)=col.names
if(export) write.csv(es_val,"~/Documents/Work/doc/diff_post_gauss_paper/data/es_O3_post_temp.csv",row.names = F,quote = F)


vs_val=matrix(NA,nrow = n_cases,ncol=n_methods+2)
vs_val[,2]=0:(n_cases-1)
vs_val[,1]=0
for(i in 1:n_cases){
  vs_val[i,(1:n_methods)+2]=paste0(round(vs_tab[i,1:n_methods],ndig)," (",round(vs_tab[i,n_methods+1:n_methods],ndig),") |")
}
colnames(vs_val)=col.names
if(export) write.csv(vs_val,"~/Documents/Work/doc/diff_post_gauss_paper/data/vs_O3_post_temp.csv",row.names = F,quote = F)


vsw_val=matrix(NA,nrow = n_cases,ncol=n_methods+2)
vsw_val[,2]=0:(n_cases-1)
vsw_val[,1]=0
for(i in 1:n_cases){
  vsw_val[i,(1:n_methods)+2]=paste0(round(vsw_tab[i,1:n_methods],ndig)," (",round(vsw_tab[i,n_methods+1:n_methods],ndig),") |")
}
colnames(vsw_val)=col.names
if(export) write.csv(vsw_val,"~/Documents/Work/doc/diff_post_gauss_paper/data/vsw_O3_post_temp.csv",row.names = F,quote = F)


#----------------------------------------------------
