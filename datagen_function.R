#Functions to define the different data going into the Radial MR estimation

#parameters taken by the function are:
#number of SNPS
#percentage that are pleiotropic
#type of pleiotropy (balanced or unbalanced)
#Note: the pleiotropic SNPs are randomly distributed across the full set of SNPs


dat_gen <- function(nSNPS, pct_pleio, type){
  
  
  X1gammas<-rnorm(30,0,1)
  X1gammas<-c(X1gammas,rep(0,30))
  X1gammas<-c(X1gammas,rep(0,30))
  X1gammas<-c(X1gammas,rnorm(30,0,1))
  X1gammas<-c(X1gammas,rnorm(30,0,1))
  X1gammas<-c(X1gammas,rep(0,30))
  X1gammas<-c(X1gammas,rnorm(30,0,1))
  
  
  X2gammas<-rep(0,30)
  X2gammas<-c(X2gammas,rnorm(30,0,1))
  X2gammas<-c(X2gammas,rep(0,30))
  X2gammas<-c(X2gammas,rnorm(30,0,1))
  X2gammas<-c(X2gammas,rep(0,30))
  X2gammas<-c(X2gammas,rnorm(30,0,1))
  X2gammas<-c(X2gammas,rnorm(30,0,1))
  
  X3gammas<-rep(0,30)
  X3gammas<-c(X3gammas,rep(0,30))
  X3gammas<-c(X3gammas,rnorm(30,0,1))
  X3gammas<-c(X3gammas,rep(0,30))
  X3gammas<-c(X3gammas,rnorm(30,0,1))
  X3gammas<-c(X3gammas,rnorm(30,0,1))
  X3gammas<-c(X3gammas,rnorm(30,0,1))
  
  
  pleiotropyvec <- as.numeric(runif(nSNPS, 0, 1) < pct_pleio)
  
  if(type == 'balanced'){
    alphavec<-(rnorm(nSNPS,0,2))*pleiotropyvec
  }
  if(type == 'unbalanced'){
    alphavec<-abs(rnorm(nSNPS,0,2))*pleiotropyvec
  }
  
  dat <- data.frame(cbind(X1gammas,X2gammas,X3gammas,alphavec))
  return(dat)
}