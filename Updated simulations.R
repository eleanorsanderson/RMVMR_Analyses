##############################################
###   Install and load relevant packages   ###
##############################################


#install.packages("remotes")
#library(remotes)
#install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
#install_github("WSpiller/RadialMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

args <- as.numeric(commandArgs(T))
set.seed((args[1]*100000))
job_id <- ((args[1]))
message("job number ", job_id)

#Set random gen.seed
set.seed(12345)

#Load libraries into R
#library(RMVMR)
#library(RadialMR)
#library(MVMR)
#library(ggplot2)

source("format_rmvmr.R")
source("strength_mvmr.R")
source("pleiotropy_rmvmr.R")
source("format_radial.R")
source("ivw_rmvmr.R")
source("ivw_radial.R")

source("datagen_function.R")

#Set number of iterations
iterations<-100

#################################
### Set simulation parameters ###
#################################
#Define effect of X1
beta1<-1

#Define effect of X2
beta2<-0.2

#Define effect of X3
beta3<-(-0.5)


#Define dataframe for the simulation output 
output <- data.frame()

#########################
###   Generate data   ###
#########################

#Set number of instruments to 210 - note dat-gen function not actually flexible to this number changing

J<-210

#Set sample size
N<-100000
k = 0


for(type in c('balanced', 'unbalanced')){
  for(pct_pleio in c(0, 0.1, 0.25, 0.5, 0.75)){

for(it in 1:iterations){
  
  #Print iteration number
  print(it)
  k = k+1

  output[k,"type"] <- type
  output[k, "pct_pleio"] <- pct_pleio
  
  
  #Generate set of instruments
  G.dat<-data.frame(rep(0,N))
  
  for(i in 1:J){
    G.dat[,i]<-rnorm(N,0,1)
  }
  
  #Generate unmeasured confounder
  U<-rnorm(N,0,1)
  
  #shared effect between exposures
  #X1X2<-rnorm(N,0,1)
  #X1X3<-rnorm(N,0,1)
  #X2X1<-rnorm(N,0,1)
  #X2X3<-rnorm(N,0,1)
  #X3X1<-rnorm(N,0,1)
  #X3X2<-rnorm(N,0,1)
  
  #Generate exposures
  X1<-1 + 1*U + rnorm(N,0,1) #+ rnorm(1,0,5)*X1X2 + rnorm(1,0,5)*X1X3 + rnorm(1,0,5)*X3X1 + rnorm(1,0,5)*X2X1
  X2<-1 + 1*U + rnorm(N,0,1) #+ rnorm(1,0,5)*X1X2 + rnorm(1,0,5)*X2X1 + rnorm(1,0,5)*X2X3 + rnorm(1,0,5)*X3X2
  X3<-1 + 1*U + rnorm(N,0,1) #+ rnorm(1,0,5)*X1X3 + rnorm(1,0,5)*X2X3 + rnorm(1,0,5)*X3X1 + rnorm(1,0,5)*X3X2
  
  #Generate outcome
  Y<-1 + rnorm(N,0,1)
  

  dat <- dat_gen(J, pct_pleio, type)
  
  X1gammas <- dat$X1gammas
  X2gammas <- dat$X2gammas
  X3gammas <- dat$X3gammas
  alphavec <- dat$alphavec
  
  #Pleiotropic effects for instruments
  for (i in c(1:J)){
    Y<-Y + alphavec[i]*G.dat[,i]
  }
  
  #X1 associations
  for (i in 1:J){
    X1<-X1 + X1gammas[i]*G.dat[,i]
  }
  
  #X2 associations
  for (i in 1:J){
    X2<-X2 + X2gammas[i]*G.dat[,i]
  }
  
  #X3 associations
  for (i in 1:J){
    X3<-X3 + X3gammas[i]*G.dat[,i]
  }
  
  #Generate outcome
  Y<- Y + beta1*X1 + beta2*X2 + beta3*X3 + 1*U + rnorm(N,0,1)
  
  #Create combined dataframe
  c.dat<-data.frame(X1,X2,X3,Y)
  c.dat<-cbind(G.dat,c.dat)
  for(i in 1:J){
    names(c.dat)[i]<-paste("rs",i,sep="")
  }
  
  #Divide sample into non-overlapping subsets for generating exposure and outcome associations
  exposure.sampleX1<-c.dat[1:25000,]
  exposure.sampleX2<-c.dat[25001:50000,]
  exposure.sampleX3<-c.dat[50001:75000,]
  outcome.sampleY<-c.dat[75001:100000,]
  
  # Estimate G-X1 associations
  
  gammaX1hat<-rep(0,J)
  segammaX1hat<-rep(0,J)
  pgammaX1hat<-rep(0,J)
  
  for(i in 1:J){
    
    gammaX1hat[i]<-summary(lm(exposure.sampleX1[,J+1]~exposure.sampleX1[,i]))$coef[2,1]
    segammaX1hat[i]<-summary(lm(exposure.sampleX1[,J+1]~exposure.sampleX1[,i]))$coef[2,2]
    pgammaX1hat[i]<-summary(lm(exposure.sampleX1[,J+1]~exposure.sampleX1[,i]))$coef[2,4]
    
  }
  
  # Estimate G-X2 associations
  
  gammaX2hat<-rep(0,J)
  segammaX2hat<-rep(0,J)
  pgammaX2hat<-rep(0,J)
  
  for(i in 1:J){
    gammaX2hat[i]<-summary(lm(exposure.sampleX2[,J+2]~exposure.sampleX2[,i]))$coef[2,1]
    segammaX2hat[i]<-summary(lm(exposure.sampleX2[,J+2]~exposure.sampleX2[,i]))$coef[2,2]
    pgammaX2hat[i]<-summary(lm(exposure.sampleX2[,J+2]~exposure.sampleX2[,i]))$coef[2,4]
    
  }
  
  # Estimate G-X3 associations
  
  gammaX3hat<-rep(0,J)
  segammaX3hat<-rep(0,J)
  pgammaX3hat<-rep(0,J)
  
  for(i in 1:J){
    gammaX3hat[i]<-summary(lm(exposure.sampleX3[,J+3]~exposure.sampleX3[,i]))$coef[2,1]
    segammaX3hat[i]<-summary(lm(exposure.sampleX3[,J+3]~exposure.sampleX3[,i]))$coef[2,2]
    pgammaX3hat[i]<-summary(lm(exposure.sampleX3[,J+3]~exposure.sampleX3[,i]))$coef[2,4]
    
  }
  
  # Estimate G-Y associations
  
  gammaYhat<-rep(0,J)
  segammaYhat<-rep(0,J)
  for(i in 1:J){
    
    gammaYhat[i]<-summary(lm(outcome.sampleY[,J+4]~outcome.sampleY[,i]))$coef[2,1]
    segammaYhat[i]<-summary(lm(outcome.sampleY[,J+4]~outcome.sampleY[,i]))$coef[2,2]
    
  }
  
  #Create dataframe containing summary data
  
  sum.data<-data.frame(names(c.dat[1:J]),gammaX1hat,gammaX2hat,gammaX3hat,
                       segammaX1hat,segammaX2hat,segammaX3hat,
                       gammaYhat,segammaYhat,pgammaX1hat,pgammaX2hat,pgammaX3hat)
  
  names(sum.data)[1]<-"SNP"
  
  
  ################################################
  ###   Radial MVMR Analyses: X1, X2, and X3   ###
  ################################################
  
  # Select instruments for any exposure
  
  X1X2X3datuni<-sum.data
  
  #Format the X1-X2 data
  X1X2X3f.data<-format_rmvmr(sum.data[,2:4], sum.data[,8], sum.data[,5:7], sum.data[,9], sum.data[,1])
  
  #Estimate conditional instrument strength
  strength_X1X2X3<- suppressWarnings(strength_mvmr(X1X2X3f.data))
  
  #Save conditional F statistic for exposure 1
  output[k,'X1f']<-strength_X1X2X3[1]
  
  #Save conditional F statistic for exposure 2
  output[k,'X2f']<-strength_X1X2X3[2]
  
  #Save conditional F statistic for exposure 3
  output[k,'X3f']<-strength_X1X2X3[3]
  
  #Estimate causal effects  using Radial MVMR
  X1X2X3res<-ivw_rmvmr(X1X2X3f.data,F)
  
  #Save effect estimate and se for exposure 1
  output[k,'X1est']<-X1X2X3res$coef[1,1]
  output[k,'X1se']<-X1X2X3res$coef[1,2]
  
  #Save effect estimate and se for exposure 2
  output[k,'X2est']<-X1X2X3res$coef[2,1]
  output[k,'X2se']<-X1X2X3res$coef[2,2]
  
  #Save effect estimate and se for exposure 3
  output[k,'X3est']<-X1X2X3res$coef[3,1]
  output[k,'X3se']<-X1X2X3res$coef[3,2] 
  
  #Estimate pleiotropic effects and detect outliers
  pleioX1X2X3<-pleiotropy_rmvmr(X1X2X3f.data,X1X2X3res)
  
  #report Q statistics
  output[k,'Qx1'] <- pleioX1X2X3$gq['Exposure_1','q_statistic']
  output[k,'Qx1_pvalue'] <- pleioX1X2X3$gq['Exposure_1','p_value']
  
  output[k,'Qx2'] <- pleioX1X2X3$gq['Exposure_2','q_statistic']
  output[k,'Qx2_pvalue'] <- pleioX1X2X3$gq['Exposure_2','p_value']
  
  output[k,'Qx3'] <- pleioX1X2X3$gq['Exposure_3','q_statistic']
  output[k,'Qx3_pvalue'] <- pleioX1X2X3$gq['Exposure_3','p_value']
  
  #############################################################
  ###   Radial MVMR Analyses: X1, X2, and X3 with pruning   ###
  #############################################################
  X1X2X3fp.data<-X1X2X3f.data
  
  outliers<-"start"
  
  while(!is.null(outliers)){
    
    X1X2X3pres<-ivw_rmvmr(X1X2X3fp.data,F)
    pleioX1X2X3p<-pleiotropy_rmvmr(X1X2X3fp.data,X1X2X3pres)
    
    ###Prune###
    
    if(min(pleioX1X2X3p$qdat$qj_p)< 0.05){
      outliers<-pleioX1X2X3p$qdat[pleioX1X2X3p$qdat$qj_p <0.05,]$snp
      outliers<-unique(outliers)
    
      
      X1X2X3fp.data<-X1X2X3fp.data[!X1X2X3fp.data$SNP %in% outliers,]
      
    }else{
      outliers<-NULL
      
    }
    
  }
  output[k,'noutliers'] <- length(X1X2X3f.data$SNP) - length(X1X2X3fp.data$SNP)
  
  output[k,'Qx1p'] <- pleioX1X2X3p$gq['Exposure_1','q_statistic']
  output[k,'Qx1p_pvalue'] <- pleioX1X2X3p$gq['Exposure_1','p_value']
  
  output[k,'Qx2p'] <- pleioX1X2X3p$gq['Exposure_2','q_statistic']
  output[k,'Qx2p_pvalue'] <- pleioX1X2X3p$gq['Exposure_2','p_value']
  
  output[k,'Qx3p'] <- pleioX1X2X3p$gq['Exposure_3','q_statistic']
  output[k,'Qx3p_pvalue'] <- pleioX1X2X3p$gq['Exposure_3','p_value']
  
  strength_X1X2X3p<- suppressWarnings(strength_mvmr(X1X2X3fp.data))
  
  #Save conditional F statistic for exposure 1
  output[k,'X1fp']<-strength_X1X2X3p[1]
  
  #Save conditional F statistic for exposure 2
  output[k,'X2fp']<-strength_X1X2X3p[2]
  
  #Save conditional F statistic for exposure 3
  output[k,'X3fp']<-strength_X1X2X3p[3]
  
  #Estimate causal effects  using Radial MVMR
  X1X2X3resp<-ivw_rmvmr(X1X2X3fp.data,F)
  
  #Save effect estimate and se for exposure 1
  output[k,'X1pest']<-X1X2X3resp$coef[1,1]
  output[k,'X1pse']<-X1X2X3resp$coef[1,2]
  
  #Save effect estimate and se for exposure 2
  output[k,'X2pest']<-X1X2X3resp$coef[2,1]
  output[k,'X2pse']<-X1X2X3resp$coef[2,2]
  
  #Save effect estimate and se for exposure 3
  output[k,'X3pest']<-X1X2X3resp$coef[3,1]
  output[k,'X3pse']<-X1X2X3resp$coef[3,2]
  
}
  }
}


save(output, file=sprintf("rmvmr_sims_%s.Rda", job_id))
