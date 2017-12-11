# Plot output from MCMC

# Parallel code number
args <- commandArgs(trailingOnly = TRUE)

# Libraries
library("truncnorm"); library("coda"); library("compiler")
library("adaptivetau"); library("plyr"); library("reshape2")
library("ggplot2")
enableJIT(1)

# Only load at beginning
source("bi_gillespie_simple_peru_code_withRecovered_allM2.R")
source("bi_mcmc_functions.R")
hh<-as.data.frame(read.csv("hhsize.csv"))
para_peru<-as.data.frame(read.csv("para_peru_fp.csv"))
fu<-as.data.frame(read.csv("fu_manip.csv"))

## Read in all unique
fp1_u <- read.csv("bi_fp1_unique_all_fits.csv")[,-1]

##********** FP1 **************************************************************************
#### Run model for this fit
fits<-fp1_u
nfits = dim(fits)[1]

inn = args[1]

# Plot output as all and ranges?  
FOIv = as.numeric(fits[inn,c(1,2)])
betav = c(as.numeric(fits[inn,c(3,4)]),0.2)
para_in = c(0.26,0.15,0.00013,0.35,0.008,0.013,0.013,0.5,0.2)
names(para_in) <- c("muA","p","phi","chi","eps","mu","muL","d","n")   
print(fits[inn,])
# initial conditions
ns = 487; nr = 213  # number f susceptible or resistant households

# The hh size distribution for this set of households
hh_vec_dr <- rMydist(nr, hh$Freq, 14)  # inputs = number, prob distribution, maxm number to sample from - here 1 to 15 people per house
hh_vec_ss <- rMydist(ns, hh$Freq, 14)  # inputs = number, prob distribution, maxm number to sample from - here 1 to 15 people per house
  
  #####******* PRE-STUDY ************************************************************************************************************************************************
  # Initiate with background infection rate etc
  # Run for 
  time_frame = 10 # years?
  
  # Summary data. First column = hh number, Second = hh size, Third = MDR (1) or SS (0) hh case
  xxr<-unlist(lapply(seq_along(hh_vec_dr),function(x){rep(x+1,hh_vec_dr[x])}))
  xxs<-unlist(lapply(seq_along(hh_vec_ss),function(x){rep(x+1,hh_vec_ss[x])}))
  summ_data<-c()
  summ_data<-as.data.frame(cbind(seq(1,(nr+ns),1),c(matrix(1,nr,1),matrix(0,ns,1)),c(xxr,xxs)))
  colnames(summ_data)<-c("hhnum","rors","hhsize") ### 1 in rors = MDR
  
  # Initial proportion in latent cases
  rr<-rnorm((nr+ns),0.5,0.1) # from Martinez. TST positive. What proportion MDR? Fit this? 
  # Split 98% vs 2% 
  summ_data$Latents = floor(rr*(summ_data$hhsize-1)) # include active case: so -1
  w<-which(summ_data$Latents < 0); if(length(w)>0){summ_data$Latents[w] = 0} # Check none negative
  ssL = rbinom(nr+ns, summ_data$Latents, 0.98) ### Do random sampling so get some latent resistance!
  ssR = summ_data$Latents - ssL
  summ_data$Ls = round(0.97*ssL,0) ### 3% to fast
  summ_data$Lr = round(0.97*ssR,0)  
  summ_data$Lfs = ssL - summ_data$Ls 
  summ_data$Lfr = ssR - summ_data$Lr
  if(sum(summ_data$Latents) != sum((summ_data$Ls + summ_data$Lr + summ_data$Lfs + summ_data$Lfr))){print(c(summ_data,"L sum error"))} # Check all latents dolled out
  # Initial conditions for active
  summ_data$As = c(matrix(0,nr,1),matrix(1,ns,1))
  summ_data$Ar = summ_data$rors # 1 if MDR, 0 o/w
  
  # Check not more people allocated than in hh 
  w<-which(rowSums(summ_data[,c("Ls","Lr","As","Ar","Lfs","Lfr")])>summ_data$hhsize)
  if(length(w>0)){c(print("ERROR with hh distribution"), w)}
  
  # timepoints for initial conditions for study
  tp <- matrix(0,nr+ns,11)
  
  for(ik in 1:(nr+ns)){
    #                   FOIv,       betav,   omegav (ws,wr,ks,kr),           iniv,                                           hhn,parain
    out<-peru_model_hh_fastprog1(FOIv, betav,c(0.8,0.64,0.26,0.4), as.numeric(summ_data[ik,c("Ls","Lr","As","Ar","Lfs","Lfr")]),as.numeric(summ_data[ik,"hhsize"]),para_in,time_frame)
    # Choose time point for each hh
    if(ik < (nr+1)){
      w<-which(out[,"AR"] > 0) # ok if more than one
      tp[ik,] = out[sample(w,1),] # Take last entry - as after that they were treated (or died) = like being found in the study
    }else{
      w<-which(out[,"AS"] > 0)
      tp[ik,] = out[sample(w,1),]
    }
  }
  colnames(tp)<-colnames(out)
  
  #####******* IN-STUDY ************************************************************************************************************************************************
  store<-c()
  # Follow-up time differs for each household
  time_frame_study_r = rMydist(nr, fu$mdr, dim(fu)[1])
  time_frame_study_s = rMydist(ns, fu$sens, dim(fu)[1])
  xxrt<-unlist(lapply(seq_along(time_frame_study_r),function(x){rep(round(fu[x,"midpoint"],0),time_frame_study_r[x])}))
  xxst<-unlist(lapply(seq_along(time_frame_study_s),function(x){rep(round(fu[x,"midpoint"],0),time_frame_study_s[x])}))
  time_framev = c(xxrt, xxst)/365 # in increasing order - need to make random
  time_framev = time_framev[sample(nr+ns)] # random permutation
  for(ik in 1:(nr+ns)){
    #                   FOIv,       betav,   omegav (ws,wr,ks,kr),           iniv,                                           hhn,parain
    out<-peru_model_hh_fastprog1(FOIv, betav,c(2,2,0.26,0.4), as.numeric(tp[ik,c("LS","LR","AS","AR","LfS","LfR")]),
                                 as.numeric(summ_data[ik,"hhsize"]),para_in, time_framev[ik])
    store<-rbind(store,cbind(ik,out))
  }
  rownames(store)<-NULL
  store<-as.data.frame(store); colnames(store)<-c("runs","time","U","LS","LR","AS","AR","LfS","LfR","C","hhn","dead")
  
  # Incidence wanted
  store$newr = 0; store$news = 0  
  for(iij in 1:(nr+ns)){w<-which(store$runs==iij);
                        if(length(w)>1){for(iip in 2:length(w)){if(store[w[iip],"AR"]>0){store[w[iip],"newr"] = max(0,store[w[iip],"AR"] - store[w[iip-1],"AR"])}}}} # if same will be 0, otherwise incidence, max removes those that were A and now not
  for(iij in 1:(nr+ns)){w<-which(store$runs==iij); 
                        if(length(w)>1){for(iip in 2:length(w)){if(store[w[iip],"AS"]>0){store[w[iip],"news"] = max(0,store[w[iip],"AS"] - store[w[iip-1],"AS"])}}}} # if same will be 0, otherwise incidence, max removes those that were A and now not
  
  # Death incidence wanted too
  store$death = 0; 
  for(iij in 1:(nr+ns)){w<-which(store$runs==iij);
                        if(length(w)>1){for(iip in 2:length(w)){if(store[w[iip],"dead"]>0){store[w[iip],"death"] = max(0,store[w[iip],"dead"] - store[w[iip-1],"dead"])}}}} # if same will be 0 - count when "new" death
  
  # If die then can't follow up - not in denominator of incidence
  store$timedead = 0;
  for(iij in 1:(nr+ns)){w<-which(store$runs==iij);
                        if(sum(store[w,"death"]==1)){ # if someone died - how much time to take off?
                          xx<-which(store[w,"death"]>0); time_of_death = store[w[xx[1]],"time"]
                          store[w[xx[1]],"timedead"] = time_framev[iij]-time_of_death # not 3 years any more but the specific time to followup
                        }
                        if(sum(store[w,"death"])>1){ # if someone died - how much time to take off?
                          xx<-which(store[w,"death"]>0);
                          for(xxj in 1:length(xx)){ # for all the deaths
                            time_of_death = store[w[xx[xxj]],"time"]
                            store[w[xx[xxj]],"timedead"] = time_framev[iij]-time_of_death # not 3 years any more but the specific time to followup
                          }}
  }
  
  dd<-ddply(store,.(runs),summarise,sums=sum(news),sumr=sum(newr),sumtdead = sum(timedead)) # sum number of new cases by hh
  
  # TB INCIDENCE
  # In numerator don't have to subtract off initial infectious case as in the 2:length(w) above
  # On demonimator need to take off the one incident case from each hh and the time off for those that died
  # in contacts of MDR index, how many any disease?
  TBIr_t = 100000 * (sum(dd[1:nr,c("sumr","sums")])) / (sum((summ_data[1:nr,"hhsize"]-1)*time_framev[1:nr]) - sum(dd[1:nr,"sumtdead"])) 
  # in contacts of S index, how many any disease?
  TBIs_t = 100000 * (sum(dd[(nr+1):(nr+ns),c("sums","sumr")])) / (sum((summ_data[(nr+1):(nr+ns),"hhsize"]-1)*time_framev[(nr+1):(nr+ns)]) - sum(dd[(nr+1):(nr+ns),"sumtdead"]))
  # in contacts of R how many s disease? and vv
  TBIr_r = round(100000 * (sum(dd[1:nr,c("sumr")])) / (sum((summ_data[1:nr,"hhsize"]-1)*time_framev[1:nr]) - sum(dd[1:nr,"sumtdead"])))
  TBIr_s = round(100000 * (sum(dd[1:nr,c("sums")])) / (sum((summ_data[1:nr,"hhsize"]-1)*time_framev[1:nr]) - sum(dd[1:nr,"sumtdead"]))) # don't need to subtract as no index
  TBIs_r = round(100000 * (sum(dd[(nr+1):(nr+ns),c("sumr")])) / (sum((summ_data[(nr+1):(nr+ns),"hhsize"])*time_framev[(nr+1):(nr+ns)]) - sum(dd[(nr+1):(nr+ns),"sumtdead"]))) # don't need to subtract as no index
  TBIs_s = round(100000 * (sum(dd[(nr+1):(nr+ns),c("sums")])) / (sum((summ_data[(nr+1):(nr+ns),"hhsize"])*time_framev[(nr+1):(nr+ns)]) - sum(dd[(nr+1):(nr+ns),"sumtdead"])))
  
  # FIT
  fitt = 1
  if(TBIs_s > 3916 && TBIs_s < 4338){fitt = 1*fitt}else{fitt = 0*fitt}
  if(TBIs_r > 13   && TBIs_r < 435) {fitt = 1*fitt}else{fitt = 0*fitt}
  if(TBIr_s > 98   && TBIr_s < 810) {fitt = 1*fitt}else{fitt = 0*fitt}
  if(TBIr_r > 1646 && TBIr_r < 2358){fitt = 1*fitt}else{fitt = 0*fitt}
  
print(c(fitt, TBIs_s,  c(3916,4338),TBIs_r,  c(13,435),TBIr_s,  c(98,810),TBIr_r,  c(1646,2358)))

# Proportion latent also of interest
#tt = which(store$time == 3)
#dd2<-ddply(store[tt,],.(runs),summarise,sumL=sum(LS,LR,LfS,LfR),total=sum(U,LS,LR,AS,AR,LfS,LfR)) # sum number of new cases by run
#dd2$latent_prop = 100 * dd2$sumL / dd2$total
#w<-which(dd2$latent_prop=="NaN") # if divide by 0 then get NaN... 
#if(length(w)>0){dd2[w,"latent_prop"]=0}
#if(mean(dd2$latent_prop)=="NaN"){print(dd2)}

### For kaplan meier?
store$time<-store$time * 365 # convert to days
# which R and S hh?
wr<-which(store$runs < (nr+1)); ws<-which(store$runs > (nr))
storer<-store[wr,]; stores<-store[ws,]
# for contacts of MDR
storer<-storer[order(storer[,"time"]),] # order by time
w<-which(storer$time==0)
total_ppl <- sum(storer[w,"hhn"])
storer$newcase <- cumsum(storer$news) + cumsum(storer$newr)
storer$pfree<-100*(total_ppl - storer$newcase)/total_ppl
kp <- cbind(1,inn,storer$time,storer$pfree)
# for contacts of DS-TB
stores<-stores[order(stores[,"time"]),] # order by time
w<-which(stores$time==0)
total_ppl <- sum(stores[w,"hhn"])
stores$newcase <- cumsum(stores$news) + cumsum(stores$newr)
stores$pfree<-100*(total_ppl - stores$newcase)/total_ppl
kp <- rbind(kp, cbind(2,inn,stores$time,stores$pfree))
colnames(kp)<-c("type","run","time","prob"); 
print(c("kp",kp))


# Output
summary_fit = c(fitt, TBIr_s, TBIr_r, TBIs_s, TBIs_r, 
                  abs(TBIr_s - 367) + abs(TBIr_r - 2253) + abs(TBIs_s - 4264) + abs(TBIs_r - 87)) #,mean(dd2$latent_prop))
#colnames(summary_fit) <- c("fit","TBIr_s", "TBIr_r", "TBIs_s", "TBIs_r","total","lat_prop")            

write.csv(store,paste("bi_fp1_store_",args[1],".csv",sep=""))


if(fitt == 1){write.csv(summary_fit,paste("bi_fp1_summary_fit_",args[1],".csv",sep=""));
              write.csv(kp,paste("bi_fp1_kp_",args[1],".csv",sep=""))} # only keep summary and kp if fit





