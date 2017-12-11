### functions for running mcmc - with burn-in

## function - distance of output from data FOR FAST PROG
dist_peru_data_fp1 <- function(parainf){
  ns = 487; nr = 213  # number f susceptible or resistant households
  
  # The hh size distribution for this set of households
  hh_vec_dr <- rMydist(nr, hh$Freq, 14)  # inputs = number, prob distribution, maxm number to sample from - here 1 to 15 people per house
  hh_vec_ss <- rMydist(ns, hh$Freq, 14)  # inputs = number, prob distribution, maxm number to sample from - here 1 to 15 people per house
  
  #####******* PRE-STUDY ************************************************************************************************************************************************
  # Initiate with background infection rate etc
  # Run for 
  time_frame = 30 # years?
  
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
  
  # If only doing minimum set
  if(length(parainf) < 5){
    FOIv = parainf[c(1,2)]
    betav = c(parainf[c(3,4)],0.2) # pf = 0.2
    para_in = c(0.26,0.15,0.00013,0.35,0.008,0.013,0.013,0.5,0.2)
    names(para_in)<-c("muA","p","phi","chi","eps","mu","muL","d","n")
  }else{
    # Input FOIv & betav
    FOIv = parainf[c(10,11)]
    betav = parainf[c(12,13,14)]
    para_in = parainf[c(1:9)]}
  
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
  
  #ff = c(fitt, TBIs_s,  c(3916,4338),TBIs_r,  c(13,435),TBIr_s,  c(98,810),TBIr_r,  c(1646,2358))
  print(c(TBIs_s,  c(3916,4338),TBIs_r,  c(13,435),TBIr_s,  c(98,810),TBIr_r,  c(1646,2358)))
  
  ### For kaplan meier?
  inn = prod(parainf[1],parainf[2],parainf[3],parainf[4],10000) # UNIQUE number? 
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
  #print(c("kp",kp))
  
  write.csv(store,paste("bi_fp1_store_",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep=""))
	print(paste("bi_fp1_store_",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep=""))  

  # Store to plot
  if (fitt == 1){
    output <- c(fitt, TBIr_s, TBIr_r, TBIs_s, TBIs_r,
                abs(TBIr_s - 367) + abs(TBIr_r - 2253) + abs(TBIs_s - 4264) + abs(TBIs_r - 87))
    write.csv(output, paste("bi_fp1_summary_fit",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep="")) # unique? 
    write.csv(kp,paste("bi_fp1_kp_",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep="")) # only keep summary and kp if fit
   
  }
  
  
  # Return
  return(fitt)
}

dist_peru_data_fp2 <- function(parainf){
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
  
  # If only doing minimum set
  if(length(parainf) < 5){
    FOIv = parainf[c(1,2)]
    betav = c(parainf[c(3,4)],0.2) # pf = 0.2
    para_in = c(0.26,0.15,0.00013,0.35,0.008,0.013,0.013,0.5,0.2)
    names(para_in)<-c("muA","p","phi","chi","eps","mu","muL","d","n")
  }else{
    # Input FOIv & betav
    FOIv = parainf[c(10,11)]
    betav = parainf[c(12,13,14)]
    para_in = parainf[c(1:9)]}
  
  # timepoints for initial conditions for study
  tp <- matrix(0,nr+ns,11)
  
  for(ik in 1:(nr+ns)){
    #                   FOIv,       betav,   omegav (ws,wr,ks,kr),           iniv,                                           hhn,parain
    out<-peru_model_hh_fastprog2(FOIv, betav,c(0.8,0.64,0.26,0.4), as.numeric(summ_data[ik,c("Ls","Lr","As","Ar","Lfs","Lfr")]),as.numeric(summ_data[ik,"hhsize"]),para_in,time_frame)
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
    out<-peru_model_hh_fastprog2(FOIv, betav,c(2,2,0.26,0.4), as.numeric(tp[ik,c("LS","LR","AS","AR","LfS","LfR")]),
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
  
  #ff = c(fitt, TBIs_s,  c(3916,4338),TBIs_r,  c(13,435),TBIr_s,  c(98,810),TBIr_r,  c(1646,2358))
  ### For kaplan meier?
  inn = prod(parainf[1],parainf[2],parainf[3],parainf[4],10000) # UNIQUE number? 
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
  #print(c("kp",kp))
  
  write.csv(store,paste("bi_fp2_store_",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep=""))
  
  # Store to plot
  if (fitt == 1){
    output <- c(fitt, TBIr_s, TBIr_r, TBIs_s, TBIs_r,
                abs(TBIr_s - 367) + abs(TBIr_r - 2253) + abs(TBIs_s - 4264) + abs(TBIs_r - 87))
    write.csv(output, paste("bi_fp2_summary_fit",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep="")) # unique? 
    write.csv(kp,paste("bi_fp2_kp_",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep="")) # only keep summary and kp if fit
    
  }
  
  # Return
  return(fitt)
}

dist_peru_data_fp3 <- function(parainf){
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
  
  # If only doing minimum set
  if(length(parainf) < 5){
    FOIv = parainf[c(1,2)]
    betav = c(parainf[c(3,4)],0.2) # pf = 0.2
    para_in = c(0.26,0.15,0.00013,0.35,0.008,0.013,0.013,0.5,0.2)
    names(para_in)<-c("muA","p","phi","chi","eps","mu","muL","d","n")
  }else{
    # Input FOIv & betav
    FOIv = parainf[c(10,11)]
    betav = parainf[c(12,13,14)]
    para_in = parainf[c(1:9)]}
  
  # timepoints for initial conditions for study
  tp <- matrix(0,nr+ns,11)
  
  for(ik in 1:(nr+ns)){
    #                   FOIv,       betav,   omegav (ws,wr,ks,kr),           iniv,                                           hhn,parain
    out<-peru_model_hh_fastprog3(FOIv, betav,c(0.8,0.64,0.26,0.4), as.numeric(summ_data[ik,c("Ls","Lr","As","Ar","Lfs","Lfr")]),as.numeric(summ_data[ik,"hhsize"]),para_in,time_frame)
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
    out<-peru_model_hh_fastprog3(FOIv, betav,c(2,2,0.26,0.4), as.numeric(tp[ik,c("LS","LR","AS","AR","LfS","LfR")]),
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
  
  #ff = c(fitt, TBIs_s,  c(3916,4338),TBIs_r,  c(13,435),TBIr_s,  c(98,810),TBIr_r,  c(1646,2358))
  
  ### For kaplan meier?
  inn = prod(parainf[1],parainf[2],parainf[3],parainf[4],10000) # UNIQUE number?
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
  #print(c("kp",kp))
  
  write.csv(store,paste("bi_fp3_store_",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep=""))
  
  # Store to plot
  if (fitt == 1){
    output <- c(fitt, TBIr_s, TBIr_r, TBIs_s, TBIs_r,
                abs(TBIr_s - 367) + abs(TBIr_r - 2253) + abs(TBIs_s - 4264) + abs(TBIs_r - 87))
    write.csv(output, paste("bi_fp3_summary_fit",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep="")) # unique? 
    write.csv(kp,paste("bi_fp3_kp_",parainf[1],parainf[2],parainf[3],parainf[4],".csv",sep="")) # only keep summary and kp if fit
  }
  
  
  # Return
  return(fitt)
}

## Evaluate the prior at each parameter 
# is the parameter in the prior? how likely is the parameter in the prior? if uniform then always 1. 
#para_peru<-read.csv("para_peru.csv",stringsAsFactors = FALSE)

prior_peru_ratio <- function(para_2){
  p1<-1; p2<-1;
  npar = length(para_2)
  #   for(i in 1:npar){
  #     p1 <-p1*dunif(as.numeric(para_1[i]),min=para_peru[which(para_peru[,"para"]==names(para_1[i])),"min"],
  #                   max=para_peru[which(para_peru[,"para"]==names(para_1[i])),"max"])
  #     print(c(i,p1,as.numeric(para_1[i]),min=para_peru[which(para_peru[,"para"]==names(para_1[i])),"min"],
  #             max=para_peru[which(para_peru[,"para"]==names(para_1[i])),"max"]))
  #   }
  for(i in 1:npar){
    min=para_peru[which(para_peru[,"para"]==names(para_2[i])),"min"]
    max=para_peru[which(para_peru[,"para"]==names(para_2[i])),"max"]
    if(min < para_2[i] && para_2[i]< max){p2 = 1*p2}else{p2 = 0*p2;print(names(para_2[i]))}
    #print(c(min,para_2[i],max))
    #  p2 <-p2*dunif(as.numeric(para_2[i]),min=para_peru[which(para_peru[,"para"]==names(para_2[i])),"min"],max=para_peru[which(para_peru[,"para"]==names(para_2[i])),"max"])
  }
  #print(c(p1,p2))
  # if all the parameters in set para1 and para2 are in the priors then the value for p1 and p2 will be the same
  # ratio of new to old
  # https://theclevermachine.wordpress.com/2012/11/04/mcmc-multivariate-distributions-block-wise-component-wise-updates/
  return(p2)#(p2/p1)
}

## function to select parameters, run model, check against last parameter
peru_mcmc <- function(target, init.theta, proposal.sd, n.iterations) {
  
  # evaluate the function "target" at "init.theta", and assign to
  # a variable called target.theta.current.
  target.theta.current <- target(init.theta)
  
  # initialise variables to store the current value of theta, the
  # vector of samples, and the number of accepted runs
  theta.current <- init.theta
  samples <- theta.current
  accepted <- 0
  
  # run MCMC for n.iteration interations
  for (i.iteration in seq_len(n.iterations)) {
    
    # draw a new theta from the (Gaussian) proposal distribution
    # and assign to a variable called theta.proposed.  
    # See "?rnorm for more information
    # Note that this step is vectorized for any arbitratry theta 
    # which will be useful when we will sample from a multivariate
    # target distribution
    theta.proposed <- rtruncnorm(n = length(theta.current),
                                 mean = as.numeric(theta.current),
                                 sd = as.numeric(proposal.sd),
                                 a = rep(0,length(theta.current)))
    
    # set the names of theta.proposed to be the same as the names of
    # theta.current
    names(theta.proposed) <- names(theta.current)
    
    # evaluate the function target at the proposed theta and
    # assign to a variable called target.theta.proposed
    # Target gives the absolute value difference to the data
    target.theta.proposed <- target(theta.proposed)
    
    # compute Metropolis-Hastings ratio (acceptance probability). Since
    # the multivariate Gaussian is symmetric, we don't need to consider
    # the proposal distribution here
    #print(theta.proposed)
    # Prior_peru_ratio gives a 0 if not in prior, 1 if in
    
    acceptance <- target.theta.proposed * prior_peru_ratio(theta.proposed)
    
    #log.acceptance <- (target.theta.proposed - target.theta.current) * prior_peru_ratio(theta.current, theta.proposed)
    #print(c(log.acceptance, target.theta.proposed - target.theta.current,
    #         prior_peru_ratio(theta.current, theta.proposed)))
    #     # draw random number number between 0 and 1 using "runif" and assign to
    #     # a variable called r.
    r <- runif(1)
    
    # test acceptance by comparing the random number to the
    # Metropolis-Hastings ratio (acceptance probability) (using
    # "exp" because we calculated the logarithm of the
    # Metropolis-Hastings ratio before)
    if (r < acceptance) { # exp(log.acceptance)
      #  if(target.theta.proposed == 1) {
      # if accepted:
      # change the current value of theta to the proposed theta
      theta.current <- theta.proposed
      
      # updated the current value of the target
      target.theta.current <- target.theta.proposed
      
      # update number of accepted proposals
      accepted <- accepted + 1
    }
    
    # add the current theta to the vector of samples
    # Note that we use `rbind` in order to deal with multivariate 
    # target. So if `theta` is a vector then `samples` is a matrix.
    samples <- rbind(samples, theta.current, deparse.level=0)
    
    # print current state of chain and acceptance rate
    # use paste() to deal with the case where `theta` is a vector
    message("iteration: ", i.iteration, ", chain:", paste(theta.current, collapse=" "),
            ", acceptance rate:", accepted / i.iteration, ", prior yes/no ", prior_peru_ratio(theta.proposed))
    
  }
  
  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
}

## function to select parameters, run model, check against last parameter
# to guess initial SD - lhs from prior rather than sampling based on current
peru_mcmc_lhs <- function(target, init.theta, proposal.sd, n.iterations) {
  
  # evaluate the function "target" at "init.theta", and assign to
  # a variable called target.theta.current.
  target.theta.current <- target(init.theta)
  
  # initialise variables to store the current value of theta, the
  # vector of samples, and the number of accepted runs
  theta.current <- init.theta
  samples <- theta.current
  accepted <- 0
  
  # generate LHS samples for a max of n.iterations (if all accepted need n.iterations)
  # Y holds all these
  X <- randomLHS(n.iterations, length(init.theta))
  Y <- matrix(0, nrow=n.iterations, ncol=length(init.theta))
  for(i in 1:length(init.theta)){
    w<-which(para_peru$para == names(parainfm)[i])
    Y[,i] <- qunif(X[,i], min=para_peru[w,"min"], max=para_peru[w,"max"])  
  }
  
  # run MCMC for n.iteration interations
  for (i.iteration in seq_len(n.iterations)) {
    
    # draw a new theta from the (Gaussian) proposal distribution
    # and assign to a variable called theta.proposed.  
    # See "?rnorm for more information
    # Note that this step is vectorized for any arbitratry theta 
    # which will be useful when we will sample from a multivariate
    # target distribution
    theta.proposed <- Y[i.iteration,]
    
    # set the names of theta.proposed to be the same as the names of
    # theta.current
    names(theta.proposed) <- names(theta.current)
    
    # evaluate the function target at the proposed theta and
    # assign to a variable called target.theta.proposed
    target.theta.proposed <- target(theta.proposed)
    
    # compute Metropolis-Hastings ratio (acceptance probability). Since
    # the multivariate Gaussian is symmetric, we don't need to consider
    # the proposal distribution here
    #print(theta.proposed)
    
    acceptance <- target.theta.proposed * prior_peru_ratio(theta.proposed)
    
    #log.acceptance <- (target.theta.proposed - target.theta.current) * prior_peru_ratio(theta.current, theta.proposed)
    #print(c(log.acceptance, target.theta.proposed - target.theta.current,
    #         prior_peru_ratio(theta.current, theta.proposed)))
    #     # draw random number number between 0 and 1 using "runif" and assign to
    #     # a variable called r.
    r <- runif(1)
    
    # test acceptance by comparing the random number to the
    # Metropolis-Hastings ratio (acceptance probability) (using
    # "exp" because we calculated the logarithm of the
    # Metropolis-Hastings ratio before)
    if (r < acceptance) { # exp(log.acceptance)
      #  if(target.theta.proposed == 1) {
      # if accepted:
      # change the current value of theta to the proposed theta
      theta.current <- theta.proposed
      
      # updated the current value of the target
      target.theta.current <- target.theta.proposed
      
      # update number of accepted proposals
      accepted <- accepted + 1
    }
    
    # add the current theta to the vector of samples
    # Note that we use `rbind` in order to deal with multivariate 
    # target. So if `theta` is a vector then `samples` is a matrix.
    samples <- rbind(samples, theta.current, deparse.level=0)
    
    # print current state of chain and acceptance rate
    # use paste() to deal with the case where `theta` is a vector
    message("iteration: ", i.iteration, ", chain:", paste(theta.current, collapse=" "),
            ", acceptance rate:", accepted / i.iteration, ", prior yes/no ", prior_peru_ratio(theta.proposed))
    
  }
  
  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
}

#sample from own distribution - decide who leaves. From http://stackoverflow.com/questions/12848736/how-to-declare-a-user-defined-distribution-in-r
rMydist <- function(n,probd,state_num) {
  ss<-sample(x = seq(1,state_num,1), size = n, prob = probd, replace=T) # Sample with given distribution 
  return(tabulate(ss,nbins=state_num)) # Return counts
}


