##*** Code needed for Peru model

# Function to calculate the dynamics over 3 years for each household with FAST PROGRESSORS SPLIT 
# fastprog 1: fitness in transmission
peru_model_hh_fastprog1 <- function(FOIv,betav, omegav, iniv,hhn, parain,time_frame){
  # FOIv = vector of force of infection. c(FOI_R, FOI_S)
  # betav = vector of beta_int and fitness c(BETA, F)
  # omegav = vector of treatment rates and treatment success c(OMEGAS, OMEGAR, KS, KR) NEEDED?? 
  # iniv = vector of initial conditions. 
  # hhn = how many in house 
  # prograte = rate of progression of fast
  
  ### Initial state
  if(length(iniv) < 6){print("Error in initial"); break}
  if(length(iniv) == 6){x0 = c(U = hhn - sum(iniv), LS = iniv[1], LR = iniv[2], AS = iniv[3], AR = iniv[4], 
                               LfS = iniv[5],LfR = iniv[6],C = 0, hhn = hhn, dead = 0)}
  
  ### States & Transitions
  # STATE ORDER IS: U LS LR AS AR LfS LfR C hhn dead
  transitions = list(c(U = -1, LS = +1),
                     c(U = -1, LR = +1),
                     c(U = -1, LfS = +1),
                     c(U = -1, LfR = +1),
                     
                     c(LfS = -1, AS = +1), 
                     c(LfR = -1, AR = +1), 
                     
                     c(LS = -1, LR = +1),
                     c(LS = -1, LfR = +1),
                     c(LS = -1, AS = +1),
                     c(LS = -1, LfS = +1),
                     
                     c(LR = -1, LS = +1),
                     c(LR = -1, LfS = +1),
                     c(LR = -1, AR = +1),
                     c(LR = -1, LfR = +1),
                     
                     c(AS = -1, AR = +1),
                     
                     c(AS = -1, LS = +1),
                     c(AR = -1, LR = +1),
                     
                     c(LS = -1,  U = +1,hhn = -1, dead = +1),
                     c(AS = -1,  U = +1,hhn = -1, dead = +1),
                     c(LR = -1,  U = +1,hhn = -1, dead = +1),
                     c(AR = -1,  U = +1,hhn = -1, dead = +1),
                     c(LfS = -1, U = +1,hhn = -1, dead = +1),
                     c(LfR = -1, U = +1,hhn = -1, dead = +1),
                     c(U = -1, C = +1,hhn = -1, dead = +1))
  
  rates <- function(x,parain,t){
    # params
    muA<-parain[["muA"]]; p<-parain[["p"]];phi<-parain[["phi"]];chi<-parain[["chi"]];eps<-parain[["eps"]];
    mu<-parain[["mu"]];muL<-parain[["muL"]];d<-parain[["d"]];n<-parain[["n"]]; 
    foi_s <- FOIv[[1]]; foi_r <- FOIv[[2]];
    beta <- betav[[1]]; f<-betav[[2]]; pf<-betav[[3]];
    omega_s <- omegav[[1]]; omega_r <- omegav[[2]]; ks <- omegav[[3]]; kr <- omegav[[4]]
    # states
    U <- x["U"]
    LR <- x["LR"]
    LS <- x["LS"]
    AS <- x["AS"]
    AR <- x["AR"]
    LfR <- x["LfR"]
    LfS <- x["LfS"]
    C <- x["C"]
    hhn <- x["hhn"]
    dead <- x["dead"]
    # output
    if(hhn > 0){return(c(U*(1-p)*(foi_s + d*beta*AS / hhn),
                         U*(1-p)*(foi_r + f*d*beta*AR / hhn),
                         U*(p)*(foi_s + d*beta*AS / hhn),
                         U*(p)*(foi_r + f*d*beta*AR / hhn),
                         LfS*pf, LfR*pf, 
                         LS*chi*(1-p)*(foi_r + f*d*beta*AR / hhn),
                         LS*chi*(p)*(foi_r + f*d*beta*AR / hhn),
                         LS*phi,
                         LS*chi*(p)*(foi_s + d*beta*AS / hhn),
                         LR*chi*(1-p)*(foi_s + d*beta*AS / hhn),
                         LR*chi*(p)*(foi_s + d*beta*AS / hhn),
                         LR*phi,
                         LR*chi*(p)*(foi_r + f*d*beta*AR / hhn),
                         AS*omega_s*eps,
                         AS*omega_s*(1-ks) + n*AS,
                         AR*omega_r*(1-kr) + n*AR,
                         LS*muL,AS*(d*muA+muL),LR*muL,AR*(d*muA+muL),LfS*muL,LfR*muL,U*muL))}
    else{
      return(rep(0,24)) # if no one in hh the noone can do anything
    }
  }
  
  parms=list(muA=parain["muA"],p=parain["p"],phi=parain["phi"],chi=parain["chi"],eps=parain["eps"],
             mu=parain["mu"],muL=parain["muL"],d=parain["d"],n=parain["n"], 
             foi_s = FOIv[1], foi_r = FOIv[2],
             beta = betav[1], f=betav[2],
             omega_s = omegav[1], omega_r = omegav[2], ks = omegav[3], kr = omegav[4])
  
  out <- ssa.adaptivetau(x0, transitions, rates, parms, tf = time_frame)
  return(out)
}


# Function to calculate the dynamics over 3 years for each household with FAST PROGRESSORS SPLIT 
#### fast prog 2 : fitness in progression
peru_model_hh_fastprog2 <- function(FOIv,betav, omegav, iniv,hhn, parain,time_frame){
  # FOIv = vector of force of infection. c(FOI_R, FOI_S)
  # betav = vector of beta_int and fitness c(BETA, F)
  # omegav = vector of treatment rates and treatment success c(OMEGAS, OMEGAR, KS, KR) NEEDED?? 
  # iniv = vector of initial conditions. 
  # hhn = how many in house 
  # prograte = rate of progression of fast
  
  ### Initial state
  if(length(iniv) < 6){print("Error in initial"); break}
  if(length(iniv) == 6){x0 = c(U = hhn - sum(iniv), LS = iniv[1], LR = iniv[2], AS = iniv[3], AR = iniv[4], 
                               LfS = iniv[5],LfR = iniv[6],C = 0, hhn = hhn, dead = 0)}
  
  ### States & Transitions
  # STATE ORDER IS: U LS LR AS AR LfS LfR C hhn dead
  transitions = list(c(U = -1, LS = +1),
                     c(U = -1, LR = +1),
                     c(U = -1, LfS = +1),
                     c(U = -1, LfR = +1),
                     
                     c(LfS = -1, AS = +1), 
                     c(LfR = -1, AR = +1), 
                     
                     c(LS = -1, LR = +1),
                     c(LS = -1, LfR = +1),
                     c(LS = -1, AS = +1),
                     c(LS = -1, LfS = +1),
                     
                     c(LR = -1, LS = +1),
                     c(LR = -1, LfS = +1),
                     c(LR = -1, AR = +1),
                     c(LR = -1, LfR = +1),
                     
                     c(AS = -1, AR = +1),
                     
                     c(AS = -1, LS = +1),
                     c(AR = -1, LR = +1),
                     
                     c(LS = -1,  U = +1,hhn = -1, dead = +1),
                     c(AS = -1,  U = +1,hhn = -1, dead = +1),
                     c(LR = -1,  U = +1,hhn = -1, dead = +1),
                     c(AR = -1,  U = +1,hhn = -1, dead = +1),
                     c(LfS = -1, U = +1,hhn = -1, dead = +1),
                     c(LfR = -1, U = +1,hhn = -1, dead = +1),
                     c(U = -1, C = +1,hhn = -1, dead = +1))
  
  rates <- function(x,parain,t){
    # params
    muA<-parain[["muA"]]; p<-parain[["p"]];phi<-parain[["phi"]];chi<-parain[["chi"]];eps<-parain[["eps"]];
    mu<-parain[["mu"]];muL<-parain[["muL"]];d<-parain[["d"]];n<-parain[["n"]]; 
    foi_s <- FOIv[[1]]; foi_r <- FOIv[[2]];
    beta <- betav[[1]]; f<-betav[[2]]; pf<-betav[[3]];
    omega_s <- omegav[[1]]; omega_r <- omegav[[2]]; ks <- omegav[[3]]; kr <- omegav[[4]]
    # states
    U <- x["U"]
    LR <- x["LR"]
    LS <- x["LS"]
    AS <- x["AS"]
    AR <- x["AR"]
    LfR <- x["LfR"]
    LfS <- x["LfS"]
    C <- x["C"]
    hhn <- x["hhn"]
    dead <- x["dead"]
    # output
    if(hhn > 0){return(c(U*(1-p)*(foi_s + d*beta*AS / hhn),
                         U*(1-p)*(foi_r + d*beta*AR / hhn),
                         U*(p)*(foi_s + d*beta*AS / hhn),
                         U*(p)*(foi_r + d*beta*AR / hhn),
                         LfS*pf, LfR*pf*f, 
                         LS*chi*(1-p)*(foi_r + d*beta*AR / hhn),
                         LS*chi*(p)*(foi_r + d*beta*AR / hhn),
                         LS*phi,
                         LS*chi*(p)*(foi_s + d*beta*AS / hhn),
                         LR*chi*(1-p)*(foi_s + d*beta*AS / hhn),
                         LR*chi*(p)*(foi_s + d*beta*AS / hhn),
                         LR*phi,
                         LR*chi*(p)*(foi_r + d*beta*AR / hhn),
                         AS*omega_s*eps,
                         AS*omega_s*(1-ks) + n*AS,
                         AR*omega_r*(1-kr) + n*AR,
                         LS*muL,AS*(d*muA+muL),LR*muL,AR*(d*muA+muL),LfS*muL,LfR*muL,U*muL))}
    else{
      return(rep(0,24)) # if no one in hh the noone can do anything
    }
  }
  
  parms=list(muA=parain["muA"],p=parain["p"],phi=parain["phi"],chi=parain["chi"],eps=parain["eps"],
             mu=parain["mu"],muL=parain["muL"],d=parain["d"],n=parain["n"], 
             foi_s = FOIv[1], foi_r = FOIv[2],
             beta = betav[1], f=betav[2],
             omega_s = omegav[1], omega_r = omegav[2], ks = omegav[3], kr = omegav[4])
  
  out <- ssa.adaptivetau(x0, transitions, rates, parms, tf = time_frame)
  return(out)
}

# Function to calculate the dynamics over 3 years for each household with FAST PROGRESSORS SPLIT 
# fastprog 3: fitness in transmission and progression
peru_model_hh_fastprog3 <- function(FOIv,betav, omegav, iniv,hhn, parain,time_frame){
  # FOIv = vector of force of infection. c(FOI_R, FOI_S)
  # betav = vector of beta_int and fitness c(BETA, F)
  # omegav = vector of treatment rates and treatment success c(OMEGAS, OMEGAR, KS, KR) NEEDED?? 
  # iniv = vector of initial conditions. 
  # hhn = how many in house 
  # prograte = rate of progression of fast
  
  ### Initial state
  if(length(iniv) < 6){print("Error in initial"); break}
  if(length(iniv) == 6){x0 = c(U = hhn - sum(iniv), LS = iniv[1], LR = iniv[2], AS = iniv[3], AR = iniv[4], 
                               LfS = iniv[5],LfR = iniv[6],C = 0, hhn = hhn, dead = 0)}
  
  ### States & Transitions
  # STATE ORDER IS: U LS LR AS AR LfS LfR C hhn dead
  transitions = list(c(U = -1, LS = +1),
                     c(U = -1, LR = +1),
                     c(U = -1, LfS = +1),
                     c(U = -1, LfR = +1),
                     
                     c(LfS = -1, AS = +1), 
                     c(LfR = -1, AR = +1), 
                     
                     c(LS = -1, LR = +1),
                     c(LS = -1, LfR = +1),
                     c(LS = -1, AS = +1),
                     c(LS = -1, LfS = +1),
                     
                     c(LR = -1, LS = +1),
                     c(LR = -1, LfS = +1),
                     c(LR = -1, AR = +1),
                     c(LR = -1, LfR = +1),
                     
                     c(AS = -1, AR = +1),
                     
                     c(AS = -1, LS = +1),
                     c(AR = -1, LR = +1),
                     
                     c(LS = -1,  U = +1,hhn = -1, dead = +1),
                     c(AS = -1,  U = +1,hhn = -1, dead = +1),
                     c(LR = -1,  U = +1,hhn = -1, dead = +1),
                     c(AR = -1,  U = +1,hhn = -1, dead = +1),
                     c(LfS = -1, U = +1,hhn = -1, dead = +1),
                     c(LfR = -1, U = +1,hhn = -1, dead = +1),
                     c(U = -1, C = +1,hhn = -1, dead = +1))
  
  rates <- function(x,parain,t){
    # params
    muA<-parain[["muA"]]; p<-parain[["p"]];phi<-parain[["phi"]];chi<-parain[["chi"]];eps<-parain[["eps"]];
    mu<-parain[["mu"]];muL<-parain[["muL"]];d<-parain[["d"]];n<-parain[["n"]]; 
    foi_s <- FOIv[[1]]; foi_r <- FOIv[[2]];
    beta <- betav[[1]]; f<-betav[[2]]; pf<-betav[[3]];
    omega_s <- omegav[[1]]; omega_r <- omegav[[2]]; ks <- omegav[[3]]; kr <- omegav[[4]]
    # states
    U <- x["U"]
    LR <- x["LR"]
    LS <- x["LS"]
    AS <- x["AS"]
    AR <- x["AR"]
    LfR <- x["LfR"]
    LfS <- x["LfS"]
    C <- x["C"]
    hhn <- x["hhn"]
    dead <- x["dead"]
    # output
    if(hhn > 0){return(c(U*(1-p)*(foi_s + d*beta*AS / hhn),
                         U*(1-p)*(foi_r + f*d*beta*AR / hhn),
                         U*(p)*(foi_s + d*beta*AS / hhn),
                         U*(p)*(foi_r + f*d*beta*AR / hhn),
                         LfS*pf, LfR*pf*f, 
                         LS*chi*(1-p)*(foi_r + f*d*beta*AR / hhn),
                         LS*chi*(p)*(foi_r + f*d*beta*AR / hhn),
                         LS*phi,
                         LS*chi*(p)*(foi_s + d*beta*AS / hhn),
                         LR*chi*(1-p)*(foi_s + d*beta*AS / hhn),
                         LR*chi*(p)*(foi_s + d*beta*AS / hhn),
                         LR*phi,
                         LR*chi*(p)*(foi_r + f*d*beta*AR / hhn),
                         AS*omega_s*eps,
                         AS*omega_s*(1-ks) + n*AS,
                         AR*omega_r*(1-kr) + n*AR,
                         LS*muL,AS*(d*muA+muL),LR*muL,AR*(d*muA+muL),LfS*muL,LfR*muL,U*muL))}
    else{
      return(rep(0,24)) # if no one in hh the noone can do anything
    }
  }
  
  parms=list(muA=parain["muA"],p=parain["p"],phi=parain["phi"],chi=parain["chi"],eps=parain["eps"],
             mu=parain["mu"],muL=parain["muL"],d=parain["d"],n=parain["n"], 
             foi_s = FOIv[1], foi_r = FOIv[2],
             beta = betav[1], f=betav[2],
             omega_s = omegav[1], omega_r = omegav[2], ks = omegav[3], kr = omegav[4])
  
  out <- ssa.adaptivetau(x0, transitions, rates, parms, tf = time_frame)
  return(out)
}





#sample from own distribution - decide who leaves. From http://stackoverflow.com/questions/12848736/how-to-declare-a-user-defined-distribution-in-r
rMydist <- function(n,probd,state_num) {
  ss<-sample(x = seq(1,state_num,1), size = n, prob = probd, replace=T) # Sample with given distribution 
  return(tabulate(ss,nbins=state_num)) # Return counts
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

