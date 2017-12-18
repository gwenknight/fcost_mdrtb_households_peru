## MCMC_PERU full algorithm to find sd *********************FOR CLUSTER*********************
# Libraries
library("truncnorm"); library("coda");library("compiler")
library("adaptivetau"); library("plyr"); library("lhs")
enableJIT(1)

# Only load at beginning
source("bi_gillespie_simple_peru_code_withRecovered_allM2.R")
source("bi_mcmc_functions.R")
hh<-as.data.frame(read.csv("hhsize.csv"))
fu<-as.data.frame(read.csv("fu_manip.csv"))
para_peru<-as.data.frame(read.csv("para_peru_fp.csv"))

# Start parameters
parainf = c(0.26,0.15,0.00013,0.35,0.008,0.013,0.013,0.5,0.2, 0.025, 0.0005, 22, 0.5)
names(parainf) <- c("muA","p","phi","chi","eps","mu","muL","d","n","fois","foir","beta","f")    
# Minimum set
parainfm = c(0.2,0.1,90,0.5) # This is for non fp model c(0.025, 0.0005, 22, 0.5)
names(parainfm) <- c("fois","foir","beta","f")    

# initial guess
sd_fits = matrix(0.01,1,length(parainfm))
print("sd_fits OK")
# Run for 10,000 with LHS sampling instead of from normal distribution for proposal 
# to get empiric guess for SD
nruns = 10000 # 10,000 should take 1.5hrs

#system.time(peru_mcmc_lhs(dist_peru_data, parainfm, sd_fits, nruns)) # to check how long for 1 = 4.46secs!

out_mcmc_lhs <- peru_mcmc_lhs(dist_peru_data_fp3, parainfm, sd_fits, nruns)
write.csv(out_mcmc_lhs,paste("bi_fp3_out_mcmc_lhs",nruns,".csv",sep=""))

m <- mcmc(out_mcmc_lhs)
print(1 - rejectionRate(m))

# what is the sd?
fits<-c();sd_fits<-c()
out_mcmc_lhs <- as.data.frame(out_mcmc_lhs)
u <- unique(out_mcmc_lhs$fois) # which new, accepted ones?
w<-c()
for (i in 1:length(u)){
  w<-c(w,which(out_mcmc_lhs$fois == u[i])[1])
}
fits <- rbind(fits, out_mcmc_lhs[w,]) # bind all fits found together
sd_fits <- numcolwise(sd)(fits) # Empiric sd can use for the future guess of the MCMC
write.csv(fits,paste("bi_fp3_fits_",nruns,".csv",sep=""))
write.csv(sd_fits,paste("bi_fp3_sdfits_",nruns,".csv",sep=""))

