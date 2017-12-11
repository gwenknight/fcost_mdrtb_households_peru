### After run MCMC_LHS what do? 
# Read in sd_fits and use in MCMC runs

# Parallel code number
args <- commandArgs(trailingOnly = TRUE)

## MCMC code
# Libraries
library("truncnorm"); library("coda"); library("compiler")
library("adaptivetau"); library("plyr")
enableJIT(1)

# Only load at beginning
source("bi_gillespie_simple_peru_code_withRecovered_allM2.R")
source("bi_mcmc_functions.R")
hh<-as.data.frame(read.csv("hhsize.csv"))
fu<-as.data.frame(read.csv("fu_manip.csv"))
para_peru<-as.data.frame(read.csv("para_peru_fp.csv"))

# Start parameters
parainf = c(0.26,0.15,0.00013,0.35,0.008,0.013,0.013,0.5,0.2, 0.025, 0.0005, 22, 0.5, 0.2)
names(parainf) <- c("muA","p","phi","chi","eps","mu","muL","d","n","fois","foir","beta","f","pf")    
# Minimum set
parainfm = c(0.2,0.1,90,0.5) #c(0.025, 0.0005, 22, 0.5)

parainfm = c(0.1447608, 0.093632, 75.73684, 0.320248) # Second set for check of fits
names(parainfm) <- c("fois","foir","beta","f")    

# From LHS - check number here
sd_fits = read.csv("bi_fp3_sdfits_10000.csv")[1,2:5]

# How many want to run? 
nruns = 5000
out_mcmc <- peru_mcmc(dist_peru_data_fp3, parainfm, sd_fits, nruns)
write.csv(out_mcmc,paste("bi_fp3_out_mcmc_",nruns,"_",args[1],".csv",sep=""))

# Output
m <- mcmc(out_mcmc)
print(1 - rejectionRate(m)) # print rejection rate

# what is the sd?
fits<-c();sd_fits<-c()
out_mcmc <- as.data.frame(out_mcmc)
u <- unique(out_mcmc$fois) # which new, accepted ones?
w<-c()
for (i in 1:length(u)){
  w<-c(w,which(out_mcmc$fois == u[i])[1])
}
fits <- rbind(fits, out_mcmc[w,]) # bind all fits found together
sd_fits <- numcolwise(sd)(fits) # Empiric sd can use for the future guess of the MCMC
write.csv(fits,paste("bi_fp3_fits_",nruns,"_",args[1],".csv",sep=""))
write.csv(sd_fits,paste("bi_fp3_sdfits_",nruns,"_",args[1],".csv",sep=""))
