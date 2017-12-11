## Plot output_parallel

library(reshape)
library(ggplot2)
library(plyr)

nfits = 5000

#### Want
# all fits = summary measures for each parameter set that gives a fit - don't write out if doesn't fit
all_fits<-c()
# kp data = only write out if parameter set fits
kp_all <- c()

for(i in 1:nfits){
  #if(file.exists(paste("store_",i,".csv",sep=""))){}
  if(file.exists(paste("bi_fp2_summary_fit_",i,".csv",sep=""))){
    st <- read.csv(paste("bi_fp2_summary_fit_",i,".csv",sep=""))[,-1]
    all_fits<-rbind(all_fits,t(st))
  }
  if(file.exists(paste("bi_fp2_kp_",i,".csv",sep=""))){
    kp <- read.csv(paste("bi_fp2_kp_",i,".csv",sep=""))[,-1]
    kp_all <- rbind(kp_all,kp)
  }
}

write.csv(all_fits,"bi_all_fits.csv")
write.csv(kp_all,"bi_kp_all.csv")

############# All fits to data
all_fits<- as.data.frame(all_fits)
colnames(all_fits) <- c("fit","TBIr_s", "TBIr_r", "TBIs_s", "TBIs_r","total")#,"lat_prop")            
all_fitsm <- melt(all_fits[1:100,c("TBIr_s", "TBIr_r", "TBIs_s", "TBIs_r")],id.vars=c()) # Perhaps only plot first 1,000? ****
# data points
datap<-as.data.frame(rbind(c(4264,3916,4338),c(87,13,435),c(344,98,810),c(2112,1646,2358)))
datap$variable<-c("TBIs_s", "TBIs_r", "TBIr_s", "TBIr_r")
colnames(datap)<-c("mean","lower95","upper95","variable")
# plot
g<-ggplot(all_fitsm, aes(x=variable, y=value))+geom_jitter()
g<- g+ geom_point(data = datap, aes(x=variable, y=mean),colour=c("blue","blue","red","red"),size=3) + theme_bw(base_size = 28)
g<- g+ geom_errorbar(data = datap, aes(x=variable,y=mean,ymin=lower95,ymax=upper95),colour=c("blue","blue","red","red"),width=0.4)
g<- g+ scale_x_discrete("Household type",
                        labels=c("DS-TB in HH\n MDR-index case","MDR-TB in HH\n MDR-index case",
                                 "DS-TB in HH\n DS-index case","MDR-TB in HH\n DS-index case")) + scale_y_continuous("TB Incidence per 100,000")  

ggsave("bi_Fit_model2.pdf",width=12,height=8)


# What else of interest? latent proportion = 40% at end
#plot(all_fits$lat_prop) # getting around 70%... 

############# Plot kaplan meier. 1 = R, 2 = S
kp_all<-as.data.frame(kp_all); colnames(kp_all)<-c("type","run","time","prob")
w<-which(kp_all$type == 2); max_run = max(kp_all$run)
kp_all[w,"run"] = kp_all[w,"run"] + max_run # need to group separately so add on max runs
kp_all[,"prob"]	= kp_all[,"prob"]/100

g<-ggplot(kp_all,aes(x=time,y=prob,col=factor(type),group=factor(run)))+geom_line()
ggsave("bi_Fit_model2_kp_plain.pdf",width=12,height=8)
g<-g+scale_colour_manual("",breaks = c("1","2"),values=c("red","green"),labels = c("Contacts of Drug\nResistant Tuberculosis","Contacts of Drug\nSusceptible Tuberculosis"))
g<-g+scale_x_continuous("Follow up Time (days)") + scale_y_continuous("Probability of remaining Free of Tuberculosis",lim=c(0.91,1))
#g<-g+ theme(legend.key.size = unit(1, "cm"))
ggsave("bi_Fit_model2_kp.pdf",width=12,height=8)

kpro<-kp_all
kpro$time = round(kp_all$time,1) # round time so as to group better. 
gg<-ddply(kpro,.(time,type),summarise,median=median(prob)) 
g<-ggplot(gg,aes(x=time,y=median,col=factor(type)))+geom_line()
g<-g+scale_colour_manual("",breaks = c("1","2"),values=c("red","green"),labels = c("Contacts of Drug\nResistant Tuberculosis","Contacts of Drug\nSusceptible Tuberculosis"))
g<-g+scale_x_continuous("Follow up Time (days)") + scale_y_continuous("Probability of remaining Free of Tuberculosis",lim=c(0.91,1))
ggsave("bi_Fit_model2_kp_round.pdf",width=12,height=8)

### parameters of fits

## Inputs
nsteps = 5000
fits1<-c()
for(i in 1:nfits){
  if(file.exists(paste("bi_fp2_fits_",nsteps,"_",i,".csv",sep=""))){
  ffn <- read.csv(paste("bi_fp2_fits_",nsteps,"_",i,".csv",sep=""))[-1,-1] # remove first row as that's the basic start guess
  fits1 <- rbind(fits1, ffn) # bind together
 }
}

## Correlations
pdf("bi_foisvsf.pdf")
plot(fits1$fois, fits1$f)
dev.off()

pdf("bi_foirvsf.pdf")
plot(fits1$foir, fits1$f)
dev.off()

pdf("bi_betavsf.pdf")
plot(fits1$beta, fits1$f)
dev.off()

library(GGally); library(ggplot2)
pdf("bi_correl_fp2.pdf")
pm <- ggpairs(all_fits)
pm
dev.off()
`
