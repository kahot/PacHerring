#Prep for BayPass

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pop_info$Population.Year)
pops17<-unique(pop_info$Population.Year[grep("17", pop_info$Population.Year)])
pws<-unique(pop_info$Population.Year[grep("PWS", pop_info$Population.Year)])
pws.ss.tb<-unique(pop_info$Population.Year[grep("PWS|SS|TB", pop_info$Population.Year)])

#2017 
AC_17<-data.frame(matrix(nrow=312774))
for (i in 1:length(pops17)){
    df<-read.table(paste0("Data/new_vcf/AC/",pops[i],"_AC.txt"))
    hist(df$V1, main=paste(pops[i]))
    colnames(df)<-c("alt","total")
    df$ref<-df$total-df$alt
    df2<-df[,c("ref","alt")]
    AC_17<-cbind(AC_17, df2)
}

AC_17<-AC_17[,-1]
write.table(AC_17, "Data/new_vcf/AC/AC17.geno",row.names = F, col.names = F)



#Load the output
source("Rscripts/baypass_utils.R")
library(ggplot2)
omega<-read.table("Data/new_vcf/AC/ph2017_mat_omega.out")

plot.omega(omega,PC=c(1,2),pop.names=paste0("Pop",1:nrow(omega)),
           main=expression("SVD of "*Omega),col=rainbow(nrow(omega)),
           pch=16,pos=2)

#XtX
xtx<-read.table("Data/new_vcf/AC/ph2017_summary_pi_xtx.out")
colnames(xtx)<-xtx[1,]
xtx<-xtx[-1,]
xtx<-data.frame(sapply(xtx, as.numeric))
colnames(xtx)[7]<-"logP"

pos<-read.table("Data/new_vcf/AC/SNP_chr.pos.txt")
colnames(pos)<-c("chr","pos")
xtx<-cbind(pos,xtx)
#some are recorded as Inf. Replace it with a large number
bigP<-which(xtx$logP==Inf)
max(xtx$logP[xtx$logP!=Inf]) #max is 15.65
#replace Inf with 20
xtx$logP[bigP]<-20

#find the top loci
#inf: 60 sites (all in chr7 or chr12 inversion areas)
xtx[bigP,] 

nrow(xtx[xtx$logP>15,]) #138
nrow(xtx[xtx$logP>10,]) #1752
nrow(xtx[xtx$logP>7,]) #2961

topPos<-xtx[xtx$logP>7,] #138


ggplot(xtx[1:100000,], aes(x=MRK, y=logP))+
    geom_point(size=0.1, color="pink")

10^-(xtx$logP)
hist(10^-(xtx$logP))




#simulate to create "pseudo-observed data sets (POD)"
omega.mat=as.matrix(omega)

bta2017.data<-geno2YN("Data/new_vcf/AC/AC17.geno")
pi.beta.coef=read.table("Data/new_vcf/AC/ph2017_summary_beta_params.out",h=T)$Mean


simu.bta<-simulate.baypass(omega.mat,nsnp=300000,beta.coef=NA,sample.size=bta2017.data$NN,
                           beta.pi=pi.beta.coef,pop.trait=0,pi.maf=0.05, suffix="btapods",
                  print.sim.params.values=FALSE,output.bayenv.format=FALSE,
                  remove.fixed.loci=FALSE,coverage=NA)

#######################################################
#get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table("Data/new_vcf/AC/anapod_mat_omega.out"))
plot(pod.omega,omega) ; abline(a=0,b=1)
fmd.dist(pod.omega,omega)
#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution from the POD analysis

pod.pi.beta.coef=read.table("anapod_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef,pi.beta.coef)  ; abline(a=0,b=1)
#######################################################
#XtX calibration
#######################################################
#get the pod XtX
pod.xtx=read.table("Data/new_vcf/AC/anapod_summary_pi_xtx.out",h=T)$M_XtX
#compute the 1% threshold
pod.thresh=quantile(pod.xtx,probs=0.999)
#add the thresh to the actual XtX plot
plot(xtx$M_XtX, pch=".", col="lightblue")
abline(h=pod.thresh,lty=2, col="red")

ggplot(xtx, aes(x=MRK, y=M_XtX))+
    geom_point(size=0.1, color="lightblue")+
    geom_hline(yintercept=pod.thresh, color="red", linetype=2)+
    theme_classic()+xlab('')+ylab("XtX")
ggsave("Output/baypass_output_xtx.pdf", width = 9, height = 7)
