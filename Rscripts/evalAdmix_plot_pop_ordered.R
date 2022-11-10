# create slurm scripts to run evalAdmix

qfiles<-list.files("Data/ngsadmix/qfiles/")
ffiles<-list.files("Data/ngsadmix/", pattern =".fopt.gz")

sink("Output/Slurmscripts/runevalAdmix.sh")
cat("#!/bin/bash -l\n\n")
cat("#SBATCH --job-name=run_evalAdmix
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error run_evalAdmix_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --mail-user=ktist@ucdavis.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL
#SBATCH -p high \n")
cat("\n")
for (i in 1: length(qfiles)){
    cat("evalAdmix -beagle ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.BEAGLE.PL.gz ")
    cat(paste0("-fname ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_k",i,".fopt.gz "))
    cat(paste0("-qname ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_k",i,".qopt "))
    cat("-P 20\n")
    #cat(paste0("-o output.corres.txt"))
    cat(paste0("mv output.corres.txt output.corres.k",i,".txt\n"))
    cat("\n")
}
sink(NULL)
#

#Create bash scripts to run locally
sink("evalAdmix_runlocal.sh")
cat("#!/bin/bash\n\n")

for (i in 1: length(qfiles)){
    cat("evalAdmix -beagle Data/ngsadmix/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.BEAGLE.PL.gz ")
    cat(paste0("-fname Data/ngsadmix/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_k",i,".fopt.gz "))
    cat(paste0("-qname Data/ngsadmix/qfiles/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_k",i,".qopt "))
    cat("-P 10 ")
    cat(paste0("-o output.corres.k",i,".txt\n"))
}
sink(NULL)




##################################
## Plot the results from evalAdmix

source("Rscripts/visFuns.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(tibble)
library(colorspace)

# create population metadata data frame and save as "Sample_metadata_892pops.csv"
pop_info <- read.table("Data/familiarize/EVOS_MasterSheet_JoeMcGirr_April2020.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
aligned <- read.table("Data/familiarize/aligned_samples.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
pop_info <- pop_info[pop_info$Sample %in% aligned$sample,]

vcf_samples <- read.table("Data/plink/plates_1_through_5_rm.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(vcf_samples)[names(vcf_samples)=="V1"] <- "Sample"
vcf_sample_info <- inner_join(vcf_samples,pop_info, by = "Sample")

write.csv(vcf_sample_info,"Data/Sample_metadata_892pops.csv", row.names = F)

###########
# Plot for various Ks
# list output files
ofiles<-list.files("Data/ngsadmix/",pattern="output.corres")

pop<-read.csv("Data/Sample_metadata_892pops.csv")

pop$Population.Year<-factor(pop$Population.Year, levels=c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))

poporder<-paste(pop$Population.Year[order(pop$Population.Year)])
pop_order<-c("TB91","TB96","TB06","TB17","PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17")
for (i in 1:length(ofiles)){
    # extract K from the file name
    oname<-ofiles[i]
    k<-as.integer(substr(oname, 16,17))
    
    #read the qopt file for k=k
    q<-read.table(paste0("Data/ngsadmix/PH_maf05_k",k,".qopt"))
    
    #order according to population and plot the NGSadmix results
    q$id<-pop$Population.Year
    q<-q[order(q$id),]
    
    
    ord<-orderInds(pop = as.vector(poporder), q = q[,1:(i+1)])
    
    xlabels<-data.frame(x=tapply(1:length(poporder),list(poporder), mean))
    xlabels$pop<-factor(rownames(xlabels), levels=pop_order)
    xlabels<-xlabels[order(xlabels$pop),]
    
    
    barplot(t(q[,1:3]),col=c(4,2,7),space=0,border=NA,xaxt="n",xlab="Population",ylab=paste0("Admixture proportions for K=",k))
    
    
    pdf(paste0("Output/ngsadmix/Admix_plot_k",k,"_small.pdf"), height = 3.5, width=8)
    barplot(t(q[,1:3])[,ord],col=c(4,2,7),space=0,border=NA,xaxt="n",xlab="Population",ylab=paste0("Admixture proportions for K=",k))
    text(xlabels$x,-0.05,xlabels$pop,xpd=T, srt=90, adj=1,cex=0.8)
    abline(v=cumsum(sapply(unique(poporder[ord]),function(x){sum(pop[ord,"Population.Year"]==x)})),col=1,lwd=1.2)
    dev.off()
    
    #Plot the correlation matrix from evalAdmix
    r<-read.table(paste0("Output/ngsadmix/",ofiles[i]))
    
    # Plot correlation of residuals
    pdf(paste0("Output/ngsadmix/evalAdmix_corplot_k",k,".pdf"), height = 12, width=15)
    plotCorRes(cor_mat = r, pop = as.vector(pop[,"Population.Year"]), ord = ord, title=paste0("Evaluation of admixture proportions with K=",k), max_z=0.1, min_z=-0.1)
    dev.off()
}



# Plot matrix within Pop_Year/Pop
#Using k=3
r<-read.table(paste0("Output/ngsadmix/output.corres.k3.txt"))
pop<-read.csv("Data/Sample_metadata_892pops.csv")
q<-read.table("Data/ngsadmix/qfiles/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_k3.qopt")
# order according to population and plot the NGSadmix reults
ord<-orderInds(pop = as.vector(pop$Population.Year), q = q)

#Create data frame of correlation estimates
r.mat<-r
colnames(r.mat)<-pop[ord,"Sample"]
rownames(r.mat)<-pop[ord,"Sample"]

#remove the lower half
r.mat[upper.tri(r.mat)] <- NA
ind.order<-rownames(r.mat)
r.mat$Sample1<-ind.order
r.df<-melt(r.mat, na.rm = T, id.vars="Sample1", variable.name='Sample2', value.name="Correlation")
r.df$Sample2<-as.character(r.df$Sample2)

# add info for filtering purpose
r.df$Pop_Yr1<-substr(r.df$Sample1,6, 10)
r.df$Pop1<-substr(r.df$Sample1,6, 7)
r.df$Pop_Yr2<-substr(r.df$Sample2,6, 10)
r.df$Pop2<-substr(r.df$Sample2,6, 7)

#set the individual in original order
r.df$Sample1<-factor(r.df$Sample1, levels=ind.order)
r.df$Sample2<-factor(r.df$Sample2, levels=ind.order)

#Select TB2017 populations
ggplot(r.df[r.df$Pop_Yr1=="TB17"&r.df$Pop_Yr2=="TB17",], 
       aes(x=Sample1, y=Sample2, fill=Correlation))+
    geom_tile(color="white")+
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, l3 = 0, p3 = .8, p4 = .6)+
    theme_minimal()+ xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 90, size=5),axis.text.y=element_text(size=5))


df<-r.df[r.df$Pop_Yr1=="TB17"&r.df$Pop_Yr2=="TB17",]
#Using scale_fill_gradientn
ggplot(df, aes(x=Sample1, y=Sample2, fill=Correlation))+
    geom_tile()+theme_bw()+
    scale_fill_gradientn(colors=c("darkblue","white", "darkred"), limits=c(-0.2,0.5), values=c(0,0.28,1))+
    theme(axis.text.x = element_text(angle = 90, size=5),axis.text.y=element_text(size=5))+
    xlab("")+ylab("")
ggsave("Output/ngsadmix/K3_corr.mat_TB17.pdf", width = 7, height = 5)


ggplot(df, aes(x=Sample1, y=Sample2, fill=Correlation))+
    geom_tile(color="white")+theme_bw()+
    scale_fill_gradient2(low="blue",mid="white",high="red", midpoint=0,name="Correlation")+
    theme(axis.text.x = element_text(angle = 90, size=6),axis.text.y=element_text(size=6))+
    xlab('')+ylab('')


#Create correlation matrix from evalAdmix results for k=3
populations<-unique(pop$Population.Year)
for (i in 1:length(populations)){
    p<-populations[i]
    df<-r.df[r.df$Pop_Yr1==p&r.df$Pop_Yr2==p,]
    
    ggplot(df, aes(x=Sample1, y=Sample2, fill=Correlation))+
        geom_tile()+theme_bw()+
        scale_fill_gradientn(colors=c("darkblue","white", "darkred"), limits=c(-0.2,0.5), values=c(0,0.28,1))+
        theme(axis.text.x = element_text(angle = 90, size=5),axis.text.y=element_text(size=5))+
        xlab("")+ylab("")+ggtitle(paste0("NGSadmix evaluation for k=3: ",p))
    ggsave(paste0("Output/ngsadmix/K3_corr.mat_",p,".pdf"), width = 7, height = 5)
}


#same plots for k=4
r<-read.table(paste0("Output/ngsadmix/output.corres.k4.txt"))
q<-read.table("Data/ngsadmix/qfiles/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_k4.qopt")
ord<-orderInds(pop = as.vector(pop$Population.Year), q = q)

r.mat<-r
colnames(r.mat)<-pop[ord,"Sample"]
rownames(r.mat)<-pop[ord,"Sample"]

#remove the lower half
r.mat[lower.tri(r.mat)] <- NA
ind.order<-rownames(r.mat)
r.mat$Sample1<-ind.order
r.df<-melt(r.mat, na.rm = T, id.vars="Sample1", variable.name='Sample2', value.name="Correlation")
r.df$Sample2<-as.character(r.df$Sample2)

# add info for filtering purpose
r.df$Pop_Yr1<-substr(r.df$Sample1,6, 10)
r.df$Pop1<-substr(r.df$Sample1,6, 7)
r.df$Pop_Yr2<-substr(r.df$Sample2,6, 10)
r.df$Pop2<-substr(r.df$Sample2,6, 7)

#set the individual in original order
r.df$Sample1<-factor(r.df$Sample1, levels=ind.order)
r.df$Sample2<-factor(r.df$Sample2, levels=ind.order)


for (i in 1:length(populations)){
    p<-populations[i]
    df<-r.df[r.df$Pop_Yr1==p&r.df$Pop_Yr2==p,]
    
    ggplot(df, aes(x=Sample1, y=Sample2, fill=Correlation))+
        geom_tile()+theme_bw()+
        scale_fill_gradientn(colors=c("darkblue","white", "darkred"), limits=c(-0.2,0.5), values=c(0,0.28,1))+
        theme(axis.text.x = element_text(angle = 90, size=5),axis.text.y=element_text(size=5))+
        xlab("")+ylab("")+ggtitle(paste0("NGSadmix evaluation for k=4: ",p))
    ggsave(paste0("Output/ngsadmix/K4_corr.mat_",p,".pdf"), width = 7, height = 5)
}


phred<-function(x){
    score<-10^-(x/10)
    return(score)
}
phred(3)    

