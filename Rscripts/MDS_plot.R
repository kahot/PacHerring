#PCA MDS results from Plink

library(ggplot2)
library(tidyverse)
library(stringr)
library(gridExtra)

cols<-c("#0072b2","#cc79a7","#009e73","#d55e00","#56b4e9","#e69f00","#f0e442")
#bla <- "#000000"
#blu <- "#0072b2"
#grb <- "#56b4e9"
#lir <- "#cc79a7"
#gre <- "#009e73"
#red <- "#d55e00"
#org <- "#e69f00"
#yel <- "#f0e442"
#gry<-  '#BBBBBB'


#eigenvec <- read.table("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/plink/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.eigenvec", header=F, stringsAsFactors = F)
#eigenval <- read.table("C:/Users/jmcgirr/Documents/Whitehead_Lab/ph/plink/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.eigenval", header=F, stringsAsFactors = F)

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[,c("Sample","pop","Year.Collected")]
colnames(pop_info)[3]<-"year"

#mds<-read.table("Output/plink/ph_maf0.05_mds_20.mds", header = T)
#mds<-merge(mds, pop_info,by.x="FID", by.y="Sample")

#Using 70% overlapping 58432SNPs 
mds<-read.table("Output/plink/maf05_ov70/ph_maf05_id_ov70_mds_10.mds", header = T)
#mds<-merge(mds, pop_info,by.x="FID", by.y="Sample")
mds<-cbind(mds,pop_info )


eigenvec <- read.table("Output/plink/maf05_ov70/ph_maf05_id_ov70_pca.eigenvec", header=F, stringsAsFactors = F)
eigenval <- read.table("Output/plink/maf05_ov70/ph_maf05_id_ov70_pca.eigenval", header=F, stringsAsFactors = F)

prop_explained <- c()
for (s in eigenval$V1) {
    #print(s / sum(eigenval$V1))
    prop_explained <- c(prop_explained,round(((s / sum(eigenval$V1))*100),2))
}

barplot(prop_explained, ylab = "% variance explained by PC", xlab = "PC",
        names.arg=c(1:length(prop_explained)))

#MDS plot
ggplot(data = mds)+
    geom_point(data = mds, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 3)+
    xlab(paste("Dim1"))+
    ylab(paste("Dim2"))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    ggtitle("70% overlaping SNPs")
ggsave("Output/plink/maf05_ov70/mds_ov70_all.pdf", width = 7.7, height = 6)    


#pca plot
pc<-cbind(eigenvec, pop_info)
ggplot(data = pc)+
    geom_point(data = pc, aes(x = V3, y = V4, fill = pop, color = pop, shape = factor(year)), size = 3)+
    xlab(paste("PCA1"))+
    ylab(paste("PCA2"))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Population")+
    ggtitle("70% overlaping SNPs")
ggsave("Output/plink/maf05_ov70/pca_ov70_all.pdf", width = 7.7, height = 6)    


######





mdfiles<-list.files("Output/plink/", pattern="ph_maf0.05_mds.+mds$")
mdsplots<-list()
notbplots<-list()
pwsplots<-list()
for (i in 1:length(mdfiles)){
    mds<-read.table(paste0("Output/plink/", mdfiles[i]), header = T)
    mds<-merge(mds, pop_info,by.x="FID", by.y="Sample")
    dim<-str_extract(mdfiles[i], "(?<=maf0.05_mds_).*(?=.mds)")
    mdsplots[[i]]<-ggplot(data = mds)+
        geom_point(data = mds, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 3)+
        xlab(paste("Dim1"))+
        ylab(paste("Dim2"))+
        theme_bw()+
        theme(panel.grid=element_blank())+
        scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
        scale_fill_manual(values=paste0(cols,"99"), guide="none")+
        scale_color_manual(values=cols, name="Population")+           
        ggtitle(paste0("Dim=",dim))
    
    mds2<-mds[mds$pop!="TB",]
    notbplots[[i]]<-ggplot()+
        geom_point(data = mds2, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 3)+
        xlab(paste("Dim1"))+
        ylab(paste("Dim2"))+
        theme_bw()+
        scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
        scale_fill_manual(values=paste0(cols,"99"), guide="none")+
        scale_color_manual(values=cols, name="Population")+  
        theme(panel.grid=element_blank())+
        ggtitle(paste0("Dim=",dim))
    
    mds3<-mds[mds$pop=="PWS",]
    pwsplots[[i]]<-ggplot()+
        geom_point(data = mds3, aes(x = C1, y = C2, color = factor(year)), size = 3, shape=16)+
        xlab(paste("Dim1"))+
        ylab(paste("Dim2"))+
        theme_bw()+
        labs(color="Year")+ 
        theme(panel.grid=element_blank())
    
}

pdf(paste0("Output/plink/MDSplots.comapre.pdf"), width = 15, height = 20)
do.call(grid.arrange, c(mdsplots, ncol=2))
dev.off()
pdf(paste0("Output/plink/MDSplots.noTB.comapre.pdf"), width = 15, height = 20)
do.call(grid.arrange, c(notbplots, ncol=2))
dev.off()
pdf(paste0("Output/plink/MDSplots.PWS.comapre.pdf"), width = 15, height = 20)
do.call(grid.arrange, c(pwsplots, ncol=2))
dev.off()


#####
### Subpopulation MDS
mdfiles<-list.files("Output/plink/Subpops/", pattern="_10.mds$")
mdsplots<-list()
for (i in 1:length(mdfiles)){
    mds<-read.table(paste0("Output/plink/Subpops/", mdfiles[i]), header = T)
    
    mds$Sample<-apply(mds[,1:2], 1, function(x) {
                            if (x[1]<1000) paste0("0",as.character(x[1]),"_",x[2])
                            else paste0(x[1],"_",x[2])})
    mds$Sample<-gsub(" ","",mds$Sample)
    mds<-merge(mds, pop_info,by="Sample")
    #dim<-str_extract(mdfiles[i], "(?<=maf0.05_mds_).*(?=.mds)")
    
    mdsplots[[i]]<-ggplot(data = mds)+
        geom_point(data = mds, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 3)+
        xlab(paste("Dim1"))+
        ylab(paste("Dim2"))+
        theme_bw()+
        theme(panel.grid=element_blank())+
        scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
        scale_fill_manual(values=paste0(cols,"99"), guide="none")+
        scale_color_manual(values=cols, name="Population")+           
        ggtitle(paste0(gsub("_maf0.05_mds_10.mds",'',mdfiles[i])))
}

pdf(paste0("Output/plink/Subpops/MDSplots.comapre.pdf"), width = 6, height = 10)
do.call(grid.arrange, c(mdsplots, ncol=1))
dev.off()

##subpopulation
mds<-read.table(paste0("Output/plink/Subpops/PWS_SS_maf0.05_mds_10.mds"), header = T)
mds$Sample<-apply(mds[,1:2], 1, function(x) {
    if (x[1]<1000) paste0("0",as.character(x[1]),"_",x[2])
    else paste0(x[1],"_",x[2])})
mds$Sample<-gsub(" ","",mds$Sample)
mds<-merge(mds, pop_info,by="Sample")
ggplot(data = mds)+
    geom_point(data = mds, aes(x = C1, y = C2, fill = factor(year), color = factor(year), shape = pop), size = 3)+
    xlab(paste("Dim1"))+
    ylab(paste("Dim2"))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_shape_manual(values=c(23,25,3,3,21), name="Population")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Year")        
    ggtitle(paste0(gsub("_maf0.05_mds_10.mds",'',mdfiles[i])))
ggsave("Output/plink/Subpops/PWS_SS_reversed.pdf", height = 4, width = 6)



## Pruned data MDS ##
mds<-read.table(paste0("Output/plink/ph_50_5_0.5_pruned_mds_10.mds"), header = T)
mds$Sample<-apply(mds[,1:2], 1, function(x) {
    if (x[1]<1000) paste0("0",as.character(x[1]),"_",x[2])
    else paste0(x[1],"_",x[2])})
mds$Sample<-gsub(" ","",mds$Sample)
mds<-merge(mds, pop_info,by="Sample")
ggplot(data = mds)+
    geom_point(data = mds, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 2.5)+
    xlab(paste("Dim1"))+
    ylab(paste("Dim2"))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Population")         
ggsave("Output/plink/All_mds_pruned.pdf", width = 6, height = 4)    

mds2<-mds[mds$pop!="TB",]
ggplot()+
    geom_point(data = mds2, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 3)+
    xlab(paste("Dim1"))+
    ylab(paste("Dim2"))+
    theme_bw()+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Population")+  
    theme(panel.grid=element_blank())
    

mds<-read.table(paste0("Output/plink/Subpops/PWS_SS_maf0.05_mds_10.mds"), header = T)
mds$Sample<-apply(mds[,1:2], 1, function(x) {
    if (x[1]<1000) paste0("0",as.character(x[1]),"_",x[2])
    else paste0(x[1],"_",x[2])})
mds$Sample<-gsub(" ","",mds$Sample)
mds<-merge(mds, pop_info,by="Sample")
ggplot(data = mds)+
    geom_point(data = mds, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 3)+
    xlab(paste("Dim1"))+
    ylab(paste("Dim2"))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(values=cols, name="Population")         
ggsave("Output/plink/Subpops/PWS_SS_mds_pruned.pdf", width = 6, height = 4)    


#color by year
ggplot()+
    geom_point(data = mds, aes(x = C1, y = C2, color = factor(year), fill=factor(year), shape=factor(pop)), size = 2.5)+
    xlab(paste("Dim1"))+
    ylab(paste("Dim2"))+
    theme_bw()+
   # labs(color="Year")+ 
    theme(panel.grid=element_blank())+
    scale_shape_manual(values=c(23,25), name="Population")+
    scale_fill_manual(values=paste0(cols[c(1,2,3,3,4)],"99"), guide="none")+
    scale_color_manual(values=cols[c(1,2,3,3,4)], name="Year")         
ggsave("Output/plink/Subpops/PWS_SS_mds_pruned_byYear.pdf", width = 6, height = 4)    

ggplot()+
    geom_point(data = mds, aes(x = C2, y = C3, color = factor(year), fill=factor(year), shape=factor(pop)), size = 2.5)+
    xlab(paste("Dim1"))+
    ylab(paste("Dim2"))+
    theme_bw()+
    # labs(color="Year")+ 
    theme(panel.grid=element_blank())+
    scale_shape_manual(values=c(23,25), name="Population")+
    scale_fill_manual(values=paste0(cols[c(1,2,3,3,4)],"99"), guide="none")+
    scale_color_manual(values=cols[c(1,2,3,3,4)], name="Year")    


library(rgl)
mycol=cols[c(1,2,3,3,4)]
mds$color<-mycol[as.numeric(factor(mds$year))]

plot3d( 
    x=mds$C1, y=mds$C2, z=mds$C3, 
    col = mds$color, 
    type = 's', 
    radius = .0007,
    xlab="Axis 1", ylab="Axis 2", zlab="Axis 3")

#PWS only
pws<-mds[mds$pop=="PWS",]
plot3d( 
    x=pws$C1, y=pws$C2, z=pws$C3, 
    col = pws$color, 
    type = 's', 
    radius = .0007,
    xlab="Axis 1", ylab="Axis 2", zlab="Axis 3")


#### PCA 
eigenvec <- read.table("Output/plink/Subpops/PWS_SS_pruned_pca.eigenvec", header=F, stringsAsFactors = F)
eigenval <- read.table("Output/plink/Subpops/PWS_SS_pruned_pca.eigenval", header=F, stringsAsFactors = F)

eigenvec$Sample<-apply(eigenvec[,1:2], 1, function(x) {
    if (x[1]<1000) paste0("0",as.character(x[1]),"_",x[2])
    else paste0(x[1],"_",x[2])})
eigenvec$Sample<-gsub(" ","",eigenvec$Sample)


vcf_sample_info <- merge(eigenvec,pop_info, by = "Sample")

pca <-data.frame(sample=vcf_sample_info$Sample, 
                 pop=vcf_sample_info$pop,
                 year=vcf_sample_info$year,
                 PC1=vcf_sample_info$V3,
                 PC2=vcf_sample_info$V4,
                 PC3=vcf_sample_info$V5,
                 PC4=vcf_sample_info$V6,
                 stringsAsFactors=FALSE)


prop_explained <- c()
for (s in eigenval$V1) {
    #print(s / sum(eigenval$V1))
    prop_explained <- c(prop_explained,round(((s / sum(eigenval$V1))*100),2))
}


ggplot()+
    geom_point(data = pca, aes(x = PC1, y = PC2, fill = pop, color = pop, shape = factor(year)), size = 3)+
    scale_fill_manual(values=cols, guide="none")+
    scale_color_manual(name = "population", values=cols)+
    scale_shape_manual(name = "year", values=c(21,22,23,23,25))+
    xlab(paste("PC1: ", prop_explained[1],"%",sep = ""))+
    ylab(paste("PC2: ", prop_explained[2],"%",sep = ""))+
    theme_bw()
ggsave("Output/plink/Subpops/PWS_SS_pca.pdf", width = 6, height = 4)


ggplot()+
    geom_point(data = pca, aes(x = PC2, y = PC3, fill = pop, color = pop, shape = factor(year)), size = 3)+
    scale_fill_manual(values=paste0(cols,"99"), guide="none")+
    scale_color_manual(name = "population", values=cols)+
    scale_shape_manual(name = "year", values=c(21,22,23,23,25))+
    xlab(paste("PC2: ", prop_explained[2],"%",sep = ""))+
    ylab(paste("PC3: ", prop_explained[3],"%",sep = ""))+
    theme_bw()
ggsave("Output/plink/Subpops/PWS_SS_pca_pc2_pc3.pdf", width = 6, height = 4)


mycol=cols[c(1,2,3,3,4)]
pca$color<-mycol[as.numeric(factor(pca$year))]

plot3d( 
    x=pca$PC1, y=pca$PC2, z=pca$PC3, 
    col = pca$color, 
    type = 's', 
    radius = .007,
    xlab="PC1", ylab="PC2", zlab="PC3")

htmlwidgets::saveWidget(rglwidget(width = 520, height = 520), 
                        file = "Output/plink/Subpops/PCA_PWA_SS_3Dscatter.html",
                        libdir = "libs",
                        selfcontained = FALSE
)



# Do this by chromosome
#create scripts to prune and run mds and pca for each chromosomes

sink("Data/vcfs/prune_vcf.sh")
cat("#!/bin/bash")
cat("\n")
cat(paste0("#SBATCH --job-name=plink_MDA_byChromosome \n"))
cat(paste0("#SBATCH --mem=8G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e plinkMDA_byChromosome.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu  \n"))
cat(paste0("#SBATCH --mail-type=ALL  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load plink")     
cat("\n\n")
for (i in 1:26){
    cat(paste0("#add variant ID names for pruning \n"))
    cat(paste0("plink --vcf /home/ktist/ph/data/split_vcf/ph_chr",i,".vcf.gz --set-missing-var-ids @:#[ph] --make-bed --out /home/ktist/ph/data/split_vcf/chr",i,"\n")) 
    cat(paste0("plink --bfile /home/ktist/ph/data/split_vcf/chr",i," --recode --tab --out /home/ktist/ph/data/split_vcf/chr",i,"\n"))
    cat(paste0("#prune the loci with 50 5 0.5 parameters \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/chr",i," --indep-pairwise 50 5 0.5 --out /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5 \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/chr",i," --extract /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5.prune.in --make-bed --out /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5_pruned \n"))
    cat(paste0("plink --bfile /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5_pruned --recode --tab --out /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5_pruned \n"))
    
    cat(paste0("#Run MDS and PCA \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5_pruned --mds-plot 10 --cluster --out /home/ktist/ph/data/plink/chr",i,"_pruned_mds_10 \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5_pruned --pca --out /home/ktist/ph/data/plink/chr",i,"_pruned_pca \n"))
    cat("\n")
}
sink(NULL)

##repeat this for noTB vcf
sink("Data/vcfs/prune_vcf_noTB.sh")
cat("#!/bin/bash -l")
cat("\n\n")
cat(paste0("#SBATCH --job-name=plink_MDA_byChromosome \n"))
cat(paste0("#SBATCH --mem=8G \n")) 
cat(paste0("#SBATCH --ntasks=1 \n")) 
cat(paste0("#SBATCH -e plinkMDA_byChromosome.err  \n"))
cat(paste0("#SBATCH --time=24:00:00  \n"))
cat(paste0("#SBATCH --mail-user=ktist@ucdavis.edu  \n"))
cat(paste0("#SBATCH --mail-type=ALL  \n"))
cat(paste0("#SBATCH -p high  \n"))
cat("\n")
cat("module load plink")     
cat("\n\n")
for (i in 1:26){
    cat(paste0("#add variant ID names for pruning \n"))
    cat(paste0("plink --vcf /home/ktist/ph/data/split_vcf/noTB_chr",i,".vcf.gz --set-missing-var-ids @:#[ph] --make-bed --out /home/ktist/ph/data/split_vcf/noTB_chr",i,"\n")) 
    cat(paste0("plink --bfile /home/ktist/ph/data/split_vcf/noTB_chr",i," --recode --tab --out /home/ktist/ph/data/split_vcf/noTB_chr",i,"\n"))
    cat(paste0("#prune the loci with 50 5 0.5 parameters \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/noTB_chr",i," --indep-pairwise 50 5 0.5 --out /home/ktist/ph/data/split_vcf/noTB_chr",i,"_50_5_0.5 \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/noTB_chr",i," --extract /home/ktist/ph/data/split_vcf/chr",i,"_50_5_0.5.prune.in --make-bed --out /home/ktist/ph/data/split_vcf/noTB_chr",i,"_50_5_0.5_pruned \n"))
    cat(paste0("plink --bfile /home/ktist/ph/data/split_vcf/noTB_chr",i,"_50_5_0.5_pruned --recode --tab --out /home/ktist/ph/data/split_vcf/noTB_chr",i,"_50_5_0.5_pruned \n"))
    
    cat(paste0("#Run MDS and PCA \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/noTB_chr",i,"_50_5_0.5_pruned --mds-plot 10 --cluster --out /home/ktist/ph/data/plink/noTB_chr",i,"_pruned_mds_10 \n"))
    cat(paste0("plink --file /home/ktist/ph/data/split_vcf/noTB_chr",i,"_50_5_0.5_pruned --pca --out /home/ktist/ph/data/plink/noTB_chr",i,"_pruned_pca \n"))
    cat("\n")
}
sink(NULL)





#chr23,24, adn 26 did not produce results (regardless of pruned or not)

fnames<-list.files("Output/plink/Subpops/Chr/", pattern= ".mds$")
mplots<-list()
for (i in 1: length(fnames)){
    mds<-read.table(paste0("Output/plink/Subpops/Chr/", fnames[i]), header = T)
    
    mds$Sample<-apply(mds[,1:2], 1, function(x) {
        if (x[1]<1000) paste0("0",as.character(x[1]),"_",x[2])
        else paste0(x[1],"_",x[2])})
    mds$Sample<-gsub(" ","",mds$Sample)
    mds<-merge(mds, pop_info,by="Sample")
    
    mplots[[i]]<-ggplot(data = mds)+
        geom_point(data = mds, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 2)+
        xlab(paste("Dim1"))+
        ylab(paste("Dim2"))+
        theme_bw()+
        theme(panel.grid=element_blank())+
        scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
        scale_fill_manual(values=paste0(cols,"99"), guide="none")+
        scale_color_manual(values=cols, name="Population")+           
        ggtitle(paste0(gsub("_pruned_mds_10.mds",'',fnames[i])))
    
    
}
pdf(paste0("Output/plink/Subpops/Chr/MDSplots.byChromosome.pdf"), width = 18, height = 30)
do.call(grid.arrange, c(mplots, ncol=3))
dev.off()

#noTB
fnames<-list.files("Output/plink/Subpops/noTB/", pattern= ".mds$")
mplots<-list()
for (i in 1: length(fnames)){
    mds<-read.table(paste0("Output/plink/Subpops/noTB/", fnames[i]), header = T)
    
    mds$Sample<-apply(mds[,1:2], 1, function(x) {
        if (x[1]<1000) paste0("0",as.character(x[1]),"_",x[2])
        else paste0(x[1],"_",x[2])})
    mds$Sample<-gsub(" ","",mds$Sample)
    mds<-merge(mds, pop_info,by="Sample")
    titl<-gsub("_pruned_mds_10.mds",'',fnames[i])
    titl<-gsub("noTB_",'',titl)
    
    mplots[[i]]<-ggplot(data = mds)+
        geom_point(data = mds, aes(x = C1, y = C2, fill = pop, color = pop, shape = factor(year)), size = 2)+
        xlab(paste("Dim1"))+
        ylab(paste("Dim2"))+
        theme_bw()+
        theme(panel.grid=element_blank())+
        scale_shape_manual(values=c(23,25,3,3,21), name="Year")+
        scale_fill_manual(values=paste0(cols,"99"), guide="none")+
        scale_color_manual(values=cols, name="Population")+           
        ggtitle(titl)
    
    
}
pdf(paste0("Output/plink/Subpops/noTB/MDSplots.byChromosome_noTB.pdf"), width = 18, height = 30)
do.call(grid.arrange, c(mplots, ncol=3))
dev.off()






#Subset vcf files by populations or years 
#1. I want to remove TB
noTB<-pop_info[pop_info$pop!="TB",]
sink("Data/popinfo/noTB.txt")
for (i in 1:nrow(noTB)){
    cat(paste0(noTB$Sample[i],"\n"))
}
sink(NULL)

pwsss<-pop_info[pop_info$pop=="PWS"|pop_info$pop=="SS",]
sink("Data/popinfo/PWS_SS.txt")
for (i in 1:nrow(pwsss)){
    cat(paste0(pwsss$Sample[i],"\n"))
}
sink(NULL)

#remove TB and CA

df<-pop_info[!(pop_info$pop=="TB"|pop_info$pop=="CA"),]
sink("Data/popinfo/NoTB_CA.txt")
for (i in 1:nrow(df)){
    cat(paste0(df$Sample[i],"\n"))
}
sink(NULL)

