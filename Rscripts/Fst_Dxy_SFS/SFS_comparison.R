#Pacific Herring SFS and resulting Fst from ANGSD/realSFS 
source("Rscripts/BaseScripts.R")

## 2D SFS (from VCF maf00) ##

# The output from ANGSD is a flatten matrix: each value is the count of sites with the corresponding joint frequency ordered as
# [0,0] [0,1] [0,2] ..

#PWS17 = 56, PWS07=46
sfs2d<-read.table("Data/new_vcf/angsd/fromVCF/2D/folded_PWS07_PWS17_maf00.sfs")
n1<-46 #pop1 # of individuals
n2<-56
sfs<-t(matrix(sfs2d, nrow=n2*2+1, ncol=n1*2+1))
sfsdf<-data.frame(lapply(data.frame(t(sfs)), unlist), stringsAsFactors=FALSE)
max(sfsdf)

#plot 2DSFS as heatmap
#create a matrix from ANGSD output
vec2mat<-function(vec, n1,n2, pop1, pop2){
    n1<-n1
    n2<-n2
    pop1<-pop1
    pop2<-pop2
    ANGSD.2D.SFS <- scan(paste(vec, sep=""), quiet=T)
    ANGSD.2D.SFS <- t(matrix(ANGSD.2D.SFS, nrow=n2*2+1, ncol=n1*2+1))
    # mask non-variant sites
    ANGSD.2D.SFS[1,1] <- 0
    ANGSD.2D.SFS[nrow(ANGSD.2D.SFS),ncol(ANGSD.2D.SFS)] <- 0
    df<-data.frame(ANGSD.2D.SFS)
    colnames(df)<-0:(ncol(df)-1)
    df$count<-0:(nrow(df)-1)
        
    return(df)
}

#Plot heatmap     

pops.info<-read.csv("Data/Sample_metadata_892pops.csv")
pops.info$yr<-''
pops.info$yr[pops.info$year==96|pops.info$year==91]<-paste0(19,pops.info$year[pops.info$year==96|pops.info$year==91])
pops.info$yr[pops.info$year==07|pops.info$year==06|pops.info$year==17]<-paste0(20,pops.info$year[pops.info$year==07|pops.info$year==06|pops.info$year==17])
pops.info$yr<-apply(pops.info["yr"], 1, function(x) {if(x==206) x=2006
                                        if (x==207) x=2007
                                        else x=x})
pops.info$yr<-as.integer(pops.info$yr)

pops<-unique(pops.info$Population.Year)

pwss<-c("PWS91","PWS96","PWS07","PWS17")
tbs<-c("TB91","TB96","TB06","PWS17")
sss<-c("SS96","SS06","SS17")
y17<-pops[grep("17",pops)]

comb1<-combn(pwss, 2)
comb1<-t(comb1)
comb2<-combn(tbs, 2)
comb2<-t(comb2)
comb3<-combn(sss, 2)
comb3<-t(comb3)
comb4<-combn(y17, 2)
comb4<-t(comb4)

#PWS
Plots<-list()
sfs.pws<-data.frame()
for (i in 1: nrow(comb1)){
    pop1<-comb1[i,1]
    pop2<-comb1[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/fromVCF/2D/folded_",pop1,"_",pop2,"_maf00.sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)
    #sfs1<-vec2mat(paste0("Data/new_vcf/angsd/fromVCF/maf05/folded_",pop1,"_",pop2,"_maf05.sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)
    
    if (pops.info$yr[pops.info$Population.Year==pop1][1]>pops.info$yr[pops.info$Population.Year==pop2][1]){
        sfs1<-data.frame(t(sfs1[,1:(ncol(sfs1)-1)]))
        colnames(sfs1)<-0:(ncol(sfs1)-1)
        sfs1$count<-0:(nrow(sfs1)-1)
        p2<-pop2
        pop2<-pop1
        pop1<-p2
    }
    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs.pws<-rbind(sfs.pws, sfsm2)
}

sfs.pws$pop1<-factor(sfs.pws$pop1, levels=c("PWS91","PWS96","PWS07","PWS17"))
sfs.pws$pop2<-factor(sfs.pws$pop2, levels=c("PWS91","PWS96","PWS07","PWS17"))

#https://stackoverflow.com/questions/49689069/heatmap-with-continuous-rainbow-colours
cols <- rev(rainbow(7)[-7]) #rainbow colors for heatmap


ggplot(sfs.pws, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=log10(value)))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="log(# of alleles)", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_PWS.png", width = 10, height = 8, dpi=300)
 
 

#TB
sfs.tb<-data.frame()
for (i in 1: nrow(comb2)){
    pop1<-comb2[i,1]
    pop2<-comb2[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/fromVCF/2D/folded_",pop1,"_",pop2,"_maf00.sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)

    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs.tb<-rbind(sfs.tb, sfsm2)
}

sfs.tb$pop1<-factor(sfs.tb$pop1, levels=c("TB91","TB96","TB06","TB17"))
sfs.tb$pop2<-factor(sfs.tb$pop2, levels=c("TB91","TB96","TB06","TB17"))

ggplot(sfs.tb, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=value))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="# of alleles", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_TB.png", width = 10, height = 8, dpi=300)

#SS
sfs.ss<-data.frame()
for (i in 1: nrow(comb3)){
    pop1<-comb3[i,1]
    pop2<-comb3[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/fromVCF/2D/folded_",pop1,"_",pop2,"_maf00.sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)
    
    if (pops.info$yr[pops.info$Population.Year==pop1][1]>pops.info$yr[pops.info$Population.Year==pop2][1]){
        sfs1<-data.frame(t(sfs1[,1:(ncol(sfs1)-1)]))
        colnames(sfs1)<-0:(ncol(sfs1)-1)
        sfs1$count<-0:(nrow(sfs1)-1)
        p2<-pop2
        pop2<-pop1
        pop1<-p2
    }
    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs.ss<-rbind(sfs.ss, sfsm2)
}
sfs.ss$pop1<-factor(sfs.ss$pop1, levels=c("SS96","SS06","SS17"))
sfs.ss$pop2<-factor(sfs.ss$pop2, levels=c("SS96","SS06","SS17"))

ggplot(sfs.ss, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=value))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="# of alleles", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_SS.png", width = 10, height = 8, dpi=300)


##2017

Plots<-list()
sfs17<-data.frame()
for (i in 1: nrow(comb4)){
    pop1<-comb4[i,1]
    pop2<-comb4[i,2]
    n1<-nrow(pops.info[pops.info$Population.Year==pop1,])
    n2<-nrow(pops.info[pops.info$Population.Year==pop2,])
    
    sfs1<-vec2mat(paste0("Data/new_vcf/angsd/fromVCF/2D/folded_",pop1,"_",pop2,"_maf00.sfs"), n1=n1, n2=n2, pop1=pop1, pop2=pop2)

    #Plot first 30
    sfs2<-sfs1[1:30,1:30]
    sfs2$count<-0:(nrow(sfs2)-1)
    sfsm2<-melt(sfs2, id.vars="count")
    sfsm2$variable<-as.integer(as.character(sfsm2$variable))
    
    #zero as white (replace with NA)
    sfsm2$value[sfsm2$value==0]<-NA
    sfsm2$pop1<-pop1
    sfsm2$pop2<-pop2
    sfs17<-rbind(sfs17, sfsm2)
}

sfs17$pop1<-factor(sfs17$pop1, levels=paste(y17))
sfs17$pop2<-factor(sfs17$pop2, levels=paste(y17))

ggplot(sfs17, aes(x=count, y=variable))+
    facet_grid(pop2~pop1)+
    geom_raster(aes(fill=log(value)))+xlab('')+ylab("")+
    scale_fill_gradientn(colors=cols, name="log(# of alleles)", na.value = "white")+
    theme_minimal()
ggsave("Output/SFS/sfs_2D_2017.png", width = 20, height = 16, dpi=300)



## Plot 1DSFS from bam (Joe's downsampled version) ##

##### Downsampled SFS from Joe's 
pops<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops$Population.Year)

#Use same Y axis
sfs1D<-data.frame()
for (i in 1:length(pops)){
    sfs <- scan(paste("Data/sfs/downsample_41/unfolded/",pops[i],"_minQ20_minMQ30_unfolded.sfs",sep = ""))
    sfs1 <- data.frame(ac=sfs)
    sfs1$count<-0:(nrow(sfs1)-1)
    #remove the invariable sites
    sfs1<-sfs1[-c(1,nrow(sfs1)),]
    sfs1$pop<-pops[i]
    sfs1D<-rbind(sfs1D, sfs1)
}
    
sfs1D$pop<-factor(sfs1D$pop, levels=c("PWS91","PWS96","PWS07","PWS17","TB91","TB96", "TB06","TB17","SS96","SS06", "SS17","BC17","WA17","CA17"))
ggplot(data=sfs1D, aes(x=count, y=ac))+
    facet_wrap(~pop, ncol=4)+
        geom_bar(stat="identity", color="gray")+xlab("Frequency bin")+ ylab("Number of alleles")+
        theme_classic()+
        scale_y_continuous(labels=scales::comma)+
        theme(strip.background = element_rect(
            color="black", fill="gray80", size=0.5, linetype="solid"))
ggsave("Output/SFS/1DSFS_all_downsampled.png", width = 12, height = 8, dpi=300)

### SFS from VCF maf00 files

#Use same Y axis
sfs1D<-data.frame()
for (i in 1:length(pops)){
    sfs <- scan(paste("Data/new_vcf/angsd/fromVCF/",pops[i],"_maf00.unfolded.sfs",sep = ""))
    sfs1 <- data.frame(ac=sfs)
    sfs1$count<-0:(nrow(sfs1)-1)
    #remove the invariable sites
    sfs1<-sfs1[-c(1,nrow(sfs1)),]
    sfs1$pop<-pops[i]
    sfs1D<-rbind(sfs1D, sfs1)
}

sfs1D$pop<-factor(sfs1D$pop, levels=c("PWS91","PWS96","PWS07","PWS17","TB91","TB96", "TB06","TB17","SS96","SS06", "SS17","BC17","WA17","CA17"))
ggplot(data=sfs1D, aes(x=count, y=ac))+
    facet_wrap(~pop, ncol=4)+
    geom_bar(stat="identity", color="gray")+xlab("Frequency bin")+ ylab("Number of alleles")+
    theme_classic()+
    scale_y_continuous(labels=scales::comma)+
    theme(strip.background = element_rect(
        color="black", fill="gray80", size=0.5, linetype="solid"))
ggsave("Output/SFS/1DSFS_all_fromVCF.png", width = 12, height = 8, dpi=300)




## Compare SFS from all samples vs. downsampled ##

# Compare Joe's downsampled SFS (capped to 100 read depth and 100M sites) vs. SFS based on all individuals and sites (capped 500 depths, >590M sites)

sfsj<-scan("Data/sfs/downsample_41/unfolded/PWS07_minQ20_minMQ30_unfolded.sfs")
sfsk<-scan("Data/new_vcf/angsd/fromBam/unfolded/PWS07_unfolded.sfs")
#plot variable sites
barplot(sfsj[-c(1,length(sfsj))])
barplot(sfsk[-c(1,length(sfsk))])

sfs<-data.frame(joe=sfsj, kaho=sfsk[1:length(sfsj)])
#standardized to Joe's'
p<-sfsj[1]/sfsk[1] #ratio of invariable sites

#0.9612
sfsj[2]/sfsk[2] #ratio of a singleton
#0.45468
sfs$k_std<-sfs$kaho*p
#remove the invariable sites
sf<-sfs[-c(1,nrow(sfs)),]
sf$count<-1:nrow(sf)

sfm<-melt(sf, id.vars = "count")
ggplot(sfm[sfm$count>=1&sfm$count<=40,], aes(x=count, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width = 0.5))

sfm2<-sfm[sfm$variable!="kaho",]
ggplot(sfm2[sfm2$count>=1&sfm2$count<=40,], aes(x=count, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))

sfsj[3]/sfsk[3] #ratio of a doubleton
#0.6204
#Shape is slightly different.




#compare folded SFS vs. unfolded SFS results for PWS07 vs. PWS17

fold<-read.delim("Data/new_vcf/angsd/Fst/fst_PWS07_PWS17_folded_50kWindow")
unf<-read.delim("Data/new_vcf/angsd/Fst/fst_PWS07_PWS17_50kWindow")

mean(fold$Nsites) #[1] 0.009639717
mean(unf$Nsites) # [1] 0.008852109
median(fold$Nsites) #0.007356
median(unf$Nsites) #0.007189

fold<-read.delim("Data/new_vcf/angsd/Fst/folded_fst_PWS07_PWS91_50kWindow")
unf<-read.delim("Data/new_vcf/angsd/Fst/fst_PWS07_PWS91_50kWindow")

mean(fold$Nsites) # 0.008773152
mean(unf$Nsites)  # 0.008241249

## folded produces higher mean Fst values ##






#### annotation of potentially selected sites #####
# top 1% outliers of PBSs (Fst/PBS from ANGSD were based on windows)
pb.out<-pbsm[order(abs(pbsm$value), decreasing = T),] #210222 windows
pb.out<-pb.out[1:1913,]

pb.out<-merge(pb.out, pbs[,c("loc","region","midPos")], by="loc", all.x=T)
table(pb.out$variable)
#PBS96 PBS07 PBS17 
#268  1084   561  

length(unique(pb.out$loc))
#1913/2102

p17<-pb.out[pb.out$variable=="PBS17",]
p07<-pb.out[pb.out$variable=="PBS07",]
p96<-pb.out[pb.out$variable=="PBS96",]

sum<-data.frame(table(pb.out$variable,pb.out$chr))
#       1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
#PBS96 15  0 15 14  4 20 20 12  7 14  8 28  5 15 15  9  4 18 15  7  8 11 26 11  1  0
#PBS07 66 38 58 54 27 59 45 39 19 87 45 65 35 57 48 59 46 47 61 38 58 28 40 35 11 12
#PBS17 23 28 31 25 26 32 21 19  9 34 28 38 24 33  5 42 32 20 37 27 15  8 40  8  8 10

chr<-aggregate(sum$Freq, by=list(sum$Var2), sum)
chr[order(chr$x),]

#extract the regions and find the annotations
#create the bed file
pops<-c("PWS96","PWS07","PWS17")
files<-c("p96","p07","p17")
for (i in 1: length(pops)){
    loc<-data.frame()
    n=1
    pbs<-get(files[i])
    while (n <=nrow(pbs)){
        x<-pbs$loc[n]
        
        if (!is.na(pbs$loc[n+1])& pbs$loc[n+1]!=(x+1)) {
            newrow=c(pbs$loc[n], pbs$chr[n],pbs$midPos[n]-25000, pbs$midPos[n]+25000 )
            loc<-rbind(loc, newrow)
            n=n+1
        }
        else if (is.na(pbs$loc[n+1])) n=n+1
        else if (pbs$loc[n+1]==(x+1)){
            k=0
            while (!is.na(pbs$loc[n+k]) & pbs$loc[n+k]==(x+k)) k=k+1
            newrow=c(pbs$loc[n], pbs$ch[n],pbs$midPos[n]-25000, pbs$midPos[n+k-1]+25000)
            loc<-rbind(loc, newrow)
            n=n+k
        }
    }        
    
    loc<-loc[,-1]
    #convert the numbers to non-scientific
    loc[,2]<-as.integer(loc[,2])
    loc[,3]<-as.integer(loc[,3])
    write.table(loc, paste0("Output/fst_pbs/PWS/", pops[i],"_OutlierPBS_loci.bed"), quote=F, sep="\t",row.names=F, col.names = F)
}





# Create a new vcf file containing the loci in p07_loc.
 ## at terminal
vcftools --gzvcf Data/vcfs/population_PWS07_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS07_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS07_pbsoutlier --recode --keep-INFO-all
vcftools --gzvcf Data/vcfs/population_PWS17_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS17_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS17_pbsoutlier --recode --keep-INFO-all
vcftools --gzvcf Data/vcfs/population_PWS96_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz --bed Output/fst_pbs/PWS/PWS96_OutlierPBS_loci.bed --out Output/fst_pbs/PWS/PWS96_pbsoutlier --recode --keep-INFO-all

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS07_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS07 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS07_PBS_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS17_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS17 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS17_PBS_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS96_pbsoutlier.recode.vcf -stats ~/Projects/PacHerring/Output/fst_pbs/PWS/PWS96 > ~/Projects/PacHerring/Output/fst_pbs/PWS/Anno.PWS96_PBS_outlier.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS07_PBS_outlier.vcf > PWS07_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS17_PBS_outlier.vcf > PWS17_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Anno.PWS96_PBS_outlier.vcf > PWS96_annotation


#snpEff gene results
pops<-c("PWS96","PWS07","PWS17")
for (f in 1:3){
    #read the annotation info
    df<-read.table(paste0("Output/fst_pbs/PWS/",pops[f],"_annotation"), header = F)
    annotations<-data.frame()
    for (i in 1: nrow(df)){
        anns<-unlist(strsplit(df$V4[i], "\\|"))
        anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
        annotations<-rbind(annotations, anns)
    }     

    colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
    Ano<-cbind(df[,1:3], annotations)
    colnames(Ano)[1:3]<-c("chr","pos","AF")
    #remove the duplicated annotations for deeper digging
    remove<-!duplicated(annotations)
    Ano2<-Ano[remove,]
    write.csv(Ano2, paste0("Output/fst_pbs/PWS/", pops[f],"_highPBS_genelist.csv"))
    
    geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
    geneids<-unique(geneids)
    geneids<-geneids[nchar(geneids)<=18]
    geneids<-unique(geneids)
    
    
    sink(paste0("Output/fst_pbs/PWS/", pops[f],"_geneid_list.txt"))
    cat(paste0(geneids,"; "))
    sink(NULL)
}



annotations<-data.frame()
for (i in 1: nrow(df)){
    anns<-unlist(strsplit(df$V4[i], "\\|"))
    anns<-anns[c(2,3,4,5,8,17,18,19,20,23)]
    annotations<-rbind(annotations, anns)
}      

colnames(annotations)<-c("Annotation","Putative_impact","Gene_name", "Gene_ID", "Transcript_biotype","Annotation2","Putative_impact2","Gene_name2", "Gene_ID2", "Transcript_biotype2")
Ano<-cbind(df[,1:3], annotations)
colnames(Ano)[1:3]<-c("chr","pos","AF")


#remove the duplicated annotations for deeper digging
remove<-!duplicated(annotations)
Ano2<-Ano[remove,]
write.csv(Ano2, "Output/fst_pbs/PWS/PWS07_highPBS_genelist.csv")            

geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
geneids<-unique(geneids)
geneids<-geneids[nchar(geneids)<=18]
geneids<-unique(geneids)


sink("Output/fst_pbs/PWS/PWS07_geneid_list.txt")
cat(paste0(geneids,"; "))
sink(NULL)



#####  from snpEff results ###
# Read the gene file

p07<-read.delim("Output/fst_pbs/PWS/PWS07.genes.txt",skip=1 )

#unique genes

#remove the no-annotated genes
p07a<-p07[grep("ENSCHA", p07$X.GeneName,invert = T),]
length(unique(p07a$X.GeneName)) #570 genes

write.table(unique(p07a$X.GeneName), "Output/fst_pbs/PWS/pws07_unique_genes.txt", quote = F, row.names = F)


## 
Ano2[Ano2$Gene_name=="FHIT",]
Ano2[Ano2$Gene_name=="styx",]

findGene<-function(gene){
    df<-Ano2[Ano2$Gene_name==gene,]
    
    return(df)
}

findGene("DDAH1")
findGene("vaspa")
findGene("hoxa2b")
Ano2[Ano2$Gene_name2=="hoxa2b",]
