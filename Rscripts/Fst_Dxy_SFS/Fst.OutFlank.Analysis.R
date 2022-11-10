# Fst Outlier analysis
#https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html
#BiocManager::install("qvalue")
#if (!("OutFLANK" %in% installed.packages())){install_github("whitlock/OutFLANK")}

library(OutFLANK)  # outflank package
library(vcfR)
library(adegenet)
library(dartR)
library(qvalue)
library(stringr)
source("Rscripts/data_subset.r")
# Italic labels
fstlab = expression(italic("F")[ST])
hetlab = expression(italic("H")[e])

#read the VCF file (all populations) 
vcf_all<-read.vcfR("Data/plink/prune/PH_pruned_all_maf05.vcf.gz")

geno_all <- extract.gt(vcf_all) # Character matrix containing the genotypes
position_all <- getPOS(vcf_all) # Positions in bp
chromosome_all <- getCHROM(vcf_all) # Chromosome information

G_all <- matrix(NA, nrow = nrow(geno_all), ncol = ncol(geno_all))

G_all[geno_all %in% c("0/0", "0|0")] <- 0
G_all[geno_all  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_all[geno_all %in% c("1/1", "1|1")] <- 2
G_all[is.na(geno_all)]<-9

G_all<-data.frame(G_all)
colnames(G_all)<-colnames(geno_all)
rownames(G_all)<-rownames(geno_all)
write.csv(G_all, gzfile("Output/Fst/MD7000/GenotypeAllPops_matrixG.csv.gz"))

#Read in the saved genotype matrix
#G_all<-read.csv("Output/Fst/MD7000/GenotypeAllPops_matrixG.csv")


#extract population IDs and position IDs.
popIDs_all<-substr(colnames(geno_all), 6,10)
posIDs_all<-rownames(geno_all)


#Subset PWS from G_all matrix
G_pws<-G_all[,grep("PWS", popIDs_all)]

popsPWS<-popIDs_all[grep("PWS", popIDs_all)]


#Calculate FST on the PWS non-pruned data
my_fst <- MakeDiploidFSTMat(t(G_pws), locusNames =posIDs_all, popNames =popsPWS)

#Data checks: Heterozygosity vs. FST
plot(my_fst$He, my_fst$FST, pch=".")

#Data checks: FST vs. FSTNoCorr
plot(my_fst$FST, my_fst$FSTNoCorr, pch=".")
abline(0,1, col="red")

## remove the pruned sites and run OutFLANK
#pruned data
vcf<-read.vcfR("Data/plink/prune/PH_pruned_all_maf05.vcf.gz")
#VCF to Outflank format (G matrix)
geno <- extract.gt(vcf) # Character matrix containing the genotypes
position <- getPOS(vcf) # Positions in bp
chromosome <- getCHROM(vcf) # Chromosome information

G_pru <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G_pru[geno %in% c("0/0", "0|0")] <- 0
G_pru[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_pru[geno %in% c("1/1", "1|1")] <- 2
G_pru[is.na(geno)]<-9

G_pru<-data.frame(G_pru)
colnames(G_pru)<-colnames(geno)
rownames(G_pru)<-rownames(geno)
#write.csv(G_all, gzfile("Output/Fst/MD7000/Pruned_AllPops_GenotypeAllPops_matrixG.csv.gz"))

#Subset PWS from pruned G_pru matrix
G_pru_pws<-G_pru[,grep("PWS", popIDs_all)]
#write.csv(G_pru_pws, gzfile("Output/Fst/MD7000/Pruned_PWS_Genotype_matrixG.csv"))


#Subset PWS pops
G91<-G_pru_pws[,grep("PWS91", colnames(G_pru_pws))]
G96<-G_pru_pws[,grep("PWS96", colnames(G_pru_pws))]
G07<-G_pru_pws[,grep("PWS07", colnames(G_pru_pws))]
G17<-G_pru_pws[,grep("PWS17", colnames(G_pru_pws))]

#write.csv(G91,"Output/Fst/MD7000/PWS91_genotypes.csv")
#write.csv(G96,"Output/Fst/MD7000/PWS96_genotypes.csv")
#write.csv(G07,"Output/Fst/MD7000/PWS07_genotypes.csv")
#write.csv(G17,"Output/Fst/MD7000/PWS17_genotypes.csv")


#find the pruned-positions
posIDs<-rownames(geno)
#reformat posIDS to match posIDs_all
posIDs<-gsub('\\[ph\\]','',posIDs)
posIDs<-gsub('\\:','\\_',posIDs)

pruned_pos<-which(posIDs_all %in% posIDs)
#write.csv(pruned_pos,"Output/Fst/MD7000/pruned_positions.csv", row.names = F)
pruned_fst<-my_fst[pruned_pos,]

#Next, run the OutFLANK() function to estimate the parameters on the neutral FST distribution.
out_trim <- OutFLANK(pruned_fst, NumberOfSamples=3, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)

#Check the fit and make sure it looks good, especially in the right tail:
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           TRUE, RightZoomFraction = 0.15, titletext = NULL)


# check the P-value histogram:
hist(out_trim$results$pvaluesRightTail)


#Using estimated neutral mean FST and df to calculate P-values for all loci
#note that it is important to run this code with the uncorrected FSTs (FSTNoCorr) and the uncorrected mean FST (FSTNoCorrbar).

P1_all <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1_all)
tail(P1_all)



P1_all$outlier<-"N"
P1_all$outlier[P1_all$OutlierFlag==TRUE]<-"Y"

ggplot(data = P1_all,aes(x = He, y = FST, colour = outlier))+
    geom_point(size=0.7, alpha)+
    scale_colour_manual(values = c("#8080804D","blue"), labels = c("Neutral SNP","Outlier SNP"), na.value = "white")+
    ggtitle("PWS all years")+
    theme_bw()+
    xlab(hetlab)+
    ylab(fstlab)+
    theme(legend.title = element_blank(),legend.text = element_text(size=11),
          plot.title = element_text(size = 13))
ggsave(paste0("Output/Fst/MD7000/Fst.vs.He_outlier_PWS_allyears.pdf"), width = 8, height=4)
#No outlier identified


# Manhattan Plot
# *For publication, we want to show the accurate estimate of FST, not the uncorrected estimate. Remember to exclude those low H loci!

P1_all$Index<-1:nrow(P1_all)
P1_all$chr<-substr(P1_all$LocusName, 1,5)
P1_all$chr<-gsub("_","",P1_all$chr)

P1_all$color<-"odd"
P1_all$color[P1_all$chr %in% evens]<-"even"
ggplot()+
    geom_point(data=P1_all[P1_all$He>0.1,], aes(x=Index, FST, color=color),size=0.5, alpha=0.4)+
    scale_color_manual(values=c("#808080","#C0C0C0"), guide="none")+
    geom_point(data=P1_all[P1_all$outlier=="Y",], aes(x=Index, FST),color="red", size=0.7)+
    xlab("Chromosome poistion")+ylab("Fst")+
    ggtitle("PWS")+
    theme_light()+theme(axis.text.x = element_blank())
ggsave(paste0("Output/Fst/MD7000/Fst_outlier_manhattan_PWS_allYears.pdf"), width = 8, height=4)


hist(P1_all$pvaluesRightTail) #looks good


#####  pairwise comparison   ###
#subset the population and run pairwise analysis

#Subset PWS from G_all matrix
pwapops<-c("PWS91", "PWS96","PWS07","PWS17")
comb<-combn(pwapops, 2)
comb<-t(comb)
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

for (i in 1:nrow(comb) ){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    g1<-G_all[,grep(paste0(pop1,"|",pop2), popIDs_all)]
    subpopIds<-popIDs_all[grep(paste0(pop1,"|",pop2), popIDs_all)]
    
    sub_fst <- MakeDiploidFSTMat(t(g1), locusNames =posIDs_all, popNames =subpopIds)
    
    #Next, run the OutFLANK() function to estimate the parameters on the neutral FST distribution.
    sub_out <- OutFLANK(sub_fst[pruned_pos,], NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1)
    
    #Check the fit and make sure it looks good, especially in the right tail:
    OutFLANKResultsPlotter(sub_out, withOutliers = TRUE,
                           NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                               FALSE, RightZoomFraction = 0.05, titletext = NULL)
    ## Zoom in on right tail
    OutFLANKResultsPlotter(sub_out , withOutliers = TRUE,
                           NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                               TRUE, RightZoomFraction = 0.15, titletext = NULL)
    
    # check the P-value histogram:
    hist(sub_out$results$pvaluesRightTail)
    #Skewed
    
    #Using estimated neutral mean FST and df to calculate P-values for all loci
    #ote that it is important to run this code with the uncorrected FSTs (FSTNoCorr) and the uncorrected mean FST (FSTNoCorrbar).
    
    P1 <- pOutlierFinderChiSqNoCorr(sub_fst, Fstbar = sub_out$FSTNoCorrbar, 
                                    dfInferred = sub_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
    
    P1$outlier<-"N"
    P1$outlier[P1$OutlierFlag==TRUE]<-"Y"
    
    ggplot(data = P1,aes(x = He, y = FST, colour = outlier))+
        geom_point(size=0.7, alpha=0.7)+
        scale_colour_manual(values = c("#8080804D","blue"), labels = c("Neutral SNP","Outlier SNP"), na.value = "white")+
        ggtitle(paste0(pop1," vs. ", pop2))+
        theme_bw()+
        xlab(hetlab)+
        ylab(fstlab)+
        theme(legend.title = element_blank(),legend.text = element_text(size=11),
              plot.title = element_text(size = 13))
    ggsave(paste0("Output/Fst/MD7000/Fst.vs.He_outlier",pop1,".",pop2,".pdf"), width = 8, height=4)
    
    #manhattan plot
    
    #Remove NA 
    P1<-P1[!is.na(P1$pvalues),]
    P1$Index<-1:nrow(P1)
    P1$chr<-substr(P1$LocusName, 1,5)
    P1$chr<-gsub("_","",P1$chr)
    
    colors<-rep(c("#808080","#C0C0C0"), times=13)
    ch<-unique(P1$chr)
    
    #count the number of sites per chromosomes
    poss<-data.frame(chr=ch)
    k=1
    for (i in 1:nrow(poss)){
        df<-P1[P1$chr==poss$chr[i],]
        poss$start[i]<-k
        poss$end[i]<-k+nrow(df)-1
        k=k+nrow(df)
    }
    
    poss$x<-poss$start+(poss$end-poss$start)/2
    P1$chr<-factor(P1$chr, levels=paste0("chr", 1:26))
    
    ggplot()+
        geom_point(data=P1[P1$He>0.1,], aes(x=Index, FST, color=chr),size=0.5, alpha=0.4)+
        scale_color_manual(values=colors, guide="none")+
        geom_point(data=P1[P1$outlier=="Y",], aes(x=Index, FST),color="red", size=0.7)+
        theme_light()+
        xlab("Chromosome poistion")+ylab("Fst")+
        scale_x_continuous(name="Chromosome", breaks=poss$x, labels=gsub("chr",'',poss$chr))+
        ggtitle(paste0(pop1," vs. ", pop2))
        
    ggsave(paste0("Output/Fst/MD7000/Fst_outlier_manhattan_",pop1,".",pop2,".png"), width = 8, height=4, dpi=300)
    
    #Output outliers
    P0<-P1[P1$outlier=="Y",]
    write.csv(P0, paste0("Output/Fst/MD7000/Outflank_outliers_",pop1,"vs",pop2,".csv"))
}


#Find how many outliers identified above are common across years
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],"vs",x[2]))
#[1] "PWS91vsPWS96" "PWS91vsPWS07" "PWS91vsPWS17" "PWS96vsPWS07" "PWS96vsPWS17" "PWS07vsPWS17"

#
nloci<-data.frame(pops=comb_pairs)
for(i in 1:6){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    df<-read.csv(paste0("Output/Fst/MD7000/Outflank_outliers_",pop1,"vs",pop2,".csv"))
    assign(comb_pairs[i], df$LocusName)
    nloci$no_outliers[i]<-nrow(df)
}

#          pops no_outliers
#1 PWS91vsPWS96         163
#2 PWS91vsPWS07          23
#3 PWS91vsPWS17          63
#4 PWS96vsPWS07         397
#5 PWS96vsPWS17         111
#6 PWS07vsPWS17          36

intersect(PWS91vsPWS96 , PWS96vsPWS07) #11
#[1] "chr1_29483547"  "chr7_21759724"  "chr10_13645866" "chr16_11053073" "chr17_10252353" "chr18_6722422"  "chr21_4513704" 
#[8] "chr21_20234426" "chr22_2379602"  "chr22_14380718" "chr26_5587775" 

intersect(PWS96vsPWS07,PWS07vsPWS17) #13
#[1] "chr1_14245753"  "chr3_10022923"  "chr4_30307656"  "chr6_29460341"  "chr8_3558924"   "chr12_26156986" "chr14_9222068" 
#[8] "chr17_1823419"  "chr17_12404988" "chr19_6807182"  "chr20_3165717"  "chr24_6128637"  "chr24_6128638"

intersect(PWS91vsPWS96,PWS07vsPWS17) #0

intersect(PWS91vsPWS96,PWS91vsPWS07 ) #1
#[1] "chr18_15568048"

intersect(PWS91vsPWS96,PWS91vsPWS17 ) #4
#[1] "chr7_7588694"   "chr7_8116563"   "chr9_5668066"   "chr11_15308218"


intersect(PWS96vsPWS07,PWS96vsPWS17) #7
#[1] "chr2_8811947"   "chr2_28913807"  "chr11_22284131" "chr12_28623155" "chr14_15982997" "chr14_15983007" "chr26_5587775" 

intersect(PWS91vsPWS07,PWS96vsPWS07) #10
# [1] "chr1_2060945"   "chr1_2159437"   "chr3_10022923"  "chr7_17491708"  "chr9_13969039"  "chr13_3549924"  "chr13_19597734"
#[8] "chr17_18383869" "chr23_23413220" "chr24_2367879" 

Reduce(intersect, list(PWS91vsPWS96 , PWS96vsPWS07,PWS07vsPWS17)) #0



#### #### #### #### #### #### #### #### #### #### #### 
#### SS populations

#Subset SS from G_all matrix
G_ss<-G_all[,grep("SS", popIDs_all)]

popsSS<-popIDs_all[grep("SS", popIDs_all)]

#Calculate FST on the PWS non-pruned data
my_fst <- MakeDiploidFSTMat(t(G_ss), locusNames =posIDs_all, popNames =popsSS)

#pruned data
pruned_fst<-my_fst[pruned_pos,]

#Next, run the OutFLANK() function to estimate the parameters on the neutral FST distribution.
out_trim <- OutFLANK(pruned_fst, NumberOfSamples=3, qthreshold = 0.05, Hmin = 0.1)

#Check the fit and make sure it looks good, especially in the right tail:
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           TRUE, RightZoomFraction = 0.15, titletext = NULL)

# check the P-value histogram:
hist(out_trim$results$pvaluesRightTail)


#Using estimated neutral mean FST and df to calculate P-values for all loci

P1_all <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                    dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)

P1_all$outlier<-"N"
P1_all$outlier[P1_all$OutlierFlag==TRUE]<-"Y"

ggplot(data = P1_all,aes(x = He, y = FST, colour = outlier))+
    geom_point(size=0.7)+
    scale_colour_manual(values = c("#8080804D","blue"), labels = c("Neutral SNP","Outlier SNP"), na.value = "white")+
    ggtitle("PWS all years")+
    theme_bw()+
    xlab(hetlab)+
    ylab(fstlab)+
    theme(legend.title = element_blank(),legend.text = element_text(size=11),
          plot.title = element_text(size = 13))
ggsave(paste0("Output/Fst/MD7000/Fst.vs.He_outlier_SS_allyears.png"), width = 8, height=4, dpi=300)

# Manhattan Plot
# *For publication, we want to show the accurate estimate of FST, not the uncorrected estimate. Remember to exclude those low H loci!

P1<-P1_all
P1<-P1[!is.na(P1$pvalues),]
P1$Index<-1:nrow(P1)
P1$chr<-substr(P1$LocusName, 1,5)
P1$chr<-gsub("_","",P1$chr)

colors<-rep(c("#808080","#C0C0C0"), times=13)
ch<-unique(P1$chr)

#count the number of sites per chromosomes
poss<-data.frame(chr=ch)
k=1
for (i in 1:nrow(poss)){
    df<-P1[P1$chr==poss$chr[i],]
    poss$start[i]<-k
    poss$end[i]<-k+nrow(df)-1
    k=k+nrow(df)
}
poss$x<-poss$start+(poss$end-poss$start)/2

ggplot()+
    geom_point(data=P1[P1$He>0.1,], aes(x=Index, FST, color=chr),size=0.5, alpha=0.4)+
    scale_color_manual(values=colors, guide="none")+
    geom_point(data=P1[P1$outlier=="Y",], aes(x=Index, FST),color="red", size=0.7)+
    theme_light()+
    xlab("Chromosome poistion")+ylab("Fst")+
    scale_x_continuous(name="Chromosome", breaks=poss$x, labels=gsub("chr",'',poss$chr))+
    ggtitle("SS")
ggsave(paste0("Output/Fst/MD7000/Fst_outlier_manhattan_SS_allYears.pdf"), width = 8, height=4)

#####  pairwise comparison   ###
#subset the population and run pairwise analysis
pops<-c("SS96","SS06","SS17")
comb<-combn(pops, 2)
comb<-t(comb)
#comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

for (j in 2:nrow(comb) ){
    pop1<-comb[j,1]
    pop2<-comb[j,2]
    
    g1<-G_all[,grep(paste0(pop1,"|",pop2), popIDs_all)]
    subpopIds<-popIDs_all[grep(paste0(pop1,"|",pop2), popIDs_all)]
    
    sub_fst <- MakeDiploidFSTMat(t(g1), locusNames =posIDs_all, popNames =subpopIds)
    
    #Next, run the OutFLANK() function to estimate the parameters on the neutral FST distribution.
    sub_out <- OutFLANK(sub_fst[pruned_pos,], NumberOfSamples=2, qthreshold = 0.05, Hmin = 0.1)
    
    #Check the fit and make sure it looks good, especially in the right tail:
    #OutFLANKResultsPlotter(sub_out, withOutliers = TRUE,
    #                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
    #                           FALSE, RightZoomFraction = 0.05, titletext = NULL)
    ## Zoom in on right tail
    #OutFLANKResultsPlotter(sub_out , withOutliers = TRUE,
    #                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
    #                           TRUE, RightZoomFraction = 0.15, titletext = NULL)
    
    # check the P-value histogram:
    #hist(sub_out$results$pvaluesRightTail)
    #Skewed
    
    #Using estimated neutral mean FST and df to calculate P-values for all loci
    
    P1 <- pOutlierFinderChiSqNoCorr(sub_fst, Fstbar = sub_out$FSTNoCorrbar, 
                                    dfInferred = sub_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
    
    
    P1$outlier<-"N"
    P1$outlier[P1$OutlierFlag==TRUE]<-"Y"
    
    #ggplot(data = P1,aes(x = He, y = FST, colour = outlier))+
    #    geom_point(size=0.7, alpha=0.7)+
    #    scale_colour_manual(values = c("#8080804D","blue"), labels = c("Neutral SNP","Outlier SNP"), na.value = "white")+
    #    ggtitle(paste0(pop1," vs. ", pop2))+
    #    theme_bw()+
    #    xlab(hetlab)+
    #    ylab(fstlab)+
    #    theme(legend.title = element_blank(),legend.text = element_text(size=11),
    #          plot.title = element_text(size = 13))
    #ggsave(paste0("Output/Fst/MD7000/Fst.vs.He_outlier",pop1,".",pop2,".png"), width = 8, height=4, dpi=300)
    
    #manhattan plot
    #Remove NA 
    P1<-P1[!is.na(P1$pvalues),]
    P1$Index<-1:nrow(P1)
    P1$chr<-substr(P1$LocusName, 1,5)
    P1$chr<-gsub("_","",P1$chr)
    
    colors<-rep(c("#808080","#C0C0C0"), times=13)
    ch<-unique(P1$chr)
    
    #count the number of sites per chromosomes
    poss<-data.frame(chr=ch)
    k=1
    for (i in 1:nrow(poss)){
        df<-P1[P1$chr==poss$chr[i],]
        poss$start[i]<-k
        poss$end[i]<-k+nrow(df)-1
        k=k+nrow(df)
    }
    
    poss$x<-poss$start+(poss$end-poss$start)/2
    P1$chr<-factor(P1$chr, levels=paste0("chr", 1:26))
    
    ggplot()+
        geom_point(data=P1[P1$He>0.1,], aes(x=Index, FST, color=chr),size=0.5, alpha=0.4)+
        scale_color_manual(values=colors, guide="none")+
        geom_point(data=P1[P1$outlier=="Y",], aes(x=Index, FST),color="red", size=0.7)+
        theme_light()+
        xlab("Chromosome poistion")+ylab("Fst")+
        scale_x_continuous(name="Chromosome", breaks=poss$x, labels=gsub("chr",'',poss$chr))+
        ggtitle(paste0(pop1," vs. ", pop2))
    ggsave(paste0("Output/Fst/MD7000/Fst_outlier_manhattan_",pop1,".",pop2,".png"), width = 8, height=4, dpi=300)
    
    #Output outliers
    P0<-P1[P1$outlier=="Y",]
    write.csv(P0, paste0("Output/Fst/MD7000/Outflank_outliers_",pop1,"vs",pop2,".csv"))
}


#Find how many outliers identified above are common across years
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],"vs",x[2]))

nloci<-data.frame(pops=comb_pairs)
for(i in 1:6){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    df<-read.csv(paste0("Output/Fst/MD7000/Outflank_outliers_",pop1,"vs",pop2,".csv"))
    assign(comb_pairs[i], df$LocusName)
    nloci$no_outliers[i]<-nrow(df)
}
nloci









## Outlier loci gene ontology


# Run GO on enriched loci
comps<-c("PWS91vsPWS96","PWS07vsPWS96")

for (i in 1:lenth(comps)){
    df<-read.csv(paste0("Output/Fst/MD7000/Outflank_outlier_results_",comps[i],".csv"), row.names=1)
    df$pos<-gsub("chr\\d+\\_", "", df$LocusName)
    write.table(df[,c("chr","pos")], paste0("Output/Fst/MD7000/",comps[i],"_snps.bed"),row.names = F, col.names = F, quote = F)
    #add 25k around Snps
    loci<-df[,c("chr","pos")]
    loci$pos<-as.integer(loci$pos)
    loci$start<-loci$pos-25000
    loci$end<-loci$pos+25000
    write.table(loci[,c("chr","start","end")], paste0("Output/Fst/MD7000/",comps[i],"_25k.bed"),row.names = F, col.names = F, quote = F)
}


# Create a new vcf file containing the loci in p07_loc.
## at terminal
vcftools --gzvcf Data/new_vcf/PWS07_maf05.vcf.gz --bed Output/Fst/MD7000/PWS07vsPWS96_25k.bed --out Output/Fst/MD7000/PWS07vsPWS96_25_outlier --recode --keep-INFO-all
vcftools --gzvcf Data/new_vcf/PWS07_maf05.vcf.gz --bed Output/Fst/MD7000/PWS91vsPWS96_25k.bed --out Output/Fst/MD7000/PWS91vsPWS96_25_outlier --recode --keep-INFO-all

#Run snpEff (run from snpEff/ directory)
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/Fst/MD7000/PWS07vsPWS96_25_outlier.recode.vcf -stats ~/Projects/PacHerring/Output/Fst/MD7000/PWS07vsPWS96 > ~/Projects/PacHerring/Output/Fst/MD7000/Anno.PWS07vsPWS96_outlier.vcf
java -Xmx8g -jar snpEff.jar Ch_v2.0.2.99 ~/Projects/PacHerring/Output/Fst/MD7000/PWS91vsPWS96_25_outlier.recode.vcf -stats ~/Projects/PacHerring/Output/Fst/MD7000/PWS07vsPWS96 > ~/Projects/PacHerring/Output/Fst/MD7000/Anno.PWS91vsPWS96_outlier.vcf

#extract the annotation information
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/Fst/MD7000/Anno.PWS07vsPWS96_outlier.vcf > Output/Fst/MD7000/PWS07vsPWS96_annotation
bcftools query -f '%CHROM %POS %INFO/AF %INFO/ANN\n' Output/Fst/MD7000/Anno.PWS91vsPWS96_outlier.vcf > Output/Fst/MD7000/PWS91vsPWS96_annotation


#snpEff gene results

for (f in 1:length(comps)){
    #read the annotation info
    df<-read.table(paste0("Output/Fst/MD7000/",comps[f],"_annotation"), header = F)
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
    write.csv(Ano2, paste0("Output/Fst/MD7000/", comps[f],"_outlier_genelist.csv"))
    
    geneids<-c(Ano2$Gene_ID, Ano2$Gene_ID2)
    geneids<-unique(geneids)
    geneids<-geneids[nchar(geneids)<=18]
    geneids<-unique(geneids)
    
    
    sink(paste0("Output/Fst/MD7000/", comps[f],"_geneid_list.txt"))
    cat(paste0(geneids,"; "))
    sink(NULL)
}



