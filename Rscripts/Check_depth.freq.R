# Check read depth and frequencies since temporal COVs show higher variance for MD7000 
source("Rscripts/BaseScripts.R")
library(gridExtra)
library(vcfR)
cols<-qualitative_hcl(8, palette="dark3")

#Compare MD2000 vs. MD7000 freq distribution
vcf_j<- read.vcfR("Data/vcfs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz")
vcf_k<-read.vcfR("Data/new_vcf/PH_DP600_7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz")
vcf_k<-read.vcfR("Data/new_vcf/PH_")


#dna<-ape::read.dna()
#gff<-

#extract depths(allelic depths) 
adj <- extract.gt(vcf_j, element = "AD")
adj1 <- masplit(adj, record = 1) #major allele
adj2 <- masplit(adj, record = 2) #minor allele
adj_1sum<-rowSums(adj1, na.rm=T)
adj_2sum<-rowSums(adj2, na.rm=T)

freqj <- adj_2sum/(adj_1sum+adj_2sum)
FJ<-data.frame(snp=rownames(adj))
FJ$freq<-freqj
FJ$chr<-gsub("_","",substr(FJ$snp,1,5))
FJ$pos<-as.integer(gsub("chr\\d+\\_","", FJ$snp))
write.csv(FJ,"Output/VCF/freq_Joe.csv")

#Same for MD7000 file
adk <- extract.gt(vcf_k, element = "AD")
adk1 <- masplit(adk, record = 1)
adk2 <- masplit(adk, record = 2)

adk_1sum<-rowSums(adk1, na.rm=T)
adk_2sum<-rowSums(adk2, na.rm=T)

freqk <- adk_2sum/(adk_1sum+adk_2sum)

FK<-data.frame(snp=rownames(adk))
FK$freq<-freqk
FK$chr<-gsub("_","",substr(FK$snp,1,5))
FK$pos<-as.integer(gsub("chr\\d+\\_","", FK$snp))
write.csv(FK,"Output/VCF/freq_newvcf_MD7000.csv")


#Combine the data
FJ$vcf<-"MD2000"
FK$vcf<-"MD7000"
F_comb<-rbind(FK, FJ[,c(1:4,7)])
#Freq<-merge(FK, FJ[,c(2,3,4,5,6)], by=c("chr","pos"),all=T)

ggplot(F_comb, aes(x=pos, y=freq, color=vcf))+
    geom_point(size=0.5, alpha=0.4)+
    facet_wrap(~chr, ncol=3)+
    theme_bw()+
    theme(axis.text.x=element_blank(), legend.title = element_blank())+xlab("")+
    scale_color_manual(values=cols[c(1,6)])
ggsave("Output/VCF/Loci_freq_comparison.pdf", width=20, height=16)



### seq depth filtering
dpk <- extract.gt(vcf_k, element = "DP")
#dpk<-matrix(as.integer(dpk),    # Convert to numeric matrix
#       ncol = ncol(dpk))

dpk<-mapply(dpk, as.integer)

dpk_sum<-rowSums(dpk, na.rm=T)
max(dpk_sum)
hist(dpk_sum)


depth<-data.frame(snp=rownames(adk))
depth$depth<-dpk_sum
hist(depth$depth)
nrow(depth[depth$depth<=5000,]) #324751
md5000_pos<-depth[depth$depth<=5000,]

plot(depth$depth, FK$freq, pch=".", xlab="Read depth", ylab="Freq",col="steelblue")

cor.test(md5000$freq, md5000_pos$depth)

snpremove<-FK$snp[!(FK$snp %in% md5000_pos$snp)] #5729

remove<-data.frame(chr=gsub("_","",substr(snpremove,1,5)), pos=as.integer(gsub("chr\\d+\\_","", snpremove)))
write.table(remove, "Output/VCF/snp.remove_MD5000.txt", sep="\t", row.names = F, col.names = F, quote = F)


### MD4000: removing at max depth at 4000
nrow(depth[depth$depth<=4000,]) #317724
md4000_pos<-depth[depth$depth<=4000,]

snpremove<-FK$snp[!(FK$snp %in% md4000_pos$snp)] #12758

remove<-data.frame(chr=gsub("_","",substr(snpremove,1,5)), pos=as.integer(gsub("chr\\d+\\_","", snpremove)))
write.table(remove, "Data/new_vcf/snp.remove_MD4000.txt", sep="\t", row.names = F, col.names = F, quote = F)

## remove loci from a vcf file using bcftools
bcftools view -T ^snp.remove_MD4000.txt PH_DP600_7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz > PH_DP600_4000_minQ20_minMQ30_NS0.5_maf05.vcf
#check the new vcf file has the right number of loci
bcftools query -f '%POS\n' PH_DP600_4000_minQ20_minMQ30_NS0.5_maf05.vcf > positions4000.txt
wc -l positions4000.txt #317724
bgzip PH_DP600_4000_minQ20_minMQ30_NS0.5_maf05.vcf
bcftools view -T ^snp.remove_MD4000.txt Timeseries_maf05.vcf.gz > Timeseries_maf05_MD4000.vcf
bgzip Timeseries_maf05_MD4000.vcf




### MD3000: removing at max depth at 4000
nrow(depth[depth$depth<=3000,]) #298926
md3000_pos<-depth[depth$depth<=3000,]

snpremove<-FK$snp[!(FK$snp %in% md3000_pos$snp)] #31556

remove<-data.frame(chr=gsub("_","",substr(snpremove,1,5)), pos=as.integer(gsub("chr\\d+\\_","", snpremove)))
write.table(remove, "Data/new_vcf/snp.remove_MD3000.txt", sep="\t", row.names = F, col.names = F, quote = F)

## remove loci from a vcf file using bcftools
bcftools view -T ^snp.remove_MD3000.txt PH_DP600_7000_minQ20_minMQ30_NS0.5_maf05.vcf.gz > PH_DP600_3000_minQ20_minMQ30_NS0.5_maf05.vcf
#check the new vcf file has the right number of loci
bcftools query -f '%POS\n' PH_DP600_3000_minQ20_minMQ30_NS0.5_maf05.vcf > positions3000.txt
wc -l positions3000.txt #298926
bgzip PH_DP600_3000_minQ20_minMQ30_NS0.5_maf05.vcf

bcftools view -T ^snp.remove_MD3000.txt Timeseries_maf05.vcf.gz > Timeseries_maf05_MD3000.vcf
bgzip Timeseries_maf05_MD3000.vcf

#calculate freq by population
pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pop_info$Population.Year)



##Examine pws07 to see why bcftools and vcftools produce different MAF 
vcf<-read.vcfR("Data/new_vcf/population/PWS07_maf05.vcf.gz")
#Same for MD7000 file
adk <- extract.gt(vcf, element = "AD")
gt<-extract.gt(vcf, element = "GT")
gt[1:2,]
adk1 <- masplit(adk, record = 1)
adk2 <- masplit(adk, record = 2) #minor allele

AF<-data.frame(snp=colnames(adk))
AF$majAD<-adk1[1,]
AF$minorAD<-adk2[1,]
AF$majCount<-0
AF$majCount[AF$majAD>0&AF$minorAD==0]<-2
AF$majCount[AF$majAD>0&AF$minorAD>0]<-1
AF$majCount[AF$majAD==0&AF$minorAD==0]<-NA

nrow(AF[!is.na(AF$majCount),]) #38
sum(AF$majCount, na.rm=T) #68
1-68/38/2
#0.1052632
#Can't use this info to get actual AF for pop
#use vcftools