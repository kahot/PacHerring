
## CDA
#library("car")
library("candisc")
#library("heplots")
library(vcfR)
library(ggplot2)

pop_info<-read.csv("Data/Sample_metadata_892pops.csv")
pop_info<-pop_info[,c("Sample","pop","Year.Collected")]
colnames(pop_info)[3]<-"year"

vcf <- read.vcfR("Data/vcfs/chr8_sex.vcf")

vcf <- read.vcfR(paste0("Data/vcfs/ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05.vcf.gz"))
gt <- extract.gt(vcf)

gt<-apply(gt,2,function(x) gsub("\\|","\\/",x))
gt<-apply(gt,2, function(y) sapply(y, function(x) {if (is.na(x)) NA
                            else {if (x=="1/1") 2
                                else if (x=="0/0") 0
                                else if (x=="0/1") 1
                            }}))

#gt2<-gt
gt<-data.frame(t(gt))
#Sample<-rownames(gt)
pop<-substr(rownames(gt),6,10)
gt<-cbind(Sample, pop, gt)

mod1<- lm(pop~., data=gt,  na.action=na.exclude)
#Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
#contrasts can be applied only to factors with 2 or more levels

#which sites have only two levels
gt %>% dplyr::mutate_all(as.factor) %>% str

ans<-apply(gt, 2, function(x) {y=unique(x)
                    if (length(y)<3) "n"
                    else NA})
#none of them

#how many NAs per locus?
ans<-apply(gt, 2, function(x) length(x[is.na(x)]))
plot(ans)
#about 300 ish


