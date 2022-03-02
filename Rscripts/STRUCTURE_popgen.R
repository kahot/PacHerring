devtools::install_github("nicholasjclark/STRUCTURE.popgen")
library(STUCTURE.popgen)

qmat<-read.table("~/programs/fastS/fastStructure/PH/ph_simple.3.varQ")
allele<-read.csv(paste0("Output/LD/gt_genetics_chr1.csv"), stringsAsFactors = T, row.names = 1)
allele<-t(allele)
ind<-rownames(allele)

qmat$ind<-ind
qmat<-qmat[,c(4,1,2,3)]

fis.STRUCTURE.popgen(qmat, allele, 3,NA.symbol=NA)

#              diff         lwr           upr        p adj
#V2-V1 -0.008187695 -0.01609120 -0.0002841914 4.024820e-02
#V3-V1  0.026522730  0.01865977  0.0343856950 5.107026e-14
#V3-V2  0.034710425  0.02680715  0.0426136975 3.175238e-14


fst.STRUCTURE.popgen(qmat, allele, 3,NA.symbol=NA)
