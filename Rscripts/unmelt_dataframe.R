tmp <- data.frame(x=gl(2,3, labels=letters[24:25]),
                  y=gl(3,1,6, labels=letters[1:3]), 
                  z=c(1,2,3,3,3,2))

#  x y z
#1 x a 1
#2 x b 2
#3 x c 3
#4 y a 3
#5 y b 3
#6 y c 2

pivot_wider(tmp, names_from = y, values_from = z)
spread(tmp, y, z)
temp<-pivot_wider(ld1, names_from = y, values_from = z)

l<-ld1[1:50,c("SNP_A","SNP_B","R2")]
common<-intersect(unique(l$SNP_A),unique(l$SNP_B)) #15
onlyA<-l$SNP_A[!(l$SNP_A %in% common)]
onlyB<-l$SNP_B[!(l$SNP_B %in% common)]

addA<-data.frame(SNP_A=)

lm<-spread(l, SNP_B, R2)


#library(DataCombine)
#lm<-cbind(rep(0,times=nrow(lm)), lm)

library(tibble)

lm<-acast(l,formula = SNP_A ~ SNP_B, value.var = 'R2',drop = F, fill=0)
lm<-data.frame(lm)
colnames(lm)<-gsub("\\.","\\:", colnames(lm))

lm<-lm[order(rownames(lm)),]
lm<-lm[,order(colnames(lm))]


for (i in 1:nrow(lm)){
#for (i in 1:11){
    rown<-rownames(lm)[i]
    coln<-colnames(lm)[i]
    if(rown==coln) next
    if(rown!=coln) {
        if (!(rown %in% colnames(lm))){
            lm<-add_column(lm, newcol=rep(0,times=nrow(lm)), .before=i)
            colnames(lm)[colnames(lm)=="newcol"]<-rown
        }
    }
}

lmcopy<-lm
lm<-lmcopy

lm<-lm[order(rownames(lm)),]
lm<-lm[,order(colnames(lm))]


#Loop through coluumns
for (i in 1:ncol(lm)){
    rown<-rownames(lm)[i]
    coln<-colnames(lm)[i]
    if (!(coln %in% rownames(lm))){
        lm<-rbind(lm,rep(0,times=ncol(lm)))
        rownames(lm)[nrow(lm)]<-coln
    }
}

i=1
write.csv(lm, paste0("Output/LD/LD_pairwise_matrix_chr",i,'.csv'))



LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
