


#Combining all chromosomes to create genome-wide

ch1<-scan(paste0("Data/new_vcf/angsd/fromBam/PWS07_unfolded_chr1.sfs"))

pws07.sfs<-data.frame(chr1=ch1)
for (i in 2:26){
    vec<-scan(paste0("Data/new_vcf/angsd/fromBam/PWS07_unfolded_chr",i,".sfs"))
    pws07.sfs[,paste0("chr",i)]<-vec
}

pws07.sfs$sum<-rowSums(pws07.sfs)
write.table(pws07.sfs$sum, "Data/new_vcf/angsd/fromBam/combined/PWS07.sfs_unfolded.txt", quote=F, row.names = F, col.names = F)
barplot(pws07.sfs$sum[-c(1,nrow(pws07.sfs))])



combineSFS2<-function(pop){
    ch1<-scan(paste0("Data/new_vcf/angsd/fromBam/",pop,"_unfolded_chr1.sfs"))
    pws.sfs<-data.frame(chr1=ch1)
    for (i in 2:26){
        vec<-scan(paste0("Data/new_vcf/angsd/fromBam/",pop,"_unfolded_chr",i,".sfs"))
        pws.sfs[,paste0("chr",i)]<-vec
    }
    pws.sfs$sum<-rowSums(pws.sfs)
    sink(paste0("Data/new_vcf/angsd/fromBam/combined/",pop,"_unfolded.sfs"))
    cat(pws.sfs$sum)
    sink(NULL)
 }

combineSFS("PWS07")
combineSFS("PWS17")


#folded
combineSFS2<-function(pop){
    ch1<-scan(paste0("Data/new_vcf/angsd/fromBam/folded/",pop,"_folded_chr1.sfs"))
    pws.sfs<-data.frame(chr1=ch1)
    for (i in 2:26){
        vec<-scan(paste0("Data/new_vcf/angsd/fromBam/folded/",pop,"_folded_chr",i,".sfs"))
        pws.sfs[,paste0("chr",i)]<-vec
    }
    pws.sfs$sum<-rowSums(pws.sfs)
    sink(paste0("Data/new_vcf/angsd/fromBam/combined/",pop,"_folded.sfs"))
    cat(pws.sfs$sum)
    sink(NULL)
}

combineSFS2("PWS07")
combineSFS2("PWS17")


## At Farm computer load this scripts in R
combineSFSfold<-function(pop){
    ch1<-scan(paste0("/home/ktist/ph/data/angsd/SFS/fromBam/folded/",pop,"_folded_chr1.sfs"))
    pws.sfs<-data.frame(chr1=ch1)
    for (i in 2:26){
        vec<-scan(paste0("/home/ktist/ph/data/angsd/SFS/fromBam/folded/",pop,"_folded_chr",i,".sfs"))
        pws.sfs[,paste0("chr",i)]<-vec
    }
    pws.sfs$sum<-rowSums(pws.sfs)
    sink(paste0("/home/ktist/ph/data/angsd/SFS/fromBam/folded/",pop,"_folded.sfs"))
    cat(pws.sfs$sum)
    sink(NULL)
}

