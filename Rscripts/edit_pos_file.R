#manipulate the position file

pos<-read.csv("~/programs/ngsLD/positions.txt", sep="\t")

pos1<-pos[pos$CHR=="chr1",]
write.table(pos1, "~/programs/ngsLD/pos_chr1.txt", row.names = F, quote = F)


pos<-read.csv("Data/new_vcf/new_positions", sep="\t", header = F)

for (i in 1:26){
    df<-pos[pos$V1==paste0("chr",i),]
    write.table(df, paste0("Data/new_vcf/positions/pos_chr",i,".txt"),sep="\t",row.names = F, quote = F, col.names = F) 
}



### 
pops_info<-read.csv("Data/Sample_metadata_892pops.csv")
pops<-unique(pops_info$Population.Year)

pwapops<-c("PWS07", "PWS17", "PWS91" ,"PWS96")
comb<-combn(pwapops, 2)
comb<-t(comb)
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

evens<-paste0("chr",seq(2,26, by=2))

pruunedPws<-read.table("Data/new_vcf/pruned/prunedpws_BEAGLE.PL", sep = "\t")
colnames(pruunedPws)<-pruunedPws[1,]
samplenames<-colnames(pruunedPws)

for (i in 2:6){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    samplenames<-colnames(pruunedPws)
    df<-pruunedPws[,c(1:3, grep(paste0(pop1,"|",pop2), samplenames))]

    write.table(df, paste0("Data/new_vcf/pruned/pruned_",pop1,"_",pop2,".BEAGLE.PL"), col.names = F,row.names = F, quote=F, sep = "\t")
}

#2. SS
sspops<-c("SS96","SS06", "SS17")
comb<-combn(sspops, 2)
comb<-t(comb)

pruSs<-read.table("Data/new_vcf/pruned/prunedss_BEAGLE.PL", sep = "\t")
colnames(pruSs)<-pruSs[1,]
samplenames<-colnames(pruSs)

for (i in 1:nrow(comb)){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    samplenames<-colnames(pruSs)
    df<-pruSs[,c(1:3, grep(paste0(pop1,"|",pop2), samplenames))]
    
    write.table(df, paste0("Data/new_vcf/pruned/pruned_",pop1,"_",pop2,".BEAGLE.PL"), col.names = F,row.names = F, quote=F, sep = "\t")
}


#2. SS
tbpops<-c("TB91","TB96","TB06", "TB17")
comb<-combn(tbpops, 2)
comb<-t(comb)

pruSs<-read.table("Data/new_vcf/pruned/prunedtb_BEAGLE.PL", sep = "\t")
colnames(pruSs)<-pruSs[1,]
samplenames<-colnames(pruSs)

for (i in 1:nrow(comb)){
    pop1<-comb[i,1]
    pop2<-comb[i,2]
    
    samplenames<-colnames(pruSs)
    df<-pruSs[,c(1:3, grep(paste0(pop1,"|",pop2), samplenames))]
    
    write.table(df, paste0("Data/new_vcf/pruned/pruned_",pop1,"_",pop2,".BEAGLE.PL"), col.names = F,row.names = F, quote=F, sep = "\t")
}



