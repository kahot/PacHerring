
df<-read.table("Data/plink/prune/ph_75_5_0.5.prune.out")
df$V1<-gsub("\\[ph\\]","",df$V1)
df$chr<-substr(df$V1, 1,5)
df$pos<-substr(df$V1, 1,5)
split<-strsplit(df$V1, ":" )

out<-data.frame(sapply(split,"[[",1))
out$pos<-sapply(split,"[[",2)

write.table(out,"Data/plink/prune/out75.txt", quote = F, row.names = F, col.names = F)

df<-read.table("Data/plink/prune/ph_75_5_0.5.prune.in", sep=":")
df<-read.table("Data/plink/prune/MD7000_maf0.05.prune.in", sep=":")
df$V2<-gsub("\\[ph\\]","",df$V2)

#out<-data.frame(sapply(split,"[[",1))
#out$pos<-sapply(split,"[[",2)

write.table(out,"Data/plink/prune/in75.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(df,"Data/plink/prune/in75MD7000.txt", sep = "\t", quote = F, row.names = F, col.names = F)
