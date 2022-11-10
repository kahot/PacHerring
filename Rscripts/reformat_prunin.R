reformat_prunin<-function(path){
    df<-read.table(path, sep=":")
    df$V2<-gsub("\\[ph+.*",'',df$V2)
    df$V2<-as.integer(df$V2)
    write.table(df, paste0(path, ".sites.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
}
