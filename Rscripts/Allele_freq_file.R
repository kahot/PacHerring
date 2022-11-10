#Create sliding window allele freq files

library(reshape2)
library(windowscanr)


# Prince William Sound
t0 <- read.table("Data/new_vcf/population/PWS91_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t1 <- read.table("Data/new_vcf/population/PWS96_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t2 <- read.table("Data/new_vcf/population/PWS07_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t3 <- read.table("Data/new_vcf/population/PWS17_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))

# Stika Sound
t1 <- read.table("Data/new_vcf/population/SS96_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t2 <- read.table("Data/new_vcf/population/SS06_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t3 <- read.table("Data/new_vcf/population/SS17_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))

# TB
t0 <- read.table("Data/new_vcf/population/TB91_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t1 <- read.table("Data/new_vcf/population/TB96_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t2 <- read.table("Data/new_vcf/population/TB06_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))
t3 <- read.table("Data/new_vcf/population/TB17_freq.frq",stringsAsFactors = FALSE,header = FALSE, skip=1, col.names = c("chr","pos","n_allele","n_sample","MajorAF","MAF"))




for (i in c(0:3)){
    df<-get(paste0("t",i))
    df$maf<-substr(df$MAF, 3,10)
    df<-df[,c(1,2,7)]
    df$maf<-as.numeric(df$maf)
    assign(paste0("t",i),df)
}

freqs <- data.frame(chr = t0$chr,
                    pos = t0$pos,
                    t0_AF = t0$maf,
                    t1_AF = t1$maf,
                    t2_AF = t2$maf,
                    t3_AF = t3$maf)

freqs$t0_transformed_freq <- asin(sqrt(freqs$t0_AF))
freqs$t1_transformed_freq <- asin(sqrt(freqs$t1_AF))
freqs$t2_transformed_freq <- asin(sqrt(freqs$t2_AF))
freqs$t3_transformed_freq <- asin(sqrt(freqs$t3_AF))

freqs$zt01 <- freqs$t1_transformed_freq - freqs$t0_transformed_freq
freqs$zt02 <- freqs$t2_transformed_freq - freqs$t0_transformed_freq
freqs$zt03 <- freqs$t3_transformed_freq - freqs$t0_transformed_freq
freqs$zt12 <- freqs$t2_transformed_freq - freqs$t1_transformed_freq
freqs$zt13 <- freqs$t3_transformed_freq - freqs$t1_transformed_freq
freqs$zt23 <- freqs$t3_transformed_freq - freqs$t2_transformed_freq


#write.table(freqs, "Data/freqs/MD7000/PWS_shifts_persite.txt")
#write.table(freqs, "Data/freqs/MD7000/SS_shifts_persite.txt")

window_size <- 50000
step_size <- 5000
freqs_win <- winScan(x = freqs, 
                     groups = "chr", 
                     position = "pos", 
                     values = c("zt01","zt02","zt03","zt12","zt13","zt23"), 
                     win_size = window_size,
                     win_step = step_size,
                     funs = c("mean"))
freqs_win <- na.omit(freqs_win)
freqs_win <- freqs_win %>% mutate(chr = as.numeric(gsub("chr", "", chr)))
freqs_win$chr_num <- factor(freqs_win$chr, levels = c(1:26))

#write.csv(freqs_win, "Data/freqs/MD7000/PWS_shifts_50k_5k.txt")
write.csv(freqs_win, "Data/freqs/MD7000/SS_shifts_50k_5k.txt")
#write.csv(freqs_win, "Data/freqs/MD7000/TB_shifts_50k_5k.txt")



freqs_win_p <- freqs_win
freqs_win_t <- freqs_win
freqs_win_s<-  freqs_win
