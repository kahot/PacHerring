#LD Decay Analysis
source("Rscripts/BaseScripts.R")

## LD Decay 
#Window size comparison (plink setting)
ld1<- read.table("Data/ngsLD/PWS07_maf05_r2.ld", header = TRUE, stringsAsFactors = FALSE)
#ld2<-read.table("Data/ngsLD/PWS07_0threshold.ld", header = TRUE, stringsAsFactors = FALSE)
ld3<-read.table("Data/ngsLD/PWS07_0.05.ld", header = TRUE, stringsAsFactors = FALSE)

ld07<-read.table("Data/plink/LD/PWS07_maf05_r2.ld", header = TRUE, stringsAsFactors = FALSE)



ld1$distance <- ld1$BP_B-ld1$BP_A
#ld2$distance <- ld2$BP_B-ld1$BP_A
ld3$distance <- ld3$BP_B-ld3$BP_A
ld07$distance<- ld07$BP_B-ld07$BP_A

plot(ld1$distance, ld1$R2, col = "white", xlim = c(0,1000), xlab = "distance (bp)" ,ylab = "R2")
smooth_ld1 <-smooth.spline(ld1$distance, ld1$R2, spar = .15)
lines(smooth_ld1, lwd = 2,lty = 2, col = red)
smooth_ld3 <-smooth.spline(ld3$distance, ld3$R2, spar = .15)
lines(smooth_ld3, lwd = 2,lty = 3, col = blu)
smooth_ld07 <-smooth.spline(ld07$distance, ld07$R2, spar = .15)
lines(smooth_ld07, lwd = 2,lty = 3, col = gre)


plot(ld2$distance, ld2$R2, col = blu, xlim = c(0,1000), xlab = "distance (bp)" ,pch=16,ylab = "R2", cex=.5)
plot(ld1$distance, ld1$R2, col = blu, xlim = c(0,10000), xlab = "distance (bp)" ,pch=16,ylab = "R2", cex=.5)
plot(ld3$distance, ld3$R2, col = blu, xlim = c(0,1000), xlab = "distance (bp)" ,pch=16, ylab = "R2", cex=.5)

hist(ld3$distance)
hist(ld1$distance)

plot(ld07$distance, ld07$R2, col = blu, xlim = c(0,10000), xlab = "distance (bp)" ,pch=16, ylab = "R2", cex=.2)



ggplot()+
    geom_smooth(data=ld1[ld1$distance<=1000,], aes(x=distance, y=R2), method = "loess", se = FALSE, span = 1/10, color = blu)+
    geom_smooth(data=ld3[ld3$distance<=1000,], aes(x=distance, y=R2), method = "loess", se = FALSE, span = 1/10, color = red)+
        theme_minimal()+
        ylim(0,1)
ggsave("Output/LD/LD_decay_curve_PWS07_differenttheresholds.pdf", width = 7, height = 5)


ggplot()+
    geom_smooth(data=ld07[ld07$distance<=1000,], aes(x=distance, y=R2), method = "loess", se = FALSE, span = 1/10, color = blu)+
     theme_minimal()+
    ylim(0,1)
ggsave("Output/LD/LD_decay_curve_PWS07_0.1threshold.pdf", width = 7, height = 5)



ggplot(fst, aes(x = midPos, y = Fst01)) + 
    geom_point(size = 1, color = gry,alpha = 0.4, shape = 1)+
    theme_minimal()+
    theme(axis.text.x=element_blank())+
    ylab("Fst\n")+ xlab("")+ 
    ggtitle(paste(pops[1],"vs. ", pops[2]))+
    geom_smooth(method = "loess", se = FALSE, span = 1/10, color = blu)+
    facet_wrap(~chr, ncol = 9)
    
ld1$id<-1:nrow(ld1)
df<-ld1[ld1$CHR_A==1,]
plot(df$id, df$R2, pch=".")

cols <- rev(rainbow(7)[-7]) #rainbow colors for heatmap

ggplot(df, aes(x=BP_A, y=BP_B, color=R2))+
    geom_point(size=0.3)+
    scale_color_gradientn(colors=cols)+
    theme_minimal()
#plot along chr1
ggplot(df, aes(x=BP_A, y=R2, color=R2))+
    geom_point(size=0.3)+
    scale_color_gradientn(colors=cols)+
    theme_minimal()



pop_names = c("PWS91","PWS96","PWS07","PWS17")
pop_line <- c(4,3,2,1)
pop <- "PWS91"
i <- 1
ld <- read.table(paste0("Data/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld"), header = TRUE, stringsAsFactors = FALSE)
ld <- read.table("Data/ngsLD/PWS07_maf05_r2.ld", header = TRUE, stringsAsFactors = FALSE)


ld$distance <- ld$BP_B-ld$BP_A

plot(ld$distance, ld$R2, col = "white", xlim = c(0,1000), xlab = "distance (bp)" ,ylab = "R2")
smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)

for (i in c(2:length(pop_names))){
    ld <- read.table(paste("Data/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
    ld$distance <- ld$BP_B-ld$BP_A
    
    #ld <- ld[ld$CHR_A == "1" & ld$CHR_B == "1" ,]
    
    smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)
}


#Create more zoomed in view to assess the cutoff
pdf("Output/LD/LDdecay_PWS_zoomedin.pdf", height = 4, width = 6)
plot(ld$distance, ld$R2, type="n", xlim = c(0,150), xlab = "distance (bp)" ,ylab = "R2")
smooth_ld <-smooth.spline(ld$distance, ld$R2, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)

for (i in c(2:length(pop_names))){
    df <- read.table(paste("Data/plink/pops/population_",pop_names[i],"_ph_filtered_snps_minDP600_maxDP2000_minQ20_minMQ30_NS0.5_maf0.05_indep_pairwise_100_10_0.8.ld",sep = ''), header = TRUE, stringsAsFactors = FALSE)
    df$distance <- df$BP_B-df$BP_A
      
    smooth_ld <-smooth.spline(df$distance, df$R2, spar = .15)
    lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)
}
dev.off()


pdf("Output/LD/LDdecay_PWS91_zoomedin.pdf", height = 4, width = 6)
plot(ld$distance, ld$R2,col="gray60", pch=".", xlim = c(0,200), xlab = "distance (bp)" ,ylab = "R2")
dev.off()
pdf("Output/LD/LDdecay_PWS17_zoomedin.pdf", height = 4, width = 6)
plot(df$distance, df$R2,col="gray60", pch=".", xlim = c(0,200), xlab = "distance (bp)" ,ylab = "R2")
dev.off()

### 
#Plot SFS

#for (pop_name in pop_names){
sfs1 <- read.table("Output/fst_pbs/folded_PWS07_TB91.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs1 <- data.frame(t(sfs1))
sfs1<-as.vector(sfs1[,1])
barplot(sfs1[-1])

sfs2 <- read.table("Output/fst_pbs/folded_PWS07_TB91_2.sfs", header = FALSE, stringsAsFactors = FALSE)
sfs2 <- data.frame(t(sfs2))
sfs2<-as.vector(sfs2[,1])
barplot(sfs2[-1])


### LD from ngsLD
#BC17 chr 1
ldbc <- read.table(paste0("Data/ngsLD/sub_ch1_BC17.ld"), header=F, sep="\t")

plot(ldbc$V3, ldbc$V4, col = "white", xlab = "distance (bp)" ,ylab = "R2", xlim=c(0,1000))
smooth_ld <-smooth.spline(ldbc$V3, ldbc$V4, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = blu)

ldbc2<-ldbc[!(is.nan(ldbc$V6)|is.infinite(ldbc$V6)),]
plot(ldbc2$V3, ldbc2$V6, col = "white", xlab = "distance (bp)" ,ylab = "D'")
smooth_ld <-smooth.spline(ldbc2$V3, ldbc2$V6, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red)

max(ld2$V7)
plot(ldbc2$V3, ldbc2$V7, col = "white", xlab = "distance (bp)" ,ylab = "R2'",xlim = c(0,200000))
smooth_ld <-smooth.spline(ldbc2$V3, ldbc2$V7, spar = .15)
lines(smooth_ld, lwd = 2,lty = pop_line[i], col = red, xlim = c(0,200000))

plot(ldbc2$V3, ldbc2$V7,col="gray60", pch=".", xlim = c(0,400000), xlab = "distance (bp)" ,ylab = "R2")

## use LDheatmap

library(snpStats)
library(LDheatmap)

sample <- read.pedfile("Data/ngsLD/CA17_maf05_ch12_7-13.Haploview.ped.gz", snps="Data/ngsLD/CA17_maf05_ch12_7-13.Haploview.info.gz")

gdat<-sample$genotypes
dist<-read.table("Data/ngsLD/CA17_maf05_ch12_7-13.Haploview.info")
dist<-dist[,-1]
LDheatmap(gdat, genetic.distances=dist, distances="physical",
          LDmeasure="r", title="Pairwise LD CA17", add.map=TRUE, add.key=TRUE,
          geneMapLocation=0.15, geneMapLabelX=NULL, geneMapLabelY=NULL,
          SNP.name=NULL, color=NULL, newpage=TRUE,
          name="ldheatmap", vp.name=NULL, pop=FALSE, flip=NULL, text=FALSE)

sample <- read.pedfile("Data/ngsLD/PWS07_maf05_ch12_7.8-7.9.Haploview.ped", snps="Data/ngsLD/PWS07_maf05_ch12_7.8-7.9.Haploview.info")

gdat<-sample$genotypes
dist<-read.table("Data/ngsLD/PWS07_maf05_ch12_7.8-7.9.Haploview.info")
dist<-dist[,-1]

LDheatmap(gdat, genetic.distances=dist, distances="physical",
          LDmeasure="r", title="Pairwise LD PWS07 chr12", add.map=TRUE, add.key=TRUE,
          geneMapLocation=0.15, geneMapLabelX=NULL, geneMapLabelY=NULL,
          SNP.name=NULL, color=NULL, newpage=TRUE,
          name="ldheatmap", vp.name=NULL, pop=FALSE, flip=NULL, text=FALSE)
