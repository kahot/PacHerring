library(reshape2)
library(ggplot2)
## Look at the 3 groups in chr8 like Petrou et al. 2021

gp8<-read.csv("Output/chr/DP7000/chr8_PCAgroups.csv", row.names = 1)
gp8$Sample[gp8$Group==1]



# Calculate Ho and He from bcftools stats output files
sfiles<-list.files("Output/Stats_window/", pattern="_statsFile")
groups<-c(1,2,3)

#rearrage into 3 groups 
group1<-data.frame()
group2<-data.frame()
group3<-data.frame()

for (i in 1: length(sfiles)){
    df<-read.table(paste0("Output/Stats_window/", sfiles[i]), sep="\t", header=F)
    df<-df[,c(3:10, 14:15)]
    colnames(df)<-c("Sample","nRefHom","nNonRefHom","nHets", "nTransitions", "nTransversions","nIndels","average depth","nMissing","window_no")
    
    df$p<-(2*df$nRefHom+df$nHets)/(rowSums(df[,c("nRefHom","nNonRefHom","nHets")])*2)
    df$q<-(2*df$nNonRefHom+df$nHets)/(rowSums(df[,c("nRefHom","nNonRefHom","nHets")])*2)
    df$He<-2*df$p*df$q
    df$Ho<-df$nHets/rowSums(df[,c("nRefHom","nNonRefHom","nHets")])

    df1<-df[df$Sample %in%gp8$Sample[gp8$Group==1], ]
    group1<-rbind(group1, df1)
    df2<-df[df$Sample %in%gp8$Sample[gp8$Group==2], ]
    group2<-rbind(group2, df2)
    df3<-df[df$Sample %in%gp8$Sample[gp8$Group==3], ]
    group3<-rbind(group3, df3)
}

write.csv(group1, "Output/chr8/het/Group1.heterozygosity.csv")
write.csv(group2, "Output/chr8/het/Group2.heterozygosity.csv")
write.csv(group3, "Output/chr8/het/Group3.heterozygosity.csv")

#create a summary for each group

Het_sum<-data.frame(group=groups)
for (i in 1:3){
    df<-get(paste0("group",i))
    Ho<-aggregate(df[,"Ho"], by=list(df$window_no), mean )
    He<-aggregate(df[,"He"], by=list(df$window_no), mean )
    het<-cbind(Ho, He$x)
    colnames(het)<-c("window_id","Ho","He")
    write.csv(het,paste0("Output/chr8/het/Het.summary_",groups[i],".csv"))
    
    Het_sum$Ho_mean[i]<-  mean(het$Ho, na.rm=T)
    Het_sum$Ho_median[i]<-median(het$Ho,na.rm=T)
    Het_sum$He_mean[i]<-  mean(het$He,na.rm=T)
    Het_sum$He_median[i]<-median(het$He,na.rm=T)
    print(i)
}

write.csv(Het_sum, "Output/chr8/het/Hetero_group_summary.csv")

Het_sum<-read.csv("Output/Stats_window/Hetero_pop_summary.csv", row.names = 1)


Hetm<-melt(Het_sum[,c(1,2,4)], id.vars="group")
Hetm2<-melt(Het_sum[,c(1,3,5)], id.vars="pop")


Hetm$pop<-factor(Hetm$pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "TB91","TB96","TB06","TB17","BC17","WA17","CA17"))
Hetm2$pop<-factor(Hetm2$pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17", "TB91","TB96","TB06","TB17","BC17","WA17","CA17"))

ggplot(Hetm,aes(x=group, y=value, color=variable))+
    geom_point()+
    theme_bw()+
    xlab('')+ylab("Heterozygosity")+
    theme(legend.title = element_blank())


#
g1<-read.csv(paste0("Output/chr8/het/Het.summary_1.csv"), row.names = 1)
g2<-read.csv(paste0("Output/chr8/het/Het.summary_2.csv"), row.names = 1)
g3<-read.csv(paste0("Output/chr8/het/Het.summary_3.csv"), row.names = 1)

g1$group<-"group1"
g2$group<-"group2"
g3$group<-"group3"

gps<-rbind(g1,g2,g3)

ggplot()+
    geom_boxplot(data=gps, aes(x=group, y=Ho, color=group, fill=group),outlier.alpha = 0.2,  alpha=0.6)+
    geom_point(data=Het_sum, aes(x=group, y=Ho_mean))+
    theme_classic()+xlab('')
ggsave("Output/chr8/het/Ho_byGroup.pdf",height =4, width = 4 )


#Distribution of groups per pop.year

#regroup by population for grouping
pops<-gsub("_statsFile",'',sfiles)

library(dplyr)
pop.sum<-gp8 %>% count(Group, yr.pop)

pop.sum$yr.pop<-factor(pop.sum$yr.pop, levels=c("PWS91","PWS96","PWS07","PWS17","SS96","SS06","SS17","BC17","WA17","CA17"))
ggplot(data=pop.sum, aes(x=yr.pop, y=n, fill=factor(Group)))+
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values=c("red", org, blu), labels=c("Group1","Group2","Group3"))+
    xlab("")+ylab("Proportion")+theme(legend.title = element_blank())
ggsave("Output/chr/DP7000/chr8_group_barplot.pdf", width = 8, height = 5.5)



