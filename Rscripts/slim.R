
#Read SliM Results




#70 genomes (35 individuals) at generation 40 
gen50<-read.table("~/Projects/SLiM/PacHerring/gen50_100ind.txt")
gen55<-read.table("~/Projects/SLiM/PacHerring/gen55_100ind.txt")
#convert to pop allele freq

gen50$AF<-gen50$V9/100
gen55$AF<-gen55$V9/100
max(gen50$AF)

pops<-c("p1","p2","p3")
time<-c(500,550,600,1000)
summary<-data.frame()
for (i in 1:length(pops)){
    results<-data.frame()
    for (j in 1:length(time)){
        df<-read.table(paste0("~/Projects/SLiM/PacHerring/",pops[i], "_g",time[j],".txt"))
        df$AF<-df$V9/200
        df$time<-time[j]
        results<-rbind(results,df)
    }
    results$pop<-pops[i]
    summary<-rbind(summary,results)
}
    
library(ggplot2)
mean<-data.frame(aggregate(summary$AF, by=list(summary$pop, summary$time), FUN=mean))

ggplot(summary, aes(x=factor(time), y=AF, color=pop))+
    geom_boxplot(position=position_dodge(width=1))+
    geom_point(data=mean, aes(x=factor(Group.2), y=x, fill=Group.1), position=position_dodge(width=1), color="gray40")+
    scale_x_discrete()

ggplot(summary[summary$pop=="p1",], aes(x=time,y=AF))+
    geom_point()   +
    geom_path(aes(x=time, y=AF, group=factor(V2)))


muts<-unique(summary$V2)
sum1<-summary[summary$V2 %in% muts[1:10],]

ggplot(sum1, aes(x=V4,y=AF, color=pop))+
    geom_point(size=1, position=position_dodge(width=1))+
    facet_wrap(~factor(time),ncol=1)
    

## Read nonWF model ran for 1000 generations (netural mutations only)

time<-c(100,500,1000)
results<-data.frame()
for (j in 1:length(time)){
    df<-read.table(paste0("~/Projects/SLiM/mutation",time[j],".txt"), skip=2)
        df$AF<-df$V9/200
        df$time<-time[j]
        results<-rbind(results,df)
}



ggplot(results, aes(x=factor(time), y=AF))+
    geom_boxplot()+
    geom_point(stat="summary", fun="mean", color="gray40")+
    scale_x_discrete()

#track the allele freq change?
muts<-unique(results$V2)
re2<-
ggplot(results, aes(x=time, y=AF))+
    geom_point()+
    facet_wrap(!V2)