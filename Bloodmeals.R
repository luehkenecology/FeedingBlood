#' ---
#' title: "Distribution of Mesionoviridae in Germany"
#' author: "Renke Luehken"
#' date: "xxxxx, 2015"
#' ---

# clear memory
rm(list = ls())

library(ggplot2)
library(reshape2)
library(vegan)
library(plyr)
library(ggdendro)

# working directory
setwd("F:/NeuAll/R/BLOODMEALS")
# setwd("~/Google Drive/R/BLOODMEALS")

# read data
data<-read.table (file = "data/New/Bloodmeals_complete_New_RL.txt",row.names=1,header=TRUE,sep="\t")
data2<-na.omit(data)
subset(data2,SPECIES=="Cx. spec.")


##############################
# correlation between the diversity of mosquitoes and the diversity of detected hosts per site
##############################

# data frame manipulation --> TRAPPINGSITE x HOST.DNA.NEW table
AS<-count(data2, c("TRAPPINGSITE", "HOST.DNA.NEW"))
AS2<-acast(AS, TRAPPINGSITE~HOST.DNA.NEW, value.var="freq")
AS2[is.na(AS2)] <- 0

# caluclation of the species richness
Diversity.mosquitoes <- diversityresult(AS2, index='richness' ,method='s', sortit=TRUE, digits=3)

# data frame manipulation --> TRAPPINGSITE x SPECIES table
AS<-count(data2, c("TRAPPINGSITE", "SPECIES"))
AS2<-acast(AS, TRAPPINGSITE~SPECIES, value.var="freq")
AS2[is.na(AS2)] <- 0

# caluclation of the species richness
Diversity.hosts <- diversityresult(AS2, index='richness' ,method='s', sortit=TRUE, digits=3)

# merge both diversity indices in one data.table
Diversity.merge<-merge(Diversity.mosquitoes,Diversity.hosts,by="row.names")

#
gg2 <- ddply(Diversity.merge, c("richness.x", "richness.y"), "nrow", .drop = T)

# spearman correlation between number of host species and the number of mosquito species
cor.test(Diversity.merge$richness.x,Diversity.merge$richness.y,method="spearman")

# plot
ggplot(data = gg2, aes(x = richness.y, y = richness.x, size = factor(nrow))) + 
  geom_point() + 
  theme_classic() +
  ylab("number of host species") +
  xlab("number of mosquito species") +
  scale_size_discrete(range = c(1, 10),guide = guide_legend(title = "trapping sites"))










#
data3<-subset(data2,HOST.GROUP=="mammal")
d_dypl<-ddply(data3,.(HOST.DNA.NEW),summarize, freq= sum(POOLSIZE))
newdata <-cbind(d_dypl[with(d_dypl, order(-freq)), ],num=seq(1,nrow(d_dypl),1))
newdata3<-merge(data3,newdata,by="HOST.DNA.NEW")
newdata3$HOST.DNA.NEW<-reorder(newdata3$HOST.DNA.NEW,newdata3$num)

data4<-subset(data2,HOST.GROUP=="bird")
d_dypl2<-ddply(data4,.(HOST.DNA.NEW),summarize, freq= sum(POOLSIZE))
newdata <-cbind(d_dypl2[with(d_dypl2, order(-freq)), ],num=seq(1,nrow(d_dypl2),1))
newdata4<-merge(data4,newdata,by="HOST.DNA.NEW")
newdata4$HOST.DNA.NEW<-reorder(newdata4$HOST.DNA.NEW,newdata4$num)

p2<-ggplot(data2, aes(x = HOST.DNA.NEW,fill=as.factor(PREVALENCE))) + 
  geom_bar(data=newdata3)+ 
  geom_bar(data=newdata4)+ 
  theme_classic()+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0,hjust=1))+
  facet_wrap(~HOST.GROUP,scales="free")+
  scale_fill_manual(guide = guide_legend(title = "Mesoni infection"),values=c("black","gray"))+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text =  element_text(face = "italic"))

png(file = "output/FeedingPreferenceMosquitoesMesoniInfection.jpg",width = 8, height=5, units = 'in', res = 1000)
plot(p2)
dev.off()

d_dypl2<-ddply(data2,.(SPECIES),summarize, freq= sum(POOLSIZE))
newdata <-cbind(d_dypl2[with(d_dypl2, order(-freq)), ],num=seq(1,nrow(d_dypl2),1))
newdata5<-merge(data2,newdata,by="SPECIES")
newdata5$SPECIES<-reorder(newdata5$SPECIES,newdata5$num)

colnames(newdata5)
p2<-ggplot(newdata5,aes(x=SPECIES,fill=HOST.GROUP.HUMAN)) +
  geom_bar(colour= "black")+
  theme_classic()+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0,hjust=1))+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text =  element_text(face = "italic"))


  png(file = "output/FeedingPreferenceMosquitoes.jpg",width = 6, height=5, units = 'in', res = 1000)
plot(p2)
dev.off()


dfdfdf<-ddply(data2,.(SPECIES, HOST.DNA.NEW),summarize,sum=sum(POOLSIZE))
ggplot(dfdfdf, aes(SPECIES, HOST.DNA.NEW)) + geom_tile(aes(fill = sum),
                                                      colour = "white") + scale_fill_gradient(low = "white",
                                                                                              high = "steelblue")
p <- ggplot(dfdfdf, aes(x=SPECIES,y=HOST.DNA.NEW))+ geom_tile(aes(fill=sum))+
  scale_fill_gradient(guide = guide_legend(title = "sum of individuals"),low="green", high="red")+
  theme_bw()+
  ylab("host species")+
  xlab("mosquito species")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0,hjust=1))+
  theme(axis.text =  element_text(face = "italic"))

  
png(file = "output/HeatMapFeedingPreference.jpg",width = 7, height=7, units = 'in', res = 1000)
plot(p)
dev.off()
  



AS<-count(data2, c("SPECIES", "HOST.DNA.NEW"))


AS<-count(data2, c("SPECIES", "HOST.GROUP.HUMAN"))

AS2<-acast(AS, SPECIES~HOST.DNA.NEW, value.var="freq")


AS2<-acast(AS, SPECIES~HOST.GROUP.HUMAN, value.var="freq")

AS2[is.na(AS2)] <- 0

#AS2<-ifelse(AS2>0,1,0)




env.norm <- decostand(AS2, method = "normalize", na.rm = TRUE, MARGIN=2)

spe.ch<-vegdist(env.norm,method="bray")
hc <- hclust(dist(spe.ch), "average")
dhc <- as.dendrogram(hc,hang=0.1)
ddata <- dendro_data(dhc, type="rectangle")

clust    <- cutree(hc,k=5)                    # find 2 clusters
clust.df <- data.frame(label=names(clust)
                       , cluster=factor(clust))
ddata[["labels"]] <- merge(ddata[["labels"]],clust.df, by="label")


p<-ggplot() + 
  geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(ddata), aes(x, y, label=label, hjust=0, color=cluster), 
            size=3,fontface=3) +
  ylab("distance")+
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme_bw()+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x =  element_text(face = "italic"),
        #text =  element_text(face = "italic"),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())


png(file = "output/FeedingPreferenceMosquitoesDendgrogram.jpg",width = 8, height=6, units = 'in', res = 1000)
plot(p)
dev.off()



dimnames(ddata$labels)[[2]]<-c("SPECIES","x","y","cluster")
OINK<-merge(data2,ddata$labels,by="SPECIES")
p2<-ggplot(OINK,aes(x=SPECIES,fill=HOST.GROUP.HUMAN)) +
  ylab("percentage")+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  geom_bar(colour="black",position="fill")  +theme_classic()+  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~cluster)

p2<-ggplot(OINK,aes(x=cluster,fill=HOST.GROUP.HUMAN)) +
  ylab("count")+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  geom_bar(colour="black")  +theme_classic()+  scale_y_continuous(expand=c(0,0))


png(file = "output/HostpreferenceDendgromCluster.jpg",width = 5, height=5, units = 'in', res = 1000)
plot(p2)
dev.off()



####################
####################
####################
# 
colnames(data2)
AAA<-ddply(data2,.(PREVALENCE,HOST.GROUP),summarize,sum=length(POOLSIZE))

AAA<-ddply(data2,.(TRAPPINGSITE,HOST.GROUP.HUMAN),summarize,sum=length(POOLSIZE))


A<-c(41,6)
B<-c(358,111)
data.table<-rbind(A,B)
chisq.test(data.table)

AAA<-ddply(data2,.(SPECIES,TRAPPINGSITE),summarize,sum=length(POOLSIZE))
sum(AAA$sum)
# dfdf
AAAA<-ddply(AAA,.(SPECIES),summarize,sum=length(sum))
AAAA[ order(AAAA[,2]), ]

table(data2$YEAR)
colnames(data2)

AS<-count(data2, c("TRAPPINGSITE", "HOST.GROUP.HUMAN"))


AS2<-acast(AS, TRAPPINGSITE~HOST.GROUP.HUMAN, value.var="freq")

AS2[is.na(AS2)] <- 0

#AS2<-ifelse(AS2>0,1,0)




env.norm <- decostand(AS2, method = "normalize", na.rm = TRUE, MARGIN=2)

spe.ch<-vegdist(env.norm,method="bray")
hc <- hclust(dist(spe.ch), "average")
dhc <- as.dendrogram(hc,hang=0.1)
ddata <- dendro_data(dhc, type="rectangle")

clust    <- cutree(hc,k=8)                    # find 2 clusters
clust.df <- data.frame(label=names(clust)
                       , cluster=factor(clust))
ddata[["labels"]] <- merge(ddata[["labels"]],clust.df, by="label")


p<-ggplot() + 
  geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(ddata), aes(x, y, label=label, hjust=0, color=cluster), 
            size=3,fontface=3) +
  ylab("distance")+
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme_bw()+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x =  element_text(face = "italic"),
        #text =  element_text(face = "italic"),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

plot(p)

png(file = "output/FeedingPreferenceMosquitoesDendgrogram.jpg",width = 8, height=6, units = 'in', res = 1000)
plot(p)
dev.off()

OINK$cluster
dfdf<-subset(OINK,cluster==3)
sum(table(dfdf$TRAPPINGSITE))

dimnames(ddata$labels)[[2]]<-c("TRAPPINGSITE","x","y","cluster")
OINK<-merge(data2,ddata$labels,by="TRAPPINGSITE")
p2<-ggplot(OINK,aes(x=SPECIES,fill=HOST.GROUP.HUMAN)) +
  ylab("percentage")+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  geom_bar(colour="black",position="fill")  +theme_classic()+  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~cluster)

p2<-ggplot(OINK,aes(x=cluster,fill=HOST.GROUP.HUMAN)) +
  ylab("count")+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  geom_bar(colour="black")  +theme_classic()+  scale_y_continuous(expand=c(0,0))
plot(p2)
