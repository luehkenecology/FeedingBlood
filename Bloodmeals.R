Test2000

#' ---
#' title: "Bloodmeals in German mosquito speces"
#' author: "Renke Luehken"
#' date: format(Sys.time(), "%d%b %Y")
#' ---

df
df
#df2

# clear memory
rm(list = ls())
library(EcoSimR)

# load libraries
library(plyr)
library(reshape2)

library(ggplot2)
library(ggdendro)

library(vegan)
library(BiodiversityR)

#library(devtools)
#install_github("cran/BiodiversityR")

# working directory
setwd("F:/NeuAll/R/BLOODMEALS")
setwd("C:/Users/RenkeLuehken/Google Drive/R/BLOODMEALS")
#setwd("~/Google Drive/R/BLOODMEALS")

# read data
data <- read.table (file = "data/blood_mesoni.txt",row.names=1,header=TRUE,sep="\t")

PathoNew<-read.table (file = "data/PathoNew.txt",row.names=1,header=TRUE,sep="\t")

data_om <- na.omit(data)

#============================================================
# glmm analysis
#============================================================
# Test case : the non-monotonic Sobol g-function
# The method of sobol requires 2 samples
# (there are 8 factors, all following the uniform distribution on [0,1])
library(sensitivity)
library(boot)
library(lhs)

A <- randomLHS(n=10000, k=1)
randomLHS(n=10000, k=1)

X1 <- cbind(Fa= randomLHS(n=100000, k=1),
            FH=randomLHS(n=100000, k=1)*0.2,
            Fh= randomLHS(n=100000, k=1),
            A=runif(100000,0.21,326),
            Cv=runif(100000,0.1,0.9))
gtz <- X1[,1] + X1[,2] + X1[,3]
gtzT <- gtz<=1
X1b<-(X1[gtzT,])

X2<-cbind(Fa= randomLHS(n=100000, k=1),
          FH=randomLHS(n=100000, k=1)*0.2,
          Fh= randomLHS(n=100000, k=1),
          A=runif(100000,0.21,326),
          Cv=runif(100000,0.1,0.9))
gtz <- X2[,1]+X2[,2]+X2[,3]
gtzT <- gtz<=1
X2b<-(X2[gtzT,])


#X2[1:10,]
#X1 <- data.frame(matrix(runif(8 * n), nrow = n))

y <- function(X) X[,1]*X[,2] * X[,3] * X[,4]

#sobol.fun

X1c <- data.frame(X1b[,c(1,2,4,5)])
X2c <- data.frame(X2b[,c(1,2,4,5)])



n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis
x <- sobol(model = sobol.fun, X1 = X1, X2 = X2, order = 2, nboot = 100)
print(x)




# sensitivity analysis
x <- sobol(model = y, X1 = X1c[1:10000,], X2 = X2c[1:10000,], order = 2,
           nboot = 10000)
print(x)
str(x)

tztzt <- data.frame(AAA = x$X[,1]*x$X[,2], BBB = x$y)
ggplot(tztzt, aes(x = AAA, y = BBB)) +
  geom_point(colour = "gray", size = 0.1) +
  theme_classic() #+
  #stat_smooth(method = loess, colour = "black")

plot(x$X[,1]*x$X[,1],
     x$y)
#plot(x)

#============================================================
# glmm analysis
#============================================================
colnames(data_om)
## table 1
ddply(data_om, .(HOST.GROUP.HUMAN, SPECIES), summarize, sum(POOLSIZE),.drop=F)
ddply(data_om, .(SPECIES), summarize, length(unique(HOST.DNA.NEW)),.drop=F)
ddply(data_om, .(HOST.GROUP.HUMAN,HOST.DNA.NEW), summarize, sum(POOLSIZE))



# data preperation
library(fifer)

# Sampling Method
XYb <- ddply(data_om,.(SAMPLING.METHOD, HOST.GROUP),
            summarize, sum1 = sum(POOLSIZE),.drop=F)
XYc <- rbind(XYb$sum1[1:2], XYb$sum1[3:4],
             XYb$sum1[5:6], XYb$sum1[7:8],
             XYb$sum1[9:10], XYb$sum1[11:12],
             XYb$sum1[13:14], XYb$sum1[15:16])
dimnames(XYc) <- list(period=c(as.character(unique(XYb$SAMPLING.METHOD))),
                      loc=c("bird", "mammal"))
chisq.test(XYc[apply(XYc, 1, function(x) !all(x==0)),])
chisq.post.hoc(XYc)


XY <- ddply(data,.(SPECIES, SAMPLING.METHOD, HOST.GROUP.HUMAN),
      summarize, sum1 = sum(POOLSIZE),.drop=F)

XYb <- subset(XY, SPECIES=="Ae. vexans")
XYb <- subset(XY, SPECIES=="Cx. pipiens pipiens")
XYb <- subset(XY, SPECIES=="Oc. cantans")

XYc <- rbind(XYb$sum1[1:3], XYb$sum1[4:6],
             XYb$sum1[7:9], XYb$sum1[10:12],
             XYb$sum1[13:15], XYb$sum1[16:18],
             XYb$sum1[19:21], XYb$sum1[22:24])
chisq.test(XYc[apply(XYc[,-1], 1, function(x) !all(x==0)),])

# period
colnames(data)
cbind(data$DATE, data$period)
XYb <- ddply(data,.(period, HOST.GROUP.HUMAN),
             summarize, sum1 = sum(POOLSIZE),.drop=F)
XYc <- rbind(XYb$sum1[1:3], XYb$sum1[4:6])
dimnames(XYc) <- list(period=c(as.character(unique(XYb$SAMPLING.METHOD))),loc=c("bird", "human","non-human mammal"))
chisq.test(XYc[apply(XYc[,-1], 1, function(x) !all(x==0)),])
chisq.post.hoc(XYc)


XY <- ddply(data_om,.(SPECIES, SAMPLING.METHOD, HOST.GROUP.HUMAN),
            summarize, sum1 = sum(POOLSIZE),.drop=F)

XYb1 <- subset(XY, SPECIES=="Ae. vexans")
XYb2 <- subset(XY, SPECIES=="Cx. pipiens pipiens")
XYb3 <- subset(XY, SPECIES=="Oc. cantans")

rbind(XYb1,XYb2,XYb3)

XYc <- rbind(XYb$sum1[1:3], XYb$sum1[4:6])
chisq.test(XYc[apply(XYc[,-1], 1, function(x) !all(x==0)),])
#Go through each row and determine if a value is zero
row_sub = apply(XYc, 1, function(row) all(row !=0 ))
##Subset as usual
rt<-XYc[row_sub,]

#
#dreDre <- subset(data, SPECIES=="Ae. vexans" | SPECIES=="Cx. pipiens pipiens" | SPECIES=="Oc. cantans")
#tzt <- na.omit(dreDre)

trap_method_1 <- ddply(data_om,.(SAMPLING.METHOD),
                       summarise,d3=length(HOST.GROUP.HUMAN))
trap_method_2 <- merge(trap_method_1, data_om, by="SAMPLING.METHOD")
trap_method_3 <- ddply(trap_method_2,.(SAMPLING.METHOD, HOST.GROUP.HUMAN),summarise,
             d4=length(HOST.GROUP.HUMAN)/d3*100)
trap_method_4 <-unique(trap_method_3[c("SAMPLING.METHOD","HOST.GROUP.HUMAN","d4")])
trap_method_5 <- merge(trap_method_4, trap_method_1, by="SAMPLING.METHOD")
trap_method_1$pos <- 101

trap_method <- ggplot(trap_method_5, aes(x = SAMPLING.METHOD,y=d4,ymax=140,fill=HOST.GROUP.HUMAN)) +
  geom_bar(colour="black",position = 'stack',stat="identity")+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  #facet_wrap(~cluster,scales="free")+
  ylab("percentage")+
  xlab("trapping method")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0,hjust=1))+
  scale_y_continuous(expand=c(0,0))+
  annotate("text", x = 1:8, y = 90, label = c(trap_method_1$d3))
























dimnames(XYc) <- list(period=c("early", "late"),loc=c("bird", "human","non-human mammal"))


XY <- ddply(data,.(period, HOST.GROUP.HUMAN),
            summarize, sum1 = sum(POOLSIZE),.drop=F)
XYb <- subset(XY, SPECIES=="Oc. cantans")
XYc <- rbind(XYb$sum1[1:3], XYb$sum1[4:6])
dimnames(XYc) <- list(period=c("early", "late"),loc=c("bird", "human","non-human mammal"))
chisq.test(XYc)


# Makes a table of observations -- similar to first example in chisq.test
M <- as.table(rbind(c(76, 32), c(48,23), c(45,34)))
dimnames(M) <- list(sex=c("Male","Female","Juv"),loc=c("Lower","Middle"))
M
chisq.test(M)
# Shows post-hoc pairwise comparisons using fdr method
chisq.post.hoc(M)
#

colnames(data)
A <- ddply(data,
           .(TRAPPINGSITE,YEAR,HOST.GROUP.HUMAN,SAMPLING.METHOD,period),
           summarize, sum1=sum(POOLSIZE))
B <- ddply(A,.(TRAPPINGSITE, YEAR),summarize,sum2=sum(sum1))
C <- merge(A,B,by=c("TRAPPINGSITE","YEAR"))
C$dep <- C$sum1/C$sum2

C2 <- subset(C, SPECIES=="Cx. pipiens pipiens")
dsub <- subset(C,HOST.GROUP.HUMAN=="human")
my.mod0 <- glm(dep ~  SAMPLING.METHOD + period, data = dsub, family = quasibinomial)
my.mod1 <- glm(dep ~  period, data = dsub, family = quasibinomial)
drop1(my.mod0,test="Chisq")

summary(my.mod)


library(lme4)
lmer(dep~1+(SPECIES|TRAPPINGSITE),data=C)
table(C$dep)

colnames(birdA)
birdA$feq <- ifelse(birdA$sum1 > 0, 1, 0)
birdA <- subset(A, HOST.GROUP.HUMAN=="bird")
TT1<- lmer(dep ~ (TRAPPINGSITE |SPECIES), family = quasibinomial, data = C)
TT1<- glmer(dep ~ 1 + (1|SPECIES:TRAPPINGSITE:YEAR),
            weights = C$sum2, family = binomial, data = C)

TT1<- glmer(feq ~ 1 + (1|TRAPPINGSITE), family = binomial, data = birdA)
TT2<- glmer(feq ~ 1 + (1|YEAR:TRAPPINGSITE:SPECIES),
            family = binomial,
            data = birdA)

pchisq(logLik(TT1) - logLik(TT2), 1)
anova(TT1, TT2)

model<-lmer(dep ~ (1|SPECIES),quasibinomial,data=C) 
summary(model)



pchisq(logLik(TT1) - logLik(TT2), 1)

vc<-VarCorr(TT2)
summary(TT)
sapply(TT, slot, "x")
summary(TT)
?lr.test()

colnames(PathoNew)
Patho2<-ddply(PathoNew,)

#
#
#
#
?weighted.mean
library(SDMTools)
Mend<-ddply(PathoNew,.(Host.class,Pathogen),summarise,
        wm = wt.mean(Prop,Number.of.positive.pools),
        wsd=wt.sd(Prop,Number.of.positive.pools))

weighted.mean(PathoNew$Prop,PathoNew$Number.of.positive.pools,na.rm=T)

library(ggplot2)
limits <- aes(ymax = wm + wsd, ymin=wm - wsd)
pc<-ggplot(Mend,aes(Pathogen,wm))+
  geom_point()+
  geom_errorbar(limits, width=0.25)+
  facet_wrap(~Host.class)+
  theme_bw()+
  ylab("percentage [%]")+
  xlab("pathogen")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust=1))
  

# save plot
png(file = "output/OUT.jpg",width = 10, height=5, units = 'in', res = 1000)
plot(pc)
dev.off()

# omit NA values (e.g. not detected host,...)
data2<-na.omit(data)

# percentage mosquitoes
data.frame(ddply(data2,.(SPECIES),summarise,freq=length(SPECIES)),
      perc=round(ddply(data2,.(SPECIES),summarise,freq=length(SPECIES))$freq/775*100,1))
46.8+12.9+12.8

ddply(data2,.(HOST.DNA.NEW,SPECIES),summarise,freq=length(SPECIES))
colnames(data2)
18/20

colnames(data2)
A<-ddply(data2,.(SPECIES,HOST.GROUP.HUMAN),summarise,freq=length(SPECIES))
B<-ddply(A,.(SPECIES),summarise,sum=sum(freq))
C<-merge(A,B,by="SPECIES")
C$perc<-C$freq/C$sum*100

C2<-C[24:32,]
ddply(C2,.(HOST.GROUP.HUMAN),summarise,mean(perc))
sum(28,4,1,8)/126*100
sum(35,6,6)/126*100

C3<-C[c(2:11,33:47),]
ddply(C3,.(HOST.GROUP.HUMAN),summarise,mean(perc))
sum(25,42,363,15,8,99,2,3,1,2,1,2,16)

data.frame(ddply(data2,.(SPECIES,HOST.GROUP.HUMAN),summarise,freq=length(SPECIES)),
           perc=round(ddply(data2,.(SPECIES,HOST.GROUP.HUMAN),summarise,freq=length(SPECIES))$freq/775*100,1))


# percentage hosts
data.frame(ddply(data2,.(HOST.DNA.NEW),summarise,freq=length(HOST.DNA.NEW)),
           perc=round(ddply(data2,.(HOST.DNA.NEW),summarise,freq=length(HOST.DNA.NEW))$freq/775*100,1))

# percentage hosts
data.frame(ddply(data2,.(HOST.GROUP.HUMAN),summarise,freq=length(HOST.GROUP.HUMAN)),
           perc=round(ddply(data2,.(HOST.GROUP.HUMAN),summarise,freq=length(HOST.GROUP.HUMAN))$freq/775*100,1))


# unique species
unique(data2$SPECIES)
unique(data2$TRAPPINGSITE)
unique(data2$HOST.DNA.NEW)

#
A6<-ddply(data2,.(HOST.DNA.NEW),summarise,freq=unique(HOST.GROUP.HUMAN))
table(A6$freq)

#hosts per mosquito species
A7<-ddply(data2,.(SPECIES),summarise,freq=unique(HOST.DNA.NEW))
species_sub<-subset(A7,! A7$SPECIES %in% c("Ae./Oc. spec.","Cx. pipiens s.l./torrentium","Cx. spec.","Oc. cantans/annulipes")) 
ddply(species_sub,.(SPECIES),summarise,freq=length(freq))



nrow(data2)

sort(unique(data2$HOST.DNA.NEW))

#============================================================
#
#============================================================
AS<-count(data2, c("SPECIES", "HOST.DNA.NEW"))
AS2<-acast(AS, SPECIES~HOST.DNA.NEW, value.var="freq")
AS2[is.na(AS2)] <- 0
AS2<-ifelse(AS2>0,1,0)

# dft<-data.frame(rbind(c(1,0,1,0,1),c(1,1,1,1,1),c(1,1,1,0,1),c(0,1,1,1,0),c(1,0,1,1,1)))

myModel <- cooc_null_model(speciesData=AS2,algo="sim2",nReps=5000)
summary(myModel)

# Cx. pipiens
data_pipiens<-subset(data2,SPECIES=="Cx. pipiens pipiens")
unique_trappingsites_pipiens<-unique(data_pipiens$HOST.DNA.NEW)
unique_trappingsites_full<-unique(data2$HOST.DNA.NEW)
matching<-match(unique_trappingsites_pipiens,
                unique_trappingsites_full)
utz<-as.character(unique_trappingsites_full[-c(matching)])


AS<-count(data_pipiens, c("TRAPPINGSITE", "HOST.DNA.NEW"))
AS2<-acast(AS, TRAPPINGSITE~HOST.DNA.NEW, value.var="freq")
AS2[is.na(AS2)] <- 0
AS2<-ifelse(AS2>0,1,0)

df = data.frame(matrix(vector(), nrow(AS2), length(utz),
                       dimnames=list(c(), utz)),
                stringsAsFactors=F)
df[is.na(df)] <- 0

AS2_2<-cbind(AS2,df)

myModel <- cooc_null_model(speciesData=AS2_2,algo="sim2",nReps=5000)
summary(myModel)

# Oc. cantans
data_cantans<-subset(data2,SPECIES=="Oc. cantans")
unique_trappingsites_cantans<-unique(data_cantans$HOST.DNA.NEW)
unique_trappingsites_full<-unique(data2$HOST.DNA.NEW)
matching<-match(unique_trappingsites_cantans,
                unique_trappingsites_full)
utz<-as.character(unique_trappingsites_full[-c(matching)])

AS<-count(data_cantans, c("TRAPPINGSITE", "HOST.DNA.NEW"))
AS2<-acast(AS, TRAPPINGSITE~HOST.DNA.NEW, value.var="freq")
AS2[is.na(AS2)] <- 0
AS2<-ifelse(AS2>0,1,0)

df = data.frame(matrix(vector(), nrow(AS2), length(utz),
                       dimnames=list(c(), utz)),
                stringsAsFactors=F)
df[is.na(df)] <- 0

AS2_2<-cbind(AS2,df)

myModel <- cooc_null_model(speciesData=AS2_2,algo="sim2",nReps=5000)
summary(myModel)

# Ae, vexabns
data_vexans<-subset(data2,SPECIES=="Ae. vexans")
unique_trappingsites_vexans<-unique(data_vexans$HOST.DNA.NEW)
unique_trappingsites_full<-unique(data2$HOST.DNA.NEW)
matching<-match(unique_trappingsites_vexans,
                unique_trappingsites_full)
utz<-as.character(unique_trappingsites_full[-c(matching)])


AS<-count(data_cantans, c("TRAPPINGSITE", "HOST.DNA.NEW"))
AS2<-acast(AS, TRAPPINGSITE~HOST.DNA.NEW, value.var="freq")
AS2[is.na(AS2)] <- 0
AS2<-ifelse(AS2>0,1,0)

df = data.frame(matrix(vector(), nrow(AS2), length(utz),
                       dimnames=list(c(), utz)),
                stringsAsFactors=F)
df[is.na(df)] <- 0

AS2_2<-cbind(AS2,df)

myModel <- cooc_null_model(speciesData=AS2_2,algo="sim2",nReps=5000)
summary(myModel)

##
##
#
#================================================================================
# Heat map per trapping site for selected
#================================================================================

# sum data per mosquito species and detected host species
data_reduce<-subset(data2,SPECIES=="Ae. vexans"|SPECIES=="Oc. cantans"|SPECIES=="Cx. pipiens pipiens")

dataHeatMap<-ddply(data_reduce,.(SPECIES, TRAPPINGSITE, HOST.DNA.NEW),summarise,sum=sum(POOLSIZE))
dataHeatMap2<-ddply(data_reduce,.(SPECIES,TRAPPINGSITE),summarise,sum=sum(POOLSIZE))

zwei<-merge(dataHeatMap,dataHeatMap2,by=c("SPECIES","TRAPPINGSITE"))

zwei$sum=zwei$sum.x/zwei$sum.y*100


# sum of mosquitoes
#d_dypl2<-ddply(data_reduce,.(TRAPPINGSITE),summarise, freq= sum(POOLSIZE))

#newdata <-cbind(d_dypl2[with(d_dypl2, order(-freq)), ],num=seq(1,nrow(d_dypl2),1))


#newdata500<-merge(zwei,newdata,by="TRAPPINGSITE")
#newdata500$SPECIES<-reorder(newdata500$SPECIES,newdata500$num)

tztz<-data.frame(no=seq(1:31),
                 HOST.DNA.NEW=sort(unique(data_reduce$HOST.DNA.NEW),decreasing=T))
tztz_2<-merge(tztz,newdata500,by="HOST.DNA.NEW")


unique_hosts_reduced<-ddply(zwei,.(SPECIES,TRAPPINGSITE),summarise,HOST.DNA.NEW=unique(HOST.DNA.NEW))
unique_trappingsites_full<-ddply(data2,.(SPECIES,TRAPPINGSITE),
                                 summarise,HOST.DNA.NEW=unique(HOST.DNA.NEW))

unique_hosts_red<-ddply(data2,.(SPECIES,TRAPPINGSITE),summarise,HOST.DNA.NEW=unique(HOST.DNA.NEW))
uniqu<-unique(data_reduce$TRAPPINGSITE)

subben<-subset(unique_hosts_red, TRAPPINGSITE %in% c(as.character(uniqu)))
dongo<-ddply(subben,.(TRAPPINGSITE),summarise,HOST.DNA.NEW=unique(HOST.DNA.NEW))
dongo2<-merge(dongo,zwei,by=c("HOST.DNA.NEW","TRAPPINGSITE"))

sum_site_host<-ddply(data2,.(TRAPPINGSITE,HOST.DNA.NEW),summarise,sum=sum(POOLSIZE))

colnames(zwei)
ggg<-merge(sum_site_host,dongo2,by=c("HOST.DNA.NEW","TRAPPINGSITE"),all=T)
ggg[is.na(ggg)] <- 0
gg2<-subset(ggg,sum.x<1)
?merge
ggg1<-ggg[,c(1:3,7)]

dt5<-merge(ggg1,zwei,by=c("HOST.DNA.NEW","TRAPPINGSITE"))
colnames(dt5)

# plot heatmap
p <- ggplot(zwei, aes(x=TRAPPINGSITE,y=HOST.DNA.NEW   ))+ geom_tile(aes(fill=sum))+
  scale_fill_gradient(guide = guide_legend(title = "percentage per trapping site"),low="green", high="red")+
  #geom_point(data=gg2, aes(x=TRAPPINGSITE, y=HOST.DNA.NEW,size=sum.y.1,z=NULL))+
  theme_bw()+
  ylab("host species")+
  xlab("trapping site")+
  ylim(rev(levels(zwei$HOST.DNA.NEW)))+
  theme(axis.text.y = element_text(vjust = 0.25,hjust=1))  +
  theme(axis.text.x = element_blank())  +
  theme(axis.text =  element_text(face = "italic"))+
  theme(strip.text = element_text(face = "italic"))+
  theme(legend.position="bottom")+
  facet_wrap(~SPECIES,scales="free_x")

# save plot
png(file = "output/SpeciesSpec.jpg",width = 10, height=7, units = 'in', res = 1000)
plot(p)
dev.off()

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

# add the number
gg2 <- ddply(Diversity.merge, c("richness.x", "richness.y"), "nrow", .drop = T)

# spearman correlation between number of host species and the number of mosquito species
cor.test(Diversity.merge$richness.x,Diversity.merge$richness.y,method="spearman")

# plot
p<-ggplot(data = gg2, aes(x = richness.y, y = richness.x, size = factor(nrow))) + 
  geom_point() + 
  theme_classic() +
  ylab("number of host species") +
  xlab("number of mosquito species") +
  scale_size_discrete(range = c(1, 10),guide = guide_legend(title = "trapping sites"))

# save plot
png(file = "output/DiversityMosquitoesHosts.jpg",width = 6, height=5, units = 'in', res = 1000)
plot(p)
dev.off()

##############################
# correlation between the diversity of hosts per mosquito species and the diversity of detected hosts per site
##############################
NEU<-ddply(data2,.(TRAPPINGSITE,SPECIES),summarise,Hosts=length(unique(HOST.DNA.NEW)))
NEU2<-data.frame(TRAPPINGSITE=rownames(Diversity.hosts),richness=Diversity.hosts[,1])

Diversity.merge.2<-merge(NEU,NEU2,by="TRAPPINGSITE")

# sum of mosquitoes
d_dypl2<-ddply(Diversity.merge.2,.(SPECIES),summarise, freq= length(Hosts))
newdata5<-merge(Diversity.merge.2,d_dypl2,by="SPECIES")
BILLS<-subset(newdata5,newdata5$freq>4)

ddply(BILLS, .(SPECIES), summarise,
      corr=(cor.test(Hosts, richness,
                     method="spearman")$p.value))

p<-ggplot(data = Diversity.merge.2, aes(x = richness, y = Hosts, colour = SPECIES)) + 
  geom_point() + 
  theme_classic() +
  facet_wrap(~SPECIES)+
  ylab("number of host species") +
  xlab("number of mosquito species") +
  scale_size_discrete(range = c(1, 10),guide = guide_legend(title = "trapping sites"))

# save plot
#png(file = "output/DiversityMosquitoesHosts.jpg",width = 6, height=5, units = 'in', res = 1000)
#plot(p)
#dev.off()
  
##############################
# host group preferences of the detected mosquito species
##############################

# sum of mosquitoes
d_dypl2<-ddply(data2,.(SPECIES),summarise, freq= sum(POOLSIZE))


newdata <-cbind(d_dypl2[with(d_dypl2, order(-freq)), ],num=seq(1,nrow(d_dypl2),1))

# merge the information of the 
newdata5<-merge(data2,newdata,by="SPECIES")

# reorder the species names
newdata5$SPECIES<-reorder(newdata5$SPECIES,newdata5$num)

# plot
p2<-ggplot(newdata5,aes(x=SPECIES)) +
  geom_bar(colour= "black")+
  theme_classic()+
  xlab("")+
  ylab("frequency")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0,hjust=1))+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  scale_y_continuous(expand=c(0,0))+
  #theme(axis.text =  element_text(face = "italic"))+
  theme(axis.text.x =  element_blank())
  


png(file = "output/FeedingPreferenceMosquitoes.jpg",width = 6, height=5, units = 'in', res = 1000)
plot(p2)
dev.off()

##############################
# Heat map on the detected hosts for each species
##############################

# sum data per mosquito species and detected host species
dataHeatMap<-ddply(data2,.(SPECIES, HOST.DNA.NEW),summarise,sum=sum(POOLSIZE))
dataHeatMap2<-ddply(data2,.(SPECIES),summarise,sum=sum(POOLSIZE))

zwei<-merge(dataHeatMap,dataHeatMap2,by="SPECIES")
zwei$sum=zwei$sum.x/zwei$sum.y*100

newdata500<-merge(zwei,newdata,by="SPECIES")
newdata500$SPECIES<-reorder(newdata500$SPECIES,newdata500$num)

tztz<-data.frame(no=seq(1:33),
                 HOST.DNA.NEW=sort(unique(newdata500$HOST.DNA.NEW),decreasing=T))
tztz_2<-merge(tztz,newdata500,by="HOST.DNA.NEW")

nonHmam<-c("Bos taurus","Canis lupus familiaris","Capra hircus",
           "Capreolus capreolus","Castor fiber","Cervus elaphus",
           "Equus caballus","Erinaceus europaeus","Felis catus",
           "Lama pacos","Lepus europaeus","Microtus guentheri",
           "Mustela/Lutrea spec.","Myocastor coypus","Myodes glareolus",
           "Myotis spec.","Oryctolagus cuniculus","Ovis aries",
           "Rattus norvegicus","Sus scrofa","Vulpes vulpes",
           "Homo sapiens",
           "Anas platyrhynchus","Corvus corone","Cyanistes caeruleus",
           "Delichon urbica","Emberiza citrinella","Erithacus rubecula",
           "Passer domesticus","Sylvia atricapilla","Sylvia communis",
           "Turdus merula","Turdus philomelos")
tztz_2$HOST.DNA.NEW <- factor(tztz_2$HOST.DNA.NEW, levels = nonHmam)
tztz_2$HOST.DNA.NEW  # notice the changed order of factor levels

# plot heatmap
p <- ggplot(tztz_2, aes(x=SPECIES,y=(HOST.DNA.NEW)))+ geom_tile(aes(fill=sum))+
  scale_fill_gradient(guide = guide_legend(title = "percentage per mosquito species"),low="green", high="red")+
  #geom_point(data=zwei, aes(x=SPECIES, y=HOST.DNA.NEW,z=NULL,size=sum.x,colour=sum))+
  theme_bw()+
  ylab("Host species")+
  xlab("Mosquito species")+
  ylim(rev(levels(tztz_2$HOST.DNA.NEW)))+
  theme(axis.text.y = element_text(vjust = 0.25,hjust=1))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust=1))  +
  theme(axis.text =  element_text(face = "italic"))+
  theme(legend.position="bottom")

# save plot
png(file = "output/HeatMapFeedingPreference.jpg",width = 8, height=7, units = 'in', res = 1000)
plot(p)
dev.off()
  
library(ggplot2)
gp1<- ggplot_gtable(ggplot_build(p+theme(legend.position="none")))
gp2<- ggplot_gtable(ggplot_build(p2+theme(legend.position="none")))
#gp3<- ggplot_gtable(ggplot_build(p3+theme(legend.position="none")))

library(grid)
maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth
#gp3$widths[2:3] <- maxWidth



tmp <- ggplot_gtable(ggplot_build(p)) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 

tmp <- ggplot_gtable(ggplot_build(p2)) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend2 <- tmp$grobs[[leg]] 

#tmp <- ggplot_gtable(ggplot_build(p3)) 
#leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
#legend3 <- tmp$grobs[[leg]] 

library(gridExtra)
#gul1<-arrangeGrob(gp1, 
#                  legend,
#                  ncol =2,
#                  widths=c(1),
#                  widths=c(1.5,0.5))
gul1<-arrangeGrob(gp2,
                  ncol =1,
                  widths=c(1))

#?arrangeGrob
gul2<-arrangeGrob(gp1,
                  legend,
                  nrow=2,
                  widths=c(1),
                  heights=c(1.9,0.1))
plot(gul2)
#gul2<-arrangeGrob(gp2,
#                  ncol=2,
#                  widths=c(1.5,0.5))

#gul3<-arrangeGrob(gp3, 
 #                 legend3,
  #                ncol=2,
   #               widths=c(1.8,0.2))



png("Fig1.png",width = 6, height=10, units = 'in', res = 600)
grid.arrange(gul1,gul2,heights=c(2/10, 8/10),nrow=2)
dev.off()

##############################
# Heat map on the detected hosts for each species
##############################
AS<-count(data2, c("SPECIES", "HOST.DNA.NEW"))


AS<-count(data2, c("SPECIES", "HOST.GROUP.HUMAN"))

AS2<-acast(AS, SPECIES~HOST.DNA.NEW, value.var="freq")


AS2<-acast(AS, SPECIES~HOST.GROUP.HUMAN, value.var="freq")

AS2[is.na(AS2)] <- 0

AS2<-ifelse(AS2>0,1,0)




env.norm <- decostand(AS2, method = "normalize", na.rm = TRUE, MARGIN=2)
spe.ch<-vegdist(env.norm,method="jaccard")
hc <- hclust(dist(spe.ch), "average")
dhc <- as.dendrogram(hc,hang=0.1)
ddata <- dendro_data(dhc, type="rectangle")

clust    <- cutree(hc,k=5)                    # find 2 clusters
clust.df <- data.frame(label=names(clust)
                       , cluster=factor(clust))
ddata[["labels"]] <- merge(ddata[["labels"]],clust.df, by="label")


?revalue
#ddata$labels$cluster<-revalue(ddata$labels$cluster,
 #                             c("1"="5","2"="4","3"="3","4"="2","5"="1"))

#ddata$labels$label<-revalue(ddata$labels$label,
#                              c("Cx. pipiens/ molestus/torrentium"="Cx. pipiens s.l./torrentium"))
p<-ggplot() + 
  geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=ddata$labels, aes(x, y,label=(label), hjust=0, color=cluster), 
            size=5,fontface=3) +
  ylab("distance")+
  coord_flip() + 
  scale_y_reverse(expand=c(0.2, 0)) + 
  theme_bw()+
  theme(text = element_text(size=20),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x =  element_text(face = "italic"),
        #text =  element_text(face = "italic"),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank(),
        legend.position = "none",
        legend.position="bottom")

png(file = "output/FeedingPreferenceMosquitoesDendgrogramMosqus.jpg",
    width = 16.5, height=8, units = 'in', res = 1000)
plot(p)
dev.off()

dimnames(ddata$labels)[[2]]<-c("SPECIES","x","y","cluster")
data2$SPECIES<-revalue(data2$SPECIES,
                            c("Cx. pipiens/ molestus/torrentium"="Cx. pipiens s.l./torrentium"))

OINK<-merge(data2,ddata$labels,by="SPECIES")

#p2<-ggplot(OINK,aes(x=SPECIES,fill=HOST.GROUP.HUMAN)) +
#  ylab("percentage")+
 # scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
 # geom_bar(colour="black",position="fill")  +theme_classic()+  scale_y_continuous(expand=c(0,0))+
 # facet_wrap(~cluster)

#colnames(OINK)
#p2<-ggplot(OINK,aes(x=cluster,fill=HOST.GROUP.HUMAN)) +
 # ylab("count")+
 # scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
 # geom_bar(colour="black")  +theme_classic()+  scale_y_continuous(expand=c(0,0))

#library('scales')

OINK7<-ddply(OINK,.(SPECIES),summarise,d3=length(HOST.GROUP.HUMAN))
OINK8<-merge(OINK,OINK7,by="SPECIES")
OINK9<-ddply(OINK8,.(SPECIES,HOST.GROUP.HUMAN,cluster),summarise,d4=length(HOST.GROUP.HUMAN)/d3*100)
OINK1<-unique(OINK9[c("SPECIES", "HOST.GROUP.HUMAN",
                      "cluster","d4")])

OINK1$cluster<-revalue(OINK1$cluster,
                             c("1"="4","2"="3","3"="1","4"="5","5"="2"))

neworder <- c("1","2","3","4","5")
OINK1 <- arrange(transform(OINK1,
                           cluster=factor(cluster,levels=neworder)),cluster)

p3<-ggplot(OINK1, aes(x = SPECIES,y=d4,fill=HOST.GROUP.HUMAN)) +
  geom_bar(colour="black",position = 'stack',stat="identity")+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  facet_wrap(~cluster,scales="free")+
  ylab("percentage")+
  xlab("mosquito species")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,face = "italic", vjust = 0,hjust=1))+
  scale_y_continuous(expand=c(0,0))

png(file = "output/HostpreferenceDendgromClusterSpecies.jpg",
    width = 8, height=8, units = 'in', res = 1000)
plot(p3)
dev.off()

















































  ####################
####################
####################
# 
colnames(data2)
AAA<-ddply(data2,.(PREVALENCE,HOST.GROUP),summarise,sum=length(POOLSIZE))

AAA<-ddply(data2,.(TRAPPINGSITE,HOST.GROUP.HUMAN),summarise,sum=length(POOLSIZE))


A<-c(41,6)
B<-c(358,111)
data.table<-rbind(A,B)
chisq.test(data.table)

AAA<-ddply(data2,.(SPECIES,TRAPPINGSITE),summarise,sum=length(POOLSIZE))
sum(AAA$sum)
# dfdf
AAAA<-ddply(AAA,.(SPECIES),summarise,sum=length(sum))
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

plot(p)

png(file = "output/FeedingPreferenceMosquitoesDendgrogramSites.jpg",width = 8, height=6, units = 'in', res = 1000)
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


png(file = "output/HostpreferenceDendgromClusterSites.jpg",width = 8, height=6, units = 'in', res = 1000)
plot(p2)
dev.off()

#==================================================
# spatial differences in host preference?
#==================================================
# read data
data<-read.table (file = "data/Bloodmeals_FINAL.txt",row.names=1,header=TRUE,sep="\t")

# omit NA values (e.g. not detected host,...)
data2<-na.omit(data)

# sum of mosquitoes
d_dypl2<-ddply(data2,.(SPECIES),summarise, freq= sum(POOLSIZE))

newdata <-cbind(d_dypl2[with(d_dypl2, order(-freq)), ],num=seq(1,nrow(d_dypl2),1))

# merge the information of the 
newdata5<-merge(data2,newdata,by="SPECIES")

# reorder the species names
newdata5$SPECIES<-reorder(newdata5$SPECIES,newdata5$num)

newdata4<-subset(newdata5,SPECIES=="Ae. vexans"|SPECIES=="Cx. pipiens pipiens"|SPECIES=="Oc. cantans")
unique(newdata4$TRAPPINGSITE)
unique(newdata4$SPECIES)


#


# plot
p24<-ggplot(newdata4,aes(x=TRAPPINGSITE,fill=HOST.GROUP.HUMAN)) +
  geom_bar(colour= "black")+
  theme_classic()+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90,face = "italic", vjust = 0,hjust=1))+
  scale_fill_manual(guide = guide_legend(title = "host group"),values=c("black","white","gray"))+
  scale_y_continuous(expand=c(0,0))+
  #theme(axis.text =  element_text(face = "italic"))+
  facet_wrap(~SPECIES,scales="free",ncol=1)

png(file = "output/NEUs.jpg",width = 6, height=12, units = 'in', res = 1000)
plot(p24)
dev.off()

#newdata4[,7]
#data3<-read.table (file = "data/Bloodmeals_FINAL.txt",row.names=1,header=TRUE,sep="\t")
#data4<-na.omit(data3)
#XY<-read.table (file = "data/XYNEW.txt",row.names=1,header=TRUE,sep="\t")
XY2<-data.frame(x = data_om$x, y = data_om$y)

unique(data$TRAPPINGSITE)
dft <- merge(XY2, data, by="TRAPPINGSITE")
sort(unique(data$TRAPPINGSITE))

library(raster)
r <- raster("data/df.tif") # CORINE land cover for Germany
isBecomes <- cbind(c(seq(1,44),48,49,50,255),
                   (c(rep(1,11), # urban
                      rep(2,1),  # acre
                      rep(2,10),  # rural
                      rep(3,17),  # near natural
                      rep(4,5),  # water bodies
                      5,5,5,5))) 
rCORINEagg <- reclassify(r , rcl=isBecomes)

plot(rCORINEagg)

VALUES_Land <- extract(rCORINEagg, XY2,
                     buffer=2500)
toS<-lapply(VALUES_Land,table)

rbind.named.fill <- function(x) {
  nam <- sapply(x, names)
  unam <- unique(unlist(nam))
  len <- sapply(x, length)
  out <- vector("list", length(len))
  for (i in seq_along(len)) {
    out[[i]] <- unname(x[[i]])[match(unam, nam[[i]])]
  }
  setNames(as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE), unam)
}
toS2<-rbind.named.fill(toS)
toS2[is.na(toS2)] <- 0

VALUES<-toS2/rowSums(toS2)*100
VALUE_FULL <- cbind(VALUES,XY2)
VALUE_FULL$IDsb <- VALUE_FULL$x + VALUE_FULL$y

dataxx<-subset(data_om, SPECIES=="Cx. pipiens pipiens")
dataxx$IDsb <- dataxx$x + dataxx$y

#WERTE <- merge(data_om,VALUE_FULL, by="TRAPPINGSITE")
trap_method_1 <- ddply(dataxx,.(IDsb),
                       summarise,d3=length(HOST.GROUP.HUMAN))
trap_method_2 <- merge(trap_method_1, dataxx, by="IDsb")
trap_method_3 <- ddply(trap_method_2,.(IDsb, HOST.GROUP.HUMAN),summarise,
                       d4=length(HOST.GROUP.HUMAN)/d3*100)
trap_method_4 <- unique(trap_method_3[c("IDsb","HOST.GROUP.HUMAN","d4")])
trap_method_5 <- merge(trap_method_4, trap_method_1, by="IDsb")

trtdd <- VALUE_FULL[!duplicated(VALUE_FULL[6:7]),]

WERTE <- merge(trap_method_4, trtdd, by = "IDsb")
human <- subset(WERTE, HOST.GROUP.HUMAN == "bird")
plot(human[,8], human$d4)

plot(human[,6], human$d4)


PCFIDL<-prcomp(toS2)$x[,1:2]
XY3<-cbind(XY2,PCFIDL)

#toS2$SUM<-rowSums(toS2)
#toS3<-toS2[,1:5]/toS2$SUM*100
#XY3<-cbind(XY2,toS3)

# percentage per trapping site and species
a = ddply(newdata4, .(TRAPPINGSITE,SPECIES), function(d) {
  data.frame(table(d$HOST.GROUP.HUMAN)/length(d$HOST.GROUP.HUMAN))
})

MERGEN4<-merge(a,XY3,by="TRAPPINGSITE")
MERGEN5<-subset(MERGEN4,SPECIES=="Ae. vexans"|SPECIES=="Cx. pipiens pipiens"|SPECIES=="Oc. cantans")
MERGEN6<-subset(MERGEN4,Var1=="bird" &
                  SPECIES=="Cx. pipiens pipiens")
plot(MERGEN6[,8],MERGEN6$Freq)
summary(glm(MERGEN6$Freq~MERGEN6[,7]))

cor.test(MERGEN6[,7],MERGEN6$Freq)

SUMME<-ddply(MERGEN5,.(SPECIES,TRAPPINGSITE,HOST.GROUP.HUMAN),
             summarise,T=length(HOST.GROUP.HUMAN),.drop = F)

MERGE1000<-merge(SUMME,XY3,by="TRAPPINGSITE")
MERGE1000$T<-ifelse(MERGE1000$T>0,1,0)

MERGE10007<-na.omit(MERGE1000)

MERGE10000<-subset(MERGE1000,SPECIES=="Cx. pipiens pipiens"&HOST.GROUP.HUMAN=="bird")
plot(MERGE10000[,7],MERGE10000[,4])
#AS<-count(data2, c("SPECIES", "HOST.DNA.NEW"))

AS<-count(MERGEN5, c("TRAPPINGSITE", "HOST.GROUP.HUMAN"))

#AS2<-acast(AS, SPECIES~HOST.DNA.NEW, value.var="freq")


AS2<-acast(AS, TRAPPINGSITE~HOST.GROUP.HUMAN, value.var="freq")

AS2[is.na(AS2)] <- 0

AS2<-ifelse(AS2>0,1,0)

env.norm <- decostand(AS2, method = "normalize", na.rm = TRUE, MARGIN=2)
spe.ch<-vegdist(env.norm,method="jaccard")
hc <- hclust(dist(spe.ch), "average")
dhc <- as.dendrogram(hc,hang=0.1)
ddata <- dendro_data(dhc, type="rectangle")

clust    <- cutree(hc,k=6)                    # find 2 clusters
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

#MERGEN5$cluster
#dfdf<-subset(MERGEN5,cluster==3)
#sum(table(dfdf$TRAPPINGSITE))

colnames(MERGEN5)

dimnames(ddata$labels)[[2]]<-c("TRAPPINGSITE","x","y","cluster")
OINK<-merge(MERGEN5,ddata$labels,by="TRAPPINGSITE")

ddply(OINK,.(cluster,SPECIES,HOST.GROUP.HUMAN),summarise,length(HOST.GROUP.HUMAN))
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


png(file = "output/DIFF.jpg",width = 14, height=6, units = 'in', res = 1000)
plot(p)
dev.off()


##############################
# Mesonivirus detected in mosquitoes linked to the blood meals detected
##############################

##### 
data3<-subset(data2,HOST.GROUP=="mammal")
d_dypl<-ddply(data3,.(HOST.DNA.NEW),summarise, freq= sum(POOLSIZE))
newdata <-cbind(d_dypl[with(d_dypl, order(-freq)), ],num=seq(1,nrow(d_dypl),1))
newdata3<-merge(data3,newdata,by="HOST.DNA.NEW")
newdata3$HOST.DNA.NEW<-reorder(newdata3$HOST.DNA.NEW,newdata3$num)

data4<-subset(data2,HOST.GROUP=="bird")
d_dypl2<-ddply(data4,.(HOST.DNA.NEW),summarise, freq= sum(POOLSIZE))
newdata <-cbind(d_dypl2[with(d_dypl2, order(-freq)), ],num=seq(1,nrow(d_dypl2),1))
newdata4<-merge(data4,newdata,by="HOST.DNA.NEW")
newdata4$HOST.DNA.NEW<-reorder(newdata4$HOST.DNA.NEW,newdata4$num)

p2<-ggplot(data2, aes(x = HOST.DNA.NEW,fill=as.factor(Mesoni))) + 
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
