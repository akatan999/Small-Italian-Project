########
library(Hmisc)
library(mgcv)
library(ggplot2)
library(reshape)
library(plyr)
library(dplyr)
library(ggpubr)

output.dir <- "~/Max/PPI/Publications/Figures"

##All individuals
BCPDPonzaeffort <- read.table("~/Max/PPI/Publications/Time and trends/BCPDPonzaeffort.txt", header = TRUE)

###Only juveniles (optional)
BCPDPonzaeffort_age <- read.table("~/Max/PPI/Publications/Time and trends/BCPDPonzaeffort_age.txt", header = TRUE)

##########Transforming at 5-days time period all
BCPDPonzaeffort$Days <- cut2(BCPDPonzaeffort$Days, seq(70,150,5))
BCPDPonzaeffort <- ddply(BCPDPonzaeffort, c("Days","Year","EURING"), summarise, mean = mean(CPUE))

##########Transforming at 5-days time period (age)
BCPDPonzaeffort_age$Days <- cut2(BCPDPonzaeffort_age$Days, seq(70,150,5))
BCPDPonzaeffort_age <- ddply(BCPDPonzaeffort_age, c("Days","Year","EURING","ETA"), summarise, mean = mean(CPUE))

##Days and EURING all
BCPDPonzaeffort[, c(1)] <- sapply(BCPDPonzaeffort[, c(1)], as.numeric)
BCPDPonzaeffort[, c(3)] <- sapply(BCPDPonzaeffort[, c(3)], as.factor)

##Days and EURING age
BCPDPonzaeffort_age[, c(1)] <- sapply(BCPDPonzaeffort_age[, c(1)], as.numeric)
BCPDPonzaeffort_age[, c(3)] <- sapply(BCPDPonzaeffort_age[, c(3)], as.factor)

#Assign names to the columns all
colnames(BCPDPonzaeffort)[4]<-'CPUE'
colnames(BCPDPonzaeffort)[1]<-'Days'

#Assign names to the columns age
colnames(BCPDPonzaeffort_age)[5]<-'CPUE'
colnames(BCPDPonzaeffort_age)[1]<-'Days'

##Save the data
write.table(BCPDPonzaeffort,file="BCPDPonzaeffort_all.txt",sep="\t",na="0",row.names=F,col.names=T)
write.table(BCPDPonzaeffort_age,file="BCPDPonzaeffort_age.txt",sep="\t",na="0",row.names=F,col.names=T)

###Eliminate Saxicola torquata & Oenanthe hispanica all
Eur <- levels(factor(BCPDPonzaeffort$EURING))
BCPDPonzaeffort <- BCPDPonzaeffort[BCPDPonzaeffort$EURING!="11390"&BCPDPonzaeffort$EURING!="11480",]

###Eliminate Saxicola torquata & Oenanthe hispanica & Fingilla for age analysis
Eur <- levels(factor(BCPDPonzaeffort_age$EURING))
BCPDPonzaeffort_age <- BCPDPonzaeffort_age[BCPDPonzaeffort_age$EURING!="11390"&BCPDPonzaeffort_age$EURING!="11480"&BCPDPonzaeffort_age$EURING!="16360",]

####CPUE models all
Eur <- levels(factor(BCPDPonzaeffort$EURING))
new.modsel <- list()
new.lam.pred <- list()
mod <- list()

for(i in 1:length(Eur)){

new.modsel[[i]] <- expand.grid(Days=(mean(BCPDPonzaeffort$Days[BCPDPonzaeffort$EURING==Eur[i]])),Year=c(2004,2005,2007:2024))

mod[[i]]<-gam((CPUE*100)~ (as.factor(Year):Days),data=BCPDPonzaeffort[BCPDPonzaeffort$Year>2003 & BCPDPonzaeffort$EURING==Eur[i],], family=tw(), method="REML")

new.lam.pred[[i]] <-predict(mod[[i]], new.modsel[[i]], type="response", se=TRUE)  
new.lam.pred[[i]] <-cbind(new.modsel[[i]],new.lam.pred[[i]])
new.lam.pred[[i]]$EURING <- Eur[i]

assign(paste("mod",  Eur[i] , sep="_"), mod[[i]])
assign(paste("new.modsel",  Eur[i] , sep="_"), new.modsel[[i]])
assign(paste("new.lam.pred",  Eur[i] , sep="_"), new.lam.pred[[i]])
  
}

###Merge the results in a dataframe
Results_mergedr <- do.call(rbind,new.lam.pred)

#Matching EURING CODE with the species 
euring <- read.csv("EURING_Code.csv")
Results_mergedr$Species <- euring$Species[match(Results_mergedr$EURING,euring$EURING)] 

##Plotting from 2005 onwards and up to 2024
Results_mergedr2 <- Results_mergedr[Results_mergedr$Year<2025&Results_mergedr$Year>2003,]
write.table(Results_mergedr2 ,file="Results_mergedr2.txt",sep="\t",na="0",row.names=F,col.names=T)
###########


####CPUE models age (juv)
BCPDPonzaeffort_age_juv <- BCPDPonzaeffort_age[BCPDPonzaeffort_age$ETA=="5",]
Eur <- levels(factor(BCPDPonzaeffort_age_juv$EURING))

new.modsel <- list()
new.lam.pred <- list()
mod <- list()

for(i in 1:length(Eur)){
  
  new.modsel[[i]] <- expand.grid(Days=(mean(BCPDPonzaeffort_age_juv$Days[BCPDPonzaeffort_age_juv$EURING==Eur[i]])),Year=c(2004,2005,2007:2024))
  
  mod[[i]]<-gam((CPUE*100)~ (as.factor(Year):Days),data=BCPDPonzaeffort_age_juv[BCPDPonzaeffort_age_juv$Year>2003 & BCPDPonzaeffort_age_juv$EURING==Eur[i],], family=tw(), method="REML")
  
  new.lam.pred[[i]] <-predict(mod[[i]], new.modsel[[i]], type="response", se=TRUE)  
  new.lam.pred[[i]] <-cbind(new.modsel[[i]],new.lam.pred[[i]])
  new.lam.pred[[i]]$EURING <- Eur[i]
  
  assign(paste("mod",  Eur[i] , sep="_"), mod[[i]])
  assign(paste("new.modsel",  Eur[i] , sep="_"), new.modsel[[i]])
  assign(paste("new.lam.pred",  Eur[i] , sep="_"), new.lam.pred[[i]])
  
}

###Merge the results in a dataframe
Results_mergedr <- do.call(rbind,new.lam.pred)

#Matching EURING CODE with the species 
euring <- read.csv("EURING_Code.csv")
Results_mergedr$Species <- euring$Species[match(Results_mergedr$EURING,euring$EURING)] 

##Plotting from 2005 onwards and up to 2024
Results_mergedr2 <- Results_mergedr[Results_mergedr$Year<2025&Results_mergedr$Year>2003,]

write.table(Results_mergedr2 ,file="Results_mergedr2_juv.txt",sep="\t",na="0",row.names=F,col.names=T)

######################################

#########Plotting
Results_mergedr2 <- read.table("~/Max/PPI/Publications/Time and trends/Results_mergedr2.txt", header = TRUE)

###Plotting the results (all species)
pdf("Ponza trends_smoothed_2025.pdf", width=10, height=8)

p <- ggplot(Results_mergedr2, aes(Year,fit))+geom_line()+facet_wrap(~Species, scales="free_y", labeller=label_wrap_gen(10)) + ggtitle("Trend Ponza") + labs(x="Year",y="CPUE (n/m*h)") + theme_bw() + theme(strip.text = element_text(face = "italic")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9))

p <- p+ geom_smooth(data =Results_mergedr2, aes(Year,fit),se=T, method="loess",show.legend = FALSE,lwd=0.7) 
p
dev.off()

###Plotting the results (all species)_linear model
pdf("Ponza trends_smoothed_2025_linear.pdf", width=10, height=8)

p <- ggplot(Results_mergedr2, aes(Year,fit))+geom_point()+facet_wrap(~Species, scales="free_y", labeller=label_wrap_gen(10)) + ggtitle("Trend Ponza") + labs(x="Year",y="CPUE (n/m*h)") + theme_bw() + theme(strip.text = element_text(face = "italic")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9))

p <- p+ geom_smooth(data =Results_mergedr2, aes(Year,fit),se=F, method="lm",show.legend = FALSE,lwd=0.7)
p + stat_regline_equation(label.x.npc =  "left", aes(label = ..eq.label..)) +
  stat_regline_equation(label.x.npc =  "top", aes(label = ..rr.label..))
dev.off()

###Extract data from the reduced species ggplot for correlation analysis (done)
pg <- ggplot_build(p)
pg1 <- pg$data[[2]]
write.table(pg1,file="pg1.txt",sep="\t",na="0",row.names=F,col.names=T)
######Needs some trick with the smoothed data to be done manually to produce smoothed data_new.txt file (done and saved)

######################################
##Create a database of the EBCC data
birds <-read.csv("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Database/RegionalIndicesForMCardinale2018_new.csv", header=T, sep=";")
birds2 <- do.call("rbind", replicate(51, birds[,c(1:3,55:57)], simplify = FALSE))
birds2$years <- rep(1966:2016, each=181)
values <- as.vector(as.matrix(birds[,4:54]))
birds2$values <- values
colnames(birds2)[7]<-'Year'
colnames(birds2)[1]<-'EURING'

####Select Central & East europe for defined species
birds3 <- birds2[birds2$CountryGroup=="East Europe",]
birds4 <- birds3[birds3$Year>2004,]
birds5 <- birds4[birds4$Species=="Sylvia communis"|birds4$Species=="Oenanthe oenanthe"|birds4$Species=="Sylvia borin"|birds4$Species=="Anthus trivialis"|birds4$Species=="Ficedula albicollis"|birds4$Species=="Jynx torquilla"|birds4$Species=="Saxicola rubetra"|birds4$Species=="Turdus philomelos"|birds4$Species=="Phoenicurus phoenicurus"|birds4$Species=="Erithacus rubecula"|birds4$Species=="Sylvia atricapilla"|birds4$Species=="Upupa epops"|birds4$Species=="Oriolus oriolus"|birds4$Species=="Luscinia megarhynchos"|birds4$Species=="Muscicapa striata"|birds4$Species=="Phylloscopus trochilus"|birds4$Species=="Phylloscopus collybita"|birds4$Species=="Hirundo rustica"|birds4$Species=="Hippolais icterina"|birds4$Species=="Ficedula hypoleuca"|birds4$Species=="Phylloscopus sibilatrix"|birds4$Species=="Merops apiaster"|birds4$Species=="Streptopelia turtur"|birds4$Species=="Acrocephalus arundinaceus"|birds4$Species=="Acrocephalus schoenobaenus"|birds4$Species=="Fringilla coelebs"|birds4$Species=="Phoenicurus ochruros"|birds4$Species=="Delichon urbica",]

####Select North europe for defined species
birds31 <- birds2[birds2$CountryGroup=="North Europe",]
birds41 <- birds31[birds31$Year>2004,]
birds51 <- birds41[birds41$Species=="Sylvia communis"|birds41$Species=="Phylloscopus trochilus"|birds41$Species=="Hippolais icterina"|birds41$Species=="Phylloscopus sibilatrix"|birds41$Species=="Oriolus oriolus"|birds41$Species=="Phylloscopus collybita"|birds41$Species=="Sylvia borin"|birds41$Species=="Sylvia atricapilla"|birds41$Species=="Phoenicurus phoenicurus"|birds41$Species=="Ficedula hypoleuca"|birds41$Species=="Erithacus rubecula"|birds41$Species=="Anthus trivialis"|birds41$Species=="Oenanthe oenanthe"|birds41$Species=="Jynx torquilla"|birds41$Species=="Saxicola rubetra"|birds41$Species=="Turdus philomelos"|birds41$Species=="Upupa epops"|birds41$Species=="Oriolus oriolus"|birds41$Species=="Luscinia megarhynchos"|birds41$Species=="Muscicapa striata"|birds41$Species=="Acrocephalus schoenobaenus"|birds41$Species=="Fringilla coelebs"|birds41$Species=="Phoenicurus ochruros"|birds41$Species=="Delichon urbica",]

####Select West europe for defined species
birds311 <- birds2[birds2$CountryGroup=="West Europe",]
birds411 <- birds311[birds311$Year>2004,]
birds511 <- birds411[birds411$Species=="Ficedula hypoleuca"|birds411$Species=="Hirundo rustica"|birds411$Species=="Muscicapa striata"|birds411$Species=="Sylvia communis"|birds411$Species=="Sylvia borin"|birds411$Species=="Sylvia atricapilla"|birds411$Species=="Phylloscopus trochilus"|birds411$Species=="Phylloscopus collybita"|birds411$Species=="Phoenicurus phoenicurus"|birds411$Species=="Hirundo rustica"|birds411$Species=="Hippolais icterina"|birds411$Species=="Ficedula albicollis"|birds411$Species=="Erithacus rubecula"|birds411$Species=="Anthus trivialis"|birds411$Species=="Oenanthe oenanthe"|birds411$Species=="Jynx torquilla"|birds411$Species=="Turdus philomelos"|birds411$Species=="Upupa epops"|birds411$Species=="Oriolus oriolus"|birds411$Species=="Luscinia megarhynchos"|birds411$Species=="Muscicapa striata"|birds411$Species=="Phylloscopus sibilatrix"|birds411$Species=="Streptopelia turtur"|birds411$Species=="Acrocephalus arundinaceus"|birds411$Species=="Acrocephalus schoenobaenus"|birds411$Species=="Fringilla coelebs"|birds411$Species=="Phoenicurus ochruros"|birds411$Species=="Delichon urbica",]

####Select west Balkan for defined species
birds61 <- birds2[birds2$CountryGroup=="SouthEast Europe",]
birds71 <- birds61[birds61$Year>2004,]
birds81 <- birds71[birds71$Species=="Lanius senator"|birds71$Species=="Sylvia borin"|birds71$Species=="Sylvia atricapilla"|birds71$Species=="Hirundo rustica"|birds71$Species=="Hippolais icterina"|birds71$Species=="Ficedula hypoleuca"|birds71$Species=="Ficedula albicollis"|birds71$Species=="Erithacus rubecula"|birds71$Species=="Anthus trivialis"|birds71$Species=="Oenanthe oenanthe"|birds71$Species=="Jynx torquilla"|birds71$Species=="Saxicola rubetra"|birds71$Species=="Turdus philomelos"|birds71$Species=="Upupa epops"|birds71$Species=="Oriolus oriolus"|birds71$Species=="Luscinia megarhynchos"|birds71$Species=="Muscicapa striata"|birds71$Species=="Phylloscopus sibilatrix"|birds71$Species=="Oenanthe hispanica"|birds71$Species=="Merops apiaster"|birds71$Species=="Streptopelia turtur"|birds71$Species=="Fringilla coelebs"|birds71$Species=="Phoenicurus ochruros"|birds71$Species=="Delichon urbica",]

####Select south for defined species
birds6 <- birds2[birds2$CountryGroup=="South Europe",]
birds7 <- birds6[birds6$Year>2004,]
birds8 <- birds7[birds7$Species=="Sylvia communis"|birds7$Species=="Sylvia cantillans"|birds7$Species=="Oenanthe hispanica"|birds7$Species=="Merops apiaster"|birds7$Species=="Streptopelia turtur"|birds7$Species=="Hirundo rustica"|birds7$Species=="Lanius senator"|birds7$Species=="Sylvia borin"|birds7$Species=="Sylvia atricapilla"|birds7$Species=="Ficedula hypoleuca"|birds7$Species=="Anthus trivialis"|birds7$Species=="Oenanthe oenanthe"|birds7$Species=="Jynx torquilla"|birds7$Species=="Saxicola rubetra"|birds7$Species=="Turdus philomelos"|birds7$Species=="Upupa epops"|birds7$Species=="Oriolus oriolus"|birds7$Species=="Luscinia megarhynchos"|birds7$Species=="Muscicapa striata"|birds7$Species=="Phylloscopus sibilatrix"|birds7$Species=="Acrocephalus schoenobaenus"|birds7$Species=="Fringilla coelebs"|birds7$Species=="Phoenicurus ochruros"|birds7$Species=="Delichon urbica",]

###Merge the results in a dataframe
Birds11 <- rbind(birds5, birds51, birds511, birds8, birds81)

##Select only complete cases
Birds12 <- Birds11[complete.cases(Birds11),]

#Birds13 <- Birds12 %>%
#  group_by(Species, CountryGroup) %>%
#  mutate(scaledValues=scale(values, center = TRUE, scale=TRUE))

##
pdf("EBCC_all regions.pdf")
#birds2b<-birds2 %>%
#  group_by(Species, CountryGroup) %>%
#  mutate(scaledValues=scale(values, center = TRUE,scale=TRUE))

birds2b <- birds2[birds2$Year<2017,]
birds2b <- birds2b[birds2$Year>2004,]
birds2b <- birds2b[birds2b$CountryGroup!="WestSouth Europe"&birds2b$CountryGroup!="East Mediterranean",]

p <- ggplot(birds2b, aes(Year,values))+geom_line(aes(color=CountryGroup))+facet_wrap(~Species, scales="free",labeller = label_wrap_gen(30)) + labs(y="Index") + theme(strip.text = element_text(face="bold", size=4.8,lineheight=10.0), axis.text.x = element_text(size=6),aspect.ratio=1) + scale_x_continuous(breaks = round(seq(min(birds2b$Year), max(birds2b$Year), by = 5),1))
p
dev.off()

####################
pdf("EBCC_trend_selected regions.pdf", width=8, height=5)
Birds12 <- Birds12[Birds12$Year<2017,]
p <- ggplot(Birds12, aes(Year,values))+geom_line(aes(color=CountryGroup))+facet_wrap(~Species, scales="free_y",labeller = label_wrap_gen(10)) + labs(y="Index") + theme(strip.text = element_text(face="bold", size=4.8,lineheight=10.0), axis.text.x = element_text(size=6),aspect.ratio=1) + theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1)) + theme_bw() + scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) 

#Results_merged_allrNew1 <- Results_merged_allrNew1[Results_merged_allrNew1$Year<2017,]
#Results_merged_allrNew1 <- Results_merged_allrNew1[Results_merged_allrNew1$Year>2004,]

#p <- p+ geom_smooth(data =Results_merged_allrNew1, aes(Year,scaledValues),se=F, method="loess",show.legend = FALSE,lwd=0.7)+ scale_x_continuous(breaks = round(seq(min(Results_merged_allrNew1$Year), max(Results_merged_allrNew1$Year), by = 5),1))

p
dev.off()

###Reduced figure
#Birds14 <- Birds13[Birds13$Species!="Phoenicurus ochruros"&Birds13$Species!="Oenanthe hispanica"&Birds13$Species!="Upupa epops",]
#Results_merged_allrNew2 <- Results_merged_allrNew1[Results_merged_allrNew1$Species!="Phoenicurus ochruros"&Results_merged_allrNew1$Species!="Oenanthe hispanica"&Results_merged_allrNew1$Species!="Upupa epops",]
  
#pdf("EBCC_trend_selected regions_reduced.pdf")
#Birds14 <- Birds14[Birds14$Year<2016,]
#p <- ggplot(Birds14, aes(Year,scaledValues))+geom_line(aes(color=CountryGroup))+facet_wrap(~Species, scales="free",labeller = label_wrap_gen(30)) + labs(y="Index") + theme(strip.text = element_text(face="bold", size=4.8,lineheight=10.0), axis.text.x = element_text(size=6),aspect.ratio=1)

#Results_merged_allrNew2 <- Results_merged_allrNew2[Results_merged_allrNew2$Year<2016,]
#Results_merged_allrNew2 <- Results_merged_allrNew2[Results_merged_allrNew2$Year>2004,]

#p <- p+ geom_smooth(data =Results_merged_allrNew2, aes(Year,scaledValues),se=F, method="loess",show.legend = FALSE,lwd=0.7, span=1)+ scale_x_continuous(breaks = round(seq(min(Results_merged_allrNew1$Year), max(Results_merged_allrNew1$Year), by = 5),1))
#p 
#dev.off()
#write.table(Birds14,file="Birds14.txt",sep="\t",na="0",row.names=F,col.names=T)
#########################

#####################
#Reduce the EBCC data 
Birds15 <- Birds12[c(2,3,7:8)]

##Load the smoothed data of the Ponza index and assign the species
datasmoothed <- read.table("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Trends/Smoothed data_new.txt", header = TRUE, sep="\t")
SpeciesList <- read.csv("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Database/Reduced species list_new.csv", header = TRUE, sep=";")
datasmoothed1 <- merge(datasmoothed,SpeciesList,by="Number")
datasmoothedF <- datasmoothed1[c(2:4)]
#datasmoothedF <- datasmoothed1[c(2,3,7)]

###Center the Ponza index values
datasmoothedFi <- datasmoothedF %>%
  group_by(Species) %>%
  mutate(scaledValuesFit=scale(Fit, center = TRUE, scale=TRUE))

Birds20 <- Birds15 %>%
       left_join(datasmoothedFi)
Birds20 <- Birds20[c(1:4,6)]

##Final matrix  
weightMatrix <- read.csv("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Trends/Weight_matrix.csv", sep=";") 

##Sensitivity matrix with changes in swallow and willow warbler only  
#weightMatrix <- read.csv("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Trends/Weight_matrix_sensitivity.csv", sep=";") 

###Need to rename East and SouthEast in birds 20 to match the weight file
levels(Birds20$CountryGroup)[levels(Birds20$CountryGroup)=="East Europe"] <- "Central & East Europe"
levels(Birds20$CountryGroup)[levels(Birds20$CountryGroup)=="SouthEast Europe"] <- "West Balkan"

databaseWeight <- merge(weightMatrix,Birds20, by=c("Species","CountryGroup"))

##Scale the EBCC index using the weight from the GIS analysis
databaseWeightF <- databaseWeight%>% 
  mutate(ScaledValues=values*Weight) %>% 
  group_by(Year,Species) %>% 
  summarise(Weight=sum(Weight), 
            ScaledValuesW=sum(ScaledValues)/sum(Weight) 
  ) 

###Center the EBCC weighed index values
databaseWeightFi <-databaseWeightF %>%
  group_by(Species) %>%
  mutate(scaledValues=scale(ScaledValuesW, center = TRUE, scale=TRUE))

databaseWeightFin <- merge(databaseWeightFi,datasmoothedFi, by=c("Species","Year"))
colnames(databaseWeightFin)[5]<-'scaledValuesW'
colnames(databaseWeightFin)[7]<-'scaledValuesFit'

###Figure final
pdf("Figure final.pdf")
p <- ggplot(databaseWeightFin, aes(Year,scaledValuesW))+geom_line(lwd=1)+ facet_wrap(~Species, scales="free",labeller = label_wrap_gen(30)) + labs(y="Index") + theme(strip.text = element_text(face="bold", size=4.8,lineheight=10.0), axis.text.x = element_text(size=6),aspect.ratio=1)+ scale_x_continuous(breaks = round(seq(min(databaseWeightFin$Year), max(databaseWeightFin$Year), by = 5),1))

p <- p + geom_line(aes(x=Year, y=scaledValuesFit), col="black", lwd=0.7, linetype=5)

#scale_linetype_manual(values = c("Solid","Dashed"),name="", labels = c("EBCC","Ponza"),guide="legend")+ theme(plot.title = element_text(lineheight=.8, face="bold"), legend.position="top") 

p
dev.off()

###Correlation analysis with smoothed data
require(broom)
#databaseWeightF$Idx <- paste(databaseWeightF$Species,databaseWeightF$CountryGroup,sep="_")
sp1 = unique(data.frame(databaseWeightFin)[,"Species"])

pd3 <- list()
for(i in 1:length(sp1)) { 
     pd1 <- subset(databaseWeightFin, Species == sp1[i])
    pd2 <- lm(pd1$scaledValuesW~pd1$scaledValuesFit)
  pd3[[i]] <- tidy(pd2)
}

###Center the Ponza index values data unsmoothed
UndatasmoothedFi <- Results_mergedr2 %>%
  group_by(Species) %>%
  mutate(scaledValuesFit=scale(fit, center = TRUE, scale=TRUE))

databaseWeightFinUNS <- merge(databaseWeightFi,UndatasmoothedFi, by=c("Species","Year"))

###Correlation analysis with unsmoothed data
#databaseWeightF$Idx <- paste(databaseWeightF$Species,databaseWeightF$CountryGroup,sep="_")
sp1 = unique(data.frame(databaseWeightFinUNS)[,"Species"])

pd3 <- list()
for(i in 1:length(sp1)) { 
  pd1 <- subset(databaseWeightFinUNS, Species == sp1[i])
  pd2 <- lm(pd1$scaledValues~pd1$scaledValuesFit)
  pd3[[i]] <- tidy(pd2)
}

###Figure final with unsmoothed Ponza data
pdf("Figure final_all.pdf")
p <- ggplot(databaseWeightFinUNS, aes(Year,scaledValues))+geom_line(lwd=1)+ facet_wrap(~Species, scales="free",labeller = label_wrap_gen(30)) + labs(y="Index") + theme(strip.text = element_text(face="bold", size=4.8,lineheight=10.0), axis.text.x = element_text(size=6),aspect.ratio=1)+ scale_x_continuous(breaks = round(seq(min(databaseWeightFinUNS$Year), max(databaseWeightFinUNS$Year), by = 5),1))

p <- p + geom_line(aes(x=Year, y=scaledValuesFit), col="black", lwd=0.7, linetype=5)

p
dev.off()

###Figure final both smoothed
pdf("Figure final_both_smoothed_new.pdf")
p <- ggplot(databaseWeightFin, aes(Year,scaledValuesFit))+geom_line(lwd=1)+ facet_wrap(~Species, scales="free_y",labeller = label_wrap_gen(10)) + labs(y="Index") + theme(strip.text = element_text(face="bold", size=5,lineheight=10.0), axis.text.x = element_text(size=6),aspect.ratio=1)+ scale_x_continuous(breaks = round(seq(min(databaseWeightFin$Year), max(databaseWeightFin$Year), by = 5),1))+
  theme_bw()

#scale_linetype_manual(values = c("Solid","Dashed"),name="", labels = c("EBCC","Ponza"),guide="legend")+ theme(plot.title = element_text(lineheight=.8, face="bold"), legend.position="top") 
p <- p+ geom_smooth(aes(Year,scaledValuesW),se=F, method="loess",show.legend = FALSE,lwd=0.7,col="black",linetype=5) 

p
dev.off()

###Extract data from the reduced species ggplot for correlation analysis (done)
pg <- ggplot_build(p)
pg1 <- pg$data[[2]]
write.table(pg1,file="pg1_EBCC_sensitivity.txt",sep="\t",na="0",row.names=F,col.names=T)
#write.table(pg1,file="pg1_Ponza.txt",sep="\t",na="0",row.names=F,col.names=T)

##Load the smoothed data of both indices
datasmoothedEBCC <- read.table("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Trends/Smoothed data_All_new.txt", header = TRUE, sep="\t")
SpeciesList <- read.csv("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Database/Reduced species list_new.csv", header = TRUE, sep=";")
datasmoothed1EBCC <- merge(datasmoothedEBCC,SpeciesList,by="Number")

###Correlation analysis with both smoothed data (final analysis)
#databaseWeightF$Idx <- paste(databaseWeightF$Species,databaseWeightF$CountryGroup,sep="_")
sp1 = unique(data.frame(datasmoothed1EBCC)[,"Species"])

pd3 <- list()
for(i in 1:length(sp1)) { 
  pd1 <- subset(datasmoothed1EBCC, Species == sp1[i])
  pd2 <- lm(pd1$Ponza~pd1$EBCC)
  pd3[[i]] <- tidy(pd2)
}

datasmoothed1EBCC$Species[datasmoothed1EBCC$Species == "Delichon urbica"] <- "Delichon urbicum"

###Figure final both smoothed
pdf("Figure final_both_smoothed_new.pdf")
p <- ggplot(datasmoothed1EBCC, aes(Year,Ponza))+geom_line(lwd=1)+ facet_wrap(~Species, scales="free_y",labeller = label_wrap_gen(10)) + labs(y="Index")+  theme_bw() + theme(strip.text.x = element_text(size=8, face="italic"), axis.text.x = element_text(size=6),aspect.ratio=1)+ scale_x_continuous(breaks = round(seq(min(datasmoothed1EBCC$Year), max(datasmoothed1EBCC$Year), by = 5),1))
#scale_linetype_manual(values = c("Solid","Dashed"),name="", labels = c("EBCC","Ponza"),guide="legend")+ theme(plot.title = element_text(lineheight=.8, face="bold"), legend.position="top") 
p <- p+ geom_smooth(aes(Year,EBCC),se=F, method="loess",show.legend = FALSE,lwd=0.7,col="black",linetype=5) 

p
dev.off()

###theme_bw() overwrite the plot
