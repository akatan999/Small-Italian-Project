###################
library(gdata)
library(mgcv)
library(PBSmapping)
library(RODBC)
library(nlme)
library(ade4)
library(tidyr)
library(dplyr)

rm(list=ls()) 

#load the data
dat3 <- read.csv("Ponza 2002_2024_0 data.csv", sep=",")
dat4 <- read.csv("Ponza 2006.csv")
effortPonza <- read.csv("Daily effort_Ponza_2002-2024.csv", sep=",")
dat1 <- read.csv("ponza2021_22.csv", sep=",")
dat2 <- read.csv("Ponza 2003_2010_corretto_3duplicates.csv", sep=";")
dat5 <- read.csv("Ponza 2023.csv", sep=",")
dat6 <- read.csv("Ponza 2024.csv", sep=",")

# join all the data in a unique data frame
Dbase_merged <- rbind(dat1,dat2,dat3,dat4,dat5,dat6)

#######Make Hispanica, Subalpine and Lanius as a species for both databases
Dbase_merged[Dbase_merged == '12651'] <- '12650'
Dbase_merged[Dbase_merged == '12652'] <- '12650'
Dbase_merged[Dbase_merged == '12655'] <- '12650'

Dbase_merged[Dbase_merged == '15231'] <- '15230'
Dbase_merged[Dbase_merged == '15232'] <- '15230'
Dbase_merged[Dbase_merged == '15233'] <- '15230'

Dbase_merged[Dbase_merged == '11481'] <- '11480'
Dbase_merged[Dbase_merged == '11482'] <- '11480'

#Select the target species
Dbase_merged_all <- Dbase_merged[Dbase_merged$EURING=="12760"|Dbase_merged$EURING=="12590"|Dbase_merged$EURING=="12750"|Dbase_merged$EURING=="13120"|Dbase_merged$EURING=="11370"|Dbase_merged$EURING=="10990"|Dbase_merged$EURING=="13080"|Dbase_merged$EURING=="13350"|Dbase_merged$EURING=="12650"|Dbase_merged$EURING=="13490"|Dbase_merged$EURING=="11220"|Dbase_merged$EURING=="13110"|Dbase_merged$EURING=="9920"|Dbase_merged$EURING=="11460"|Dbase_merged$EURING=="12000"|Dbase_merged$EURING=="12770"|Dbase_merged$EURING=="11040"|Dbase_merged$EURING=="10090"|Dbase_merged$EURING=="15080"|Dbase_merged$EURING=="15230"|Dbase_merged$EURING=="13480"|Dbase_merged$EURING=="8400"|Dbase_merged$EURING=="8460"|Dbase_merged$EURING=="6870"|Dbase_merged$EURING=="8480"|Dbase_merged$EURING=="11480"|Dbase_merged$EURING=="15231"|Dbase_merged$EURING=="15232"|Dbase_merged$EURING=="15233"|Dbase_merged$EURING=="11210"|Dbase_merged$EURING=="11390"|Dbase_merged$EURING=="12530"|Dbase_merged$EURING=="16360"|Dbase_merged$EURING=="12430"|Dbase_merged$EURING=="10010"|Dbase_merged$EURING=="12651"|Dbase_merged$EURING=="12652"|Dbase_merged$EURING=="12655"|Dbase_merged$EURING=="11481"|Dbase_merged$EURING=="11482",]

#Optional for extended species (added 4 new species for age analysis)
Dbase_merged_age <- Dbase_merged[Dbase_merged$EURING=="12760"|Dbase_merged$EURING=="12590"|Dbase_merged$EURING=="12750"|Dbase_merged$EURING=="13120"|Dbase_merged$EURING=="11370"|Dbase_merged$EURING=="10990"|Dbase_merged$EURING=="13080"|Dbase_merged$EURING=="13350"|Dbase_merged$EURING=="12650"|Dbase_merged$EURING=="13490"|Dbase_merged$EURING=="11220"|Dbase_merged$EURING=="13110"|Dbase_merged$EURING=="9920"|Dbase_merged$EURING=="11460"|Dbase_merged$EURING=="12000"|Dbase_merged$EURING=="12770"|Dbase_merged$EURING=="11040"|Dbase_merged$EURING=="10090"|Dbase_merged$EURING=="15080"|Dbase_merged$EURING=="15230"|Dbase_merged$EURING=="13480"|Dbase_merged$EURING=="8400"|Dbase_merged$EURING=="8460"|Dbase_merged$EURING=="6870"|Dbase_merged$EURING=="8480"|Dbase_merged$EURING=="11480"|Dbase_merged$EURING=="15231"|Dbase_merged$EURING=="15232"|Dbase_merged$EURING=="15233"|Dbase_merged$EURING=="11210"|Dbase_merged$EURING=="11390"|Dbase_merged$EURING=="12530"|Dbase_merged$EURING=="16360"|Dbase_merged$EURING=="12430"|Dbase_merged$EURING=="10010"|Dbase_merged$EURING=="12651"|Dbase_merged$EURING=="12652"|Dbase_merged$EURING=="12655"|Dbase_merged$EURING=="11481"|Dbase_merged$EURING=="11482"|Dbase_merged$EURING=="9810"|Dbase_merged$EURING=="7780"|Dbase_merged$EURING=="10110"|Dbase_merged$EURING=="7390",]
###########

#Matching EURING CODE with the species
euring <- read.csv("EURING_Code.csv")
Dbase_merged_all$Species <- euring$Species[match(Dbase_merged_all$EURING,euring$EURING)]
sppList_paper = unique(Dbase_merged_all$Species)

#Matching RING CODE with the species (optional)
ring <- read.csv("Ring.csv")
Dbase_merged_ring = Dbase_merged_all
Dbase_merged_ring$Ring <- ring$RING[match(Dbase_merged_ring$Species,ring$Species)] 
Dbase_merged_ring = Dbase_merged_ring[Dbase_merged_ring$CATCH>0,]
colnames(Dbase_merged_ring)[2]<-'DATE'

###Optional for ring analysis
BCPDPonzaeffortring <- merge(Dbase_merged_ring, effortPonza, by="DATE")

#Select the target species where it is possible to distinguish the age class (optional)
Dbase_merged_age <- Dbase_merged_age[Dbase_merged_age$EURING=="12750"|Dbase_merged_age$EURING=="11370"|Dbase_merged_age$EURING=="10990"|Dbase_merged_age$EURING=="12650"|Dbase_merged_age$EURING=="13490"|Dbase_merged_age$EURING=="11220"|Dbase_merged_age$EURING=="13110"|Dbase_merged_age$EURING=="11460"|Dbase_merged_age$EURING=="12000"|Dbase_merged_age$EURING=="12770"|Dbase_merged_age$EURING=="11040"|Dbase_merged_age$EURING=="10090"|Dbase_merged_age$EURING=="15080"|Dbase_merged_age$EURING=="15230"|Dbase_merged_age$EURING=="13480"|Dbase_merged_age$EURING=="8400"|Dbase_merged_age$EURING=="8460"|Dbase_merged_age$EURING=="6870"|Dbase_merged_age$EURING=="8480"|Dbase_merged_age$EURING=="11480"|Dbase_merged_age$EURING=="11210"|Dbase_merged_age$EURING=="11390"|Dbase_merged_age$EURING=="16360",]

Dbase_merged_age = Dbase_merged_age[Dbase_merged_age$ETA=="5"|Dbase_merged_age$ETA=="6",]

##Save the merged data
colnames(Dbase_merged_all)[2]<-'DATE'
Dbase_merged_all <- merge(Dbase_merged_all, effortPonza, by="DATE")
write.table(Dbase_merged_all,file="Database_individuals_all.txt",sep="\t",na="0",row.names=F,col.names=T)

write.table(Dbase_merged_age,file="Database_individuals_age.txt",sep="\t",na="0",row.names=F,col.names=T)

##Calculate n of birds per species caught per day 
BCPD <- Dbase_merged_all %>%                           # Specify data frame
  group_by(EURING,DATA) %>%                        # Specify group indicator
  summarise_at(vars(CATCH),                        # Specify column
               list(CATCH = sum),na.rm = FALSE)    # Specify function

##Calculate n of birds per species caught per day only for age (optional)
BCPD_age <- Dbase_merged_age %>%                   # Specify data frame
  group_by(EURING,DATA,ETA) %>%                        # Specify group indicator
  summarise_at(vars(CATCH),                        # Specify column
               list(CATCH = sum),na.rm = FALSE)    # Specify function

##Calculate n of birds per per day 
BCPD_all <- Dbase_merged_all %>%                           # Specify data frame
  group_by(DATA) %>%                        # Specify group indicator
  summarise_at(vars(CATCH),                        # Specify column
               list(CATCH = sum),na.rm = FALSE)    # Specify function

BCPD_age = BCPD_age[complete.cases(BCPD_age),]

##Merge effort with BCPD_all
colnames(BCPD_all)[1]<-'DATE'
BCPD_all = merge(BCPD_all,effortPonza, by="DATE")
BCPD_all$cpue = BCPD_all$CATCH/BCPD_all$Effort

###Calculate the julian day, age analysis
BCPD_all$Julian <- do.call(paste, list(BCPD_all$MONTH, BCPD_all$Day, BCPD_all$Year))
BCPD_all$Julian <- as.Date(BCPD_all$Julian, format=c("%m %d %Y"))
BCPD_all$Days <- as.numeric(format(BCPD_all$Julian, "%j"))
write.table(BCPD_all ,file="BCPD_all.txt",sep="\t",na="0",row.names=F,col.names=T)

################
##Calculate n of birds per species caught per day and hour (optional)
#BCPD <- aggregate(Dbase_merged$CATCH,list(Dbase_merged$EURING,Dbase_merged$DATA,Dbase_merged$ORA),sum)

##Calculate n of birds per species caught per day and sex
#BCPD <- aggregate(Dbase_merged$CATCH,list(Dbase_merged$EURING,Dbase_merged$DATA, Dbase_merged$SESSO),sum)

#Assign names to the columns in case of hour and sex as extra variables
#colnames(BCPD)[3]<-'Hour'
#colnames(BCPD)[4]<-'Catch'
#colnames(BCPD)[3]<-'Sex'

##Save the merged data with hour
#write.table(BCPD,file="Database_selected_hour.txt",sep="\t",na="0",row.names=F,col.names=T)

###Optional for hour analysis
#BCPDPonzaeffortHour <- BCPDPonzaeffort[c(2,3,5,7,8,10,12)]

###Add daylight optional
#daylight <- read.csv("Daylight.csv", sep=";")
#daylight <- daylight[c(2,5)]

#Merge catch data with daylight
#BCPDPonzaeffortHour <- merge(BCPDPonzaeffortHour, daylight, by="Days")
#BCPDPonzaeffortHour <- BCPDPonzaeffortHour[c(2:11,1,12)]

##Scale for the daylight (optional)
#BCPDPonzaeffortHour$CPUE1 <- BCPDPonzaeffortHour$CPUE*(0.470/BCPDPonzaeffortHour$Decimals) 

#write.table(BCPDPonzaeffort,file="BCPDPonzaeffort_sex.txt",sep="\t",na="0",row.names=F,col.names=T)
#write.table(BCPDPonzaeffortHour,file="BCPDPonzaeffort_hour.txt",sep="\t",na="0",row.names=F,col.names=T)
####################

#Assign names to the columns
colnames(BCPD)[1]<-'EURING'
colnames(BCPD)[2]<-'DATE'
#colnames(BCPD)[3]<-'Station'
colnames(BCPD)[3]<-'Catch'

#Assign names to the columns age (optional)
colnames(BCPD_age)[1]<-'EURING'
colnames(BCPD_age)[2]<-'DATE'
colnames(BCPD_age)[4]<-'Catch'

#Merge catch data with effort data (all)
BCPDPonzaeffort <- merge(BCPD, effortPonza, by="DATE")

#Merge catch data with effort data (age)
BCPDPonzaeffort_age <- merge(BCPD_age, effortPonza, by="DATE")

#Calcolate CPUE (catch per unit of effort, all)
BCPDPonzaeffort$CPUE <- BCPDPonzaeffort$Catch/BCPDPonzaeffort$Effort*100

#Calcolate CPUE (catch per unit of effort, age)
BCPDPonzaeffort_age$CPUE <- BCPDPonzaeffort_age$Catch/BCPDPonzaeffort_age$Effort*100

###Calculate the julian day
require(lubridate)
BCPDPonzaeffort$Julian <- do.call(paste, list(BCPDPonzaeffort$MONTH, BCPDPonzaeffort$Day, BCPDPonzaeffort$Year))
BCPDPonzaeffort$Julian <- as.Date(BCPDPonzaeffort$Julian, format=c("%m %d %Y"))
BCPDPonzaeffort$Days <- as.numeric(format(BCPDPonzaeffort$Julian, "%j"))

###Calculate the julian day, age analysis
BCPDPonzaeffort_age$Julian <- do.call(paste, list(BCPDPonzaeffort_age$MONTH, BCPDPonzaeffort_age$Day, BCPDPonzaeffort_age$Year))
BCPDPonzaeffort_age$Julian <- as.Date(BCPDPonzaeffort_age$Julian, format=c("%m %d %Y"))
BCPDPonzaeffort_age$Days <- as.numeric(format(BCPDPonzaeffort_age$Julian, "%j"))

##Save the year catch file
#write.table(Total ,file="Total.txt",sep="\t",na="0",row.names=F,col.names=T)
write.table(BCPDPonzaeffort ,file="BCPDPonzaeffort.txt",sep="\t",na="0",row.names=F,col.names=T)
write.table(BCPDPonzaeffort_age ,file="BCPDPonzaeffort_age.txt",sep="\t",na="0",row.names=F,col.names=T)



