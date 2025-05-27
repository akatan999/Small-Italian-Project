#rm(list=ls()) 
#dev.off(

# Updated 2019-11-05

### **************************************************** #
### Ponza bird timing data     
### **************************************************** #
 
# Methods description 
# Two differemt methods were used to study timing of migration. Method 1 calculated the central tendency of the catch per unit effort for each species each year. The central tendency model, described in Colebrook 1979, (Colebrook, J.M., 1979. Continuous Plankton Records: Seasonal Cycles of Phytoplankton and Copepods in the North Atlantic Ocean and the North Sea. Mar. Biol. 51, 23–32.). The central tendency model looks at the cumulative value of catch over the season (compensated for effort) and defines peak timing as the day at the 50% of the cumulative catch is reached. Similarly, the tails of the distribution is defined as the point when 10 % and 90 % are is reached. The length of the migration period is defined as the number of days between p10 and p90. 
# This model has a drawback in that it gives only accurate peak timing estimates if the data covers the whole migration season. If the start is missed, our peak estimate will become biased towards later dates, and the opposite if the end of the migration is missed. We took this into account by first modelling the average migration period for each species by fitting a normal distribution function to the pooled data for all years for each species, and then only analysed species where we had the whole migration period, as defined by the dates where 10 % and 90 % of the probability distribution of CPUE. 

# Method 2 used another method to define timing and subset data. It used a 1-week (7 days) the moving average of the CPUE, and defines peak timing for each year and species as the day with the highest CPUE of the smoothed data. The start and the end are defined as the days when the CPUE value is down to 10 % of the peak value. If the season starts after the migration onset and/or ends before the migration termination, the tails of the distributions cannot be defined. In the cases where both the start, peak and end are defined, it is possible to robustly calculate the migration period. In cases when any of these are missing,  it is still possible to look at the inter-annual variation of the individual metrics.          
### **************************************************** #
### LIBRARIES
library(plyr)
library(caTools)
library(mgcv)
library(ggplot2)
library(dplyr)
library(gdata)

### **************************************************** #
### FUNCTIONS
# 
length.real <- function(x) { length(x[is.na(x) == FALSE]) }
length.not0 <- function(x) { length(x[x != 0 & is.na(x) == FALSE]) }
lmtrend <- function(x) { b = 1:length(x); a <- lm(scale(as.numeric(x)) ~ b); return(a$coef[2]) }

# Function for central tendency (peak timing)
# Choose in the if statement how many data points in minimum you want for the calculation; in this case we set to a minimum of 5 individuals and estimate 10,50 and 90 percentile
ct <- function(x) {  
a1 <- subset(x, CPUE > 0)  
  if(nrow(a1) < 5) { p50 = NA; p10 = NA; p90 = NA; days = NA }
else {
df <- data.frame(Days = min(a1[,"Days"]):max(a1[,"Days"]))

a2 = merge(a1, df, by = "Days", all.y = TRUE)
a2[,"CPUE"][is.na(a2[,"CPUE"])] <- 0
a3 <- cumsum(a2[,"CPUE"])
a4 <- a3/a3[length(a3)]
perc10 <- a4 > .1
perc50 <- a4 > .5
perc90 <- a4 > .9

p10 <- a2[,"Days"][match("TRUE", perc10)]
p50 <- a2[,"Days"][match("TRUE", perc50)]
p90 <- a2[,"Days"][match("TRUE", perc90)]
days <- p90-p10
}

return(data.frame(p10 = p10, p50 = p50, p90 = p90, days = days))
}

## Function for peak timing based on moving average
## This function is looking at smoothed catch per effort over the season (moving average, 7 days) and defines the peak timing as the day of the max of the smoothed data. 
## It further defines the beginning and the end of the season as the days of less than 10 % of peak CPUE (i.e. lim <.1). 

peak2 <- function(x) {
a1 <- subset(x, CPUE > -1) 
  if(length.not0(a1[,"CPUE"]) < 4) { peak = NA; pstart = NA; pend = NA; ndays = NA; prepeak = NA; postpeak = NA }
else {

#a1 <- subset(dat2, Species == "12590" & Yr == "2006")

df <- data.frame(Days = min(a1[,"Days"]):max(a1[,"Days"]))
a2 = merge(a1, df, by = "Days", all.y = TRUE)
a2[,"CPUE"][is.na(a2[,"CPUE"])] <- 0

a3 <- runmean(a2[,"CPUE"], k = 7)

peak <- match(max(a3), a3)
peakf <- a3/max(a3)
lim <- .1
pstart <- match("TRUE", peakf > lim)
if (any(peakf[1:pstart] < lim) == FALSE) { pstart <- NA }
pend <- 1 + length(a3) - match("TRUE", rev(peakf) > lim)
if (any(peakf[length(peakf):pend] < lim) == FALSE) { pend <- NA }
peak <- a2[,"Days"][peak]
pstart <- a2[,"Days"][pstart]
pend <- a2[,"Days"][pend]

if (length(pstart) > 1) { pstart <- NA }
if (length(pend) > 1) { pend <- NA }
ndays = pend-pstart
prepeak = peak-pstart
postpeak = pend-peak

}
return(data.frame(peak = peak, pstart = pstart, pend = pend, ndays = ndays, prepeak = prepeak, postpeak = postpeak))
}

# ******************************************************* #
### WORKING DIRECTORY

# Jonas
#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/data/work - synced/research/max collaboration birds/")

# Max
setwd("~/Max/PPI/Publications")
output.dir <- "~/Max/PPI/Publications"

# ******************************************************* #
### Read data Jonas
#dat <- read.delim("BCPDPonzaeffort_sex.txt")
#species <- read.delim("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Database/Species_list_Cardinale2.txt"); species <- species[,1:3]

### Read data Max
#dat <- read.delim("BCPDPonzaeffort.txt")
dat <- read.delim("~/Max/PPI/Publications/Time and trends/BCPDPonzaeffort.txt")

#dat <- dat[,c(1:8,12:14)]

species <- read.delim("~/Max/PPI/Publications/Time and trends/Species_list_Cardinale.txt"); species <- species[,1:2]

#ebcc <- read.csv("C:/Users/massi/OneDrive/Dokument/Max files for backup/Documents/Bird migration/Publications/European birds trends/Fenology and EBCC papers/Database/RegionalData_Cardinale2018.csv") 

# Aggregate subspecies
dat[dat[,"EURING"] %in% c("15231", "15232", "15233"),"EURING"] <- "15230"
dat[dat[,"EURING"] %in% c("12651", "12652", "12655"),"EURING"] <- "12650"
dat$Species <- species[match(dat[,"EURING"], species[,"EURING"]),1]

# Define Sex on some species
#dat$SexSpec = species[match(dat[,"EURING"], species[,"EURING"]),"Sex"]
#dat$Spec_sex = ifelse(dat$SexSpec == "MF", paste(dat$Species, ", Sex = ", dat$Sex, sep = ""), as.character(dat$Species))

# Aggreate into catch per effort by sex
#dat2 <- aggregate(list(CPUE = dat[,"CPUE"], Catch = dat[,"Catch"]), list(Species = dat[,"Spec_sex"], Yr = dat[,"Year"], Days = dat[,"Days"]), sum)

dat2 <- aggregate(list(CPUE = dat[,"CPUE"], Catch = dat[,"Catch"]), list(Species = dat[,"Species"], Yr = dat[,"Year"], Days = dat[,"Days"]), sum)
# ******************************************************* #

### TRENDS IN EBCC DATA
#ebcc2 <- subset(ebcc, CountryGroup %in% c("East Europe", "North Europe", "South Europe"))[,c(2,3,40:53)]
#ebcc2$trend <- apply(ebcc2[,4:ncol(ebcc2)], 1, lmtrend)
#ebcc3 <- reshape(ebcc2[,c(1,2, ncol(ebcc2))], idvar = 1, timevar = 2, v.names = 3, direction = "wide")[,c(1,4:6)]
#ebcc3$trend <- ebcc3[,]

#Save ebcc data 
#write.table(ebcc2,file="ebcc2.txt",sep="\t",na="0",row.names=F,col.names=T)
# ******************************************************* #

### PART ONE, MODEL MIGRATION PERIOD
colnames(species)[1]<-'Species'
sp1 = unique(species[,"Species"])

period = 50:160

df <- matrix(rep(NA, length(sp1)*5), nrow = length(sp1), ncol = 5)
for(i in 1:length(sp1)) { # 
  
datax <- subset(dat2, Species == sp1[i])

if(nrow(datax) > 10) {
mod1 <- gam(log(CPUE+1)~s(Days, k=3),data=datax,
            family=gaussian(link = "log"))
new.modsel <- expand.grid(Days=period)
new.lam.pred <-predict(mod1, new.modsel, type="response", se=TRUE)

c0 <- exp(as.vector(unlist(new.lam.pred$fit)))-1
c1 <- cumsum(c0)
c2 = data.frame(Days = period, Cumsum = c1/c1[length(c1)])
p10 <- c2[,2] > .15; p10 <- match("TRUE", p10)
p50 <- c2[,2] > .50; p50 <- match("TRUE", p50)
p90 <- c2[,2] > .85; p90 <- match("TRUE", p90)

limits <- c(period[1]+p10, period[1]+p90, period[1]+p50)

df[i,1] = sp1[i]
df[i,2] = limits[1]
df[i,3] = limits[2]
df[i,4] = limits[3]
#df[i,5] = df[i,3]-df[i,2]
}
}

colnames(df) = c("Species", "p10", "p90", "p50", "90%length")
df <- as.data.frame(df)
df$Species <- sp1

# ******************************************************* #

# First and last day per year
datx1 <- aggregate(dat[,"Days"], list(Year = dat[,"Year"]), function(x) min(x, na.rm = TRUE))
datx2 <- aggregate(dat[,"Days"], list(Year = dat[,"Year"]), function(x) max(x, na.rm = TRUE))
datx3 <- cbind(datx1, datx2[,"x"]); colnames(datx3) <- c("Year", "first", "last")

# ******************************************************* #

# Pick out species where we cover the full migration period
speciesperiod <- matrix(NA, nrow = length(sp1)*20, ncol = 4)
speciesperiod[,2] <- rep(2002:2021, length(sp1))
colnames(speciesperiod) <- c("Species", "Yr", "okstart", "okend")
speciesperiod <- as.data.frame(speciesperiod)
speciesperiod[,"Species"] <- rep(sp1, each = 20)

years <- 2002:2024
for(i in 1:length(years)) {

#i = 10
d1 <- subset(datx3, Year == years[i])  
for(j in 1:length(sp1)) {
  
#j = 4  
d2 <- df[df[,"Species"] == sp1[j],][1,]  
if(d1[,"first"] < d2[2]) { speciesperiod[speciesperiod[,"Species"] == sp1[j] & speciesperiod[,"Yr"] == years[i],"okstart"]<- "ok"} else { speciesperiod[speciesperiod[,"Species"] == sp1[j] & speciesperiod[,"Yr"] == years[i],"okstart"] <- "no" }
if(d1[,"last"] > d2[3]) { speciesperiod[speciesperiod[,"Species"] == sp1[j] & speciesperiod[,"Yr"] == years[i],"okend"] <- "ok"} else { speciesperiod[speciesperiod[,"Species"] == sp1[j] & speciesperiod[,"Yr"] == years[i],"okend"] <- "no" }
  
}}

# ******************************************************* #

# Subset data for analysis based on table 
spperiod.ok <- subset(speciesperiod, okstart == "ok" & okend == "ok") 
newdat <- dat2[paste(dat2[,"Species"], dat2[,"Yr"]) %in%
paste(spperiod.ok[,"Species"], spperiod.ok[,"Yr"]),] 

# ******************************************************* #
### PART TWO, CALCULATE TIMING DATA

# Central tendency, quartiles, etc. 
## See ct function in the beginning of the document 

# Calculate peak timing
timing <- ddply(newdat, .(Yr, Species), ct) # Based on central tendency
timing <- timing[is.na(timing[,3]) == FALSE,]

#Save timing data for the first index
write.table(timing,file="timing.txt",sep="\t",na="0",row.names=F,col.names=T)

# PLOTS
pdf(file = "birds_period.pdf")
tab <- table(timing[,1], timing[,2])
sp <- colnames(tab)[colSums(tab) > 1]

# Comparing first capture and peak timing
#Yr = 2012
for(i in 1:length(sp)) {
pd1 <- subset(dat2, Species == sp[i] & Yr == Yr)
df <- data.frame(Days = min(pd1[,"Days"]):max(pd1[,"Days"]))
pd2 <- merge(pd1, df, by = "Days", all.y = TRUE)
x = barplot(pd2[,"CPUE"], names.arg = pd2[,"Days"], xaxt = "n", las = 2, xlab = "Julian day", ylab = "CPUE (n/m of net)", main = paste("", sp[i]))
axis(side = 1, at = c(1, 500, 1000), labels = c("70", "100", "150"))

pd3 <- subset(timing, Species == sp[i] & Yr == Yr)
m1 <- match(pd2[,"Days"], pd3[,4:3])
m2 <- x[m1 == 1]; m2 <- m2[is.na(m2) == FALSE]  
m3 <- x[m1 == 2]; m3 <- m3[is.na(m3) == FALSE]  
abline(v = m2, col ="red", lwd = 2)
abline(v = m3, col ="blue", lwd = 2)
}

dev.off()

pdf(file = "birds_timing.pdf")

# Comparing two measure of timing
for(i in 1:length(sp)) {

pd1 <- subset(timing, Species == sp[i])
pd2 <- merge(data.frame(2002:2021), pd1, by = 1, all.x = TRUE)
plot(pd2[,1], pd2[,"p50"], type = "o", col = "blue", pch = 19, xlab = "Year", ylab = "Julian day", ylim = c(min(pd1[,"p10"])*.9, max(pd1[,"p90"])*1.2), main = paste("", sp[i]))
abline(lsfit(pd2[,1], pd2[,"p50"]), lty = 2, col = "blue")  
    
# p10
points(pd2[,1], pd2[,"p10"], type = "o", col = "red", pch = 19)
abline(lsfit(pd2[,1], pd2[,"p10"]), lty = 2, col = "red")  

# p90
points(pd2[,1], pd2[,"p90"], type = "o", col = "red", pch = 19)
abline(lsfit(pd2[,1], pd2[,"p90"]), lty = 2, col = "red")  

legend("top", pch = 19, lty = 1, col = c("red", "blue"), bty = "n", c("Percentiles", "Peak timing"))

}

dev.off()

# ******************************************************* #
## New function for defining peak based on moving average
## See peak2 function above including descriptions. 
## This method does not rely on a pre-defined migration period, it just uses the data and the empirical distribution of CPUE, which I think is better! 

require(grDevices)

# Run 
timing <- ddply(dat2, .(Yr, Species), peak2)
names(timing) [3] <- "p50"
names(timing) [4] <- "p10"
names(timing) [5] <- "p90"
names(timing) [6] <- "days"
timing <- timing[,c(1:6)]

#levels(timing$Species)[11] <- "Jynx torquilla"
#levels(timing$Species)
#levels(timing$Species)[13] <- "Luscinia megarhynchos"
#levels(timing$Species)[30] <- "Turdus philomelos"

#Save timing data for the second index
write.table(timing,file="timing.txt",sep="\t",na="0",row.names=F,col.names=T)
#write.csv(timing,file="timing.csv")
#write.table(timing,file="timing_sex2.txt",sep="\t",na="0",row.names=F,col.names=T)

###############If you want to restart from the timing file already saved#######
timing <- read.table("timing.txt", header=TRUE)
sppList <- (unique(timing$Species))

names(timing) [3] <- "peak"
names(timing) [4] <- "pstart"
names(timing) [5] <- "pend"
names(timing) [6] <- "days"

# PLOT New timing measure
pdf(file = "birds_timing1_p1.pdf", width=12, height = 6)

par(mfrow=c(2,4), pty="s")

#for(i in 1:length(sp1)) {  #:length(sp)
for(i in 1:8) {  #:length(sp)
  pd1 <- subset(timing, Species == sp[i])
slope.start = NA; slope.peak = NA; slope.end = NA 
 if(nrow(pd1) > 5) { 
  pd2 <- merge(data.frame(2003:2021), pd1, by = 1, all.x = TRUE)
#quartz(); 

# Peak 
plot(pd2[,1], pd2[,"peak"], type = "o", col = "black", pch = 19, xlab = "Year", ylab = "Julian day", ylim = c(.7*min(pd1[,"peak"], na.rm = TRUE), 1.3*max(pd1[,"peak"], na.rm = TRUE)), main = sp[i])

if(length.real(pd2[,"peak"]) > 4) { abline(lsfit(pd2[,1], pd2[,"peak"]), lty = 2, col = "black") 
   slope.peak <- lm(pd2[,"peak"] ~ pd2[,1])$coef[2] }
   
# Start of migration
points(pd2[,1], pd2[,"pstart"], type = "o", col = "red", pch = 19)
if (length.real(pd2[,"pstart"]) > 5)
 { abline(lsfit(pd2[,1], pd2[,"pstart"]), lty = 2, col = "red")  
 slope.start <- lm(pd2[,"pstart"] ~ pd2[,1])$coef[2]
}

# End of migration 
points(pd2[,1], pd2[,"pend"], type = "o", col = "red", pch = 19)
if (length.real(pd2[,"pend"]) > 5)
 { abline(lsfit(pd2[,1], pd2[,"pend"]), lty = 2, col = "red")
  slope.end <- lm(pd2[,"pend"] ~ pd2[,1])$coef[2]	  }

}
if(i == 1) { slopes <- data.frame() }
slopes <- rbind(slopes, data.frame(Species = sp[i], Slope.start = slope.start, Slope.peak = slope.peak, Slope.end = slope.end))

legend("top", pch = 19, lty = 1, col = c("red", "black"), bty = "n", c("Percentiles", "Peak timing"))
}

dev.off()

############
###Figures for the paper
#####
#timing <- read.csv("Timing_final.csv", sep=";")
timing <- read.table("Timing.txt", header=T)
timing[timing==0] <- NA
timing <- timing[timing$Species != "Delichon urbicum",]
sppList <- (unique(timing$Species))

pdf(file = "birds_timing_final1.pdf")
ggplot(filter(timing,Species %in% sppList[1:15])) +
  #ggplot(data=timing, aes(x=Yr,y=p50,ymin=p10,ymax=p90), col=1) +  
  geom_linerange(data=filter(timing,Species %in% sppList[1:15]), aes(x=Yr,ymin=p50,ymax=p90), col=1, na.rm = TRUE)+geom_linerange(data=filter(timing,Species %in% sppList[1:15]), aes(x=Yr,ymin=p10,ymax=p50), col=1, na.rm = TRUE)+geom_point(data=filter(timing,Species %in% sppList[1:15]), aes(x=Yr,y=p50), col=1, na.rm = TRUE)+
  geom_smooth(aes(Yr,p10),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5, na.rm = TRUE) +
  geom_smooth(aes(Yr,p50),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5, colour="black", na.rm = TRUE) + 
  geom_smooth(aes(Yr,p90),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5, na.rm = TRUE) +
  facet_wrap(~Species, scale="free_y", labeller=label_wrap_gen(10)) +
  theme_bw() + theme(strip.text = element_text(face = "italic")) +
  xlab("Year") + ylab("Time span (Days)") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))
dev.off()

pdf(file = "birds_timing_fina2.pdf")
ggplot(filter(timing,Species %in% sppList[16:31])) +
  #ggplot(data=timing, aes(x=Yr,y=p50,ymin=p10,ymax=p90), col=1) +  
  geom_linerange(data=filter(timing,Species %in% sppList[16:31]), aes(x=Yr,ymin=p50,ymax=p90), col=1, na.rm = TRUE)+geom_linerange(data=filter(timing,Species %in% sppList[16:31]), aes(x=Yr,ymin=p10,ymax=p50), col=1, na.rm = TRUE)+geom_point(data=filter(timing,Species %in% sppList[16:31]), aes(x=Yr,y=p50), col=1, na.rm = TRUE)+
  geom_smooth(aes(Yr,p10),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5, na.rm = TRUE) +
  geom_smooth(aes(Yr,p50),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5, colour="black", na.rm = TRUE) + 
  geom_smooth(aes(Yr,p90),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5, na.rm = TRUE) +
  facet_wrap(~Species, scale="free_y", labeller=label_wrap_gen(10)) +
  theme_bw() + theme(strip.text = element_text(face = "italic")) +
  xlab("Year") + ylab("Time span (Days)")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))
dev.off()

pdf(file = "birds_timing_final_all species.pdf")
ggplot(filter(timing,Species %in% sppList[1:30])) +
  #ggplot(data=timing, aes(x=Yr,y=p50,ymin=p10,ymax=p90), col=1) +  
  geom_linerange(data=filter(timing,Species %in% sppList[1:30]), aes(x=Yr,ymin=p50,ymax=p90), col=1, na.rm = TRUE, linetype=4)+geom_linerange(data=filter(timing,Species %in% sppList[1:30]), aes(x=Yr,ymin=p10,ymax=p50), col=1, na.rm = TRUE, linetype=4)+geom_point(data=filter(timing,Species %in% sppList[1:30]), aes(x=Yr,y=p50), col=1, na.rm = TRUE)+
  geom_smooth(aes(Yr,p10),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=4, na.rm = TRUE) +
  geom_smooth(aes(Yr,p50),se=F, method="loess",show.legend = FALSE,lwd=0.7, linetype=5, colour="red", na.rm = TRUE) + 
  geom_smooth(aes(Yr,p90),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=4, na.rm = TRUE) +
  facet_wrap(~Species, scale="free_y", labeller=label_wrap_gen(10)) +
  theme_bw() + theme(strip.text = element_text(face = "italic")) +
  xlab("Year") + ylab("Julian day")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9))
dev.off()

####################
pdf(file = "birds_timing_final3.pdf")
ggplot(filter(timing,Species %in% sppList[1:15])) +
  geom_line(aes(Yr,days), col=1, linetype=5) +
  geom_smooth(aes(Yr,days),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5) +
  geom_point(aes(Yr,days)) + facet_wrap(~Species, scale="free_y", labeller=label_wrap_gen(10)) +
  theme_bw() +
  xlab("Year") + ylab("Time span(days)")
dev.off()

########################
pdf(file = "birds_timing_final4.pdf")
ggplot(filter(timing,Species %in% sppList[16:31])) +
  geom_line(aes(Yr,days), col=1, linetype=5) +
  geom_smooth(aes(Yr,days),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5) +
  geom_point(aes(Yr,days)) + facet_wrap(~Species, scale="free_y", labeller=label_wrap_gen(10)) +
  theme_bw() +
  xlab("Year") + ylab("Time span (days)")
dev.off()
################################################

############

###Figures for the paper (by sex)
timingsex <- read.delim("timing_sex2.txt")
pdf(file = "birds_timing_final_sex2.pdf")
ggplot(timingsex2M, aes(x=Yr,y=p50,ymin=p10,ymax=p90), col=1) +  geom_pointrange(data=timingsex, aes(x=Yr,y=p50,ymin=p10,ymax=p90), col=1)+ geom_smooth(aes(Yr,p10),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5) +geom_smooth(aes(Yr,p50),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5, colour="black") + geom_smooth(aes(Yr,p90),se=F, method="lm",show.legend = FALSE,lwd=0.7, linetype=5) + facet_wrap(~Species+Sex, scale="free_y", labeller=label_wrap_gen(30)) + theme_bw() + xlab("Year") + ylab("Time span (Days)")
dev.off()

# Time span (90 %)
for(i in 1:length(sp)) {
  pd1 <- subset(timing, Species == sp[i])
  pd2 <- merge(data.frame(2002:2019), pd1, by = 1, all.x = TRUE)
  plot(pd2[,1], pd2[,"days"], type = "o", col = "blue", pch = 19, xlab = "Year", ylab = "Julian day", main = paste("Species =", sp[i]))
  abline(lsfit(pd2[,1], pd2[,"days"]), lty = 2, col = "blue")  
  axis(4)  
}

dev.off()
###################################
###Tables
dat2 <- aggregate(Catch~Year+Species,data = dat, FUN="sum")

library(tidyr)
dat3 <- spread(dat2,Year,Catch)
dateffort<-aggregate(Effort~Year,data = dat, FUN="mean")

##Save the data
write.table(dat3,file="TableStatistcs.txt",sep="\t",na="0",row.names=T,col.names=T)
##Save the data
write.table(dateffort,file="Tableeffort.txt",sep="\t",na="0",row.names=T,col.names=T)

#####
timing <- read.csv("Timing_final.csv",sep=";")
timing$p50[timing$p50==0] <- NA
timing$p10[timing$p10==0] <- NA
timing$p90[timing$p90==0] <- NA
timing$days[timing$days==0] <- NA
sp1 <- (unique(timing$Species))

require(broom)
pd3 <- list()
pd4 <- list()
pd5 <- list()
pd6 <- list()

i=14

#for(i in 1:length(sp1)) { 

pd1 <- subset(timing, Species == sp1[i])
pd2 <- lm(pd1$p10~pd1$Yr)
pd3[[i]] <- tidy(pd2)
pd2 <- lm(pd1$p50~pd1$Yr)
pd4[[i]] <- tidy(pd2)
pd2 <- lm(pd1$p90~pd1$Yr)
pd5[[i]] <- tidy(pd2)
pd2 <- lm(pd1$days~pd1$Yr)
pd6[[i]] <- tidy(pd2)

}

DataFrame <- do.call(rbind.data.frame, pd3)
DataFrame1 <- do.call(rbind.data.frame, pd4)
DataFrame2 <- do.call(rbind.data.frame, pd5)
DataFrame3 <- do.call(rbind.data.frame, pd6)

TableStat <- cbind(DataFrame,DataFrame1,DataFrame2,DataFrame3)

##Save the data
write.table(TableStat,file="TableStat_quantile.txt",sep="\t",na="0",row.names=T,col.names=T)


#####
##By sex, males
timingsex <- read.csv("Timing_sex.csv",sep=";")

timingsexM <- timingsex[timingsex$Sex=="2",]
sp1 <- (unique(timingsexM$Species))
timingsexM$p50[timingsexM$p50==0] <- NA
timingsexM$p10[timingsexM$p10==0] <- NA
timingsexM$p90[timingsexM$p90==0] <- NA

pd3 <- list()
pd4 <- list()
pd5 <- list()
pd6 <- list()

#for(i in 1:length(sp1)) { 

  i= 7
  
  pd1 <- subset(timingsexM, Species == sp1[i])
  pd2 <- lm(pd1$p10~pd1$Year)
  pd3[[i]] <- tidy(pd2)
  pd2 <- lm(pd1$p50~pd1$Year)
  pd4[[i]] <- tidy(pd2)
  pd2 <- lm(pd1$p90~pd1$Year)
  pd5[[i]] <- tidy(pd2)
  pd2 <- lm(pd1$days~pd1$Year)
  pd6[[i]] <- tidy(pd2)

  }

summary(pd2)


DataFrame <- do.call(rbind.data.frame, pd3)
DataFrame1 <- do.call(rbind.data.frame, pd4)
DataFrame2 <- do.call(rbind.data.frame, pd5)
DataFrame3 <- do.call(rbind.data.frame, pd6)

TableStat <- cbind(DataFrame,DataFrame1,DataFrame2,DataFrame3)

##Save the data
write.table(TableStat,file="TableStatsexM.txt",sep="\t",na="0",row.names=T,col.names=T)

##By sex, females
timingsexF <- timingsex[timingsex$Sex=="2",]
sp1 <- (unique(timingsexF$Species))
timing$p50[timing$p50==0] <- NA
timing$p10[timing$p10==0] <- NA
timing$p90[timing$p90==0] <- NA

pd3 <- list()
pd4 <- list()
pd5 <- list()
pd6 <- list()

for(i in 1:length(sp1)) { 
  pd1 <- subset(timingsexF, Species == sp1[i])
  pd2 <- lm(pd1$Yr~pd1$p10)
  pd3[[i]] <- tidy(pd2)
  pd2 <- lm(pd1$Yr~pd1$p50)
  pd4[[i]] <- tidy(pd2)
  pd2 <- lm(pd1$Yr~pd1$p90)
  pd5[[i]] <- tidy(pd2)
  pd2 <- lm(pd1$Yr~pd1$days)
  pd6[[i]] <- tidy(pd2)
}

DataFrame <- do.call(rbind.data.frame, pd3)
DataFrame1 <- do.call(rbind.data.frame, pd4)
DataFrame2 <- do.call(rbind.data.frame, pd5)
DataFrame3 <- do.call(rbind.data.frame, pd6)

TableStat <- cbind(DataFrame,DataFrame1,DataFrame2,DataFrame3)

##Save the data
write.table(TableStat,file="TableStatsexF.txt",sep="\t",na="0",row.names=T,col.names=T)

###Quantile regression
library(quantreg)
require(glmnet)
library(MCMCpack)
library(arm)

timing <- read.table("timing.txt", header =TRUE)
sp1 <- (unique(timing$Species))
timing$p50[timing$p50==0] <- NA
timing$p10[timing$p10==0] <- NA
timing$p90[timing$p90==0] <- NA
timing$days[timing$days==0] <- NA

i=3

pd1 <- subset(timing, Species == sp1[i])

model1 <- rq(pd1$p10~pd1$Yr, tau=.5)
model2 <- rq(pd1$p50~pd1$Yr, tau=.5)
model3 <- rq(pd1$p90~pd1$Yr, tau=.9)
model4 <- rq(pd1$days~pd1$Yr, tau=.5)

summary(model1, se="ker")
summary(model2, se="ker")
summary(model3, se="ker")
summary(model4, se="ker")

plot(pd1$p10~pd1$Yr, data = pd1, pch = 16, main = "mpg ~ wt")
lm(pd1$p10~pd1$Yr, data = pd1, col = "red", lty = 2)
abline(rq(pd1$p10~pd1$Yr, data = pd1), col = "blue", lty = 2)

lm <- lm(pd1$p50~pd1$Yr)
summary(lm)

###Other regressions
pd2 <- bayesglm(pd1$Yr~pd1$p10)
pd3 <- summary(pd2)
pd4 <- bayesglm(pd1$Yr~pd1$p50)
pd5 <- summary(pd4)
pd6 <- bayesglm(pd1$Yr~pd1$p90)
pd7 <- summary(pd6)

model <- MCMCregress(pd1$Yr~pd1$p10, data=pd1)







