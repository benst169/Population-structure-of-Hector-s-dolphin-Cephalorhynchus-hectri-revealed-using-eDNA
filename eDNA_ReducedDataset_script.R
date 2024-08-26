#Workspace setup----
rm(list=ls(all=TRUE))# Clear workspace
library(dplyr)
library(raster)
library(sf)
library(terra)
library(mapplots)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Data----
ceph.seq92 <- read.csv("CephSequences_n92.csv") #List of sequences names (incl. length and sequence) that were matched as Cephalorhynchus
site.seq92 <- read.csv("SiteSequence_92.csv") #Matrix of all sequences (summer 2022, n=46), sites and number of reads
hap.match<-read.csv("Haplotype_matching.csv") 
edna.meta <- read.csv("eDNA.metadata.csv",na.strings = c("NA",""))

##Data management ----
### Retain only Cephalorhynchus ASVs
ceph.seq92$V1<-gsub(">","",ceph.seq92$V1)
site.seq92<-merge(ceph.seq92, site.seq92, by.x = "V1", by.y = "sequence_names", all.x = TRUE)
site.seq92[,4:98]<-lapply(site.seq92[,4:98],as.integer)

### Remove alternative sampling method 
site.seq92 <- site.seq92[,-c(7,25,44,53)]

### Manipulation into three forms of the same data: dat1 = Number of reads in every sequence (rows) against samples (columns); dat2 = Number of reads in every sample (rows) against sequence (columns); dat3 = as dat2 but proportions
dat1<-site.seq92
dat1[is.na(dat1)] <- 0
dat1<-merge(dat1,hap.match[,c(1,8:10)],by = "V1",all.x=TRUE)
dat1<-dat1[,c(1:3,95:97,4:94)]

#Change dataframe orientation
dat2<-t(dat1[,c(6:97)]); colnames(dat2)<-dat2[1,]; dat2<-as.data.frame(dat2[-1,])
dat2[] <- lapply(dat2, as.numeric)
dat2<-dat2[, order(names(dat2))]
dat2$Total.Reads<-apply(dat2[1:114],1,FUN=sum)

#Convert read count data to proportion of reads in a sample
dat3<-dat2[,1:114]/dat2$Total.Reads
dat3[dat3==0]<-NaN


##Error rate----
#Isolate all previously identified likely error ASVs, then their occurrence
pot.errors<-dat3[,11:59]

vec.list<-NA
for(i in 1:ncol(pot.errors)){
  for (j in 1:nrow(pot.errors)) {
    if(is.na(pot.errors[j,i])==FALSE){
      c <- pot.errors[j,i]
      vec.list<-append(vec.list,c)
    }
  }
}
vec.list<-vec.list[-1]

## Calculate the error rate DT
n <- length(vec.list); alpha <- 0.01; mean.vec<-mean(vec.list); sd.vec <- sd(vec.list); t.value<-qt((1-(alpha/2)),n-1)
DT.mutated <- mean.vec +t.value*sd.vec

##Remove probable errors based on the DT
dat3 <- as.data.frame(t(dat3[,1:114]))

for(i in 1:ncol(dat3)){
  for(j in 1:nrow(dat3)){
    if(is.na(dat3[j,i])==FALSE & dat3[j,i] <= DT.mutated){
      dat3[j,i] = 0
    }
  }
}

dat3[is.na(dat3)] <- 0
dat3<-subset(dat3,apply(dat3,1,sum) > 0) 

dat3$NoSamples <- NA
for(i in 1:nrow(dat3)){
  dat3[i,92] <- sum(dat3[i,1:91]>0)
}

#Cleanup workspace
rm(i,j,c,mean.vec,n,t.value,vec.list,alpha,sd.vec,pot.errors,ceph.seq92,site.seq92)

##Final Haplotype Criteria----
HapCriteria <- as.data.frame(as.vector(row.names(dat3)))
colnames(HapCriteria)<-c("HapNames")
# Criteria 1: Does the ASV occur across multiple samples
HapCriteria$MultipleSamples <- ifelse(dat3$NoSamples > 1,1,0)

# Criteria 2: Does the ASV occur in the highest proportion at least once in a sample
HighestProp<-rep(NA,91)
for(i in 1:91){
  max.row <- which.max(dat3[,i])
  HighestProp[i]<- row.names(dat3)[max.row]
}

HighestProp<-as.data.frame(unique(HighestProp))
colnames(HighestProp)<-"HapNames"
HighestProp$HighestInSample<-1
HapCriteria <- merge(HapCriteria,HighestProp, by = "HapNames",all.x = TRUE)
HapCriteria[is.na(HapCriteria)] <- 0
HapCriteria$NoCriteriaFit<-apply(HapCriteria[,2:3],1,sum)


## Final Dataframes----
#If ASVs match one criteria than they are considered, if match two than more confident they are real
FinalConfHaplotypes <- as.data.frame(HapCriteria$HapNames[which(HapCriteria$NoCriteriaFit == 2)])
FinalAllHaplotypes <- as.data.frame(HapCriteria$HapNames[which(HapCriteria$NoCriteriaFit >= 1)])
colnames(FinalConfHaplotypes)<-"names"
colnames(FinalAllHaplotypes)<-"names"
dat3$names<-rownames(dat3)


FinalConfHaplotypes <- merge(FinalConfHaplotypes,dat3[,c(1:91,93)], by= "names", all.x = TRUE, all.y = FALSE)#All sequences we are confident are true haplotypes
FinalAllHaplotypes <- merge(FinalAllHaplotypes,dat3[,c(1:91,93)], by= "names", all.x = TRUE, all.y = FALSE) #All sequences including those we are less confident about
FinalAllHaplotypes[9:10,1] <- c("eDNA1","eDNA2")

rownames(FinalConfHaplotypes)<-FinalConfHaplotypes$names
FinalConfHaplotypes<-as.data.frame(t(FinalConfHaplotypes[,2:92]))
FinalConfHaplotypes[FinalConfHaplotypes==0]<-NA  

rownames(FinalAllHaplotypes)<-FinalAllHaplotypes$names
FinalAllHaplotypes<-as.data.frame(t(FinalAllHaplotypes[,2:92]))
FinalAllHaplotypes[FinalAllHaplotypes==0]<-NA  

#Cleanup workspace
rm(hap.match,HapCriteria,HighestProp,max.row,i)




#Data summaries----
##Easy functions ----
se<-function(x){
  sd(x)/sqrt(length(x))
}

count2<- function(x){
  length(which(x == 2))
}

count3<- function(x){
  length(which(x == 3))
}

count1<- function(x){
  length(which(x == 1))
}

count0<- function(x){
  length(which(x == 0))
}

count1plus<- function(x){
  length(which(x > 0))
}

#Results -----
##Summary statistics----
HaplotypeSummaryReads <- dat2

#Remove erroneous sequences
for(i in 1:114){
  for(j in 1:nrow(HaplotypeSummaryReads)){
    error <- HaplotypeSummaryReads[j,115]*DT.mutated
    if(HaplotypeSummaryReads[j,i]-error < 0){
      HaplotypeSummaryReads[j,i] <- 0
    }else{
      HaplotypeSummaryReads[j,i]<- HaplotypeSummaryReads[j,i]
    }
  }
}

HaplotypeSummaryReads <- HaplotypeSummaryReads[,c(1,3:8,10,60,89)]
colnames(HaplotypeSummaryReads) <- colnames(FinalAllHaplotypes)

HaplotypeSummaryReads$TotalConfReads<-apply(HaplotypeSummaryReads[,c(1:6,8)],1,sum)
HaplotypeSummaryReads$TotalConfHaps<-apply(HaplotypeSummaryReads[,c(1:6,8)],1,count1plus)

HaplotypeSummaryReads$SampleID<-rownames(HaplotypeSummaryReads)
HaplotypeSummaryReads$Site <- edna.meta[match(HaplotypeSummaryReads$SampleID,edna.meta$EppendorfID),1]
HaplotypeSummaryReads$Season <- edna.meta[match(HaplotypeSummaryReads$SampleID,edna.meta$EppendorfID),2]
HaplotypeSummaryReads$SampleType<-edna.meta[match(HaplotypeSummaryReads$SampleID,edna.meta$EppendorfID),10]
HaplotypeSummaryReads$SampleType<-ifelse(is.na(HaplotypeSummaryReads$SampleType) == TRUE,"Lab control",HaplotypeSummaryReads$SampleType)


###Table 1: Summary of sampling methods and haplotypes using only haplotypes we are confident in ----
SampleTypeSum<-aggregate(TotalConfReads~SampleType,data = HaplotypeSummaryReads[,-c(7,9,10)],FUN = mean)
SampleTypeSE<-aggregate(TotalConfReads~SampleType,data = HaplotypeSummaryReads[,-c(7,9,10)],FUN = se)
SampleTypeN<-aggregate(TotalConfHaps~SampleType,data = HaplotypeSummaryReads[,-c(7,9,10)],FUN = length)
SampleTypeN0<-aggregate(TotalConfHaps~SampleType,data = HaplotypeSummaryReads[,-c(7,9,10)],FUN = count0)
SampleTypeN1<-aggregate(TotalConfHaps~SampleType,data = HaplotypeSummaryReads[,-c(7,9,10)],FUN = count1)
SampleTypeN2<-aggregate(TotalConfHaps~SampleType,data = HaplotypeSummaryReads[,-c(7,9,10)],FUN = count2)
SampleTypeN3<-aggregate(TotalConfHaps~SampleType,data = HaplotypeSummaryReads[,-c(7,9,10)],FUN = count3)

SampleTypeSum$SE<-SampleTypeSE$TotalConfReads
SampleTypeSum$N<-SampleTypeN$TotalConfHaps
SampleTypeSum$N0<-SampleTypeN0$TotalConfHaps
SampleTypeSum$N1<-SampleTypeN1$TotalConfHaps
SampleTypeSum$N2<-SampleTypeN2$TotalConfHaps
SampleTypeSum$N3<-SampleTypeN3$TotalConfHaps

#write.csv(SampleTypeSum,"C:/Users/steph/Desktop/Publishing/eDNA/Tables/SampleTypeSummary_reduced.csv")

#Extract the count and associated frequency for each Haplotype for each site
colnames(HaplotypeSummaryReads)[2]<-"Hap_Cplus"
HapSites<-aggregate(cbind(Hap_A,Hap_Cplus,Hap_D,Hap_H,Hap_I,Hap_K,Hap_S)~Site,data=HaplotypeSummaryReads,FUN=count1plus)
HapSitesFreq<-HapSites
for(i in 1:nrow(HapSitesFreq)){
  HapSitesFreq[i,2:8]<-HapSitesFreq[i,2:8]/sum(HapSites[i,2:8])
}


#write.csv(HapSites,"HaplotypeFrequenciesSites_reduced.csv")

rm(SampleTypeSE,SampleTypeN,SampleTypeN1,SampleTypeN0,SampleTypeN2,SampleTypeN3,error,i,j,DT.mutated,HapSites,HapSitesFreq,SampleTypeSum)



#Figures----
#Figure 1. Map of study area----
nzcoast<-shapefile("./Shapefiles/nzcoast.shp")
south.island <-crop(nzcoast, extent(166,174.5,-47,-40))

study.extent<-shapefile("./Shapefiles/StudyExtent.shp")
study.extent<-study.extent[study.extent$id != 3,]
study.extent<-spTransform(study.extent,crs(nzcoast))

#inset.limits<- c(0.1309,0.35,0.71,0.9602)
inset.limits<- c(0.05,0.4,0.6,0.96)

tiff("C:/Users/steph/Desktop/Publishing/eDNA/Figures/Figure 1. Map of study area.tiff", width = 200, height = 200, res=300, units = "mm")
par(mar=c(0.1,0.1,0.1,0.1))
plot(south.island,col="gray",lwd = 0.05)
#plot(bathy,add=TRUE)
plot(study.extent,add=TRUE,col = "blue")
#lines(x=c(173.5,173.5), y=c(-46,-43.4),lwd=1.2);lines(x=c(170.2,173.5),y=c(-43.4,-43.4),lwd=1.2);lines(x=c(170.2,173.5),y=c(-46,-46),lwd=1.2);lines(x=c(170.2,170.2),y=c(-43.4,-46),lwd=1.2) #Lines for original study area map
lines(x=c(174.5,174.5), y=c(-47,-40.5),lwd=1.2);lines(x=c(166.035,174.5),y=c(-40.5,-40.5),lwd=1.2);lines(x=c(166.035,174.5),y=c(-47,-47),lwd=1.2);lines(x=c(166.035,166.035),y=c(-40.5,-47),lwd=1.2)

text(x=172,y=-43.73,labels = "Banks Peninsula",cex = 0.9);text(x=170.9,y=-44.4,labels = "Timaru",cex = 0.9);text(x=170.2,y=-45.78,labels = "Dunedin",cex=0.9);text(x=173.25,y=-42.35,labels = "Kaikōura",cex = 0.8);text(x=173.5,y=-41.65,labels = "Cloudy/Clifford \nBay",cex = 0.8);text(x=170.55,y=-45.1,labels = "Oamaru",cex = 0.8);text(x=170.45,y=-45.35,labels = "Moeraki",cex = 0.8);text(x=169.35,y=-46.3,labels = "Kaka Point",cex = 0.8);text(x=168.85,y=-46.46,labels = "Toetoe Bay",cex = 0.8);text(x=167.6,y=-46.05,labels = "Te Waewae Bay",cex = 0.8)
scalebar(92.6,xy=c(172.5,-46.95),type="bar",divs = 4, lonlat = TRUE, label = NA,cex = 1); text(x=174,y=-46.92,"50nm",cex=1)

text(173.8,-43.73,"n=28")
text(171.9,-44.4,"n=9")
text(171.35,-45.78,"n=11")


par(fig=inset.limits, new=TRUE, mar=c(0,0,0,0) )
plot(nzcoast,col="gray",lwd=0.001)
lines(x=c(174.5,174.5), y=c(-47,-40.5),lwd=0.8);lines(x=c(166.04,174.5),y=c(-40.5,-40.5),lwd=0.8);lines(x=c(166.04,174.5),y=c(-47,-47),lwd=0.8);lines(x=c(166.04,166.04),y=c(-40.5,-47),lwd=0.8)
north(xy = c(165,-38),type=2,cex=1.5)
axis(1,at = seq(165,180,by=5),labels = c("165.00º","170.00º","175.00º","180.00º"),cex.axis = 0.8)
axis(4, at = seq(-50,-35, by = 5),labels = c("-50.00º","-45.00º","-40.00º","-35.00º"),las=2,cex.axis = 0.8)


box(lwd=0.9,col="black",lty = "solid")
dev.off()

rm(south.island,study.extent,inset.limits)

##Figure 2. eDNA framework info----
#Using C+ and a couple of samples with unknowns as examples:
x.locs = 1:5
unknowns = dat2[c(9,14,15,22,39),60:114]/dat2[c(9,14,15,22,39),115]
unknowns[unknowns==0]<-NA

tiff("C:/Users/steph/Desktop/Publishing/eDNA/Figures/Reduced Data/Figure 2. Example of Haplotype, likely errors and unknown sequences.tiff", width= 120, height = 80, res=300, units="mm")
par(mgp = c(2.3,0.4,0), mai = c(1.2,0.9,0.1,0.1))
plot(x.locs, FinalConfHaplotypes[c(9,14,15,22,39),2],col = eDNA.conf.colours[2],ylim = c(0,1.1),axes = F, pch = 19, xlab = "", ylab = "Proportion of ASV \nin sample",yaxs = "i",cex = 0.8,, cex.lab = 0.7)
axis(2,tck = -0.01, cex.axis = 0.6)
axis(1, at = seq(1,5, by=1), labels = rownames(FinalConfHaplotypes[c(9,14,15,22,39),]),tck = -0.01, las = 2,cex.axis=0.5)
abline(h=0)
abline(h=DT.mutated, lty = "dotted")
for(i in 1:32){
  points(x.locs,pot.errors[c(9,14,15,22,39),i],pch = 21, bg =  "orange", col = "black",cex= 0.7)
  points(x.locs,unknowns[,i],pch = 21, bg= "gray", col = "black",cex= 0.7)
}
par(xpd = NA)
points(x = 0.5, y=-0.6, pch = 19, col = eDNA.conf.colours[2])
text(x = 1.4, y=-0.605, "Previously identified haplotype",cex= 0.6)
points(x = 2.4, y=-0.6, pch = 21, col = "black",bg="orange")
text(x = 3.2, y=-0.605, "Likely erroneous sequence",cex= 0.6)
points(x = 4.1, y=-0.6, pch = 21, bg = "gray",col= "black")
text(x = 4.7, y=-0.605, "Unknown sequence",cex= 0.6)
dev.off()

##Figure 3. Map of haplotype locations three panels----
##Data----
BP.area<-crop(nzcoast, extent(c(172.5,173.8,-44.37,-43.52)))
TIM.area<-crop(nzcoast, extent(c(171, 172.2, -44.6, -44.2)))
OTA.area<-crop(nzcoast, extent(c(170.38, 174, -45.9, -45.6))) 

FinalConfHaplotypes$SampleID<- rownames(FinalConfHaplotypes)
FinalConfHaplotypes[is.na(FinalConfHaplotypes)] = 0
NoHaplotypes <- rep(NA, length = nrow(FinalConfHaplotypes))
for(i in 1:length(NoHaplotypes)){
  NoHaplotypes[i] <- sum(FinalConfHaplotypes[i,1:7]>0)
}

mapdata<-merge(FinalConfHaplotypes, edna.meta[,c(5,6,19)],by.x = "SampleID", by.y = "EppendorfID",all.x = TRUE, all.y = FALSE)
mapdata<-mapdata[complete.cases(mapdata),]
mapdata$Site<-NA
mapdata[1:40,"Site"]<-"BP"
mapdata[41:54,"Site"]<-"DUN"
mapdata[55:87,"Site"]<-"TIM"
mapdata$Site<-as.factor(mapdata$Site)

no.occur <- function(x){
  length(which(x > 0))
}

site.presence<-aggregate(.~Site,data=mapdata[,c(11,2:8)],FUN = no.occur)
site.presence<-as.data.frame(t(site.presence[,-1]))
colnames(site.presence)<-c("BP","OTA","TIM")

hapA <- subset(mapdata[,c(2,9,10)], mapdata$Hap_A > 0);hapCplus <- subset(mapdata[,c(3,9,10)], mapdata$`Hap_C+` > 0);hapD <- subset(mapdata[,c(4,9,10)], mapdata$Hap_D > 0);hapH <- subset(mapdata[,c(5,9,10)], mapdata$Hap_H > 0);hapI <- subset(mapdata[,c(6,9,10)], mapdata$Hap_I > 0);hapK <- subset(mapdata[,c(7,9,10)], mapdata$Hap_K > 0);hapS <- subset(mapdata[,c(8,9,10)], mapdata$Hap_S > 0)

hapA<-SpatialPointsDataFrame(coords = hapA[,c(3,2)],data = hapA,proj4string = CRS("+init=EPSG:4326"));hapCplus<-SpatialPointsDataFrame(coords = hapCplus[,c(3,2)],data = hapCplus,proj4string = CRS("+init=EPSG:4326"));hapD<-SpatialPointsDataFrame(coords = hapD[,c(3,2)],data = hapD,proj4string = CRS("+init=EPSG:4326"));hapH<-SpatialPointsDataFrame(coords = hapH[,c(3,2)],data = hapH,proj4string = CRS("+init=EPSG:4326"));hapI<-SpatialPointsDataFrame(coords = hapI[,c(3,2)],data = hapI,proj4string = CRS("+init=EPSG:4326"));hapK<-SpatialPointsDataFrame(coords = hapK[,c(3,2)],data = hapK,proj4string = CRS("+init=EPSG:4326"));hapS<-SpatialPointsDataFrame(coords = hapS[,c(3,2)],data = hapS,proj4string = CRS("+init=EPSG:4326"))

a<-rgb(255/255,208/255,176/255); c<-rgb(201/255,48/255,35/255); d<-rgb(240/255,208/245,250/255); e<-rgb(248/255,236/255,54/255); g<-rgb(255/255,173/255,59/255); h<-rgb(16/255,161/255,216/255); i<-rgb(176/255,84/255,210/255); j<-rgb(45/255,98/255,179/255); k<-rgb(15/255,15/255,15/255); l<-rgb(94/255,193/255,0/255); m<-rgb(26/255,111/255,39/255); n<-"darkorange"; o<-rgb(134/255,153/255,168/255); p<-rgb(202/255,142/255,112/255); q<-rgb(139/255,89/255,88/255); r<-rgb(222/255,238/255,254/255); s<-rgb(13/255,113/255,119/255); v<-rgb(254/255,253/255,174/255); e1<-rgb(98/255,232/255,216/255); e2<-rgb(222/255,80/255,230/255)

eDNA.conf.colours<-c(a,c,d,h,i,k,s)
eDNA.all.colours<-c(a,c,d,h,i,k,m,s,e1,e2)

rm(a,c,d,e,g,h,i,j,k,l,m,n,o,p,q,r,s,v,e1,e2)

tiff("C:/Users/steph/Desktop/Publishing/eDNA/Figures/Reduced Data/Fig. 3. Sample locations and haplotype detections_reduced.tiff",width = 120, height = 200,res=300, unit="mm")
par(mfrow=c(3,1),mai=c(0.01,0.01,0.01,0.2))
plot(BP.area, col = "gray")
polygon(x=c(172.5,172.587,172.587,172.5),y = c(-43.5,-43.5,-43.9,-43.9),border ="white", col = "white")
lines(x = c(172.587,172.587),y=c(-43.519,-43.835))


points(hapA,bg=eDNA.conf.colours[1],pch=21,cex=2); points(hapCplus,bg=eDNA.conf.colours[2],pch=21,cex=2); points(hapD,bg=eDNA.conf.colours[3],pch=21,cex=2); points(hapH,bg=eDNA.conf.colours[4],pch=21,cex=2); points(hapI,bg=eDNA.conf.colours[5],pch=21,cex=2)
scalebar(1.852,xy=c(172.615,-43.81),type="line",lonlat = TRUE,label = NA); text(x=172.62,y=-43.821,"1nm"); north(xy = c(172.615,-43.58),type=2,cex=1.5)

add.pie(x=173.192, y= -43.60, radius=0.076, site.presence$BP,labels="",clockwise=TRUE,col=eDNA.conf.colours,cex=0.5,xpd=NA)
text(x=173.192, y= -43.60, labels = " BP \n(n=32)",pos=1,offset=4.3, cex = 1.2)

plot(TIM.area,col = "gray")
points(hapCplus,bg=eDNA.conf.colours[2],pch=21,cex=2); points(hapS,bg=eDNA.conf.colours[7],pch=21,cex=2.5); points(hapH,bg=eDNA.conf.colours[4],pch=21,cex=1.5); points(hapA,bg=eDNA.conf.colours[1],pch=21,cex=2); points(hapK,bg=eDNA.conf.colours[6],pch=21,cex=2)
scalebar(1.852,xy=c(171.029,-44.575),type="line",lonlat = TRUE, label = NA); text(x=171.039,y=-44.585,"1nm")

add.pie(x=171.638, y= -44.4, radius=0.079, site.presence$TIM,labels="",clockwise=TRUE,col=eDNA.conf.colours,cex=0.5)

text(x=171.638, y= -44.4, labels = "  TIM \n(n=32)",pos=1,offset=4.3, cex = 1.2)


plot(OTA.area,col="gray")
points(hapS,bg=eDNA.conf.colours[7],pch=21,cex=2,xpd = NA);points(hapCplus,bg=eDNA.conf.colours[2],pch=21,cex=1.25)
scalebar(1.852,xy=c(170.400,-45.884),type="line",lonlat = TRUE, label = NA); text(x=170.412,y=-45.892,"1nm")
legend(y=-45.603,x=170.385,legend=c("A","C+","D","H","I","K","S"),fill =eDNA.conf.colours,bty="o",bg="white",cex=1.3)

add.pie(x=170.87, y= -45.72, radius=0.0595, site.presence$OTA,labels="",clockwise=TRUE,col=eDNA.conf.colours,cex=0.5)

text(x=170.87, y= -45.72, labels = " DUN \n(n=14)",pos=1,offset=4.3, cex=1.2)

dev.off()



##Figure 4: Comparison of previous sampling methods and eDNA detections while being less confident on haplotypes and more confident -----
###Data setup for less confident haplotypes

FinalAllHaplotypes$SampleID<- rownames(FinalAllHaplotypes)
FinalAllHaplotypes[is.na(FinalAllHaplotypes)] = 0
NoHaplotypes <- rep(NA, length = nrow(FinalAllHaplotypes))
for(i in 1:length(NoHaplotypes)){
  NoHaplotypes[i] <- sum(FinalAllHaplotypes[i,1:10]>0)
}

mapdata1<-merge(FinalAllHaplotypes, edna.meta[,c(5,6,19)],by.x = "SampleID", by.y = "EppendorfID",all.x = TRUE, all.y = FALSE)
mapdata1<-mapdata1[complete.cases(mapdata1),]
mapdata1$Site<-NA
mapdata1[1:40,14]<-"BP"
mapdata1[41:54,14]<-"OTA"
mapdata1[55:87,14]<-"TIM"
mapdata1$Site<-as.factor(mapdata1$Site)

no.occur<-function(x){
  length(which(x > 0))
}
site.presence1<-aggregate(.~Site,data=mapdata1[,c(14,2:11)],FUN = no.occur)
site.presence1<-as.data.frame(t(site.presence1[,-1]))
colnames(site.presence1)<-c("BP","OTA","TIM")

hapM <- subset(mapdata1[,c(8,9,10)], mapdata1$Hap_M > 0);eDNA1 <- subset(mapdata1[,c(10,9,10)], mapdata1$eDNA1 > 0);eDNA2 <- subset(mapdata1[,c(11,9,10)], mapdata1$eDNA2 > 0)

hapM<-SpatialPointsDataFrame(coords = hapM[,c(3,2)],data = hapM,proj4string = CRS("+init=EPSG:4326"));eDNA1<-SpatialPointsDataFrame(coords = eDNA1[,c(3,2)],data = eDNA1,proj4string = CRS("+init=EPSG:4326"));eDNA2<-SpatialPointsDataFrame(coords = eDNA2[,c(3,2)],data = eDNA2,proj4string = CRS("+init=EPSG:4326"))


tiff("C:/Users/steph/Desktop/Publishing/eDNA/Figures/Reduced Data/Figure 4. Not confident haplotype distribution_reduced.tiff",width = 120, height = 200,res=300, unit="mm")
par(mfrow=c(3,1),mai=c(0.01,0.01,0.01,0.2))
plot(BP.area, col = "gray")
polygon(x=c(172.5,172.587,172.587,172.5),y = c(-43.5,-43.5,-43.9,-43.9),border ="white", col = "white")
lines(x = c(172.587,172.587),y=c(-43.519,-43.835))


points(hapA,bg=eDNA.all.colours[1],pch=21,cex=2); points(hapCplus,bg=eDNA.all.colours[2],pch=21,cex=2); points(hapD,bg=eDNA.all.colours[3],pch=21,cex=2); points(hapH,bg=eDNA.all.colours[4],pch=21,cex=2); points(hapI,bg=eDNA.all.colours[5],pch=21,cex=2); points(hapM,bg=eDNA.all.colours[7],pch=21,cex=2); points(eDNA2,bg=eDNA.all.colours[10],pch=21,cex=2)
scalebar(1.852,xy=c(172.615,-43.81),type="line",lonlat = TRUE,label = NA); text(x=172.62,y=-43.821,"1nm"); north(xy = c(172.615,-43.58),type=2,cex=1.5)

add.pie(x=173.192, y= -43.60, radius=0.076, site.presence1$BP,labels="",clockwise=TRUE,col=eDNA.all.colours,cex=0.5,xpd=NA)
text(x=173.192, y= -43.60, labels = " BP \n(n=34)",pos=1,offset=4.3, cex = 1.2)

plot(TIM.area,col = "gray")
points(hapCplus,bg=eDNA.all.colours[2],pch=21,cex=2); points(hapA,bg=eDNA.all.colours[1],pch=21,cex=2); points(hapS,bg=eDNA.all.colours[8],pch=21,cex=2.5); points(hapH,bg=eDNA.all.colours[4],pch=21,cex=1.5); points(hapK,bg=eDNA.all.colours[6],pch=21,cex=2)
scalebar(1.852,xy=c(171.029,-44.575),type="line",lonlat = TRUE, label = NA); text(x=171.039,y=-44.585,"1nm")

add.pie(x=171.638, y= -44.4, radius=0.079, site.presence1$TIM,labels="",clockwise=TRUE,col=eDNA.all.colours,cex=0.5)

text(x=171.638, y= -44.4, labels = "  TIM \n(n=32)",pos=1,offset=4.3, cex = 1.2)

plot(OTA.area,col="gray")
points(hapS,bg=eDNA.all.colours[8],pch=21,cex=2,xpd = NA);points(hapCplus,bg=eDNA.all.colours[2],pch=21,cex=1.25); points(eDNA1,bg=eDNA.all.colours[9],pch=21,cex=2)
scalebar(1.852,xy=c(170.400,-45.884),type="line",lonlat = TRUE, label = NA); text(x=170.412,y=-45.892,"1nm")
legend(y=-45.604,x=170.385,legend=c("A","C+","D","H","I","K","M","S","eDNA1","eDNA2"),fill =eDNA.all.colours,bty="o",bg="white",cex=1)

add.pie(x=170.87, y= -45.72, radius=0.0595, site.presence1$OTA,labels="",clockwise=TRUE,col=eDNA.all.colours,cex=0.5)

text(x=170.87, y= -45.72, labels = " DUN \n(n=15)",pos=1,offset=4.3, cex=1.2)
dev.off()

rm(hapA,hapCplus,hapD,hapH,hapI,hapK,hapM,hapS,eDNA1,eDNA2,OTA.area,TIM.area,BP.area,NoHaplotypes,no.occur,site.presence,site.presence1,mapdata1,eDNA.all.colours,eDNA.conf.colours,nzcoast)



#Appendices----

###Appendix 2: Summary of Area based sampling results----

SiteYearSUM<-aggregate(TotalConfReads~Site+Season, data= HaplotypeSummaryReads[,-c(7,9,10)],FUN = mean)
SiteYearSE<-aggregate(TotalConfReads~Site+Season, data= HaplotypeSummaryReads[,-c(7,9,10)],FUN = se)
SiteYearN<-aggregate(TotalConfHaps~Site+Season, data= HaplotypeSummaryReads[,-c(7,9,10)],FUN = length)
SiteYearN0<-aggregate(TotalConfHaps~Site+Season, data= HaplotypeSummaryReads[,-c(7,9,10)],FUN = count0)
SiteYearN1<-aggregate(TotalConfHaps~Site+Season, data= HaplotypeSummaryReads[,-c(7,9,10)],FUN = count1)
SiteYearN2<-aggregate(TotalConfHaps~Site+Season, data= HaplotypeSummaryReads[,-c(7,9,10)],FUN = count2)
SiteYearN3<-aggregate(TotalConfHaps~Site+Season, data= HaplotypeSummaryReads[,-c(7,9,10)],FUN = count3)

names(SiteYearSUM)<- c("Site","Season","Mean.reads")
SiteYearSUM$SE.reads<-SiteYearSE$TotalConfReads
SiteYearSUM$N.reads<-SiteYearN$TotalConfHaps
SiteYearSUM$N0.reads<-SiteYearN0$TotalConfHaps
SiteYearSUM$N1.reads<-SiteYearN1$TotalConfHaps
SiteYearSUM$N2.reads<-SiteYearN2$TotalConfHaps
SiteYearSUM$N3.reads<-SiteYearN3$TotalConfHaps


#write.csv(SampleTypeSum,"C:/Users/steph/Desktop/Publishing/eDNA/Tables/SiteYearSummary.csv")



##Appendices X. Discovery curve----
discovery<-mapdata[,1:11]
for(i in 2:8){
  for(j in 1:nrow(discovery)){
    if(discovery[j,i] > 0){
      discovery[j,i] <- names(discovery[i])
    }
  }
}

discovery$Date<-as.Date(edna.meta[match(discovery$SampleID,edna.meta$EppendorfID),3], format = "%d/%m/%Y")
discovery<-discovery[order(discovery$Date),]

discovery.data<-as.data.frame(rep(0:87));names(discovery.data) = "No.Samples"
discovery.data$all.haps<-c(1,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7)
                           
discovery.data$bp.haps<- c(1,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
discovery.data$tim.haps<-c(0,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
discovery.data$dun.haps<-c(1,1,1,1,1,1,2,2,2,2,2,2,2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)


log.reg<-nls(all.haps~SSlogis(No.Samples,a,b,c),data=discovery.data)
log.reg1<-nls(bp.haps~SSlogis(No.Samples,a,b,c),data=discovery.data)
log.reg2<-nls(tim.haps~SSlogis(No.Samples,a,b,c),data=discovery.data)
log.reg3<-nls(dun.haps~SSlogis(No.Samples,a,b,c),data=discovery.data)

tiff("C:/Users/steph/Desktop/Publishing/eDNA/Figures/Reduced Data/Appendices X. Discovery curve with logistic curves.tiff", width= 120, height=80, res=300, units="mm")
par(mai = c(0.7,0.7,0.2,0.2), mgp = c(1.3,0.2,0))
plot(all.haps~No.Samples, data = discovery.data,axes=FALSE, xlab = "Number of samples", ylab = "Total number of \nunique haplotypes", bty="n", ylim=c(0,10),xaxs = "i",yaxs="i",pch=20,cex.lab=0.6,cex=0.4)
axis(1, at=seq(0,140, by=10),cex.axis = 0.5,tck =-0.01)
axis(2, at=seq(0,10, by=1),cex.axis = 0.5,tck=-0.01)
lines(discovery.data[,1], predict(log.reg, discovery.data[,1]), lwd = 1.5)
points(discovery.data$bp.haps~discovery.data$No.Samples,col = rgb(201/255,48/255,35/255),pch=20,cex=0.4)
lines(discovery.data[1:41,1],predict(log.reg1,discovery.data[1:41,1]),col=rgb(201/255,48/255,35/255),lwd=1.5)
points(discovery.data$tim.haps~discovery.data$No.Samples,col = "#E69F00",pch=20,cex=0.4)
lines(discovery.data[1:33,1],predict(log.reg2,discovery.data[1:33,1]),col="#E69F00",lwd=1.5)
points(discovery.data$dun.haps~discovery.data$No.Samples,col = rgb(13/255,113/255,119/255),pch=20,cex=0.4)
lines(discovery.data[1:13,1],predict(log.reg3,discovery.data[1:13,1]),col=rgb(13/255,113/255,119/255),lwd=1.5)
legend(x=5,y=9.9,legend=c("All","Banks Peninsula","Timaru","Dunedin"), fill = c("black",rgb(201/255,48/255,35/255),"#E69F00",rgb(13/255,113/255,119/255)), border = c("black",rgb(201/255,48/255,35/255),"#E69F00",rgb(13/255,113/255,119/255)),cex=0.5,box.lwd = 0.7)
dev.off()

tiff("C:/Users/steph/Desktop/Publishing/eDNA/Figures/Reduced Data/Appendices X. Discovery curve.tiff", width= 120, height=80, res=300, units="mm")
par(mai = c(0.7,0.7,0.2,0.2), mgp = c(1.3,0.2,0))
plot(all.haps~No.Samples, data = discovery.data,axes=FALSE, xlab = "Number of samples", ylab = "Total number of \nunique haplotypes", bty="n", ylim=c(0,10),xaxs = "i",yaxs="i",type="l",cex.lab=0.6,cex=0.4)
axis(1, at=seq(0,140, by=10),cex.axis = 0.5,tck =-0.01)
axis(2, at=seq(0,10, by=1),cex.axis = 0.5,tck=-0.01)
points(y=discovery.data$bp.haps, x=discovery.data$No.Samples+0.2,col = rgb(201/255,48/255,35/255),type="l",cex=0.4)
points(y=discovery.data$tim.haps,x=discovery.data$No.Samples+0.4,col = "#E69F00",type="l",cex=0.4)
points(y=discovery.data$dun.haps,x=discovery.data$No.Samples,col = rgb(13/255,113/255,119/255),type="l",cex=0.4)
legend(x=5,y=9.9,legend=c("All","Banks Peninsula","Timaru","Dunedin"), fill = c("black",rgb(201/255,48/255,35/255),"#E69F00",rgb(13/255,113/255,119/255)), border = c("black",rgb(201/255,48/255,35/255),"#E69F00",rgb(13/255,113/255,119/255)),cex=0.5,box.lwd = 0.7)
dev.off()

##Appendices X. Off target amplification----
seq.scores<-read.csv("AllSeq_100plus.csv")
seq.scores<-subset(seq.scores, seq.scores$Cephalorhynchus == 0)
seq.scores$seqname<-gsub(">","",seq.scores$seqname)
offtarget.sites<-merge(seq.scores[,c(5,9)],all.asv, by.x = "seqname", by.y = "sequence_names",all.x=TRUE, all.y=FALSE)
offtarget.sites[is.na(offtarget.sites)]<-0

offtarget.sites<-aggregate(. ~ Species.1,offtarget.sites[,2:97],FUN =sum)

offtarget.sites$Total.Reads<-apply(offtarget.sites[2:96],1,FUN=sum)
offtarget.sites$NoSamples<-NA
for(i in 1:nrow(offtarget.sites)){
  offtarget.sites[i,98]<- sum(offtarget.sites[i,2:96]>0)
}

offtarget.sites$BP<-NA
offtarget.sites$TIM<-NA
offtarget.sites$DUN<-NA

for(i in 1:17){
  offtarget.sites[i,99]<-ifelse(sum(offtarget.sites[i,c(3:45)])>0,1,0) 
  offtarget.sites[i,100]<-ifelse(sum(offtarget.sites[i,c(64:96)])>0,1,0)   
  offtarget.sites[i,101]<-ifelse(sum(offtarget.sites[i,c(46:60)])>0,1,0) 
}

offtarget.sites$Average.Reads<- apply(offtarget.sites[2:96],1,FUN=mean)
offtarget.sites$sd.Reads<- apply(offtarget.sites[2:96],1,FUN=sd)

offtarget.sites<-offtarget.sites[-1,c(1,97:103)]

seq.scores<-seq.scores[,c(9,12,14)]
seq.scores$Query.Cover.1 <-as.numeric(gsub("%","",seq.scores$Query.Cover.1))
seq.scores$Per..Ident.1 <-as.numeric(gsub("%","",seq.scores$Per..Ident.1))
seq.scores<-aggregate(. ~ Species.1, seq.scores, mean)
offtarget.sites<-merge(offtarget.sites,seq.scores, by = "Species.1")

#write.csv(offtarget.sites,"C:/Users/steph/Desktop/Publishing/eDNA/Figures/Reduced Data/OffTargetAmplification.csv",row.names=FALSE)


##Appendices X----

fieldnotes<-read.csv("Fieldnotes.csv")

edna.meta$Sample.ranking = as.factor(edna.meta$Sample.ranking)
edna.meta$SampleType = as.factor(edna.meta$SampleType)
edna.meta$Site = as.factor(edna.meta$Site)
edna.meta$Prox2dolph = as.factor(edna.meta$Prox2dolph)

edna.meta$Depth = fieldnotes[match(edna.meta$EncounterID,fieldnotes$EncounterID),9]
edna.meta$Turbidity = fieldnotes[match(edna.meta$EncounterID,fieldnotes$EncounterID),10]
edna.meta$SST = fieldnotes[match(edna.meta$EncounterID,fieldnotes$EncounterID),11]
edna.meta$BF = fieldnotes[match(edna.meta$EncounterID,fieldnotes$EncounterID),7]
edna.meta$Swell = fieldnotes[match(edna.meta$EncounterID,fieldnotes$EncounterID),8]
edna.meta$CollectDateTime = as.POSIXct(as.character(paste(edna.meta$Date.collected, edna.meta$Time.collected)), format="%d/%m/%Y %I:%M:%S %p")

mapdata$Date<-edna.meta[match(mapdata$SampleID,edna.meta$EppendorfID),2]

edna.meta$FilterDateTime = NA
for(i in 1:nrow(edna.meta)){
  if(i < 65){
    edna.meta[i,"FilterDateTime"] = edna.meta[i,"Time.Date.filtered"]
  }else{
    edna.meta[i,"FilterDateTime"] = paste(edna.meta[i,"Date.collected"],edna.meta[i,"Time.Date.filtered"],sep = " ")
  }
  
}
edna.meta$FilterDateTime = as.POSIXct(edna.meta$FilterDateTime,format = "%d/%m/%Y %H:%M")
edna.meta$Time2Filter = as.numeric(edna.meta$FilterDateTime-edna.meta$CollectDateTime)


HaplotypeSummaryReads$Ceph.detect<-ifelse(HaplotypeSummaryReads$TotalConfReads>0,1,0)

edna.meta$TotalCephReads <- HaplotypeSummaryReads[match(edna.meta$EppendorfID,rownames(HaplotypeSummaryReads)),"TotalConfReads"]
edna.meta$TotalCephHaps <- HaplotypeSummaryReads[match(edna.meta$EppendorfID,rownames(HaplotypeSummaryReads)),"TotalConfHaps"]
edna.meta$CephDetected <-  HaplotypeSummaryReads[match(edna.meta$EppendorfID,rownames(HaplotypeSummaryReads)),"Ceph.detect"]

edna.meta <- subset(edna.meta,edna.meta$Season == 2023)

rm(fieldnotes)


##MODELS: Cephalorhynchus detected----
m1.1 <- glm(CephDetected~Site,data=edna.meta,family = binomial(link = logit))
m1.2 <- glm(CephDetected~Turbidity,data=edna.meta,family = binomial(link = logit))
m1.3 <- glm(CephDetected~Depth,data=edna.meta,family = binomial(link = logit)) 
m1.4 <- glm(CephDetected~SST,data=edna.meta,family = binomial(link = logit))
m2.1<-glm(CephDetected~Sample.ranking,data=edna.meta,family=binomial(link = logit))
m2.2<-glm(CephDetected~tot.group.size,data=edna.meta,family=binomial(link = logit))
m2.3<-glm(CephDetected~boat.group.size,data=edna.meta,family=binomial(link = logit))
m2.4<-glm(CephDetected~Prox2dolph,data=edna.meta,family=binomial(link = logit))
m3.1<-glm(CephDetected~Amount.collected,data=edna.meta,family = binomial(link = logit))
m3.2<-glm(CephDetected~Amount.filtered,data=edna.meta,family = binomial(link = logit))
m3.3<-glm(CephDetected~Time2Filter,data=edna.meta,family = binomial(link = logit))

summary(m1.1)
summary(m1.2)
summary(m1.3)
summary(m1.4)
summary(m2.1)
summary(m2.2)
summary(m2.3)
summary(m2.4)
summary(m3.1)
summary(m3.2)
summary(m3.3)



##Plotting----
###Plot data----
#Standardised site to total number of samples collected
site.mat<-matrix(c(length(which(edna.meta$CephDetected == 0 & edna.meta$Site =="Dunedin"))/length(which(edna.meta$Site =="Dunedin")),length(which(edna.meta$CephDetected == 1 & edna.meta$Site =="Dunedin"))/length(which(edna.meta$Site =="Dunedin")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Site =="Timaru"))/length(which(edna.meta$Site =="Timaru")),length(which(edna.meta$CephDetected == 1 & edna.meta$Site =="Timaru"))/length(which(edna.meta$Site =="Timaru")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Site =="Banks Peninsula"))/length(which(edna.meta$Site =="Banks Peninsula")),length(which(edna.meta$CephDetected == 1 & edna.meta$Site =="Banks Peninsula"))/length(which(edna.meta$Site =="Banks Peninsula"))),nrow=2,ncol=3)

rank.mat<-matrix(c(length(which(edna.meta$CephDetected == 0 & edna.meta$Sample.ranking =="Poor"))/length(which(edna.meta$Sample.ranking =="Poor")),length(which(edna.meta$CephDetected == 1 & edna.meta$Sample.ranking =="Poor"))/length(which(edna.meta$Sample.ranking =="Poor")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Sample.ranking =="Ok"))/length(which(edna.meta$Sample.ranking =="Ok")),length(which(edna.meta$CephDetected == 1 & edna.meta$Sample.ranking =="Ok"))/length(which(edna.meta$Sample.ranking =="Ok")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Sample.ranking =="Good"))/length(which(edna.meta$Sample.ranking =="Good")),length(which(edna.meta$CephDetected == 1 & edna.meta$Sample.ranking =="Good"))/length(which(edna.meta$Sample.ranking =="Good")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Sample.ranking =="Excellent"))/length(which(edna.meta$Sample.ranking =="Excellent")),length(which(edna.meta$CephDetected == 1 & edna.meta$Sample.ranking =="Excellent"))/length(which(edna.meta$Sample.ranking =="Excellent"))),nrow=2,ncol=4)

prox.mat<-matrix(c(length(which(edna.meta$CephDetected == 0 & edna.meta$Prox2dolph =="<1m"))/length(which(edna.meta$Prox2dolph =="<1m")),length(which(edna.meta$CephDetected == 1 & edna.meta$Prox2dolph =="<1m"))/length(which(edna.meta$Prox2dolph =="<1m")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Prox2dolph =="1m"))/length(which(edna.meta$Prox2dolph =="1m")),length(which(edna.meta$CephDetected == 1 & edna.meta$Prox2dolph =="1m"))/length(which(edna.meta$Prox2dolph =="1m")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Prox2dolph =="2m"))/length(which(edna.meta$Prox2dolph =="2m")),length(which(edna.meta$CephDetected == 1 & edna.meta$Prox2dolph =="2m"))/length(which(edna.meta$Prox2dolph =="2m")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Prox2dolph =="3m"))/length(which(edna.meta$Prox2dolph =="3m")),length(which(edna.meta$CephDetected == 1 & edna.meta$Prox2dolph =="3m"))/length(which(edna.meta$Prox2dolph =="3m")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Prox2dolph =="4m"))/length(which(edna.meta$Prox2dolph =="4m")),length(which(edna.meta$CephDetected == 1 & edna.meta$Prox2dolph =="4m"))/length(which(edna.meta$Prox2dolph =="4m")),
                   length(which(edna.meta$CephDetected == 0 & edna.meta$Prox2dolph ==">5m"))/length(which(edna.meta$Prox2dolph ==">5m")),length(which(edna.meta$CephDetected == 1 & edna.meta$Prox2dolph ==">5m"))/length(which(edna.meta$Prox2dolph ==">5m"))), nrow=2,ncol=6)


tiff("C:/Users/steph/Desktop/Publishing/eDNA/Figures/Appendix X. Field sampling figures.tiff",height = 200, width=160, res = 300, unit="mm")
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,0),4,3,byrow = FALSE))
par(mai = c(0.5,0.5,0.2,0.1),mgp = c(2,0.5,0))
#Environment
bp<-barplot(site.mat,beside=TRUE,ylab = "Proportion of Samples \n(#/total collected)",xlab = "Site",ylim=c(0,1), cex.axis = 0.9, cex.lab = 0.9)
axis(1,at=c(2,5,8), labels = c("DUN","TIM","BP"),cex.axis = 0.8)
abline(h=0)
boxplot(Turbidity~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Secchi disk depth (m)", xlab = "Cephalorhynchus detected")
boxplot(Depth~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Water depth (m)", xlab = "Cephalorhynchus detected")
boxplot(SST~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Sea surface temperature (ºC)", xlab = "Cephalorhynchus detected")
#Collection
rbp<-barplot(rank.mat,beside=TRUE,ylab = "Proportion of Samples \n(#/total collected)",xlab = "Sample ranking",ylim=c(0,1), cex.axis = 0.9, cex.lab = 0.9)
axis(1,at=c(2,5,8,11), labels = c("Poor","Ok","Good","Excellent"),cex.axis = 0.8)
abline(h=0)
rbp<-barplot(prox.mat,beside=TRUE,ylab = "Proportion of Samples \n(#/total collected)",xlab = "Proximity to dolphins",ylim=c(0,1), cex.axis = 0.9, cex.lab = 0.9)
axis(1,at=c(2,5,8,11,14,17), labels = c("<1m","1m","2m","3m","4m",">5m"),cex.axis = 0.8)
abline(h=0)
boxplot(tot.group.size~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Total group size", xlab = "Cephalorhynchus detected")
boxplot(boat.group.size~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Group size \nat vessel", xlab = "Cephalorhynchus detected")
#Filtration
boxplot(Time2Filter~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Time from collection \nto filtration (hours)", xlab = "Cephalorhynchus detected")
boxplot(Amount.filtered~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Amount filtered (L)", xlab = "Cephalorhynchus detected")
boxplot(Amount.collected~CephDetected,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Amount Collected (L)", xlab = "Cephalorhynchus detected")

text(x = -4.5, y = 18.9, "A)", xpd=NA,cex = 1.2)
text(x = -1.5, y = 18.9, "B)",xpd=NA,cex = 1.2)
text(x = 1.5, y = 18.9, "C)",xpd=NA,cex = 1.2)
dev.off()




##MODELS: Total Reads ----
m1.1 <- glm(TotalCephReads~Site,data=edna.meta,family = gaussian)
m1.2 <- glm(TotalCephReads~Turbidity,data=edna.meta,family = gaussian)
m1.3 <- glm(TotalCephReads~Depth,data=edna.meta,family = gaussian) 
m1.4 <- glm(TotalCephReads~SST,data=edna.meta,family = gaussian)
m2.1<-glm(TotalCephReads~Sample.ranking,data=edna.meta,family=gaussian)
m2.2<-glm(TotalCephReads~tot.group.size,data=edna.meta,family=gaussian)
m2.3<-glm(TotalCephReads~boat.group.size,data=edna.meta,family=gaussian)
m2.4<-glm(TotalCephReads~Prox2dolph,data=edna.meta,family=gaussian)
m3.1<-glm(TotalCephReads~Amount.collected,data=edna.meta,family = gaussian)
m3.2<-glm(TotalCephReads~Amount.filtered,data=edna.meta,family = gaussian)
m3.3<-glm(TotalCephReads~Time2Filter,data=edna.meta,family = gaussian)

summary(m1.1)
summary(m1.2)
summary(m1.3)
summary(m1.4)
summary(m2.1)
summary(m2.2)
summary(m2.3)
summary(m2.4)
summary(m3.1)
summary(m3.2)
summary(m3.3)


#For exploring the relationships.
boxplot(TotalCephReads~Site,data=edna.meta,ylab = "Total number of reads",ylim = c(0,100000),xlab = "Site", cex.axis = 0.9, cex.lab = 0.9)
text(x=2,y=90000,"*",cex = 1.5)
plot(TotalCephReads~Turbidity,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Secchi disk depth (m)", xlab = "Cephalorhynchus detected",ylim = c(0,10000))
plot(TotalCephReads~Depth,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Water depth (m)", xlab = "Cephalorhynchus detected")
plot(SST~TotalCephReads,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Sea surface temperature (ºC)", xlab = "Cephalorhynchus detected")
#Collection
rbp<-barplot(rank.mat,beside=TRUE,ylab = "Proportion of Samples \n(#/total collected)",xlab = "Sample ranking",ylim=c(0,1), cex.axis = 0.9, cex.lab = 0.9)
axis(1,at=c(2,5,8,11), labels = c("Poor","Ok","Good","Excellent"),cex.axis = 0.8)
abline(h=0)
rbp<-barplot(prox.mat,beside=TRUE,ylab = "Proportion of Samples \n(#/total collected)",xlab = "Proximity to dolphins",ylim=c(0,1), cex.axis = 0.9, cex.lab = 0.9)
axis(1,at=c(2,5,8,11,14,17), labels = c("<1m","1m","2m","3m","4m",">5m"),cex.axis = 0.8)
abline(h=0)
plot(tot.group.size~TotalCephReads,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Total group size", xlab = "Cephalorhynchus detected")
plot(boat.group.size~TotalCephReads,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Group size \nat vessel", xlab = "Cephalorhynchus detected")
#Filtration
plot(Time2Filter~TotalCephReads,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Time from collection \nto filtration (hours)", xlab = "Cephalorhynchus detected")
plot(Amount.filtered~TotalCephReads,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Amount filtered (L)", xlab = "Cephalorhynchus detected")
plot(Amount.collected~TotalCephReads,data=edna.meta, cex.axis = 0.9, cex.lab = 0.9,ylab = "Amount Collected (L)", xlab = "Cephalorhynchus detected")

text(x = -4.5, y = 18.9, "A)", xpd=NA,cex = 1.2)
text(x = -1.5, y = 18.9, "B)",xpd=NA,cex = 1.2)
text(x = 1.5, y = 18.9, "C)",xpd=NA,cex = 1.2)
dev.off()
