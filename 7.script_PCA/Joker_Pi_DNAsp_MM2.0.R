###############################################
### STARTING by giving the path to the data ###
###############################################

library("zoo")
library("grDevices")
library(fields)
library(plyr) 
library(data.table)
library(missMethyl)
library(edgeR)
library("DESeq2")
library("Rcpp")
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library("vsn")
library("gplots")
library("airway")
library("pasilla")
library("Biobase")
library("BiocParallel")
library("ggplot2")
library("reshape")
register(MulticoreParam(3))
library(goseq)
library(GO.db)
library(ape)
library(Biostrings)


```

#####################
### LOAD THE DATA ###
#####################

```
rm(list=ls())
## Give the working folder
PATH_STAT = "~/Desktop/Papers/DnaSP/Dadi/"
REF = "~/Documents/Solution/Transcriptomics/Coverage/SMEL_PerGene/SMEL_V3.2016_11_01.CDS.uniq.fa"
species = "Eggplant"
Colspecies = "#A164D0"
ChrLOC=read.table("/Users/stephaniearnoux/Documents/Solution/references/SMEL_V3_gene_Loc.tab")
setwd(PATH_STAT)

## Give sample details 
pop2 = "S. melongena" 
pop3 = "S. insanum"

Species= "Eggplant"
title = "MM_DNAsp"

CropCol_Light = "#A164D0"
CropCol_Dark = "#5D04A1"

## Settle the file names
StatCropNS = "MM_DNAsp_CROP_NS_dadi.noMiss.VCF.Clean.out" 
StatCropSyn = "MM_DNAsp_CROP_SYN_dadi.noMiss.VCF.Clean.out" 

StatWildNS = "MM_DNAsp_WILD_NS_dadi.noMiss.VCF.Clean.out"
StatWildSyn = "MM_DNAsp_WILD_SYN_dadi.noMiss.VCF.Clean.out"

################### NOW DONT TOUCH BUT RUN IT THROUGH
## Reads the files
Stat_CropNS = read.table(paste(StatCropNS,sep=""),header = TRUE)
colnames(Stat_CropNS)[1] <- "Gene"
Stat_CropSYN = read.table(paste(StatCropSyn,sep=""),header = TRUE)
colnames(Stat_CropSYN)[1] <- "Gene"

Stat_WildNS = read.table(paste(StatWildNS,sep=""),header = TRUE)
colnames(Stat_WildNS)[1] <- "Gene"
Stat_WildSYN = read.table(paste(StatWildSyn,sep=""),header = TRUE)
colnames(Stat_WildSYN)[1] <- "Gene"



```
#####################################
### Sanity Check and data reading ###
#####################################

```
dtCropNS=data.frame(Gene=as.character(Stat_CropNS$Gene),PiCNS=as.numeric(as.character(Stat_CropNS$Pi)))
dtCropSyn=data.frame(Gene=as.character(Stat_CropSYN$Gene),PiCSyn=as.numeric(as.character(Stat_CropSYN$Pi)))
CropPiNS_and_PiSyn = merge(dtCropNS,dtCropSyn, by="Gene")
CropPiNS_and_PiSyn$PiNpiS =  CropPiNS_and_PiSyn$PiCNS/CropPiNS_and_PiSyn$PiCSyn

dtWildNS=data.frame(Gene=as.character(Stat_WildNS$Gene),PiWNS=as.numeric(as.character(Stat_WildNS$Pi)))
dtWildSyn=data.frame(Gene=as.character(Stat_WildSYN$Gene),PiWSyn=as.numeric(as.character(Stat_WildSYN$Pi)))
WildPiNS_and_PiSyn = merge(dtWildNS,dtWildSyn, by="Gene")
WildPiNS_and_PiSyn$PiNpiS =  WildPiNS_and_PiSyn$PiWNS/WildPiNS_and_PiSyn$PiWSyn

Both_Pi = merge(CropPiNS_and_PiSyn,WildPiNS_and_PiSyn, by="Gene")
dtChrLOC = data.table(ChrLOC, key = "V1")

## Get the locations 
pdtCrop <- data.table(CropPiNS_and_PiSyn, key = "Gene")
dtCrop <- pdtCrop[dtChrLOC]
chr <- substr(dtCrop$Gene,7,8)
dtCrop$V2 = as.numeric(as.character(dtCrop$V2))
pFdtCrop <- cbind(dtCrop,chr=as.numeric(as.character(chr)))
#FdtCrop <- na.omit(pFdtCrop)

pdtWild = data.table(WildPiNS_and_PiSyn, key = "Gene") 
dtWild <- pdtWild[dtChrLOC]
chr <- substr(dtWild$Gene,7,8)
dtWild$V2 = as.numeric(as.character(dtWild$V2))
pFdtWild <- cbind(dtWild,chr=as.numeric(as.character(chr)))
#FdtWild <- na.omit(pFdtWild)


FdtCrop = pFdtCrop[order(pFdtCrop$V2),]
FdtCrop = FdtCrop[order(FdtCrop$chr),] 
FdtWild = pFdtWild[order(pFdtWild$V2),]
FdtWild = FdtWild[order(FdtWild$chr),] 


```
########################
### Delta Pi / Taj D ###
########################

```

delta.pi.dist=data.frame(as.character(FdtWild$Gene),as.numeric(FdtWild$chr), as.numeric(FdtCrop$V2),  as.numeric(FdtWild$PiNpiS),as.numeric(FdtCrop$PiNpiS))
colnames(delta.pi.dist) = c("GENE","CHROM","POS", "PiNPiS_Wild","PiNPiS_Crop"); dim(delta.pi.dist); head(delta.pi.dist)
delta.pi.dist = na.omit(delta.pi.dist)

```
#######################
### GET THE EVAL GO ###
#######################

```

d = delta.pi.dist
d$ITAG <- as.character(d$GENE) #just in case...sometime factors act wierd in  

# Need length of each ITAG, because goseq adjusts for this use Biostrings to calculate this

itagSeqs <- readDNAStringSet(file = REF) 
itagLength <- nchar(itagSeqs) #length of each ITAG
names(itagLength) <- names(itagSeqs)

# head(names(itagLength))
names(itagLength) <- substr(names(itagLength),1,20) # not needed if you use seqinR
# head(d$ITAG)
d$ITAG[!d$ITAG %in% names(itagLength)] #Looks OK


colnames(d)

##Filter
#maxPiWild = quantile(d$PiNPiS_Wild, 0.99, na.rm=TRUE)[[1]]
#maxPiCrop = quantile(d$PiNPiS_Crop, 0.57, na.rm=TRUE)[[1]]
maxPiWild = quantile(d$PiNPiS_Wild, 0.9,na.rm=TRUE)[[1]]
maxPiCrop = quantile(d$PiNPiS_Crop, 0.9,na.rm=TRUE)[[1]] 
which((d$PiNPiS_Crop) >= maxPiCrop)
# 11220 14263 15632 // #1074
which((d$PiNPiS_Wild) >= maxPiWild)
# 2523  9125  9966 11944 12120 12515 17304 // #96

##Filter
d$PiNPiS_Crop[(d$PiNPiS_Crop) >= maxPiCrop] <- NA
d$PiNPiS_Wild[(d$PiNPiS_Wild) >= maxPiWild] <- NA
d=na.omit(d)


which(((d$PiNPiS_Crop) <= 0) & ((d$PiNPiS_Wild) <= 0))
# 350

##Filter
d$PiNPiS_Crop[((d$PiNPiS_Crop) <= 0) & ((d$PiNPiS_Wild) <= 0)]<- NA
d=na.omit(d)
##### WE CREATE AN EVALGO FUNCTION THAT PERFORM THE ENRICHMENT TEST #####

gene.names=delta.pi.dist$GENE
PiC=delta.pi.dist$PiNPiS_Crop
PiW=delta.pi.dist$PiNPiS_Wild
PiC.thresh=quantile(d$PiNPiS_Crop, 0.9,na.rm=TRUE)[[1]]
#PiC.threshmin=0.05
PiW.thresh=quantile(d$PiNPiS_Wild, 0.9,na.rm=TRUE)[[1]]
#PiW.threshmin=0.05
ilength=itagLength

#add GO: header if needed


#get length list to match gene names
ilength <- ilength[names(ilength) %in% gene.names]


#filter genes based on criterion
up <- as.integer(PiC > PiC.thresh ) 
# up <- as.integer(PiC > PiC.thresh & PiW < PiW.thresh & PiW > PiW.threshmin)
names(up) <- gene.names #upSelected genes

down <- as.integer(PiW > PiW.thresh) 
# down <- as.integer(PiW > PiW.thresh & PiC < PiC.thresh & PiC > PiC.threshmin) 
names(down) <- gene.names #downregulated genes


length(up[up==1])  #A 196 188 11
length(down[down==1]) #B 374 370 25

########
##### LIST OF A and B Genes !!!
########
OUTGO.A<-paste(title,"_Pi_SelectedC_list.txt",sep="")
write.table(names(up[up==1]), OUTGO.A, row.names = FALSE, col.names = "GENE")
OUTGO.B<-paste(title,"_Pi_SelectedW_list.txt",sep="")
write.table(names(down[down==1]), OUTGO.B, row.names = FALSE, col.names = "GENE")

########
#### FIGURES 
########
imageoutput=file.path("~/Desktop/Papers",paste(title,"Hist_PiNPiS.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 13, units = 'cm',res = 300)
hist(WildPiNS_and_PiSyn$PiNpiS,nclass=80, col='dimgray',main='',xlab=expression(paste(pi,"N/",pi,"S",sep="")))
hist(CropPiNS_and_PiSyn$PiNpiS,nclass=80,add=TRUE, col=Colspecies,main='',xlab='')
title(main= paste(Species),lwd=1.8,lty=5,cex.main=1.4)
dev.off()

imageoutput=file.path("~/Desktop/Papers",paste(title,"Hist_PiNPiS_zoom.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 13, units = 'cm',res = 300)
hist(WildPiNS_and_PiSyn$PiNpiS,ylim = c(0,100), xlim=c(0,20),nclass=1500, col='dimgray',main='',xlab=expression(paste(pi,"N/",pi,"S",sep="")))
hist(CropPiNS_and_PiSyn$PiNpiS,ylim = c(0,100), xlim=c(0,20),nclass=2500,add=TRUE, col=Colspecies,main='',xlab='')
abline(v=PiW.thresh,col="gray20",lwd=1.2,lty=5)
abline(v=PiC.thresh,col=Colspecies,lwd=1.2,lty=5)
title(main= paste(Species),lwd=1.8,lty=5,cex.main=1.4)
dev.off()

imageoutput=file.path("~/Desktop/Papers",paste(title,"Pi_GO_filt.wall.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 10, units = 'cm',res = 300)
par(mar=c(4.2,5.1,3.1,2.1))
plot(d$PiNPiS_Crop~d$PiNPiS_Wild, cex.axis=1,pch=20,cex.lab=1.4,xlab =expression(paste(pi,"N/",pi,"S Wild",sep="")),ylab = expression(paste(pi,"N/",pi,"S Crop",sep="")), col = ifelse((d$PiNPiS_Crop > PiC.thresh | d$PiNPiS_Wild > PiW.thresh  ),'dimgray',Colspecies),las=1)
#xaxt="n", yaxt="n"
#xlim=c(0,0.07), ylim=c(0,0.07), yaxt="n", xlim=c(0,0.1),xaxt="n",
#axis(1, at=c(0,0.02,0.04,0.06,0.08,0.1), cex.axis=0.8,las=1) #xaxis
#axis(2, at=c(0,0.02,0.04,0.06,0.08,0.1), cex.axis=0.8,las=1) #0.8
#axis(1, at=c(0,0.1,0.2,0.3,0.4,0.5), cex.axis=0.8,las=1) #xaxis
#axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5), cex.axis=0.8,las=1) #0.8

abline(h=PiC.thresh,col="gray20",lwd=1.2,lty=5)
abline(v=PiW.thresh,col="gray20",lwd=1.2,lty=5)

title(main= paste(Species),lwd=1.8,lty=5,cex.main=1.4)
#title(main= "Gene nucleotide diversity changes \n due to domestication",lwd=2,lty=5,cex.main=1.4)
dev.off()

imageoutput=file.path("~/Desktop/Papers",paste(title,"Pi_GO_filt.wall_zoom.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 10, units = 'cm',res = 300)
par(mar=c(4.2,5.1,3.1,2.1))
plot(d$PiNPiS_Crop~d$PiNPiS_Wild, cex.axis=1,pch=20,cex.lab=1.4,xlab =expression(paste(pi,"N/",pi,"S Wild",sep="")),ylab = expression(paste(pi,"N/",pi,"S Crop",sep="")), col = ifelse((d$PiNPiS_Crop > PiC.thresh | d$PiNPiS_Wild > PiW.thresh  ),'dimgray',Colspecies),las=1, xlim=c(0,10),ylim=c(0,10))
#xaxt="n", yaxt="n"
#xlim=c(0,0.07), ylim=c(0,0.07), yaxt="n", xlim=c(0,0.1),xaxt="n",
#axis(1, at=c(0,0.02,0.04,0.06,0.08,0.1), cex.axis=0.8,las=1) #xaxis
#axis(2, at=c(0,0.02,0.04,0.06,0.08,0.1), cex.axis=0.8,las=1) #0.8
#axis(1, at=c(0,0.1,0.2,0.3,0.4,0.5), cex.axis=0.8,las=1) #xaxis
#axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5), cex.axis=0.8,las=1) #0.8

abline(h=PiC.thresh,col="gray20",lwd=1.2,lty=5)
abline(v=PiW.thresh,col="gray20",lwd=1.2,lty=5)

title(main= paste(Species),lwd=1.8,lty=5,cex.main=1.4)
#title(main= "Gene nucleotide diversity changes \n due to domestication",lwd=2,lty=5,cex.main=1.4)

dev.off()

ks.test(d$PiNPiS_Crop,d$PiNPiS_Wild, alternative = "greater") #
ks.test(d$PiNPiS_Crop,d$PiNPiS_Wild, alternative = "less") # SAME PiN/PiS


################################################################################################################
################################################################################################################
################################################################################################################

### Sanety check
#listofx = c(FdtCrop, FdtWild)
numchroms=length(unique(pFdtWild$chr))
chrlist=unique(pFdtWild$chr)[-1]
chrlist2=chrlist
x = pFdtWild
BadGENES=NULL
for (i in 1:length(chrlist)){
  Carp = x[chr==i,]
  for (j in 1:(nrow(Carp)-2)){
    if((Carp[j,]$V2 > Carp[j+1,]$V2) & (Carp[j+1,]$V2 < Carp[j+2,]$V2)){
      BadGENES=rbind(BadGENES, Carp[j+1,])
    }else if((Carp[j,]$V2 < Carp[j+1,]$V2) & (Carp[j+1,]$V2 > Carp[j+2,]$V2)){
      BadGENES=rbind(BadGENES, Carp[j+1,])
    }else{
      BadGENES=BadGENES
    }
  }
  for (j in nrow(Carp)){
    if(Carp[j,]$V2 < Carp[j-1,]$V2){
      BadGENES=rbind(BadGENES, Carp[j,])
    }else{
      BadGENES=BadGENES
    }
  }
}
#NULL
head(BadGENES) ; nrow(BadGENES) #217

FdtCrop <- pFdtCrop[!Gene%in%BadGENES$Gene,]
FdtCrop = FdtCrop[order(FdtCrop$V2),]
FdtCrop = FdtCrop[order(FdtCrop$chr),] 
FdtWild <- pFdtWild[!Gene%in%BadGENES$Gene,]
FdtWild = FdtWild[order(FdtWild$V2),]
FdtWild = FdtWild[order(FdtWild$chr),] 

###############################
###############################
###  Stats Pi & Tajima's D  ###
###############################
###############################
## Pi per chr
aggregate(FdtWild$PiWNS~FdtWild$chr, FUN=mean)[,2]
aggregate(FdtWild$PiWSyn~FdtWild$chr, FUN=mean)[,2]
aggregate(FdtCrop$PiCNS~FdtCrop$chr, FUN=mean)[,2]
aggregate(FdtCrop$PiCSyn~FdtCrop$chr, FUN=mean)[,2]
#Global Wild PiN
mean(FdtWild$PiWNS, na.rm=TRUE) # [1]  0.000519216
max(FdtWild$PiWNS, na.rm=TRUE) # [1] 0.2727273
min(FdtWild$PiWNS, na.rm=TRUE) # [1] 0
median(FdtWild$PiWNS, na.rm=TRUE)# [1] 0
#Global Wild PiSyn
mean(FdtWild$PiWSyn, na.rm=TRUE) # [1] 0.000347476
max(FdtWild$PiWSyn, na.rm=TRUE) # [1] 0.2727273
min(FdtWild$PiWSyn, na.rm=TRUE) # [1] 0
median(FdtWild$PiWSyn, na.rm=TRUE)# [1] 0
# Global Crop PiN
mean(FdtCrop$PiCNS, na.rm=TRUE) # [1] 0.0004339988
max(FdtCrop$PiCNS, na.rm=TRUE) # [1] 0.5274725
min(FdtCrop$PiCNS, na.rm=TRUE) # [1] 0
median(FdtCrop$PiCNS, na.rm=TRUE) # [1] 0
# Global Crop PiS
mean(FdtCrop$PiCSyn, na.rm=TRUE) # [1] 0.0002697723
max(FdtCrop$PiCSyn, na.rm=TRUE) # [1] 0.1758242
min(FdtCrop$PiCSyn, na.rm=TRUE) # [1] 0
median(FdtCrop$PiCSyn, na.rm=TRUE) # [1] 0
#t.test(FdtCrop$Pi,FdtWild$Pi)

##################################
##################################
### PiNS NUCLEOTIDE DIVERSITY  ###
##################################
##################################
NSFdtCrop = FdtCrop[!(FdtCrop$PiCNS=='NA'),]
NSFdtWild = FdtWild[!(FdtWild$PiWNS=='NA'),]
Rolling.crop = rollmean(NSFdtCrop$PiCNS,50,align="left",na.rm=TRUE)
data.pi.crop=data.frame(cbind(as.numeric(NSFdtCrop$chr), as.numeric(NSFdtCrop$V2), Rolling.crop))
colnames(data.pi.crop) = c("CHROM","POS","PI"); dim(data.pi.crop); head(data.pi.crop)

Rolling.wild = rollmean(NSFdtWild$PiWNS,50,align="left")
data.pi.wild=data.frame(cbind(as.numeric(NSFdtWild$chr), as.numeric(NSFdtWild$V2), Rolling.wild))
colnames(data.pi.wild) = c("CHROM","POS","PI"); dim(data.pi.wild); head(data.pi.wild)

###############################

CROP=data.pi.crop
WILD=data.pi.wild

###############################
################
## TICKS CROP ##

#attach(CROP)
CROP$pos = NA

# Ticks Creation
ticks=NULL
lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(CROP,CHROM==chrlist[i-1])$POS, 1)
    CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS+lastbase
  }
  ticks=c(ticks, ((((max(CROP[CROP$CHROM==chrlist[i], ]$pos))-lastbase)/2)+lastbase))
}

################
## TICKS WILD ##

#attach(WILD)
WILD$pos = NA

lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(WILD,CHROM==chrlist[i-1])$POS, 1)
    WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS+lastbase
  }
}

###############################

imageoutput=file.path("~/Desktop/Papers",paste(title,"PiNS_GenomeWide.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 23, units = 'cm',res = 300)
layout(matrix(c(1,1,1,1,1,2), nrow = 1, ncol = 6, byrow = TRUE))
par(oma = c(2,2,2,0),mar=c(4.1,4.1,1,0))
with(CROP, plot(pos,rep(0.1,nrow(CROP)),xlim=c(0,max(CROP[CROP$CHROM==chrlist[12], ]$pos)),ylim=c(0,0.0125),main="",bty="n",ylab= expression(paste("Mean nucleotide diversity (", pi,"N)",sep="")),xlab="Chromosomes",xaxt="n",yaxt="n",col="white",cex.lab=1.2))
axis(1, at=ticks, lab=chrlist2, cex.axis=1) #0.9
axis(2, at=c(0,0.0025,0.005,0.0075,0.01,0.0125), cex.axis=1,line = -1.5,las=1) 
#########
#i in even chromosomes and j in uneven chromosmes...
#for (i in c(8))
#{
#  rect(min(WILD[WILD$CHROM==chrlist[i],]$pos)-min(WILD[WILD$CHROM==chrlist[i],]$POS),0.2785,max(WILD[WILD$CHROM==chrlist[i],]$pos),0.28,col="dimgrey",border = NA)
#}
#for (j in c(9))
#{
#  rect(min(WILD[WILD$CHROM==chrlist[j],]$pos)-min(WILD[WILD$CHROM==chrlist[j],]$POS),0.2785,max(WILD[WILD$CHROM==chrlist[j],]$pos),0.28,col="darkgrey",border = NA)
#}

################
## Lines WILD ##

for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[((2*i)-1)],]$pos,WILD[WILD$CHROM==chrlist[((2*i)-1)],]$PI,type="l",col="darkgrey",lwd=1.3)
}
for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[2*i],]$pos,WILD[WILD$CHROM==chrlist[2*i],]$PI,type="l",col="dimgrey",lwd=1.3)
}

################
## Lines Crop ##

for (i in (1:12))
{
  i = (2*i)-1
  points(CROP[CROP$CHROM==chrlist[i],]$pos,CROP[CROP$CHROM==chrlist[i],]$PI,type="l",col=CropCol_Light,lwd=1.3)
} #pch=20
for (i in (1:12))
{
  points(CROP[CROP$CHROM==chrlist[2*i],]$pos,CROP[CROP$CHROM==chrlist[2*i],]$PI,type="l",col=CropCol_Dark,lwd=1.3)
}

#################################
### Pi NUC. DIV. DISTRIBUTION ###
#################################

par(mar=c(4.2,0,1,2.1))
plot(x=density(CROP[!CROP$CHROM==0,]$PI)$y,y=density(CROP[!CROP$CHROM==0,]$PI)$x, col="white",ylim=c(0,0.0125),ylab="",main = "",xlab="Density", xaxt="n",yaxt="n",cex.lab=1,col.axis="white",bty="n")
polygon(x=density(WILD[!WILD$CHROM==0,]$PI)$y,y=density(WILD[!WILD$CHROM==0,]$PI)$x, col=adjustcolor( "darkgrey", alpha.f = 0.8), border="dimgrey",lwd=1.3)
polygon(x=density(CROP[!CROP$CHROM==0,]$PI)$y,y=density(CROP[!CROP$CHROM==0,]$PI)$x, col=adjustcolor( CropCol_Light, alpha.f = 0.8), border=CropCol_Dark,lwd=1.3)

legend("topright", inset=0.05, c(pop2,pop3), fill=c(adjustcolor( CropCol_Light, alpha.f = 0.8),adjustcolor( "darkgrey", alpha.f = 0.8)), border=c(CropCol_Dark,"dimgrey"),bty="n" ,cex=1)

title(paste(Species,sep=""),cex.main=1.8, outer=TRUE)

dev.off()

##################################
##################################
### PiNS NUCLEOTIDE DIVERSITY  ###
##################################
##################################
SFdtCrop = FdtCrop[!(FdtCrop$PiCSyn=='NA'),]
SFdtWild = FdtWild[!(FdtWild$PiWSyn=='NA'),]
Rolling.crop = rollmean(SFdtCrop$PiCSyn,50,align="left",na.rm=TRUE)
data.pi.crop=data.frame(cbind(as.numeric(SFdtCrop$chr), as.numeric(SFdtCrop$V2), Rolling.crop))
colnames(data.pi.crop) = c("CHROM","POS","PI"); dim(data.pi.crop); head(data.pi.crop)

Rolling.wild = rollmean(SFdtWild$PiWSyn,50,align="left")
data.pi.wild=data.frame(cbind(as.numeric(SFdtWild$chr), as.numeric(SFdtWild$V2), Rolling.wild))
colnames(data.pi.wild) = c("CHROM","POS","PI"); dim(data.pi.wild); head(data.pi.wild)

###############################

CROP=data.pi.crop
WILD=data.pi.wild

###############################
################
## TICKS CROP ##

#attach(CROP)
CROP$pos = NA

# Ticks Creation
ticks=NULL
lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(CROP,CHROM==chrlist[i-1])$POS, 1)
    CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS+lastbase
  }
  ticks=c(ticks, ((((max(CROP[CROP$CHROM==chrlist[i], ]$pos))-lastbase)/2)+lastbase))
}

################
## TICKS WILD ##

#attach(WILD)
WILD$pos = NA

lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(WILD,CHROM==chrlist[i-1])$POS, 1)
    WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS+lastbase
  }
}

###############################

imageoutput=file.path("~/Desktop/Papers",paste(title,"PiSyn_GenomeWide.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 23, units = 'cm',res = 300)
layout(matrix(c(1,1,1,1,1,2), nrow = 1, ncol = 6, byrow = TRUE))
par(oma = c(2,2,2,0),mar=c(4.1,4.1,1,0))
with(CROP, plot(pos,rep(0.1,nrow(CROP)),xlim=c(0,max(CROP[CROP$CHROM==chrlist[12], ]$pos)),ylim=c(0,0.0125),main="",bty="n",ylab= expression(paste("Mean nucleotide diversity (", pi,"S)",sep="")),xlab="Chromosomes",xaxt="n",yaxt="n",col="white",cex.lab=1.2))
axis(1, at=ticks, lab=chrlist2, cex.axis=1) #0.9
axis(2, at=c(0,0.0025,0.005,0.0075,0.01,0.0125), cex.axis=1,line = -1.5,las=1) 
#########
#i in even chromosomes and j in uneven chromosmes...
#for (i in c(8))
#{
#  rect(min(WILD[WILD$CHROM==chrlist[i],]$pos)-min(WILD[WILD$CHROM==chrlist[i],]$POS),0.2785,max(WILD[WILD$CHROM==chrlist[i],]$pos),0.28,col="dimgrey",border = NA)
#}
#for (j in c(9))
#{
#  rect(min(WILD[WILD$CHROM==chrlist[j],]$pos)-min(WILD[WILD$CHROM==chrlist[j],]$POS),0.2785,max(WILD[WILD$CHROM==chrlist[j],]$pos),0.28,col="darkgrey",border = NA)
#}

################
## Lines WILD ##

for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[((2*i)-1)],]$pos,WILD[WILD$CHROM==chrlist[((2*i)-1)],]$PI,type="l",col="darkgrey",lwd=1.3)
}
for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[2*i],]$pos,WILD[WILD$CHROM==chrlist[2*i],]$PI,type="l",col="dimgrey",lwd=1.3)
}

################
## Lines Crop ##

for (i in (1:12))
{
  i = (2*i)-1
  points(CROP[CROP$CHROM==chrlist[i],]$pos,CROP[CROP$CHROM==chrlist[i],]$PI,type="l",col=CropCol_Light,lwd=1.3)
} #pch=20
for (i in (1:12))
{
  points(CROP[CROP$CHROM==chrlist[2*i],]$pos,CROP[CROP$CHROM==chrlist[2*i],]$PI,type="l",col=CropCol_Dark,lwd=1.3)
}

#################################
### Pi NUC. DIV. DISTRIBUTION ###
#################################

par(mar=c(4.2,0,1,2.1))
plot(x=density(CROP[!CROP$CHROM==0,]$PI)$y,y=density(CROP[!CROP$CHROM==0,]$PI)$x, col="white",ylim=c(0,0.0125),ylab="",main = "",xlab="Density", xaxt="n",yaxt="n",cex.lab=1,col.axis="white",bty="n")
polygon(x=density(WILD[!WILD$CHROM==0,]$PI)$y,y=density(WILD[!WILD$CHROM==0,]$PI)$x, col=adjustcolor( "darkgrey", alpha.f = 0.8), border="dimgrey",lwd=1.3)
polygon(x=density(CROP[!CROP$CHROM==0,]$PI)$y,y=density(CROP[!CROP$CHROM==0,]$PI)$x, col=adjustcolor( CropCol_Light, alpha.f = 0.8), border=CropCol_Dark,lwd=1.3)

legend("topright", inset=0.05, c(pop2,pop3), fill=c(adjustcolor( CropCol_Light, alpha.f = 0.8),adjustcolor( "darkgrey", alpha.f = 0.8)), border=c(CropCol_Dark,"dimgrey"),bty="n" ,cex=1)

title(paste(Species,sep=""),cex.main=1.8, outer=TRUE)

dev.off()