rm(list = ls())
### After creating a file from plink .ped into PCA.
#install.packages("scatterplot3d")
library(wesanderson)
library(ggplot2)
library(plyr)
library(reshape2)
library("scatterplot3d")


MM_C = "S. melongena" 
MM_W = "S. insanum"
PM_C = "C. annuum" 
PM_W = "C. annuum \n var.glabriusculum"
LA_C = "S. lycopersicum" 
LA_W = "S. pimpinellifolium"
MM_col.list <- c("#A164D0","dimgrey")
PM_col.list <- c("#FF8C00","dimgrey")
LA_col.list <- c("firebrick2","dimgrey")

#Set the environmnent with the path 
WrkingPath = "~/Desktop/PCA_Pruned_LD_Paper"
setwd(WrkingPath)

#####################
###  PCA 2D & 3D  ###
#####################
MMpattern_file = "MM_CSF"
PMpattern_file = "PM_CSF"
LApattern_file = "LA_Nochr9_CSF"

#Import the data
MMtab <- read.table(paste(MMpattern_file,".eigenvec", sep = ""))
MMpc.percent <- read.table(paste(MMpattern_file,".eigenval", sep = ""))
colnames(MMtab) <- c("pop","sample.id","EV1","EV2","EV3")
MMtab <- subset(MMtab,!(MMtab$sample.id %in% c("MM1678","MM0710","MM1789","MM1838","MM1900")))

#Import the data
PMtab <- read.table(paste(PMpattern_file,".eigenvec", sep = ""))
PMpc.percent <- read.table(paste(PMpattern_file,".eigenval", sep = ""))
colnames(PMtab) <- c("pop","sample.id","EV1","EV2","EV3")
PMtab <- subset(PMtab,!(PMtab$sample.id %in% c("PM0647","PM0828")))

LApattern_file = "LA_CSF"
#Import the data
LAtab <- read.table(paste(LApattern_file,".eigenvec", sep = ""))
LApc.percent <- read.table(paste(LApattern_file,".eigenval", sep = ""))
colnames(LAtab) <- c("pop","sample.id","EV1","EV2","EV3")


##########################################
##########################################
###########     Plot    ##################
##########################################
##########################################
max_x <- max(0.8)
min_x <- min(-0.4)
max_y <- max(1)
min_y <- min(-1)
max_z <- max(0.5)
min_z <- min(-1)

## 3D ##
pch.list =c(18,20)
imageoutput=file.path("~/Desktop/PCA_LA3D_avC2.tiff", sep = "")
tiff(file=imageoutput,height = 12, width = 14, units = 'cm',res = 300)
colors <- LA_col.list[as.numeric(LAtab$pop)]
s3d <- scatterplot3d(LAtab$EV1, LAtab$EV2, LAtab$EV3,
                     main = "Tomato",
                     xlab= paste("C1:",format(LApc.percent[1,],digits=3),"%",sep = " "),
                     ylab= paste("C2:",format(LApc.percent[2,],digits=3),"%",sep = " "),
                     zlab= paste("C3:",format(LApc.percent[3,],digits=3),"%",sep = " "),
                     angle=55, 
                     pch = pch.list[as.numeric(LAtab$pop)],
                     color=colors,
                     type= "h",
                     xlim=c(min_x,max_x),
                     ylim=c(min_y,max_y),
                     zlim=c(min_z,max_z),
                     lab=c(4,4),
                     lab.z=c(4),
                     grid=TRUE, box=FALSE, cex.lab = 1.3, cex.axis = 1.2, cex.symbols = 1.4)
legend("topright", legend = c(LA_C,LA_W),col = LA_col.list, pch = pch.list, bty="n", cex = 1.2) 
dev.off()

colors <- MM_col.list[as.numeric(MMtab$pop)]
imageoutput=file.path("~/Desktop/PCA_MM3D_avC2.tiff", sep = "")
tiff(file=imageoutput,height = 12, width = 14, units = 'cm',res = 300)
s3d <- scatterplot3d(MMtab$EV1, MMtab$EV2, MMtab$EV3,
                     main= paste("Eggplant") ,
                     xlab= paste("C1:",format(MMpc.percent[1,],digits=3),"%",sep = " "),
                     ylab= paste("C2:",format(MMpc.percent[2,],digits=3),"%",sep = " "),
                     zlab= paste("C3:",format(MMpc.percent[3,],digits=3),"%",sep = " "),
                     angle=55, 
                     pch = pch.list[as.numeric(MMtab$pop)],
                     color=colors,
                     type= "h",
                     xlim=c(min_x,max_x),
                     ylim=c(min_y,max_y),
                     zlim=c(min_z,max_z),
                     lab=c(4,4),
                     lab.z=c(4),
                     grid=TRUE, box=FALSE, cex.lab = 1.3, cex.axis = 1.2, cex.symbols = 1.4)
legend("topright", legend = c(MM_C,MM_W),col = MM_col.list, pch = pch.list, bty="n", cex = 1.2)  
dev.off()


colors <- PM_col.list[as.numeric(PMtab$pop)]
imageoutput=file.path("~/Desktop/PCA_PM3D_avC2.tiff", sep = "")
tiff(file=imageoutput,height = 12, width = 14, units = 'cm',res = 300)
s3d <- scatterplot3d(PMtab$EV1, PMtab$EV2, PMtab$EV3,
                     main= paste("Pepper") ,
                     xlab= paste("C1:",format(PMpc.percent[1,],digits=3),"%",sep = " "),
                     ylab= paste("C2:",format(PMpc.percent[2,],digits=3),"%",sep = " "),
                     zlab= paste("C3:",format(PMpc.percent[3,],digits=3),"%",sep = " "),
                     angle=55, 
                     pch = pch.list[as.numeric(PMtab$pop)],
                     color=colors,
                     type= "h",
                     xlim=c(min_x,max_x),
                     ylim=c(min_y,max_y),
                     zlim=c(min_z,max_z),
                     lab=c(4,4),
                     lab.z=c(4),
                     grid=TRUE, box=FALSE, cex.lab = 1.3, cex.axis = 1.2, cex.symbols = 1.4)
legend("topleft", legend = c(PM_C,PM_W) ,col = PM_col.list, pch = pch.list, bty="n", cex = 1.2) 
dev.off()

