## plot the AIC
rm(list = ls()) 
setwd("/path/to/Results/from/DADI_inferrence/")
list.files()
SpeZ = "Species_**"
AIC = read.table("File_**_sumstat.txt")
AIC = unique(AIC)
names(AIC) = c("Run_ID","Loglkd","AIC")
AIC$Models = sub("(.*)_2017_.*", "\\1", AIC$Run_ID)

#AIC$Models = ordered(AIC$Models, levels=c("SI_C", "IM_C", "SI_BcCw","IM_BcCw", "IM_C_2M","IM_BcCw_2M","SI_C_2N", "IM_C_2N","SI_BcCw_2N","IM_BcCw_2N","IM_C_2N_2M","IM_BcCw_2N_2M"))
#AIC$Models = ordered(AIC$Models, levels=c("SI_C", "IM_C", "SI_BcCw","IM_BcCw", "IM_C_2M","IM_BcCw_2M" ))

png(filename = paste(SpeZ,'LOG.png',sep = "_"),width = 1760, height = 780, res=100)
boxplot(AIC$Loglkd ~ AIC$Models, outpch = NA, border = "steelblue", par(mar = c(10, 5, 4, 2)+ 0.1),las = 2) #, ylim=c(-10000,0)) #ylab = "Log likelihood",las=2)
stripchart(Loglkd ~ Models, data = AIC, 
           vertical = TRUE, method = "jitter", 
           pch = 21, bg = "turquoise4", col = "darkblue", 
           add = TRUE) 
mtext("Models", side = 1, line = 8, cex = 1.5, font = 3)
mtext("Log likelihood", side = 2, line = 3.8, cex = 1.5, font = 3)
title(paste('dadi Log likelihood model convergence in',SpeZ,spe=" "))
dev.off()


results=c()
#for (Mod in  levels(AIC$Models)){
for (Mod in  AIC$Models){
  x = AIC[AIC$Models==Mod,]
  print(x)
  for (i in x$Loglkd){
    imin = i - 2
    imax = i + 2
    for (j in x$Loglkd){
      if (j != i & imin <= j & j <= imax) {
        res = c(j,i)
        Nresults = rbind(c(Mod,res))
        #print(Nresults)
      } else { Nresults = NULL}
      }
  results=rbind(results,Nresults)
  }
}
#resultsF =as.data.frame(results)
resultsF =as.data.frame(unique(results))
names(resultsF) = c("Models","LOGi","LOGj")
Res_log = aggregate(LOGj ~ Models+LOGi, data = resultsF, FUN = function(x){NROW(x)})

Res_log = Res_log[order(Res_log$Models),]

preRes_AIC = aggregate(AIC ~ Models, data = AIC, FUN = function(x){min(x)})
Res_AIC = c()
for (i in preRes_AIC$AIC){
  res = AIC[AIC$AIC == i,]
  Res_AIC=rbind(Res_AIC,res)
}
preResBAD_AIC = aggregate(AIC ~ Models, data = AIC, FUN = function(x){max(x)})
ReBADs_AIC = c()
for (i in preResBAD_AIC$AIC){
  res = AIC[AIC$AIC == i,]
  ReBADs_AIC=rbind(ReBADs_AIC,res)
}
#Res_AIC$Models = ordered(Res_AIC$Models, levels=c("SI_C", "IM_C", "SI_BcCw","IM_BcCw", "IM_C_2M","IM_BcCw_2M","SI_C_2N", "IM_C_2N","SI_BcCw_2N","IM_BcCw_2N","IM_C_2N_2M","IM_BcCw_2N_2M"))
#Res_AIC$Models = ordered(Res_AIC$Models, levels=c("SI_C", "IM_C", "SI_BcCw","IM_BcCw", "IM_C_2M","IM_BcCw_2M" ))
#boxplot( Res_AIC$AIC ~  levels(droplevels(Res_AIC$Models)), border = "steelblue", ylab = "AIC", xlab = "Models" )

png(filename = paste(SpeZ,'AIC.png',sep = "_"),width = 1760, height = 780, res=100)
boxplot( Res_AIC$AIC ~  Res_AIC$Models, border = "steelblue", par(mar = c(12, 5, 4, 2)+ 0.1),las = 2) #ylab = "AIC", xlab = "Models" ,
mtext("Models", side = 1, line = 8, cex = 1.5, font = 3)
mtext("AIC", side = 2, line = 3.8, cex = 1.5, font = 3)
title(paste('dadi model AIC comparison in', SpeZ, sep=" "))
dev.off()

Res_AIC
write.table( Res_AIC$Run_ID, paste("Best_AIC_Models_",SpeZ,".txt",sep=""), col.names = FALSE, quote=FALSE, row.names = FALSE)

write.table( Res_AIC$Run_ID, paste("Worst_AIC_Models_",SpeZ,".txt",sep=""), col.names = FALSE, quote=FALSE, row.names = FALSE)

# Run_ID    Loglkd       AIC     Models
# 8      IM_BcCw_2017_9_26_184248 -4531.191  9072.381    IM_BcCw
# 38  IM_BcCw_2M_2017_9_26_193624 -6026.614 12067.229 IM_BcCw_2M

library(dplyr)
options(dplyr.print_max = 1e9)
set.seed(1)
Num = AIC %>% 
  group_by(Models) %>%
  summarise(no_rows = length(Models))
#survey <- group_by(AIC, Models)
#summarize(survey, no_rows = length(Models))

Res_AIC
Res_AIC[,c(3,4)]
Num


