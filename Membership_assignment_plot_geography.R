#This code is based on
#R package [assignPOP v1.1] tutorial
#Alex Chen
#2016-12-26


rm(list = ls())

library(vcfR)
library(poppr)
library(adegenet)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(tidyr)

#read input the VCF file
rubi.VCF.RfGeo <- read.vcfR("data/poparea_91.recode.vcf") #RfGeo dataset
rubi.VCF.FullSNP <- read.vcfR("data/SNP.Anno1.vcf") #FullSNP dataset

#read input metadata file
pop.data <- read.csv2("data/samples_info.csv", sep = ";", header = TRUE)

#converting to genlight object
gl.rubi.RfGeo <- vcfR2genlight(rubi.VCF.RfGeo)
gl.rubi.FullSNP <- vcfR2genlight(rubi.VCF.FullSNP)

#diploid data
ploidy(gl.rubi.RfGeo) <- 2
ploidy(gl.rubi.FullSNP) <- 2

pop(gl.rubi.RfGeo) <- pop.data$geography #assign population to genlight
pop(gl.rubi.FullSNP) <- pop.data$geography

#Discriminant analysis of principal components
set.seed(42)
#pnw.dapc.RfGeo <- dapc(gl.rubi.RfGeo, n.pca = 40, n.da = 10)
temp <- optim.a.score(pnw.dapc.RfGeo) # Optimal number of PCA 7
pnw.dapc.RfSex <- dapc(gl.rubi.RfGeo, n.pca = 7, n.da = 10)

set.seed(42)
#pnw.dapc.FullSNP <- dapc(gl.rubi.FullSNP, n.pca = 95, n.da = 10)
temp <- optim.a.score(pnw.dapc.FullSNP) # Optimal number of PCA 28
pnw.dapc.FullSNP <- dapc(gl.rubi.FullSNP, n.pca = 28, n.da = 10)

#########RfGeo####################
# Extract the data from DAPC list
dapc.results.RfGeo <- as.data.frame(pnw.dapc.RfGeo$posterior)
dapc.results.RfGeo$pop <- pop(gl.rubi.RfGeo)
dapc.results.RfGeo$indNames <- rownames(dapc.results.RfGeo)

#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns.
dapc.results.RfGeo <- pivot_longer(dapc.results.RfGeo, -c(pop, indNames))

#Create a color vector
clr <-c( "#e5191d",  "#3a7db9",  "#54ad49",  "#999999")

#change colnames of the dapc.results.RfGeo
colnames(dapc.results.RfGeo) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#Create a level for plotting
dapc.results.RfGeo$Assigned_Pop <- factor(dapc.results.RfGeo$Assigned_Pop, levels=c("RA", "VE", "MD", "NN"))
dapc.results.RfGeo$Original_Pop <- factor(dapc.results.RfGeo$Original_Pop, levels=c("RA", "VE", "MD", "NN"))

####################FullSNP##############################################
# Extract the data from DAPC list
dapc.results.FullSNP <- as.data.frame(pnw.dapc.FullSNP$posterior)
dapc.results.FullSNP$pop <- pop(gl.rubi.FullSNP)
dapc.results.FullSNP$indNames <- rownames(dapc.results.FullSNP)

#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns.
dapc.results.FullSNP <- pivot_longer(dapc.results.FullSNP, -c(pop, indNames))

#Create a color vector
clr <-c( "#e5191d",  "#3a7db9",  "#54ad49",  "#999999")

#change colnames of the dapc.results
colnames(dapc.results.FullSNP) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#Create a level for plotting
dapc.results.FullSNP$Assigned_Pop <- factor(dapc.results.FullSNP$Assigned_Pop, levels=c("RA", "VE", "MD", "NN"))
dapc.results.FullSNP$Original_Pop <- factor(dapc.results.FullSNP$Original_Pop, levels=c("RA", "VE", "MD", "NN"))



###### Membership plot for FullSNP dataset  ############
p_pop_full <- ggplot(dapc.results.FullSNP, aes(x= Sample, y=Posterior_membership_probability, fill=Assigned_Pop))

a <-  p_pop_full+ geom_bar(stat='identity')+
  xlab("Individuals")+ ylab("Posterior membership probability")+
  scale_fill_manual(values = clr)+
  facet_grid(~Original_Pop, scales = "free")+
  theme(axis.text.y = element_text(size=14), axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), axis.title = element_text(size=20), axis.text.x.top = element_text(size=20))+
  theme(strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))+
  theme(legend.text = element_text(size=14), legend.title = element_text(size=14))+
  labs(fill ="Assigned Pop")

a

ggsave("fig/membership_plot_pop_full_19428.pdf", dpi=300, width =7, height=5)

###### Membership plot for RfGeo ############

p_rfGeo <- ggplot(dapc.results.RfGeo, aes(x= Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
b <- p_rfGeo + geom_bar(stat='identity')+
  xlab("Individuals")+ ylab("Posterior membership probability")+
  scale_fill_manual(values = clr)+
  facet_grid(~Original_Pop, scales = "free")+
  theme(axis.text.y = element_text(size=14), axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), axis.title = element_text(size=20), axis.text.x.top = element_text(size=20))+
  theme(strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))+
  theme(legend.text = element_text(size=14), legend.title = element_text(size=14))+
  labs(fill ="Assigned Pop")
b

ggsave("fig/membership_plot_rfGeo_91loci.pdf", dpi=300, width =7, height=5)


#Find posterior population assignment probability for full dataset
amic <- filter(dapc.results.FullSNP, dapc.results.FullSNP$Original_Pop == dapc.results.FullSNP$Assigned_Pop) #Filter the assigned population probability that is correctly assigned

mean(amic$Posterior_membership_probability)

write.table(amic, "data/assigned_PMP_FullSNP.txt", quote = FALSE, row.names=FALSE)

#Find posterior population assignment probability for RfGeo dataset

bmic <- filter(dapc.results.RfGeo, dapc.results.RfGeo$Original_Pop == dapc.results.RfGeo$Assigned_Pop) #Filter the assigned population probability that is correctly assigned

mean(bmic$Posterior_membership_probability)

write.table(bmic, "data/assigned_PMP_Rfgeo.txt", quote = FALSE, row.names=FALSE)
