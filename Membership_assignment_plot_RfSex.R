rm(list = ls())

library(vcfR)
library(poppr)
library(adegenet)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(tidyr)

#read input the VCF file
rubi.VCF.RfSex <- read.vcfR("data/snp_14.recode.vcf") #RfSex dataset
rubi.VCF.FullSNP <- read.vcfR("data/SNP.Anno1.vcf") #FullSNP dataset

#read input metadata file
pop.data <- read.csv2("data/samples_info.csv", sep = ";", header = TRUE)

#converting to genlight object
gl.rubi.RfSex <- vcfR2genlight(rubi.VCF.RfSex)
gl.rubi.FullSNP <- vcfR2genlight(rubi.VCF.FullSNP)

#diploid data
ploidy(gl.rubi.RfSex) <- 2
ploidy(gl.rubi.FullSNP) <- 2

pop(gl.rubi.RfSex) <- pop.data$sex #assign population to genlight
pop(gl.rubi.FullSNP) <- pop.data$sex

#Discriminant analysis of principal components
set.seed(42)
pnw.dapc.RfSex <- dapc(gl.rubi.RfSex, n.pca = 40, n.da = 10)
set.seed(42)
pnw.dapc.FullSNP <- dapc(gl.rubi.FullSNP, n.pca = 40, n.da = 10)


#########RfSex####################
# Extract the data from DAPC list
dapc.results.RfSex <- as.data.frame(pnw.dapc.RfSex$posterior)
dapc.results.RfSex$pop <- pop(gl.rubi.RfSex)
dapc.results.RfSex$indNames <- rownames(dapc.results.RfSex)

#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns.
dapc.results.RfSex <- pivot_longer(dapc.results.RfSex, -c(pop, indNames))

#change colnames of the dapc.results.RfSex
colnames(dapc.results.RfSex) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#Create a level for plotting
dapc.results.RfSex$Assigned_Pop <- factor(dapc.results.RfSex$Assigned_Pop, levels=c("M", "F"))
dapc.results.RfSex$Original_Pop <- factor(dapc.results.RfSex$Original_Pop, levels=c("M", "F"))

####################FullSNP##############################################
# Extract the data from DAPC list
dapc.results.FullSNP <- as.data.frame(pnw.dapc.FullSNP$posterior)
dapc.results.FullSNP$pop <- pop(gl.rubi.FullSNP)
dapc.results.FullSNP$indNames <- rownames(dapc.results.FullSNP)

#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns.
dapc.results.FullSNP <- pivot_longer(dapc.results.FullSNP, -c(pop, indNames))

#change colnames of the dapc.results
colnames(dapc.results.FullSNP) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#Create a level for plotting
dapc.results.FullSNP$Assigned_Pop <- factor(dapc.results.FullSNP$Assigned_Pop, levels=c("M", "F"))
dapc.results.FullSNP$Original_Pop <- factor(dapc.results.FullSNP$Original_Pop, levels=c("M", "F"))


clr <-c(  "#00bfc4", "#f8766d" )


p_sex_full <- ggplot(dapc.results.FullSNP, aes(x= Sample, y=Posterior_membership_probability, fill=Assigned_Pop))

c <- p_sex_full+ geom_bar(stat='identity')+
  xlab("Individuals")+ ylab("Posterior membership probability")+
  scale_fill_manual(values = clr)+
  facet_grid(~Original_Pop, scales = "free")+
  theme(axis.text.y = element_text(size=14), axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), axis.title = element_text(size=20), axis.text.x.top = element_text(size=20))+
  theme(strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="lightyellow", colour="black",size=1))+
  theme(legend.text = element_text(size=14), legend.title = element_text(size=14))+
  labs(fill ="Assigned Sex")
c
ggsave("fig/membership_plot_sex_full_19428.pdf", dpi=300, width =7, height=5)



p_rfSex <- ggplot(dapc.results.RfSex, aes(x= Sample, y=Posterior_membership_probability, fill=Assigned_Pop))

d <- p_rfSex+ geom_bar(stat='identity')+
  xlab("Individuals")+ ylab("Posterior membership probability")+
  scale_fill_manual(values = clr)+
  facet_grid(~Original_Pop, scales = "free")+
  theme(axis.text.y = element_text(size=14), axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), axis.title = element_text(size=20), axis.text.x.top = element_text(size=20))+
  theme(strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="lightyellow", colour="black",size=1))+
  theme(legend.text = element_text(size=14), legend.title = element_text(size=14))+
  labs(fill ="Assigned Sex")
d

ggsave("fig/membership_plot_rfSex_14loci.pdf", dpi=300, width =7, height=5)


#Find posterior population assignment probability for full dataset
cmic <- filter(dapc.results.FullSNP, dapc.results.FullSNP$Original_Pop == dapc.results.FullSNP$Assigned_Pop) #Filter the assigned population probability that is correctly assigned

mean(cmic$Posterior_membership_probability)

write.table(cmic, "data/assigned_PMP_FullSNP_sex.txt", quote = FALSE, row.names=FALSE)

#Find posterior population assignment probability for RfGeo dataset

dmic <- filter(dapc.results.RfSex, dapc.results.RfSex$Original_Pop == dapc.results.RfSex$Assigned_Pop) #Filter the assigned population probability that is correctly assigned

mean(dmic$Posterior_membership_probability)

write.table(dmic, "data/assigned_PMP_Rfgeo.txt", quote = FALSE, row.names=FALSE)
