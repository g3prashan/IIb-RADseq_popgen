#The code is derived from https://pubmed.ncbi.nlm.nih.gov/27543860/ 
#Benestan et al. 2016

rm(list = ls())

BiocManager::install("qvalue")

library(qvalue)
library(vcfR)
library(tidyverse)
library(ggrepel)

vcf <- read.vcfR("data/SNP.Anno1.vcf", verbose=FALSE)

### Let's make a character vector of locus names with Chromosome ID and Position:

locinames <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep="_")
SNPb <- as.data.frame(locinames)

bayescan=read.table("data/geo_bayescanLongRunOD10_fst.txt") #read output file from bayescan 

#Merge the names of the outliers with the results from the bayescan dataframe.

bayescan <- cbind(SNPb, bayescan)

#Rename columns.
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST")
head(bayescan)

write.table(bayescan, "data/lsal_bayescan_results_p10_geo.txt", quote=FALSE, sep="\t", row.names=FALSE) #save the attached dataframe for file for future use if needed

detach(bayescan) #only use if the bayescan is attached previously. 

attach(bayescan)

class(bayescan$Q_VALUE) #numeric is correct
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) #if not numeric change to niumeric

bayescan[bayescan$Q_VALUE <= 0.0001,"Q_VALUE"] = 0.0001  #if vaklue is 0 change to  0.0001

#Round the values.

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4))
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4))
bayescan$ALPHA <- (round(bayescan$ALPHA, 4))
bayescan$FST <- (round(bayescan$FST, 6))

#Add a column for the type of selection grouping based on a Q-VALUE < 0.05. You can also choose a Q-VALUE < 0.01 if you want to be more conservative.

bayescan$SELECTION <- ifelse(bayescan$ALPHA >= 0 & bayescan$Q_VALUE <= 0.05,"diversifying",
                             ifelse(bayescan$ALPHA >= 0 & bayescan$Q_VALUE > 0.05,"neutral", "balancing"))
bayescan$SELECTION <- factor(bayescan$SELECTION)
levels(bayescan$SELECTION)

#Save the results of the SNPs potentially under positive (divergent) and balancing selection (qvalue < 0.05).

positive <- bayescan[bayescan$SELECTION=="diversifying",]
neutral <- bayescan[bayescan$SELECTION=="neutral",]
balancing <- bayescan[bayescan$SELECTION=="balancing",]

#Check the number of SNPs belonging to each category.
xtabs(data=bayescan, ~SELECTION)

#Write the results of the SNPs potentially under selection (qvalue < 0.05).

write.table(neutral, "data/neutral_pro10.txt", row.names=F, quote=F)
write.table(balancing, "data/balancingpro10.txt", row.names=F, quote=F)
write.table(positive, "data/positive_pro10.txt", row.names=F, quote=F)

#Transformation Log of the Q value in order to create the ggplot graph.

range(bayescan$Q_VALUE)

bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) #log transformation
head(bayescan)

#Create a title for the ggplot graph.
x_title="Log(q-value)"
y_title="Fixation index (Fst)"

#Make the ggplot graph.

ggplot(bayescan,aes(x=LOG10_Q, y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+
  scale_fill_manual(name="Selection",values=c("black","red","blue"))+
  labs(x=x_title)+
  labs(y=y_title)+
  theme_bw()+
  geom_label_repel(aes(label = ifelse(LOG10_Q > 1.3, SNP, '')),
                   box.padding=0.8, point.padding=0.1, segment.color="grey50", size=6, nudge_y =0 )+
  theme(axis.title = element_text(size=24), legend.text = element_text(size=16), legend.title = element_text(size=16),
        axis.text = element_text(size=18) )



ggsave("fig/bayescan_prOdds=10_geo.pdf", dpi=600, width=7, height=5)
dev.off()
