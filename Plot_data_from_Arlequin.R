rm(list = ls())

#Load the library tidyverse
library(tidyverse)
library(vcfR) #to get the data and locus
library(ggrepel)
library(splines) #needed for quantile lines

### Let's make a character vector of locus names with Chromosome ID and Position:

vcf <- read.vcfR("data/SNP.Anno1.vcf", verbose=FALSE)
locinames <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep="_") #Extract Contig name and SNP position and bind with _
SNPb <- as.data.frame(locinames)

pop4_arleq <- read.delim("data/geo_fdist2_ObsOut.txt", header = TRUE) #Read output data from Arlequin
pop4_arleq <- pop4_arleq[-6]  #remove 6th column
pop4_arleq <- cbind(SNPb, pop4_arleq)
header <- c("SNP", "Locus","heterozygosity", "fst", "pval", "fstquantile") #specify name of the column
names(pop4_arleq) <- header #assign names to column
pval_adjusted <- p.adjust(pop4_arleq$pval, method = "BH") #adjust the pvalue usingBH correction
pop4_arleq  <- cbind(pop4_arleq, pval_adjusted) #add adjusted Pvalue to dataframe


#####################Now we plot the data #########################################

p <- ggplot(data = pop4_arleq, mapping=aes(x=heterozygosity, y=fst))
p
 p+
  geom_point(position="jitter", aes(colour=cut(pval_adjusted, c(-Inf, 0.05, Inf))))+
  geom_label_repel(aes(label = ifelse(pval_adjusted < 0.01, Locus, '')),
                   box.padding=unit(0.7, "lines"), point.padding=unit(0.1, "lines"), segment.color="grey50",
                                     size=6)+
  labs( x = "Heterozygosity", y = "Fixation index (Fst)")+
  theme_bw()+
  scale_color_manual(name = "Significance:",
                     values = c("(-Inf,0.05]" = "red",
                                "(0.05, Inf]" = "black"),  labels = c("p <= 0.05", "p > 0.05") )+
  theme(legend.text=element_text(size=24), legend.title = element_text(size=24))+
  theme(axis.title.x=element_text(size=30))+
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))+
  theme(axis.title.y=element_text(size=30))+
    stat_quantile(formula = y ~ bs(x, 11), quantiles = 0.95, position ="identity", method = "rq",
              show.legend = F, inherit.aes = TRUE, linetype= "dashed", color="green", size=0.8)+
  stat_quantile(formula = y ~ bs(x, 11), quantiles = 0.5, position ="identity", method = "rq",
                show.legend = F, inherit.aes = TRUE, linetype= "dashed", color="yellow", size=0.8)

ggsave("fig/Arlequin_loci_scan_geo2.pdf", dpi=300, width=11, height=8)

################################Plotting finished#######################################################################


#Get the SNPs that are under directional selection###################

pdiversifying <- pop4_arleq[pop4_arleq$pval_adjusted < 0.05,]

write.csv(pdiversifying, "data/geo_diversifyingsnps_fromArlequin.txt")

sessionInfo()

