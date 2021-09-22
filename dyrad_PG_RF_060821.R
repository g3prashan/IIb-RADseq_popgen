# Introductory tutorial for conducting Random Forest on data with a categorical response variable
# (i.e. Random Forest using classification trees)



# The code relies on the package randomForest (Liaw and Wiener 2002)
rm(list = ls())
library(randomForest)
library(dplyr)
library(missForest) #to impute the missing values
library(caret)

# Import the example data set, which comprises 95 individuals genotyped at 19428 biallelic loci, where 0=homozygote 1, 1=heterozygote, 2=homozygote 2
# The objective is to identify loci associated with geography

# class_data <- read.csv("new_data/genotype_data_19428.txt", header = TRUE) #This is the simulated data from Table S8
# head(class_data[1:100])
#
# imputed<-rfImpute(class_data[,-c(1:4)], as.factor(df.train$geography), iter=10, ntree=2000)
# imputedfeatures<-imputed[,-1]
# class_data_imputed <- cbind(class_data[,1:4], imputedfeatures)
# write.csv(class_data_imputed,file="new_data/genotype_data_19428_imputed.csv",row.names=FALSE)

#A vector file for creating test and train  dataset
class_data_imputed <- read.csv("new_data/genotype_data_19428_imputed.csv")

#A vector file for creating test and train  dataset
test <- c(1,5, 9,14, 16,21, 25,29, 32, 38, 42,47, 50, 55, 58, 63, 65,69, 73,79, 82,89, 90,95)

df.train <- class_data_imputed[-test,] #create training dataset
df.validate <- class_data_imputed[test,] #create test dataset


# We will use these sample sizes in conjunction with the strata option (see below)

###########################################################################################################################################
###########################################################################################################################################

# Now run Random Forest analysis.  we need to conduct a classification RF

# First, we need to optimize mtry by running different values of mtry at different values of ntree.

# We will run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p, where p is the number of loci
# We will initially run each of these mtry values at ntree=1000 to 10000 (by increments of 1000).
# We are looking for a plateau where the out-of-bag error rate (OOB-ER) stops decreasing with larger values of ntree
# Once we reach the plateau, we will choose the mtry value that minimizes the OOB-ER.

###########################################################################################################

# This plot shows that mtry=19380 p is the best in terms of OOB-ER , and that the OOB-ER has reached a plateau.
# Therefore, we will use mtry=p for our Random Forest analyses.
# Note that this plot differs from Figure 3 in the manuscript and thus demonstrates that optimal parameter values will vary based on each data set.

###########################################################################################################################################
###########################################################################################################################################
# Now begin the full Random Forest analyses

# Recall that even though we optimized mtry, we must now run a larger number of trees in order to achieve convergence of importance values between forests.
# As a starting point, we will grow 25,000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# since the OOB-ER reached a plateau for a given mtry value after about 400 trees.

# As a starting point, we will grow 150,000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# since the OOB-ER reached a plateau for a given mtry value after about 400 trees.




rf_all_1 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=19380, ntree=400000)


save(rf_all_1,file="new_data/rf_all_1.Rdata")


rf_all_2 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=19380, ntree=400000 )
save(rf_all_2,file="new_data/rf_all_2.Rdata")


#Check correlation of locus importance values between forests
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")


##### 400000 trees gave correlation of 0.8776236 for importance values


cor(importance_rf_all_1,importance_rf_all_2) # A correlation of ???? for locus importance values between forests is extremely good, so we'll use ??? trees for the remaining forests

rf_all_3 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=1942, ntree=400000 )
save(rf_all_3,file="rf_all_3.Rdata")
importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")



############################################################################################################################################
############################################################################################################################################

# The predictive ability of classification trees is measured by the out-of-bag error rate.  An error rate is calculated for each tree within a forest.
# We will use the error rate from the last tree in the forest, which takes all previous trees into account and thus represents the error rate after the model stabilizes/converges

rf_all_1_err.rate <- rf_all_1$err.rate[400000]
rf_all_2_err.rate <- rf_all_2$err.rate[400000]
rf_all_3_err.rate <- rf_all_3$err.rate[400000]

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")

# Export importance values for future reference
write.csv(importance_rf_all,file="new_data/rf_importance_values_all_loci_classification_19428.csv",row.names=FALSE)

############################################################################################################################################
############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of disease resistance.
# For each subset, we will use mtry=p since that is the optimal setting that we previously found.

##### Best 2%

names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes
genotypes_2perc<-df.train[,colnames(df.train) %in% names_best_2perc_unique]

# Now conduct RF on this subset
rf_2perc_1 = randomForest(x = genotypes_2perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=400000 )
 save(rf_2perc_1,file="rf_2perc_1.Rdata")

rf_2perc_2 = randomForest(x = genotypes_2perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=400000 )
 save(rf_2perc_2,file="rf_2perc_2.Rdata")

rf_2perc_3 = randomForest(x = genotypes_2perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=400000 )
 save(rf_2perc_3,file="rf_2perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_2perc_1_err.rate <- rf_2perc_1$err.rate[400000]
rf_2perc_2_err.rate <- rf_2perc_2$err.rate[400000]
rf_2perc_3_err.rate <- rf_2perc_3$err.rate[400000]

rm(rf_2perc_1,rf_2perc_2,rf_2perc_3) # remove the objects to save memory in R

##### Best 3%

names_best_3perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.97))]
names_best_3perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.97))]
names_best_3perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.97))]
names_best_3perc_unique<-unique(c(names_best_3perc_1,names_best_3perc_2,names_best_3perc_3))

# Extract genotypes
genotypes_3perc<-df.train[,colnames(df.train) %in% names_best_3perc_unique]

# Now conduct RF on this subset
rf_3perc_1 = randomForest(x = genotypes_3perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=400000 )
 save(rf_3perc_1,file="rf_3perc_1.Rdata")

rf_3perc_2 = randomForest(x = genotypes_3perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=400000 )
 save(rf_3perc_2,file="rf_3perc_2.Rdata")

rf_3perc_3 = randomForest(x = genotypes_3perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=400000 )
 save(rf_3perc_3,file="rf_3perc_3.Rdata")


# Extract and save the out of bag error rate for this subset of loci
rf_3perc_1_err.rate <- rf_3perc_1$err.rate[400000]
rf_3perc_2_err.rate <- rf_3perc_2$err.rate[400000]
rf_3perc_3_err.rate <- rf_3perc_3$err.rate[400000]

rm(rf_3perc_1,rf_3perc_2,rf_3perc_3) # remove the objects to save memory in R

# ##### Best 4%
#
# names_best_4perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.96))]
# names_best_4perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.96))]
# names_best_4perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.96))]
# names_best_4perc_unique<-unique(c(names_best_4perc_1,names_best_4perc_2,names_best_4perc_3))
#
# # Extract genotypes
# genotypes_4perc<-df.train[,colnames(df.train) %in% names_best_4perc_unique]
# rownames(genotypes_4perc)
#
# # Now conduct RF on this subset
# rf_4perc_1 = randomForest(x = genotypes_4perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000 )
# # <<<<<<< Updated upstream
#  save(rf_4perc_1,file="rf_4perc_1.Rdata")
#
# rf_4perc_2 = randomForest(x = genotypes_4perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000 )
#  save(rf_4perc_2,file="rf_4perc_2.Rdata")
#
# rf_4perc_3 = randomForest(x = genotypes_4perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000 )
#  save(rf_4perc_3,file="rf_4perc_3.Rdata")
#
# # Extract and save the out of bag error rate for this subset of loci
# rf_4perc_1_err.rate <- rf_4perc_1$err.rate[25000]
# rf_4perc_2_err.rate <- rf_4perc_2$err.rate[25000]
# rf_4perc_3_err.rate <- rf_4perc_3$err.rate[25000]
#
# rm(rf_4perc_1,rf_4perc_2,rf_4perc_3) # remove the objects to save memory in R

##### Best 5%

names_best_5perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.95))]
names_best_5perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.95))]
names_best_5perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.95))]
names_best_5perc_unique<-unique(c(names_best_5perc_1,names_best_5perc_2,names_best_5perc_3))

# Extract genotypes
genotypes_5perc<-df.train[,colnames(df.train) %in% names_best_5perc_unique]

# Now conduct RF on this subset
rf_5perc_1 = randomForest(x = genotypes_5perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=400000 )
# <<<<<<< Updated upstream
 save(rf_5perc_1,file="rf_5perc_1.Rdata")

rf_5perc_2 = randomForest(x = genotypes_5perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=400000 )
 save(rf_5perc_2,file="rf_5perc_2.Rdata")

rf_5perc_3 = randomForest(x = genotypes_5perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=400000 )
 save(rf_5perc_3,file="rf_5perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_5perc_1_err.rate <- rf_5perc_1$err.rate[400000]
rf_5perc_2_err.rate <- rf_5perc_2$err.rate[400000]
rf_5perc_3_err.rate <- rf_5perc_3$err.rate[400000]

rm(rf_5perc_1,rf_5perc_2,rf_5perc_3) # remove the objects to save memory in R

# ##### Best 10%
#
# names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
# names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
# names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
# names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))
#
# # Extract genotypes
# genotypes_10perc<-df.train[,colnames(df.train) %in% names_best_10perc_unique]
#
# # Now conduct RF on this subset
# rf_10perc_1 = randomForest(x = genotypes_10perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000 )
#
#  save(rf_10perc_1,file="rf_10perc_1.Rdata")
#
# rf_10perc_2 = randomForest(x = genotypes_10perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000 )
#  save(rf_10perc_2,file="rf_10perc_2.Rdata")
#
# rf_10perc_3 = randomForest(x = genotypes_10perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000 )
#  save(rf_10perc_3,file="rf_10perc_3.Rdata")
#
# # Extract and save the out of bag error rate for this subset of loci
# rf_10perc_1_err.rate <- rf_10perc_1$err.rate[25000]
# rf_10perc_2_err.rate <- rf_10perc_2$err.rate[25000]
# rf_10perc_3_err.rate <- rf_10perc_3$err.rate[25000]
#
# rm(rf_10perc_1,rf_10perc_2,rf_10perc_3) # remove the objects to save memory in R
#
# ##### Best 20%
#
# names_best_20perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.80))]
# names_best_20perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.80))]
# names_best_20perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.80))]
# names_best_20perc_unique<-unique(c(names_best_20perc_1,names_best_20perc_2,names_best_20perc_3))
#
# # Extract genotypes
# genotypes_20perc<-df.train[,colnames(df.train) %in% names_best_20perc_unique]
#
# # Now conduct RF on this subset
# rf_20perc_1 = randomForest(x = genotypes_20perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000 )
#  save(rf_20perc_1,file="rf_20perc_1.Rdata")
#
# rf_20perc_2 = randomForest(x = genotypes_20perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000 )
#  save(rf_20perc_2,file="rf_20perc_2.Rdata")
#
# rf_20perc_3 = randomForest(x = genotypes_20perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000 )
# save(rf_20perc_3,file="rf_20perc_3.Rdata")
#
# # Extract and save the out of bag error rate for this subset of loci
# rf_20perc_1_err.rate <- rf_20perc_1$err.rate[25000]
# rf_20perc_2_err.rate <- rf_20perc_2$err.rate[25000]
# rf_20perc_3_err.rate <- rf_20perc_3$err.rate[25000]
#
# rm(rf_20perc_1,rf_20perc_2,rf_20perc_3) # remove the objects to save memory in R
#
# ##### Best 30%
#
# names_best_30perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.70))]
# names_best_30perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.70))]
# names_best_30perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.70))]
# names_best_30perc_unique<-unique(c(names_best_30perc_1,names_best_30perc_2,names_best_30perc_3))
#
# # Extract genotypes
# genotypes_30perc<-df.train[,colnames(df.train) %in% names_best_30perc_unique]
#
# # Now conduct RF on this subset
# rf_30perc_1 = randomForest(x = genotypes_30perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000 )
# save(rf_30perc_1,file="new_data/rf_30perc_1.Rdata")
#
# rf_30perc_2 = randomForest(x = genotypes_30perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000 )
# save(rf_30perc_2,file="new_data/rf_30perc_2.Rdata")
#
# rf_30perc_3 = randomForest(x = genotypes_30perc, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000 )
# save(rf_30perc_3,file="new_data/rf_30perc_3.Rdata")
#
# # Extract and save the out of bag error rate for this subset of loci
# rf_30perc_1_err.rate <- rf_30perc_1$err.rate[25000]
# rf_30perc_2_err.rate <- rf_30perc_2$err.rate[25000]
# rf_30perc_3_err.rate <- rf_30perc_3$err.rate[25000]
#
# rm(rf_30perc_1,rf_30perc_2,rf_30perc_3) # remove the objects to save memory in R


# # Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach
#
# All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
#                               cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
#                               cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
#                               cbind(rf_4perc_1_err.rate,rf_4perc_2_err.rate,rf_4perc_3_err.rate),
#                               cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate),
#                               cbind(rf_10perc_1_err.rate,rf_10perc_2_err.rate,rf_10perc_3_err.rate),
#                               cbind(rf_20perc_1_err.rate,rf_20perc_2_err.rate,rf_20perc_3_err.rate),
#                               cbind(rf_30perc_1_err.rate,rf_30perc_2_err.rate,rf_30perc_3_err.rate))
#
# # Plot error rates for the various subsets
# All_initial_err.rate<-data.frame(All_initial_err.rate)
# All_initial_err.rate$Number_loci<-c(1942,length(names_best_2perc_unique),length(names_best_3perc_unique),length(names_best_4perc_unique),length(names_best_5perc_unique),length(names_best_10perc_unique),length(names_best_20perc_unique),length(names_best_30perc_unique))
# rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
# All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)
#
# # Write error rates to file for future reference
# write.csv(All_initial_err.rate,file="new_data/All_initial_err_rate_classification_tutorial.csv")


#############COMBINATION####################
# Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach

All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
                              cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
                              cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
                              cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate))

# Plot error rates for the various subsets
All_initial_err.rate<-data.frame(All_initial_err.rate)
All_initial_err.rate$Number_loci<-c(19428,length(names_best_2perc_unique),length(names_best_3perc_unique),
                                    length(names_best_5perc_unique))
rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best5%")
All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_err.rate,file="new_data/All_initial_err_rate_classification_tutorial.csv")


######################COMBINATION###########################################


# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5)


#optimize the number of trees before backward purging








# Based on this table and plot, the best 2% of loci have the lowest error rate

#################### Backward purging approach
names_purging <- names_best_2perc_unique

genotypes_purging<-df.train[,colnames(df.train) %in% names_purging]


#####To check the correlation between runs
# rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
# rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
# importance_rf_purging_1<-data.frame(importance(rf_purging_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
# colnames(importance_rf_purging_1)<-c("importance")
# importance_rf_purging_2<-data.frame(importance(rf_purging_2,type=1))
# colnames(importance_rf_purging_2)<-c("importance")
# cor(importance_rf_purging_1,importance_rf_purging_2) #

######25000 trees gave 94.34% correlation so we will continue with this.


rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
save(rf_purging_1,file="new_data/rf_purging_1.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
save(rf_purging_2,file="new_data/rf_purging_2.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
save(rf_purging_3,file="new_data/rf_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[length(names_purging)]]<-names_purging
error_rate_best<-data.frame(V1=1:length(names_purging),V2=1:length(names_purging),V3=1:length(names_purging))
rownames(error_rate_best)<-1:length(names_purging)
error_rate_best[length(names_purging),] <- c(rf_purging_1$err.rate[25000],rf_purging_2$err.rate[25000],rf_purging_3$err.rate[25000])

start.time<-proc.time()

for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1<-data.frame(importance(rf_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging)
  colnames(imp_purging_1)<-"Mean_Decrease_Accuracy1"
  rownames(imp_purging_2)<-colnames(genotypes_purging)
  colnames(imp_purging_2)<-"Mean_Decrease_Accuracy2"
  rownames(imp_purging_3)<-colnames(genotypes_purging)
  colnames(imp_purging_3)<-"Mean_Decrease_Accuracy3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]]<-names_keep
  genotypes_purging<-df.train[,colnames(df.train) %in% names_keep]
  rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000)
  rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000)
  rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000)
  error_rate_best[length(names_purging)-i,] <- c(rf_purging_1$err.rate[25000],rf_purging_2$err.rate[25000],rf_purging_3$err.rate[25000])
}

stop.time<-proc.time()
run.time<-stop.time -start.time
print(run.time)

error_rate_best$Average<-apply(error_rate_best,1,mean)
write.csv(error_rate_best, file="new_data/Backward_purging_OOB-ER_classification_2_perc_400ktrees-25000_trees.csv") # Save the error rates

# error_rate_best <- read.csv("new_data/Backward_purging_OOB-ER_classification_2_perc.csv")

# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus
plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1,cex.axis=1,pch=16)
# abline(v=c(112, 135), h=  c(0.087, 0.092),col="red")
dev.print(pdf, "error_rate_loci_400000trees_25k_purge.pdf")
dev.off()

# Which group of loci yields the lowest error rate?
which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) #.............. loci have the lowest OOB-ER

error_rate_best$Average[c(206)]

# Export the names of the predictor loci
write.csv(names_all_iterations[[211]],file="new_data/Predictor_loci_classification_211.csv")


class_data_211 <- cbind(class_data_imputed[, 1:4], class_data_imputed[, names_all_iterations[[211]]])

write.csv(class_data_211,file="new_data/211_loci_full_sample_400k.csv")

class_data_206 <- cbind(class_data_imputed[, 1:4], class_data_imputed[, names_all_iterations[[206]]])

write.csv(class_data_206,file="new_data/206_loci_full_sample_400k.csv")



#when making model agaiun from saved data, remember to divide the validate and training data first


