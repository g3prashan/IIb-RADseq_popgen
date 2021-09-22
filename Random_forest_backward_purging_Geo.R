# Conducting Random Forest on SNP data with a categorical response variable
# (i.e. Random Forest using classification trees)

# The code relies on the package randomForest (Liaw and Wiener 2002)
#Additionally this cold follows the random forest tutorial by Breiuc et. al, 2018. https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12773

rm(list = ls()) #clear the memory
library(randomForest)
library(dplyr)
library(caret)

# Import the example data set, which comprises 95 individuals genotyped at 19428 biallelic loci, where 0=homozygote 1, 1=heterozygote, 2=homozygote 2
# The objective is to identify loci associated with geography

full_snpdata <- read.csv("final_analysis_pop/genotype_data_19428.txt", header = TRUE) #This is the full SNP data from VCF file
head(full_snpdata[1:100])

####PRE-PROCESSING OF DATA###############################################################################################################
###########################################################################################################################################
#Divide the training and test data

set.seed(1234)
trainIndex <- createDataPartition(full_snpdata$geography, p = .8,
                                  list = FALSE,
                                  times = 1)
write.csv(trainIndex, "trainindex.csv")

snp.train <- full_snpdata[trainIndex,] #create training dataset
snp.test <- full_snpdata[-trainIndex,] #create test dataset
dim(df.test)


#We will impute the missing values in the training dataset using rfImpute function using geography as dependent variable#######

imputed<-rfImpute(snp.train[, -c(1:4)], as.factor(df.train$geography), iter=10, ntree=2000) #imputation of mising values in training data
imputedfeatures<-imputed[,-1] #removal of first column
df.train <- cbind(snp.train[,1:4], imputedfeatures) #binding training data with sample features
write.csv(df.train,file="final_analysis_pop/training_data_19428_imputed.csv",row.names=FALSE) #write a imputed dataset file
write.csv(df.test, file ="final_analysis_pop/test_data_19428.csv", row.names = FALSE) #write a training dataset file for testing purpose

#read the file for training data that has been imputed
df.train <- read.csv("final_analysis_pop/training_data_19428_imputed.csv")

###########################################################################################################################################
###########################################################################################################################################

# Tuning of mtry value for the best OOB Error

set.seed(123)

for( i in c(35, 70, 140, 1400, 3000, 5000, 10000, 19380, 19432)){
  tuneRF(df.train[, -c(1:4)], as.factor(df.train$geography), mtryStart = i)
}

# This shows that mtry=19380 is the best in terms of OOB-ER

###########################################################################################################################################
###########################################################################################################################################

# Now begin the full Random Forest analyses

# As a starting point, we will grow 150,000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# since the OOB-ER reached a plateau for a given mtry value after about 400 trees.
# 150000 trees gave correlation of 51%
#We had to increase the number of trees all the way to 400000 to get a correlation of 87.8%

rf_all_1 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=19380, ntree=400000)
save(rf_all_1,file="final_analysis_pop/rf_all_1.Rdata")

rf_all_2 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=19380, ntree=400000 )
save(rf_all_2,file="final_analysis_pop/rf_all_2.Rdata")

#Check correlation of locus importance values between forests
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2)
# A correlation of 87.8% for locus importance values between forests is good
##### 400000 trees gave correlation of 0.8776236 for importance values

rf_all_3 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=19380, ntree=400000 )
save(rf_all_3,file="final_analysis_pop/rf_all_3.Rdata")
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
write.csv(importance_rf_all,file="final_analysis_pop/rf_importance_values_all_loci_classification_19428.csv",row.names=FALSE)

############################################################################################################################################
############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of population differentiation.
# For each subset, we will use mtry=19380 since that is the optimal setting that we previously found.

##### Best 2%

names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes
genotypes_2perc<-df.train[,colnames(df.train) %in% names_best_2perc_unique] #select the loci from training dataset which are best 2%

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

#Due to computational limitation we only checked for until 5% of the best loci....

#############COMBINATION####################
# Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach

All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
                              cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
                              cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
                              cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate))

#Combine error rate of all random forest runs above

# Plot error rates for the various subsets
All_initial_err.rate<-data.frame(All_initial_err.rate)
All_initial_err.rate$Number_loci<-c(19428,length(names_best_2perc_unique),length(names_best_3perc_unique),
                                    length(names_best_5perc_unique))
rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best5%")
All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_err.rate,file="final_analysis_pop/All_initial_err_rate_classification.csv")
######################COMBINATION###########################################


All_initial_err.rate <- read.csv("final_analysis_pop/All_initial_err_rate_classification.csv")

# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5)
# Based on this table and plot, the best 2% of loci have the lowest error rate

#################### Backward purging approach###########################

names_purging <- names_best_2perc_unique
genotypes_purging<-df.train[,colnames(df.train) %in% names_purging]


##############optimize the number of trees before backward purging####################

#####To check the correlation between runs
# rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
# rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
# importance_rf_purging_1<-data.frame(importance(rf_purging_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
# colnames(importance_rf_purging_1)<-c("importance")
# importance_rf_purging_2<-data.frame(importance(rf_purging_2,type=1))
# colnames(importance_rf_purging_2)<-c("importance")
# cor(importance_rf_purging_1,importance_rf_purging_2) #
######25000 trees gave 94.34% correlation so we will continue with this#################
##############optimize the number of trees before backward purging####################



############################################################################################################
############################ We continue with 25000 trees hereon ###########################################

rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
save(rf_purging_1,file="final_analysis_pop/rf_purging_1.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
save(rf_purging_2,file="final_analysis_pop/rf_purging_2.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(df.train$geography), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000 )
save(rf_purging_3,file="final_analysis_pop/rf_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[length(names_purging)]]<-names_purging
error_rate_best<-data.frame(V1=1:length(names_purging),V2=1:length(names_purging),V3=1:length(names_purging))
rownames(error_rate_best)<-1:length(names_purging)
error_rate_best[length(names_purging),] <- c(rf_purging_1$err.rate[25000],rf_purging_2$err.rate[25000],rf_purging_3$err.rate[25000])


################Backward Purging #############################
#Donot stop and run otherwise it will present a error#########

start.time<-proc.time()

for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1<-data.frame(importance(rf_purging_1,type=1)) #type = 1 means Mean decrease in accuracy.
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

error_rate_best$Average <- apply(error_rate_best,1,mean)

write.csv(error_rate_best, file="final_analysis_pop/Backward_purging_OOB-ER_classification_2_perc_400ktrees-25000_trees_purging.csv") # Save the error rates

# error_rate_best <- read.csv("final_analysis_pop/Backward_purging_OOB-ER_classification_2_perc.csv")

# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus

plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1,cex.axis=1,pch=16)

#Print using R graphics
dev.print(pdf, "error_rate_loci_400000trees_25k_purge.pdf")
dev.off()

# Which group of loci yields the lowest error rate?

which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) #.............. loci have the lowest OOB-ER

#91 loci gave the best error rate

error_rate_best$Average[c(91)]

# Export the names of the predictor loci
write.csv(names_all_iterations[[91]],file="final_analysis_pop/Predictor_loci_classification_91.csv")


train_91 <- cbind(df.train[, 1:4], df.train[, names_all_iterations[[91]]])
test_91 <- cbind(df.test[, 1:4], df.test[, names_all_iterations[[91]]])

write.csv(train_91,file="final_analysis_pop/91_loci_full_sample_train.csv")
write.csv(test_91,file="final_analysis_pop/91_loci_full_sample_test.csv")

#train data with no imputed values
train_91_NA <- cbind(snp.train[, 1:4], snp.train[, names_all_iterations[[91]]])
write.csv(train_91_NA,file="final_analysis_pop/91_loci_full_sample_train_NA.csv")


#dataset with all samples and NA values #Need this for UMAP
all_samples_91 <- cbind(full_snpdata[, 1:4], full_snpdata[, names_all_iterations[[91]]])
write.csv(all_samples_91,file="final_analysis_pop/all_samples_91_loci_with_NA.csv")












