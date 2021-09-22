#Conducting Random Forest on data with a categorical response variable Sex

# The code relies on the package randomForest (Liaw and Wiener 2002)
rm(list = ls())

library(randomForest)
library(dplyr)

# Import the example data set, which comprises 95 individuals genotyped at 19428 biallelic loci, where 0=homozygote 1, 1=heterozygote, 2=homozygote 2
# The objective is to identify loci associated with sex differentiation

sex_data <- read.csv("new_data/genotype_data_19428.txt", header = TRUE)


####PRE-PROCESSING OF DATA###############################################################################################################
###########################################################################################################################################

#imputing the NA values

imputed<-rfImpute(sex_data[,-c(1:4)], as.factor(sex_data$sex), iter=10, ntree=2000)
imputedfeatures<-imputed[,-1]
sex_data_imputed <- cbind(sex_data[,1:4], imputedfeatures)
write.csv(sex_data_imputed,file="sex_data/sex_genotype_data_19428_imputed.csv",row.names=FALSE)

sex_data_imputed <- read.csv("sex_data/sex_genotype_data_19428_imputed.csv")
head(sex_data_imputed[1:6])

#A vector file for creating test and train  dataset

test_sex <- c(4, 8, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93 )

df.train <- sex_data_imputed[-test_sex,] #create training dataset
df.validate <- sex_data_imputed[test_sex,] #create test dataset


###########################################################################################################################################
###########################################################################################################################################

# Tuning of mtry value for the best OOB Error

for( i in c(35, 70, 140, 1943, 3886, 7772, 10000, 19380, 19432)){
  tuneRF(x = df.train[,5:19432], y = as.factor(df.train$sex), importance=TRUE, mtryStart = i)
}

# This shows that mtry=1943is the best in terms of OOB-ER

###########################################################################################################################################
###########################################################################################################################################
# Now begin the full Random Forest analyses

#As a starting point, we will grow 25000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# 40000 trees gave correlation of 87.06 %

rf_all_1 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=1943, ntree=40000))
save(rf_all_1,file="sex_data/rf_all_1_sex.Rdata")

rf_all_2 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=1943, ntree=40000 )
save(rf_all_2,file="sex_data/rf_all_2_sex.Rdata")

#Check correlation of locus importance values between forests
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2) # A correlation of 87% for locus importance values between forests is extremely good, so we'll use 25,000 trees for the remaining forests

rf_all_3 = randomForest(x = df.train[,5:19432], y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=3886, ntree=25000)
save(rf_all_3,file="sex_data/rf_all_3_sex.Rdata")

importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")

############################################################################################################################################
############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of sex differentiation.
# For each subset, we will use mtry=1943 since that is the optimal setting that we previously found.

rf_all_1_err.rate <- rf_all_1$err.rate[40000]
rf_all_2_err.rate <- rf_all_2$err.rate[40000]
rf_all_3_err.rate <- rf_all_3$err.rate[40000]

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")

# Export importance values for future reference
write.csv(importance_rf_all,file="sex_data/rf_importance_values_all_loci_classification_19428_sex.csv",row.names=FALSE)

############################################################################################################################################
############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of sex difference.
# For each subset, we will use mtry=p since that is the optimal setting that we previously found.

##### Best 2%

names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes
genotypes_2perc<-df.train[,colnames(df.train) %in% names_best_2perc_unique]

# Now conduct RF on this subset
rf_2perc_1 = randomForest(x = genotypes_2perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000 )
save(rf_2perc_1,file="sex_data/rf_2perc_1_sex.Rdata")

rf_2perc_2 = randomForest(x = genotypes_2perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000 )
save(rf_2perc_2,file="sex_data/rf_2perc_2.Rdata")

rf_2perc_3 = randomForest(x = genotypes_2perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=25000 )
save(rf_2perc_3,file="sex_data/rf_2perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_2perc_1_err.rate <- rf_2perc_1$err.rate[25000]
rf_2perc_2_err.rate <- rf_2perc_2$err.rate[25000]
rf_2perc_3_err.rate <- rf_2perc_3$err.rate[25000]

rm(rf_2perc_1,rf_2perc_2,rf_2perc_3) # remove the objects to save memory in R


##### Best 3%

names_best_3perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.97))]
names_best_3perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.97))]
names_best_3perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.97))]
names_best_3perc_unique<-unique(c(names_best_3perc_1,names_best_3perc_2,names_best_3perc_3))


# Extract genotypes
genotypes_3perc<-df.train[,colnames(df.train) %in% names_best_3perc_unique]

# Now conduct RF on this subset
rf_3perc_1 = randomForest(x = genotypes_3perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000  )
save(rf_3perc_1,file="sex_data/rf_3perc_1_sex.Rdata")

rf_3perc_2 = randomForest(x = genotypes_3perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000  )
save(rf_3perc_2,file="sex_data/rf_3perc_2_sex.Rdata")

rf_3perc_3 = randomForest(x = genotypes_3perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=25000  )
save(rf_3perc_3,file="sex_data/rf_3perc_3_sex.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_3perc_1_err.rate <- rf_3perc_1$err.rate[25000]
rf_3perc_2_err.rate <- rf_3perc_2$err.rate[25000]
rf_3perc_3_err.rate <- rf_3perc_3$err.rate[25000]

rm(rf_3perc_1,rf_3perc_2,rf_3perc_3) # remove the objects to save memory in R

##### Best 4%

names_best_4perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.96))]
names_best_4perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.96))]
names_best_4perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.96))]
names_best_4perc_unique<-unique(c(names_best_4perc_1,names_best_4perc_2,names_best_4perc_3))

# Extract genotypes
genotypes_4perc<-df.train[,colnames(df.train) %in% names_best_4perc_unique]
rownames(genotypes_4perc)

# Now conduct RF on this subset
rf_4perc_1 = randomForest(x = genotypes_4perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000  )
save(rf_4perc_1,file="sex_data/rf_4perc_1_sex.Rdata")

rf_4perc_2 = randomForest(x = genotypes_4perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000  )
save(rf_4perc_2,file="sex_data/rf_4perc_2_sex.Rdata")

rf_4perc_3 = randomForest(x = genotypes_4perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=25000  )
save(rf_4perc_3,file="sex_data/rf_4perc_3_sex.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_4perc_1_err.rate <- rf_4perc_1$err.rate[25000]
rf_4perc_2_err.rate <- rf_4perc_2$err.rate[25000]
rf_4perc_3_err.rate <- rf_4perc_3$err.rate[25000]

rm(rf_4perc_1,rf_4perc_2,rf_4perc_3) # remove the objects to save memory in R

##### Best 5%

names_best_5perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.95))]
names_best_5perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.95))]
names_best_5perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.95))]
names_best_5perc_unique<-unique(c(names_best_5perc_1,names_best_5perc_2,names_best_5perc_3))

# Extract genotypes
genotypes_5perc<-df.train[,colnames(df.train) %in% names_best_5perc_unique]

# Now conduct RF on this subset
rf_5perc_1 = randomForest(x = genotypes_5perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000  )
save(rf_5perc_1,file="sex_data/rf_5perc_1_sex.Rdata")

rf_5perc_2 = randomForest(x = genotypes_5perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000  )
save(rf_5perc_2,file="sex_data/rf_5perc_2_sex.Rdata")

rf_5perc_3 = randomForest(x = genotypes_5perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=25000  )
save(rf_5perc_3,file="sex_data/rf_5perc_3_sex.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_5perc_1_err.rate <- rf_5perc_1$err.rate[25000]
rf_5perc_2_err.rate <- rf_5perc_2$err.rate[25000]
rf_5perc_3_err.rate <- rf_5perc_3$err.rate[25000]

rm(rf_5perc_1,rf_5perc_2,rf_5perc_3) # remove the objects to save memory in R

##### Best 10%

names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))

# Extract genotypes
genotypes_10perc<-df.train[,colnames(df.train) %in% names_best_10perc_unique]

# Now conduct RF on this subset
rf_10perc_1 = randomForest(x = genotypes_10perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000  )
save(rf_10perc_1,file="sex_data/rf_10perc_1_sex.Rdata")

rf_10perc_2 = randomForest(x = genotypes_10perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000  )
save(rf_10perc_2,file="sex_data/rf_10perc_2_sex.Rdata")

rf_10perc_3 = randomForest(x = genotypes_10perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=25000  )
save(rf_10perc_3,file="sex_data/rf_10perc_3_sex.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_10perc_1_err.rate <- rf_10perc_1$err.rate[25000]
rf_10perc_2_err.rate <- rf_10perc_2$err.rate[25000]
rf_10perc_3_err.rate <- rf_10perc_3$err.rate[25000]

rm(rf_10perc_1,rf_10perc_2,rf_10perc_3) # remove the objects to save memory in R

##### Best 20%

names_best_20perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.80))]
names_best_20perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.80))]
names_best_20perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.80))]
names_best_20perc_unique<-unique(c(names_best_20perc_1,names_best_20perc_2,names_best_20perc_3))

# Extract genotypes
genotypes_20perc<-df.train[,colnames(df.train) %in% names_best_20perc_unique]

# Now conduct RF on this subset
rf_20perc_1 = randomForest(x = genotypes_20perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000  )
save(rf_20perc_1,file="sex_data/rf_20perc_1_sex.Rdata")

rf_20perc_2 = randomForest(x = genotypes_20perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000  )
save(rf_20perc_2,file="sex_data/rf_20perc_2_sex.Rdata")

rf_20perc_3 = randomForest(x = genotypes_20perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=25000  )
save(rf_20perc_3,file="sex_data/rf_20perc_3_Sex.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_20perc_1_err.rate <- rf_20perc_1$err.rate[25000]
rf_20perc_2_err.rate <- rf_20perc_2$err.rate[25000]
rf_20perc_3_err.rate <- rf_20perc_3$err.rate[25000]

rm(rf_20perc_1,rf_20perc_2,rf_20perc_3) # remove the objects to save memory in R

##### Best 30%

names_best_30perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.70))]
names_best_30perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.70))]
names_best_30perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.70))]
names_best_30perc_unique<-unique(c(names_best_30perc_1,names_best_30perc_2,names_best_30perc_3))

# Extract genotypes
genotypes_30perc<-df.train[,colnames(df.train) %in% names_best_30perc_unique]

# Now conduct RF on this subset
rf_30perc_1 = randomForest(x = genotypes_30perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000  )
save(rf_30perc_1,file="sex_data/rf_30perc_1_Sex.Rdata")

rf_30perc_2 = randomForest(x = genotypes_30perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000  )
save(rf_30perc_2,file="sex_data/rf_30perc_2_sex.Rdata")

rf_30perc_3 = randomForest(x = genotypes_30perc, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=25000  )
save(rf_30perc_3,file="sex_data/rf_30perc_3_sex.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_30perc_1_err.rate <- rf_30perc_1$err.rate[25000]
rf_30perc_2_err.rate <- rf_30perc_2$err.rate[25000]
rf_30perc_3_err.rate <- rf_30perc_3$err.rate[25000]

rm(rf_30perc_1,rf_30perc_2,rf_30perc_3) # remove the objects to save memory in R


# Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach

All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
                              cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
                              cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
                              cbind(rf_4perc_1_err.rate,rf_4perc_2_err.rate,rf_4perc_3_err.rate),
                              cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate),
                              cbind(rf_10perc_1_err.rate,rf_10perc_2_err.rate,rf_10perc_3_err.rate),
                              cbind(rf_20perc_1_err.rate,rf_20perc_2_err.rate,rf_20perc_3_err.rate),
                              cbind(rf_30perc_1_err.rate,rf_30perc_2_err.rate,rf_30perc_3_err.rate))

# Plot error rates for the various subsets
All_initial_err.rate<-data.frame(All_initial_err.rate)
All_initial_err.rate$Number_loci<-c(1942,length(names_best_2perc_unique),length(names_best_3perc_unique),length(names_best_4perc_unique),length(names_best_5perc_unique),length(names_best_10perc_unique),length(names_best_20perc_unique),length(names_best_30perc_unique))
rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_err.rate,file="sex_data/All_initial_err_rate_classification_sex.csv")


# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5)

# Based on this table and plot, the best 2% of loci have the lowest error rate
# I'll run backward purging RF with the best 2% loci


#################### Backward purging approach###########################

names_purging <- names_best_2perc_unique

genotypes_purging<-df.train[,colnames(df.train) %in% names_purging]

############################################################################################################
############################ We continue with 25000 trees hereon ###########################################


rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000  )
save(rf_purging_1,file="sex_data/rf_purging_1_sex.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000  )
save(rf_purging_2,file="sex_data/rf_purging_2_sex.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=25000  )
save(rf_purging_3,file="sex_data/rf_purging_3_sex.Rdata")

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

  rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000  )
  rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000  )
  rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(df.train$sex), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=25000  )

  error_rate_best[length(names_purging)-i,] <- c(rf_purging_1$err.rate[25000],rf_purging_2$err.rate[25000],rf_purging_3$err.rate[25000])
}


error_rate_best$Average<-apply(error_rate_best,1,mean)
write.csv(error_rate_best, file="sex_data/Backward_purging_OOB-ER_classification_2_perc_sex.csv") # Save the error rates

# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus
plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1,cex.axis=1,pch=16)

#Print using R graphics
dev.print(pdf, "error_rate_loci_25k_purge_sex.pdf")
dev.off()




# Which group of loci yields the lowest error rate?
which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) #115, 116, 123, 124, 132 loci have the lowest OOB-ER

error_rate_best$Average[]

# Export the names of the predictor loci
write.csv(names_all_iterations[[4]],file="sex_data/Predictor_loci_classification_sex_4.csv")
write.csv(names_all_iterations[[5]],file="sex_data/Predictor_loci_classification_sex_5.csv")
write.csv(names_all_iterations[[6]],file="sex_data/Predictor_loci_classification_sex_6.csv")
write.csv(names_all_iterations[[7]],file="sex_data/Predictor_loci_classification_sex_7.csv")
write.csv(names_all_iterations[[8]],file="sex_data/Predictor_loci_classification_sex_8.csv")
write.csv(names_all_iterations[[9]],file="sex_data/Predictor_loci_classification_sex_9.csv")
write.csv(names_all_iterations[[10]],file="sex_data/Predictor_loci_classification_sex_10.csv")
write.csv(names_all_iterations[[12]],file="sex_data/Predictor_loci_classification_sex_12.csv")
write.csv(names_all_iterations[[13]],file="sex_data/Predictor_loci_classification_sex_13.csv")
write.csv(names_all_iterations[[14]],file="sex_data/Predictor_loci_classification_sex_14.csv")



#14 loci sex SNPs with imputed values train

#14 løoci train with sex
df.train_14 <- cbind(df.train[, 1:4], df.train[, names_all_iterations[[14]]])
sex_data_14 <- cbind(sex_data_imputed[, 1:4], sex_data_imputed[, names_all_iterations[[14]]])

write.csv(sex_data_14,file="sex_data/14_loci_full_sample.csv")
write.csv(df.train_14,file="sex_data/14_loci_train.csv")

#14 loci sex SNPs with imputed values train
df.test_14 <- cbind(df.test[, 1:4], df.test[, names_all_iterations[[14]]])
sex_data_14 <- cbind(sex_data_imputed[, 1:4], sex_data_imputed[, names_all_iterations[[14]]])

write.csv(sex_data_14,file="sex_data/14_loci_full_sample.csv")
write.csv(df.test_14,file="sex_data/14_loci_test.csv")


# 14 loci sex with NA
sex_data_14 <- cbind(sex_data[, 1:4], sex_data[, names_all_iterations[[14]]])
write.csv(sex_data_14,file="sex_data/14_loci_full_sample_NA.csv")
