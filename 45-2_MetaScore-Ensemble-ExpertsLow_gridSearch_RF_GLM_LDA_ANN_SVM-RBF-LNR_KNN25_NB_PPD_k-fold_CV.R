# Ex of execution: Rscript 45-2_ExpertCommitteeLowRank_gridSearch_RF_GLM_LDA_ANN_SVM-RBF-LNR_KNN25_NB_PPD_k-fold_CV.R /gpfs/home/yuj114/group/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HADDOCK-Val+interfacialPSSM/Undersampled/ /gpfs/home/yuj114/group/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HADDOCK-Val+interfacialPSSM/ /gpfs/home/yuj114/group/YJ_DOCK_PP/RFBased_ScoringFunction/RF_Result/HADDOCK-Val+interfacialPSSM/10-foldCV/ 10 RF 500 2
# in Hammer server

library(class)
#library(gplots)
#library(affy)
#library(samr)
library(multtest)
library(e1071)     #naiveBayes, SVM classification method
#library(subselect) #Genetic Algorithm
#library(survival)
#library(genalg)
#library(limma)
library(randomForest)
library(boot)
library(rpart)
library(MASS)
library(nnet)
library(reshape)
require(plyr)

#-----------------

#TestData_DIR = '/gpfs/home/yuj114/group/YJ_PRI_2nd/processedData/2014Before/PR_PSSMn3merFreq_40_25/'
#TrainData_DIR = '/gpfs/home/yuj114/group/YJ_PRI_2nd/processedData/2014Before/PR_PSSMn3merFreq_40_25/Undersampled/Selected/'
#outputDIR = '/gpfs/home/yuj114/group/YJ_PRI_2nd/RF_Result/2014Before/5-foldCV/PR_PSSMn3merFreq_40_25/'

Args <- commandArgs(T)
TrainData_DIR = Args[1]
TestData_DIR = Args[2]
outputDIR = Args[3] #paste('/gpfs/home/yuj114/group/YJ_PRI_2nd/', classifier, '_Result/2014Before/', kfold, '-foldCV/', strsplit(TestData_DIR, '/')[[1]][9], '/', sep = '')
kfold = as.integer(Args[4])
classifier = Args[5]
trees = as.integer(Args[6])
rf_vars = as.integer(Args[7])
rm(Args)


if (! file.exists(outputDIR)){
	dir.create(outputDIR, recursive = T)
}

#-----------------

file_list_train <- list.files(TrainData_DIR, pattern = '.csv')
file_list_test <- list.files(TestData_DIR, pattern = '.csv')

set.seed(111)
test_idx <- sample(rep(1:kfold, length(file_list_train) / (kfold - 1), length.out = length(file_list_train)))
print (test_idx)
#--- read training data into memory

PSSMwindows_Train <- NULL

for (k in 1:length(file_list_train)){

	filetrain = paste(TrainData_DIR, file_list_train[k], sep = '')

	if(!file.exists(filetrain)){
		print (paste(filetrain, 'does not exist'))
		stop("Training file does not exist");
	}

	Train <- read.csv(filetrain, comment.char = '#')

	PSSMwindows_Train = c(PSSMwindows_Train, list(Train)) #-PSSMwindows_Train is a list of Data frames
}

rm(Train)
summary(PSSMwindows_Train)

for (i in 1:kfold) {

	print (paste("fold = ", i, "begins"))

	test_list <- c()
	PSSMwindows_Test <- NULL

	for (j in which(test_idx == i)) {
		test_list <- rbind(test_list, paste(i, file_list_test[j], sep = ','))

		filetest = paste(TestData_DIR, file_list_test[j], sep = '')
		print (j)
		print (filetest)
		if(!file.exists(filetest)){
			print (paste(filetest, 'does not exist'))
			stop("Testing file does not exist");
		}

		Test <- read.csv(filetest, comment.char = '#')

		if (nrow(Test) == 0) {
			print ("None data. Pass!")
			next
		}
		PSSMwindows_Test = c(PSSMwindows_Test, list(Test)) #-PSSMwindows_Test is a list of Data frames
	}

	write.csv(test_list, file = paste(outputDIR, 'TestSet_', i, '.txt' , sep = ''))

	rm(Test)
	summary(PSSMwindows_Test)

#	Complex		1
#	Class		2
#	iRMSD		3
#	HDScore		4
#	Evdw		5
#	Eelec		6
#	Edesolv		7
#	BSA		8
#	IC_avg		9
#	IC_max		10
#	IC_min		11
#	P1_IC_avg	12
#	P1_IC_max	13
#	P1_IC_min	14
#	P2_IC_avg	15
#	P2_IC_max	16
#	P2_IC_min	17
#	Dist10		18
#	Dist9		19
#	Dist8		20
#	Dist7		21
#	Dist6		22
#	Dist5		23
#	Dist4		24
#	Dist3		25
#	Dist2		26
#	Dist1		27
#	NumLink		28
#	NumIntRes	29
#	LinkDensity	30
#	normHDScore	31
#	normEvdw	32
#	normEelec	33
#	normEdesolv	34
#	normBSA		35
#	IP_sum		36
#	IP_avg		37
#	IPwPSSM_sum	38
#	IPwPSSM_avg	39
#	distIPwPSSM_sum	40
#	distIPwPSSM_avg	41
#	CX_mean		42
#	CX_SD		43
#	HydPhoPair_mean	44
#	rASA_mean	45
#	ssH_ratio	46
#	ssG_ratio	47
#	ssI_ratio	48
#	ssE_ratio	49
#	ssB_ratio	50
#	ssT_ratio	51
#	ssC_ratio	52

	test <- rbind.fill(PSSMwindows_Test)
	atomResNumPair <- test$Complex
	rm(PSSMwindows_Test)
#	test <- as.data.frame(sapply(test, as.numeric)) #result in '?' converted to 1

	test <- test[, c(-1, -3, -4, -5, -8, -9, -10, -11, -12, -13, -14, -15)]                                #All-HighRanks

#	test <- test[, -c(1, 3, 4, 31)]                          #All-HDScore-NormHDScore
#	test <- test[, -c(1, 3, 4, 9:17, 31)]                    #All-HDScore-NormHDScore-PSSMs

#	test <- test[, -c(1, 3, 4:8, 31:35)]                     #All-HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms
#	test <- test[, -c(1, 3, 9:11)]                           #All-PSSMIC-All
#	test <- test[, -c(1, 3, 9:11, 18:30, 36:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+PSSMIC-P1P2
#	test <- test[, -c(1, 3, 9:27, 30, 36:52)]                #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+NumLink+IntfRes
#	test <- test[, -c(1, 3, 9:30, 36:45)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+SSCategory
#	test <- test[, -c(1, 3, 9:30, 36:44, 46:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+rASA
#	test <- test[, -c(1, 3, 9:30, 36:43, 45:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+HydPho
#	test <- test[, -c(1, 3, 9:30, 36:41, 44:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+CX
#	test <- test[, -c(1, 3, 9:29, 36:41)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+LinkDensity
#	test <- test[, -c(1, 3, 9:17, 28:30, 36:41)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+Dist-Top10
#	test <- test[, -c(1, 3, 9:30, 36:39)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+distIPwPSSM_SA
#	test <- test[, -c(1, 3, 9:30, 36:37, 40:41)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+IPwPSSM_SA
#	test <- test[, -c(1, 3, 9:30, 38:41)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+IP_SA
#	test <- test[, -c(1, 3, 12:30, 36:41)]                   #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+InterfacialPSSMIC-Para3-All
#	test <- test[, -c(1, 3, 9:30, 36:41)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms
#	test <- test[, -c(1, 3, 4:8, 9:11, 18:30, 36:41)]        #NormHDScore-NormHDTerms+InterfacialPSSMIC-Para3-P1P2
#	test <- test[, -c(1, 3, 4:8, 9:29, 36:41)]               #NormHDScore-NormHDTerms+LinkDensity
#	test <- test[, -c(1, 3, 4:17, 28:30, 36:41)]             #NormHDScore-NormHDTerms+Dist-Top10
#	test <- test[, -c(1, 3, 4:8, 12:30, 36:41)]              #NormHDScore-NormHDTerms+InterfacialPSSMIC-Para3-All
#	test <- test[, -c(1, 3, 4:30, 36:41)]                    #NormHDScore-NormHDTerms
#	test <- test[, -c(1, 3, 9:11, 18:30, 31:35, 36:41)]      #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-P1P2
#	test <- test[, -c(1, 3, 9:29, 31:35, 36:41)]             #HDScore_HADDOCK-Terms+LinkDensity
#	test <- test[, -c(1, 3, 9:17, 28:30, 31:35, 36:41)]      #HDScore_HADDOCK-Terms+Dist-Top10
#	test <- test[, -c(1, 3, 12:30, 31:35, 36:41)]            #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-All
#	test <- test[, -c(1, 3, 9:30, 31:35, 36:41)]             #HDScore_HADDOCK-Terms
#	test <- test[, -c(1, 3, 28, 30, 31:35, 36:41)]           #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10_10.0AforNeg 
#	test <- test[, -c(1, 3, 4, 28, 30, 31:35, 36:41)]        #HADDOCK-Terms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10_10.0AforNeg
#	test <- test[, -c(1, 3, 12:17, 28, 30, 31:35, 36:41)]    #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-All+Dist-Top10_10.0AforNeg
#	test <- test[, -c(1, 3, 4, 12:17, 28, 30, 31:35, 36:41)] #HADDOCK-Terms+InterfacialPSSMIC-Para3-All+Dist-Top10_10.0AforNeg

	#--train RF

	train <- rbind.fill(PSSMwindows_Train[-which(test_idx == i)]) #
#	final_train <- as.data.frame(sapply(train, as.numeric)) # result in NA from '?'
#	final_train <- final_train[, c(-1, -3)]

	final_train <- train[, c(-1, -3, -4, -5, -8, -9, -10, -11, -12, -13, -14, -15)]                                #All-HighRanks

#	final_train <- train[, -c(1, 3, 4, 31)]                          #All-HDScore-NormHDScore

#	final_train <- train[, -c(1, 3, 4:8, 31:35)]                     #All-HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms
#	final_train <- train[, -c(1, 3, 9:11)]                           #All-PSSMIC-All
#	final_train <- train[, -c(1, 3, 9:11, 18:30, 36:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+PSSMIC-P1P2
#	final_train <- train[, -c(1, 3, 9:27, 30, 36:52)]                #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+NumLink+IntfRes
#	final_train <- train[, -c(1, 3, 9:30, 36:45)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+SSCategory
#	final_train <- train[, -c(1, 3, 9:30, 36:44, 46:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+rASA
#	final_train <- train[, -c(1, 3, 9:30, 36:43, 45:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+HydPho
#	final_train <- train[, -c(1, 3, 9:30, 36:41, 44:52)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+CX
#	final_train <- train[, -c(1, 3, 9:29, 36:41)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+LinkDensity
#	final_train <- train[, -c(1, 3, 9:17, 28:30, 36:41)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+Dist-Top10
#	final_train <- train[, -c(1, 3, 9:30, 36:39)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+distIPwPSSM_SA
#	final_train <- train[, -c(1, 3, 9:30, 36:37, 40:41)]             #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+IPwPSSM_SA
#	final_train <- train[, -c(1, 3, 9:30, 38:41)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+IP_SA
#	final_train <- train[, -c(1, 3, 12:30, 36:41)]                   #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms+InterfacialPSSMIC-Para3-All
#	final_train <- train[, -c(1, 3, 9:30, 36:41)]                    #HDScore_HADDOCK-Terms+NormHDScore-NormHDTerms
#	final_train <- train[, -c(1, 3, 4:8, 9:11, 18:30, 36:41)]        #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-P1P2
#	final_train <- train[, -c(1, 3, 4:8, 9:29, 36:41)]               #HDScore_HADDOCK-Terms+LinkDensity
#	final_train <- train[, -c(1, 3, 4:17, 28:30, 36:41)]             #NormHDScore-NormHDTerms+Dist-Top10
#	final_train <- train[, -c(1, 3, 4:8, 12:30, 36:41)]              #NormHDScore-NormHDTerms+InterfacialPSSMIC-Para3-All
#	final_train <- train[, -c(1, 3, 4:30, 36:41)]                    #NormHDScore-NormHDTerms
#	final_train <- train[, -c(1, 3, 9:11, 18:30, 31:35, 36:41)]      #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-P1P2
#	final_train <- train[, -c(1, 3, 9:29, 31:35, 36:41)]             #HDScore_HADDOCK-Terms+LinkDensity
#	final_train <- train[, -c(1, 3, 9:17, 28:30, 31:35, 36:41)]      #HDScore_HADDOCK-Terms+Dist-Top10
#	final_train <- train[, -c(1, 3, 12:30, 31:35, 36:41)]            #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-All
#	final_train <- train[, -c(1, 3, 9:30, 31:35, 36:41)]             #HDScore_HADDOCK-Terms
#	final_train <- train[, -c(1, 3, 28, 30, 31:35, 36:41)]           #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10_10.0AforNeg 
#	final_train <- train[, -c(1, 3, 4, 28, 30, 31:35, 36:41)]        #HADDOCK-Terms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10_10.0AforNeg
#	final_train <- train[, -c(1, 3, 12:17, 28, 30, 31:35, 36:41)]    #HDScore_HADDOCK-Terms+InterfacialPSSMIC-Para3-All+Dist-Top10_10.0AforNeg
#	final_train <- train[, -c(1, 3, 4, 12:17, 28, 30, 31:35, 36:41)] #HADDOCK-Terms+InterfacialPSSMIC-Para3-All+Dist-Top10_10.0AforNeg

	rm(train)
	colnames(final_train)[1] <- "Class"
	colnames(test)[1] <- "Class"

	final_train$Class = as.factor(final_train$Class)

	PRI.clsfier <- NULL

	if (classifier == 'GLM') {
		PRI.clsfier <- glm(Class ~ ., data = final_train, family = binomial(logit))
	} else if (classifier == 'RF') {
		set.seed(111)
		PRI.clsfier <- randomForest(Class ~ ., data = final_train, na.action = na.omit, ntree = trees, mtry = rf_vars)
	} else if (classifier == 'LDA') {
		PRI.clsfier <- lda(Class ~ ., data = final_train)
	} else if (classifier == 'ANN') {
		PRI.clsfier <- nnet(Class ~ ., data = final_train, size = 1)
	} else if (classifier == 'SVM-RBF') {
		PRI.clsfier <- svm(Class ~ ., data = final_train, cost = 1.0, gamma = 0.01, kernel = "radial", probability = TRUE)
	#	PRI.clsfier <- svm(Class ~ ., data = final_train, cost = 1.0, gamma = 0.01, kernel = "linear", probability = TRUE)
	} else if (classifier == 'SVM-LNR') {
	#	PRI.clsfier <- svm(Class ~ ., data = final_train, cost = 1.0, gamma = 0.01, kernel = "radial", probability = TRUE)
		PRI.clsfier <- svm(Class ~ ., data = final_train, cost = 1.0, gamma = 0.01, kernel = "linear", probability = TRUE)
	} else if (classifier == 'KNN') {
		print("KNN will start shortly.")
	} else if (classifier == 'NB') {
	#	PRI.clsfier <- naiveBayes(fmla, data = class_imputed.lum[, 1:(length(colnames(class_imputed.lum)) - 1)])
		PRI.clsfier <- naiveBayes(Class ~ ., data = final_train)
	} else {
		stop("No classification model designated.")
	}

	colnames_train = colnames(final_train)

	#------------------------- test

	#--set the test column names the same as training data
	colnames(test) = colnames_train

	#--
	test$Class <- as.factor(test$Class)
	tclass <- test$Class

	#--

	PRI.pred <- NULL

	if (classifier == 'GLM') {
		rm(final_train)
		temp <- predict(PRI.clsfier, test[, -1], type = 'response')
		PRI.pred <- cbind(1 - temp, temp)
	} else if (classifier == 'RF') {
		rm(final_train)
		PRI.pred <- predict(PRI.clsfier, test[, -1], type = 'prob')
	} else if (classifier == 'LDA') {
		rm(final_train)
		PRI.pred <- predict(PRI.clsfier, test[, -1], type = 'prob')$posterior
	} else if (classifier == 'ANN') {
		rm(final_train)
		temp <- predict(PRI.clsfier, test[, -1], type = 'raw')
		PRI.pred <- cbind(1 - temp, temp)
	} else if (classifier == 'SVM-RBF') {
		rm(final_train)
		temp <- predict(PRI.clsfier, test[, -1], type = 'prob', probability = TRUE, decision.values = TRUE)
		PRI.pred <- cbind(attr(temp, "probabilities")[, 2], attr(temp, "probabilities")[, 1])
	} else if (classifier == 'SVM-LNR') {
		rm(final_train)
		temp <- predict(PRI.clsfier, test[, -1], type = 'prob', probability = TRUE, decision.values = TRUE)
		PRI.pred <- cbind(attr(temp, "probabilities")[, 2], attr(temp, "probabilities")[, 1])
	} else if (classifier == 'KNN') {
		knn_isolet <- class::knn(final_train[, -1], test[, -1], final_train[, 1], prob = TRUE, k = 25)
		prob <- attr(knn_isolet, "prob")
		#prob <- 2 * ifelse(knn_isolet == "-1", 1 - prob, prob) - 1
		prob <- ifelse(knn_isolet == "0", 1 - prob, prob)
		PRI.pred <- cbind(1 - prob, prob)
		rm(knn_isolet)
		rm(final_train)
	} else if (classifier == 'NB') {
		rm(final_train)
		PRI.pred <- predict(PRI.clsfier, test[, -1], type = 'raw')
	} else {
		stop("No classification model designated.")
	}

	#PRI.pred <- predict(PRI.clsfier, test[, -1], type = 'prob')
	#table(observed = test$class, predicted = BM3.pred
	rm(test)

	#--
	#write prediction result file
	write.csv(cbind(as.character(atomResNumPair), PRI.pred, as.numeric(as.character(tclass))), file = paste(outputDIR, 'PredictedScores_TestSet_', i, '.csv' , sep = ''))
	print(sprintf("output file: %s", paste('PredictedScores_TestSet_', i, '.csv' , sep = '')))
	rm(PRI.pred)
	rm(PRI.clsfier)
}
