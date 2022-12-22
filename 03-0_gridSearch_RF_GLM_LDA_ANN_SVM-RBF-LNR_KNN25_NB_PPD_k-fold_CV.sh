#!/bin/bash
# Yong Jung
# May. 3rd 2017
# in Hammer server
#./03-0_gridSearch_RF_GLM_LDA_ANN_SVM-RBF-LNR_KNN25_NB_PPD_k-fold_CV.sh HDScore_HADDOCK-Terms 10.0 2 10 5 RF 

#HDScore-HDTerms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormHDScore-NormHDTerms_14.0AforNeg

feat=$1 #'HADDOCK-Val+interfacialPSSM'
negcut=$2 #'10.0'
underSmplRto=$3
fol=$4 #'10'
feanum=$5 #'5'
classifier=$6 #'RF'

module load R

#for trees in {250,500,750,1000}; do
#for trees in {10..200..10}; do
for trees in {10..500..10}; do
#	for mtry in {1,2,$(echo "scale=2;sqrt($feanum)"| bc -l),3,4,5,6,7}; do
#	for mtry in {1,2,3,4,5,6,7}; do
	for mtry in {1,4,7,10,13,16,19,22,25,28}; do
		echo ""
		echo "ntrees = $trees, mtry = $mtry"
#		echo "Rscript 03-0_gridSearch_RF_GLM_LDA_ANN_SVM-RBF-LNR_KNN25_NB_PPD_k-fold_CV.R /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_$feat/Undersampled_Neg${underSmplRto}x/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_$feat/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/RF_Result/$feat/$trees\_$mtry\_$fol-foldCV/ $fol $classifier $trees $mtry"
		echo "Rscript 03-0_gridSearch_RF_GLM_LDA_ANN_SVM-RBF-LNR_KNN25_NB_PPD_k-fold_CV.R /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_${negcut}AforNeg/Undersampled_Neg${underSmplRto}x/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_${negcut}AforNeg/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/RF_Result/$feat\_${negcut}AforNeg_Undersampled_Neg${underSmplRto}x/$trees\_$mtry\_$fol-foldCV/ $fol $classifier $trees $mtry"
		
#		Rscript 03-0_gridSearch_RF_GLM_LDA_ANN_SVM-RBF-LNR_KNN25_NB_PPD_k-fold_CV.R /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_$feat/Undersampled_Neg${underSmplRto}x/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_$feat/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/RF_Result/$feat/$trees\_$mtry\_$fol-foldCV/ $fol $classifier $trees $mtry
		Rscript /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/03-0_gridSearch_RF_GLM_LDA_ANN_SVM-RBF-LNR_KNN25_NB_PPD_k-fold_CV.R /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios\_${negcut}AforNeg/Undersampled_Neg${underSmplRto}x/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios\_${negcut}AforNeg/ /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/RF_Result/$feat\_${negcut}AforNeg_Undersampled_Neg${underSmplRto}x/$trees\_$mtry\_$fol-foldCV/ $fol $classifier $trees $mtry

done
done