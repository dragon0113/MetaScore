# Ex of execution: python 44_ExpertCommittee_DataConstruction_BM4+BM5+BM4UnderSmpl.py
# in Hammer server
#from difflib import SequenceMatcher
import glob
import sys
import os
import math
import time
#import random
#import time 

# Decoy path: "/home/yjung/GraphKernel_scoringFun/data/BM4_dimers/bound/pdbFLs_refined/1AK4.pdb"
# PSSM path: "/home/yjung/GraphKernel_scoringFun/data/BM4_dimers/bound/PSSM"
# Interfacial Residue path: "/home/yjung/GraphKernel_scoringFun/data/BM4_dimers/bound/pdbFLs_refined/Atom_distance_15/CA-CA_Dist_8"

#/gpfs/home/yuj114/group/YJ_DOCK_PP/i-rmsd/BM5
#HAD_Score_list = glob.glob('/gpfs/home/yuj114/group/YJ_DOCK_PP/haddockScore/BM5/*.haddockScore')
#IRMSD_list = glob.glob('/gpfs/home/yuj114/group/YJ_DOCK_PP/i-rmsd/BM5/*.irmsd')
#PSSM_list = glob.glob('/gpfs/home/yuj114/group/YJ_DOCK_PP/BM5_dimers/unbound/PSSM/*.protein1.ResNumPSSM')

#HD_list = glob.glob('/gpfs/home/yuj114/group/YJ_DOCK_PP/haddockScore/BM5/it0/*.haddockScore')
BM4_Data_list = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/*.csv')

BM4_UnderSmpl_Data_list = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/Undersampled_Neg1x/*.csv')

BM5_Data_list = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/Test-water_processed_HDScore-HDTerms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/*.csv')


PDBs = []

for hd in BM4_Data_list + BM5_Data_list:
	PDBs.append(hd.split('/')[-1].split('.')[0].split('_')[0])

PDBs.sort()

PDBs_Final = list(set(PDBs))
PDBs_Final.sort()

outputBM4 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_ExpertCommittee/'
outputBM5 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/Test-water_processed_ExpertCommittee/'
outputBM4_UnderSmpl = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_ExpertCommittee/Undersampled_Neg1x/'

os.system('rm -rf ' + outputBM4)
os.system('mkdir ' + outputBM4)
os.system('mkdir ' + outputBM4_UnderSmpl)

os.system('rm -rf ' + outputBM5)
os.system('mkdir ' + outputBM5)

for hd in BM4_Data_list:
	fw1 = open(outputBM4 + hd.split('/')[-1], 'w')
	dic_ExpCom = {}
	print >>fw1, "Complex,Class,iRMSD,HADDOCK, NormHADDOCK,PYDOCK,NormPYDOCK,DFIRE2,NormDFIRE2,PISA,NormPISA,DFIRE,NormDFIRE,MJ3H,NormMJ3H,SWARMDOCK,NormSWARMDOCK,TOBI,NormTOBI,SIPPER,NormSIPPER,ISCORE,NormISCORE"
	
	for sc in ['HD', 'PYDOCK' , 'DFIRE2', 'PISA', 'DFIRE', 'MJ3H', 'SWARMDOCK', 'TOBI', 'SIPPER', 'ISCORE']:
		dic_Score = {}
		
		inputBM4 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_' + sc + 'Score-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/'
		
		fr = open(inputBM4 + hd.split('/')[-1], 'r')
		lines = fr.read().split('\n')
		fr.close()
		
		for i in range(1, len(lines) - 1):
			
			if lines[i].split(',')[0].find('_') == -1:
				continue
			
			if sc == 'HD':
				dic_ExpCom[lines[i].split(',')[0].replace('.pdb', '')] = lines[i].split(',')[0] + ',' + lines[i].split(',')[1] + ',' + lines[i].split(',')[2] + ',' + lines[i].split(',')[3] + ',' + lines[i].split(',')[30]
			
			else:
				dic_ExpCom[lines[i].split(',')[0].replace('.pdb', '')] += ',' + lines[i].split(',')[3] + ',' + lines[i].split(',')[30]
	
	for d in dic_ExpCom.keys():
		if dic_ExpCom[d].count(',') < 15:
			continue
		
		print >>fw1, dic_ExpCom[d]
	
	fw1.close()


for hd in BM5_Data_list:
	fw1 = open(outputBM5 + hd.split('/')[-1], 'w')
	dic_ExpCom = {}
	print >>fw1, "Complex,Class,iRMSD,HADDOCK, NormHADDOCK,PYDOCK,NormPYDOCK,DFIRE2,NormDFIRE2,PISA,NormPISA,DFIRE,NormDFIRE,MJ3H,NormMJ3H,SWARMDOCK,NormSWARMDOCK,TOBI,NormTOBI,SIPPER,NormSIPPER,ISCORE,NormISCORE"
	
	for sc in ['HD', 'PYDOCK' , 'DFIRE2', 'PISA', 'DFIRE', 'MJ3H', 'SWARMDOCK', 'TOBI', 'SIPPER', 'ISCORE']:
		dic_Score = {}
		
		inputBM5 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/Test-water_processed_' + sc + 'Score-HDTerms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/'
		
		fr = open(inputBM5 + hd.split('/')[-1], 'r')
		lines = fr.read().split('\n')
		fr.close()
		
		for i in range(1, len(lines) - 1):
			
			if lines[i].split(',')[0].find('_') == -1:
				continue
			
			if sc == 'HD':
				dic_ExpCom[lines[i].split(',')[0].replace('.pdb', '')] = lines[i].split(',')[0] + ',' + lines[i].split(',')[1] + ',' + lines[i].split(',')[2] + ',' + lines[i].split(',')[3] + ',' + lines[i].split(',')[30]
			
			else:
				dic_ExpCom[lines[i].split(',')[0].replace('.pdb', '')] += ',' + lines[i].split(',')[3] + ',' + lines[i].split(',')[30]
	
	for d in dic_ExpCom.keys():
		if dic_ExpCom[d].count(',') < 15:
			continue
		
		print >>fw1, dic_ExpCom[d]
	
	fw1.close()


for hd in BM4_UnderSmpl_Data_list:
	fw1 = open(outputBM4_UnderSmpl + hd.split('/')[-1], 'w')
	dic_ExpCom = {}
	print >>fw1, "Complex,Class,iRMSD,HADDOCK, NormHADDOCK,PYDOCK,NormPYDOCK,DFIRE2,NormDFIRE2,PISA,NormPISA,DFIRE,NormDFIRE,MJ3H,NormMJ3H,SWARMDOCK,NormSWARMDOCK,TOBI,NormTOBI,SIPPER,NormSIPPER,ISCORE,NormISCORE"
	
	for sc in ['HD', 'PYDOCK' , 'DFIRE2', 'PISA', 'DFIRE', 'MJ3H', 'SWARMDOCK', 'TOBI', 'SIPPER', 'ISCORE']:
		dic_Score = {}
		
		inputBM4_UnderSmpl = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_' + sc + 'Score-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/Undersampled_Neg1x/'
		
		fr = open(inputBM4_UnderSmpl + hd.split('/')[-1], 'r')
		lines = fr.read().split('\n')
		fr.close()
		
		for i in range(1, len(lines) - 1):
			
			if lines[i].split(',')[0].find('_') == -1:
				continue
			
			if sc == 'HD':
				dic_ExpCom[lines[i].split(',')[0].replace('.pdb', '')] = lines[i].split(',')[0] + ',' + lines[i].split(',')[1] + ',' + lines[i].split(',')[2] + ',' + lines[i].split(',')[3] + ',' + lines[i].split(',')[30]
			
			else:
				dic_ExpCom[lines[i].split(',')[0].replace('.pdb', '')] += ',' + lines[i].split(',')[3] + ',' + lines[i].split(',')[30]
	
	for d in dic_ExpCom.keys():
		if dic_ExpCom[d].count(',') < 15:
			continue
		
		print >>fw1, dic_ExpCom[d]
	
	fw1.close()


