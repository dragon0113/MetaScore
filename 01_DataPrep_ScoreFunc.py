#PSSM_IC A & B, HADDOCK score added as features

# Ex of execution: python 01_DataPrep_ScoreFunc_HDScore-HDTerms+InterfPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LD+NormHDScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA.py
# in Hammer server

#from difflib import SequenceMatcher
import glob
import sys
import os
import math
import numpy as np

#from math import sqrt
#import random
#import time 

# Decoy path: "/home/yjung/GraphKernel_scoringFun/data/BM4_dimers/bound/pdbFLs_refined/1AK4.pdb"
# PSSM path: "/home/yjung/GraphKernel_scoringFun/data/BM4_dimers/bound/PSSM"
# Interfacial Residue path: "/home/yjung/GraphKernel_scoringFun/data/BM4_dimers/bound/pdbFLs_refined/Atom_distance_15/CA-CA_Dist_8"

def sigmoid(x):
	return 1 / (1 + math.exp(-x))

# Amino Acid name encoding from three letters to one letter
ama = {}
ama['ALA'] = 'A'
ama['ARG'] = 'R'
ama['ASN'] = 'N'
ama['ASP'] = 'D'
ama['CYS'] = 'C'
ama['GLU'] = 'E'
ama['GLN'] = 'Q'
ama['GLY'] = 'G'
ama['HIS'] = 'H'
ama['ILE'] = 'I'
ama['LEU'] = 'L'
ama['LYS'] = 'K'
ama['MET'] = 'M'
ama['PHE'] = 'F'
ama['PRO'] = 'P'
ama['SER'] = 'S'
ama['THR'] = 'T'
ama['TRP'] = 'W'
ama['TYR'] = 'Y'
ama['VAL'] = 'V'

# Protein-Protein interaction propensity from InterEvol database
ppIP = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/PP_InteractionPropensity.txt', 'r')
ppIP_lines = ppIP.read().split('\n')
ppIP.close()

ppIP_dic = {}

for i in range(0, len(ppIP_lines) - 1):
	
	ppIP_dic[ppIP_lines[i].split('\t')[0] + '_' + ppIP_lines[i].split('\t')[1]] = float(ppIP_lines[i].split('\t')[3])
	ppIP_dic[ppIP_lines[i].split('\t')[1] + '_' + ppIP_lines[i].split('\t')[0]] = float(ppIP_lines[i].split('\t')[3])

# Min, Max of HADDOCK features


PSSM_list = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/AddResNumPSSM_seqFLs_bound/*.protein1.ResNumPSSM')

PSSM_PDBs = []

for pssm in PSSM_list:
	PSSM_PDBs.append(pssm.split('/')[-1].replace(".protein1.ResNumPSSM", ""))

PSSM_PDBs.sort()

HAD_Score_list = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/haddockScore/water/*.haddockScore')

# HADDOCK score for each decoy

dic_HADDOCK = {}

HD_PDBs = []

dic_HADDOCK_Evdw_List = {}
dic_HADDOCK_Eelec_List = {}
dic_HADDOCK_Edesolv_List = {}
dic_HADDOCK_BSA_List = {}
dic_HADDOCK_Score_List = {}

dic_HADDOCK_Evdw_Mean = {}
dic_HADDOCK_Eelec_Mean = {}
dic_HADDOCK_Edesolv_Mean = {}
dic_HADDOCK_BSA_Mean = {}
dic_HADDOCK_Score_Mean = {}

dic_HADDOCK_Evdw_SD = {}
dic_HADDOCK_Eelec_SD = {}
dic_HADDOCK_Edesolv_SD = {}
dic_HADDOCK_BSA_SD = {}
dic_HADDOCK_Score_SD = {}

for hd in HAD_Score_list:
	
	HD_PDBs.append(hd.split('/')[-1].split('.')[0])
	
	f_HAD = open(hd, 'r')
	HAD = f_HAD.read().split('\n')
	f_HAD.close()
	
	HDScore_list = []
	
	dic_HADDOCK_Score_List[hd.split('/')[-1].split('.')[0]] = []
	
	for i in range(0, len(HAD) - 1):
		dic_HADDOCK[HAD[i].split()[0].split('.')[0]] = float(HAD[i].split()[1])
		
		dic_HADDOCK_Score_List[hd.split('/')[-1].split('.')[0]].append(float(HAD[i].split()[1]))
		
		HDScore_list.append(float(HAD[i].split()[1]))
	
	f_bound = open("/gpfs/group/vuh14/legacy/YJ_DOCK_PP/BM4_dimers/bound/pdbFLs_refined/" + hd.split('/')[-1].split('.')[0] + '.pdb', 'r')
	bound = f_bound.read().split('\n')
	f_bound.close()
	
	vdw = 0
	elec = 0
	AIR = 0
	desolv = 0
	
	for i in range(0, len(bound) - 1):
		if bound[i].find('ATOM') == 0:
			break
		
		if bound[i].find('REMARK energies:') == 0:
			vdw = float(bound[i].split(', ')[5])
			elec = float(bound[i].split(', ')[6])
			AIR = float(bound[i].split(', ')[7])
		
		if bound[i].find('REMARK Desolvation energy:') == 0:
			desolv = float(bound[i].split(': ')[1])
	
#	dic_HADDOCK[hd.split('/')[-1].split('.')[0]] = min(HDScore_list) - 1
	dic_HADDOCK[hd.split('/')[-1].split('.')[0]] = 1.0 * vdw + 0.2 * elec + 0.1 * AIR + 1 * desolv
	
	dic_HADDOCK_Score_List[hd.split('/')[-1].split('.')[0]].append(dic_HADDOCK[hd.split('/')[-1].split('.')[0]])

PDBs_Final = list(set(PSSM_PDBs).intersection(set(HD_PDBs)))
PDBs_Final.sort()

pn = 1
tn = 1
ori_tn = 1

for p in PDBs_Final:
	dic_HADDOCK_Evdw = ''
	dic_HADDOCK_Eelec = ''
	dic_HADDOCK_Edesolv = ''
	dic_HADDOCK_BSA = ''
	
	dn = 1
	
	ori_tn += 1
	
	print p + '        ' + 'p = ' + str(pn) + ', d = ' + str(dn) + ', processed_total = ' + str(tn) + ', original_total = ' + str(ori_tn)
	
	# Interfacial Residue Extraction for each decoy
	f_interf = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/pdbFLs_refined/Atom_distance_15/CA-CA_dist_8/' + p + '.contacts', 'r')
	interf = f_interf.read().split('\n')
	f_interf.close()
	
	## Discard decoys which have less than 10 interfacial pairs
#	if len(interf) - 1 < 10:
#		continue
	
	tn += 1
	
	# Class assignment
	c = '1'
	
	# HADDOCK values for each decoy
	HADDOCK_vals = os.popen('bash /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/extract_haddock_terms.sh /gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/pdbFLs_refined/' + p + '.pdb').read()
	
	dic_HADDOCK_Evdw = HADDOCK_vals.split('\n')[1].split('\t')[1]
	dic_HADDOCK_Eelec = HADDOCK_vals.split('\n')[1].split('\t')[2]
	dic_HADDOCK_Edesolv = HADDOCK_vals.split('\n')[1].split('\t')[3]
	dic_HADDOCK_BSA = HADDOCK_vals.split('\n')[1].split('\t')[4]
	
	dic_HADDOCK_Evdw_List[p] = []
	dic_HADDOCK_Eelec_List[p] = []
	dic_HADDOCK_Edesolv_List[p] = []
	dic_HADDOCK_BSA_List[p] = []
	
	dic_HADDOCK_Evdw_List[p].append(float(dic_HADDOCK_Evdw))
	dic_HADDOCK_Eelec_List[p].append(float(dic_HADDOCK_Eelec))
	dic_HADDOCK_Edesolv_List[p].append(float(dic_HADDOCK_Edesolv))
	dic_HADDOCK_BSA_List[p].append(float(dic_HADDOCK_BSA))

pn_U = 1
tn_U = 1
ori_tn_U = 1

for p_U in PDBs_Final:
	# Each decoy
	decoy_list_U = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/preproc_decoys/' + p_U + '_*w.pdb')
	
	decoy_PDBs_U = []
	
	for decoy_U in decoy_list_U:
		decoy_PDBs_U.append(decoy_U.split('/')[-1])
	
	decoy_PDBs_U.sort()
	
	dic_HADDOCK_Evdw_U = {}
	dic_HADDOCK_Eelec_U = {}
	dic_HADDOCK_Edesolv_U = {}
	dic_HADDOCK_BSA_U = {}
	
	dn_U = 1
	for d in decoy_PDBs_U:
		ori_tn_U += 1
		
		if d.find('w.pdb') == -1:
			continue
		
		print p_U + ' - ' + d + '        ' + 'p = ' + str(pn_U) + ', d = ' + str(dn_U) + ', processed_total = ' + str(tn_U) + ', original_total = ' + str(ori_tn_U)
		
		# Interfacial Residue Extraction for each decoy
		f_interf_U = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/preproc_decoys/Atom_distance_15/CA-CA_dist_8/' + d.replace('.pdb', '.contacts'), 'r')
		interf_U = f_interf_U.read().split('\n')
		f_interf_U.close()
		
		## Discard decoys which have less than 10 interfacial pairs
#		if len(interf_U) - 1 < 10:
#			continue
		
		tn_U += 1
		
		# HADDOCK values for each decoy
		HADDOCK_vals_U = os.popen('bash /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/extract_haddock_terms.sh /gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/' + d).read()
		
		dic_HADDOCK_Evdw_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[1]
		dic_HADDOCK_Eelec_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[2]
		dic_HADDOCK_Edesolv_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[3]
		dic_HADDOCK_BSA_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[4]
		
		dic_HADDOCK_Evdw_List[d.split('_')[0]].append(float(HADDOCK_vals_U.split('\n')[1].split('\t')[1]))
		dic_HADDOCK_Eelec_List[d.split('_')[0]].append(float(HADDOCK_vals_U.split('\n')[1].split('\t')[2]))
		dic_HADDOCK_Edesolv_List[d.split('_')[0]].append(float(HADDOCK_vals_U.split('\n')[1].split('\t')[3]))
		dic_HADDOCK_BSA_List[d.split('_')[0]].append(float(HADDOCK_vals_U.split('\n')[1].split('\t')[4]))
		
		dn_U += 1
	
	pn_U += 1

dic_HADDOCK_Evdw_Max = {}
dic_HADDOCK_Eelec_Max = {}
dic_HADDOCK_Edesolv_Max = {}
dic_HADDOCK_BSA_Max = {}
dic_HADDOCK_Score_Max = {}

dic_HADDOCK_Evdw_Min = {}
dic_HADDOCK_Eelec_Min = {}
dic_HADDOCK_Edesolv_Min = {}
dic_HADDOCK_BSA_Min = {}
dic_HADDOCK_Score_Min = {}

for k in dic_HADDOCK_Score_List.keys():
	
	elements_Evdw_List = np.array(dic_HADDOCK_Evdw_List[k])
	elements_Eelec_List = np.array(dic_HADDOCK_Eelec_List[k])
	elements_Edesolv_List = np.array(dic_HADDOCK_Edesolv_List[k])
	elements_BSA_List = np.array(dic_HADDOCK_BSA_List[k])
	elements_Score_List = np.array(dic_HADDOCK_Score_List[k])
	
	dic_HADDOCK_Evdw_Mean[k] = np.mean(elements_Evdw_List, axis = 0)
	dic_HADDOCK_Eelec_Mean[k] = np.mean(elements_Eelec_List, axis = 0)
	dic_HADDOCK_Edesolv_Mean[k] = np.mean(elements_Edesolv_List, axis = 0)
	dic_HADDOCK_BSA_Mean[k] = np.mean(elements_BSA_List, axis = 0)
	dic_HADDOCK_Score_Mean[k] = np.mean(elements_Score_List, axis = 0)
	
	dic_HADDOCK_Evdw_SD[k] = np.std(elements_Evdw_List, axis = 0)
	dic_HADDOCK_Eelec_SD[k] = np.std(elements_Eelec_List, axis = 0)
	dic_HADDOCK_Edesolv_SD[k] = np.std(elements_Edesolv_List, axis = 0)
	dic_HADDOCK_BSA_SD[k] = np.std(elements_BSA_List, axis = 0)
	dic_HADDOCK_Score_SD[k] = np.std(elements_Score_List, axis = 0)
	
#	elements_Evdw_List
	num = 0
	
	new_hd_list = []
	
	for hd in elements_Evdw_List:
		if hd > dic_HADDOCK_Evdw_Mean[k] + 2 * dic_HADDOCK_Evdw_SD[k]:
			num += 1
	
		else:
			new_hd_list.append(hd)
	
	if num != 0:
		print k + '\t' + str(num), "Evdw"
	
	dic_HADDOCK_Evdw_Min[k] = min(new_hd_list)
	dic_HADDOCK_Evdw_Max[k] = max(new_hd_list)
	
#	elements_Eelec_List
	num = 0
	
	new_hd_list = []
	
	for hd in elements_Eelec_List:
		if hd > dic_HADDOCK_Eelec_Mean[k] + 2 * dic_HADDOCK_Eelec_SD[k]:
			num += 1
	
		else:
			new_hd_list.append(hd)
	
	if num != 0:
		print k + '\t' + str(num), "Eelec"
	
	dic_HADDOCK_Eelec_Min[k] = min(new_hd_list)
	dic_HADDOCK_Eelec_Max[k] = max(new_hd_list)
	
	
#	elements_Edesolv_List
	num = 0
	
	new_hd_list = []
	
	for hd in elements_Edesolv_List:
		if hd > dic_HADDOCK_Edesolv_Mean[k] + 2 * dic_HADDOCK_Edesolv_SD[k]:
			num += 1
		
		else:
			new_hd_list.append(hd)
	
	if num != 0:
		print k + '\t' + str(num), "Edesolv"
	
	dic_HADDOCK_Edesolv_Min[k] = min(new_hd_list)
	dic_HADDOCK_Edesolv_Max[k] = max(new_hd_list)
	
	
#	elements_BSA_List
	num = 0
	
	new_hd_list = []
	
	for hd in elements_BSA_List:
		if hd > dic_HADDOCK_BSA_Mean[k] + 2 * dic_HADDOCK_BSA_SD[k]:
			num += 1
		
		else:
			new_hd_list.append(hd)
	
	if num != 0:
		print k + '\t' + str(num), "BSA"
	
	dic_HADDOCK_BSA_Min[k] = min(new_hd_list)
	dic_HADDOCK_BSA_Max[k] = max(new_hd_list)
	
	
#	elements_Score_List
	num = 0
	
	new_hd_list = []
	
	for hd in elements_Score_List:
		if hd > dic_HADDOCK_Score_Mean[k] + 2 * dic_HADDOCK_Score_SD[k]:
			num += 1
		
		else:
			new_hd_list.append(hd)
	
	if num != 0:
		print k + '\t' + str(num), "Score"
	
	dic_HADDOCK_Score_Min[k] = min(new_hd_list)
	dic_HADDOCK_Score_Max[k] = max(new_hd_list)


PSSM_list = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/AddResNumPSSM_seqFLs_bound/*.protein1.ResNumPSSM')

PSSM_PDBs = []

for pssm in PSSM_list:
	PSSM_PDBs.append(pssm.split('/')[-1].replace(".protein1.ResNumPSSM", ""))

PSSM_PDBs.sort()

HAD_Score_list = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/haddockScore/water/*.haddockScore')

# HADDOCK score for each decoy

dic_HADDOCK = {}

HD_PDBs = []

for hd in HAD_Score_list:
	
	HD_PDBs.append(hd.split('/')[-1].split('.')[0])
	
	f_HAD = open(hd, 'r')
	HAD = f_HAD.read().split('\n')
	f_HAD.close()
	
	HDScore_list = []
	
	for i in range(0, len(HAD) - 1):
		dic_HADDOCK[HAD[i].split()[0].split('.')[0]] = float(HAD[i].split()[1])
		HDScore_list.append(float(HAD[i].split()[1]))
	
	f_bound = open("/gpfs/group/vuh14/legacy/YJ_DOCK_PP/BM4_dimers/bound/pdbFLs_refined/" + hd.split('/')[-1].split('.')[0] + '.pdb', 'r')
	bound = f_bound.read().split('\n')
	f_bound.close()
	
	vdw = 0
	elec = 0
	AIR = 0
	desolv = 0
	
	for i in range(0, len(bound) - 1):
		if bound[i].find('ATOM') == 0:
			break
		
		if bound[i].find('REMARK energies:') == 0:
			vdw = float(bound[i].split(', ')[5])
			elec = float(bound[i].split(', ')[6])
			AIR = float(bound[i].split(', ')[7])
		
		if bound[i].find('REMARK Desolvation energy:') == 0:
			desolv = float(bound[i].split(': ')[1])
	
#	dic_HADDOCK[hd.split('/')[-1].split('.')[0]] = min(HDScore_list) - 1
	dic_HADDOCK[hd.split('/')[-1].split('.')[0]] = 1.0 * vdw + 0.2 * elec + 0.1 * AIR + 1 * desolv

PDBs_Final = list(set(PSSM_PDBs).intersection(set(HD_PDBs)))
PDBs_Final.sort()

outputDIR2 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/Bound_HDScore-HDTerms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormHDScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA/'

os.system('rm -rf ' + outputDIR2)

os.system('mkdir ' + outputDIR2)

fw2 = open(outputDIR2 + 'Training.csv', 'w')

pn = 1
tn = 1
ori_tn = 1

for p in PDBs_Final:

	# PSSM for two unbound proteins
	f_PSSM_1 = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/AddResNumPSSM_seqFLs_bound/' + p + '.protein1.ResNumPSSM', 'r')
	PSSM1 = f_PSSM_1.read().split('\n')
	f_PSSM_1.close()
	
	f_PSSM_2 = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/AddResNumPSSM_seqFLs_bound/' + p + '.protein2.ResNumPSSM', 'r')
	PSSM2 = f_PSSM_2.read().split('\n')
	f_PSSM_2.close()
	
	dic_PSSM1_IC = {}
	dic_PSSM1_AA = {}
	dic_PSSM1_Value = {}
	
	PSSM_AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	for i in range(4, len(PSSM1) - 7):
		if PSSM1[i] == '':
			continue
		
		dic_PSSM1_IC[PSSM1[i].split()[0]] = PSSM1[i].split()[-2]
		dic_PSSM1_AA[PSSM1[i].split()[0]] = PSSM1[i].split()[2]
		
		for q in range(3, 23):
			dic_PSSM1_Value[PSSM1[i].split()[0] + '_' + PSSM_AAs[q - 3]] = round(sigmoid(float(PSSM1[i].split()[q])), 3)
	
	dic_PSSM2_IC = {}
	dic_PSSM2_AA = {}
	dic_PSSM2_Value = {}
	
	for i in range(4, len(PSSM2) - 7):
		if PSSM2[i] == '':
			continue
		
		dic_PSSM2_IC[PSSM2[i].split()[0]] = PSSM2[i].split()[-2]
		dic_PSSM2_AA[PSSM2[i].split()[0]] = PSSM2[i].split()[2]
		
		for q in range(3, 23):
			dic_PSSM2_Value[PSSM2[i].split()[0] + '_' + PSSM_AAs[q - 3]] = round(sigmoid(float(PSSM2[i].split()[q])), 3)
	
	dic_HADDOCK_Evdw = ''
	dic_HADDOCK_Eelec = ''
	dic_HADDOCK_Edesolv = ''
	dic_HADDOCK_BSA = ''
	
	dn = 1
	
	ori_tn += 1
	
	print p + '        ' + 'p = ' + str(pn) + ', d = ' + str(dn) + ', processed_total = ' + str(tn) + ', original_total = ' + str(ori_tn)

	# Interfacial Residue Extraction for each decoy
	f_interf = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/pdbFLs_refined/Atom_distance_15/CA-CA_dist_8/' + p + '.contacts', 'r')
	interf = f_interf.read().split('\n')
	f_interf.close()

	## Discard decoys which have less than 10 interfacial pairs
	if len(interf) - 1 < 10:
		continue

	tn += 1

	# Class assignment
	c = '1'

	# HADDOCK values for each decoy
	HADDOCK_vals = os.popen('bash /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/extract_haddock_terms.sh /gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/pdbFLs_refined/' + p + '.pdb').read()

	dic_HADDOCK_Evdw = HADDOCK_vals.split('\n')[1].split('\t')[1]
	dic_HADDOCK_Eelec = HADDOCK_vals.split('\n')[1].split('\t')[2]
	dic_HADDOCK_Edesolv = HADDOCK_vals.split('\n')[1].split('\t')[3]
	dic_HADDOCK_BSA = HADDOCK_vals.split('\n')[1].split('\t')[4]

	#		HADDOCK_vals = os.popen('bash extract_haddock_terms.sh /data/benchmark/docking-benchmark4/runs-cmrestraints/1A2K/run1/structures/it1/water/1A2K_51w.pdb').read()
	#
	#		HADDOCK_Evdw = HADDOCK_vals.split('\n')[1].split('\t')[1]
	#		HADDOCK_Eelec = HADDOCK_vals.split('\n')[1].split('\t')[2]
	#		HADDOCK_Edesolv = HADDOCK_vals.split('\n')[1].split('\t')[3]
	#		HADDOCK_BSA = HADDOCK_vals.split('\n')[1].split('\t')[4]
	#
	#		print HADDOCK_Evdw, HADDOCK_Eelec, HADDOCK_Edesolv, HADDOCK_BSA

	dic_Inter_P1 = []
	dic_Inter_P2 = []
	Inter_Dist_List = []
	Inter_ResA_List = []
	Inter_ResB_List = []

	IP_sum = 0
	IP_avg = 0

	IPwPSSM_sum = 0
	IPwPSSM_avg = 0

	dist_IPwPSSM_sum = 0
	dist_IPwPSSM_avg = 0

	num_error = 0

	for i in range(0, len(interf) - 1):
		dic_Inter_P1.append(interf[i].split('\t')[2])
		dic_Inter_P2.append(interf[i].split('\t')[5])
	#		dic_Inter_Dist[interf[i].split('\t')[2] + '-' + interf[i].split('\t')[5]] = interf[i].split('\t')[6]
		Inter_Dist_List.append(float(interf[i].split('\t')[6]))
		Inter_ResA_List.append(interf[i].split('\t')[1] + '_' + interf[i].split('\t')[2])
		Inter_ResB_List.append(interf[i].split('\t')[4] + '_' + interf[i].split('\t')[5])
		
		try:
			IP_sum += ppIP_dic[dic_PSSM1_AA[interf[i].split('\t')[2]] + '_' + dic_PSSM2_AA[interf[i].split('\t')[5]]]

		except KeyError:
			num_error += 1
		
		for aa1 in PSSM_AAs:
			for aa2 in PSSM_AAs:
				try:
					IPwPSSM_sum += ppIP_dic[aa1 + '_' + aa2] * dic_PSSM1_Value[interf[i].split('\t')[2] + '_' + aa1] * dic_PSSM2_Value[interf[i].split('\t')[5] + '_' + aa2]

				except KeyError:
					continue
		
				try:
					dist_IPwPSSM_sum += (ppIP_dic[aa1 + '_' + aa2] * dic_PSSM1_Value[interf[i].split('\t')[2] + '_' + aa1] * dic_PSSM2_Value[interf[i].split('\t')[5] + '_' + aa2]) / float(interf[i].split('\t')[6])

				except KeyError:
					continue
	
	IP_avg = IP_sum / float(len(interf) - num_error)
	IPwPSSM_avg = IPwPSSM_sum / float(len(interf) - num_error)
	dist_IPwPSSM_avg = dist_IPwPSSM_sum / float(len(interf) - num_error)

#	redund_dic_Inter_P1 = dic_Inter_P1
#	redund_dic_Inter_P2 = dic_Inter_P2
	
	Inter_ResA_List = list(set(Inter_ResA_List))
	Inter_ResB_List = list(set(Inter_ResB_List))
	
	Inter_Res_List = list(set(Inter_ResA_List + Inter_ResB_List))
	Inter_Dist_List.sort(reverse = True)
	
	dic_Inter_P1 = list(set(dic_Inter_P1))
	dic_Inter_P2 = list(set(dic_Inter_P2))
	
#	redund_IC_list_P1 = []
#	redund_IC_list_P2 = []
	
	IC_list_P1 = []
	IC_list_P2 = []
	
	for int_p1 in dic_Inter_P1:
		try:
			IC_list_P1.append(float(dic_PSSM1_IC[int_p1]))
		
		except KeyError:
			continue
		
#	for r_int_p1 in redund_dic_Inter_P1:
#		try:
#			redund_IC_list_P1.append(float(dic_PSSM1_IC[r_int_p1]))
#
#		except KeyError:
#			continue
	
	for int_p2 in dic_Inter_P2:
		try:
			IC_list_P2.append(float(dic_PSSM2_IC[int_p2]))
		
		except KeyError:
			continue
	
	IC_list_All = IC_list_P1 + IC_list_P2
	IC_list_All.sort(reverse = True)
	
#	for r_int_p2 in redund_dic_Inter_P2:
#		try:
#			redund_IC_list_P2.append(float(dic_PSSM2_IC[r_int_p2]))
#
#		except KeyError:
#			continue
	
	final_Score = 0
	if float(dic_HADDOCK[p]) > dic_HADDOCK_Score_Mean[p] + 2 * dic_HADDOCK_Score_SD[p]:
		final_Score = dic_HADDOCK_Score_Max[p]
	
	else:
		final_Score = float(dic_HADDOCK[p])
	
	final_norm_Score = (final_Score - dic_HADDOCK_Score_Min[p]) / (dic_HADDOCK_Score_Max[p] - dic_HADDOCK_Score_Min[p])
	
	final_Evdw = 0
	if float(dic_HADDOCK_Evdw) > dic_HADDOCK_Evdw_Mean[p] + 2 * dic_HADDOCK_Evdw_SD[p]:
		final_Evdw= dic_HADDOCK_Evdw_Max[p]
	
	else:
		final_Evdw = float(dic_HADDOCK_Evdw)
	
	final_norm_Evdw = (final_Evdw - dic_HADDOCK_Evdw_Min[p]) / (dic_HADDOCK_Evdw_Max[p] - dic_HADDOCK_Evdw_Min[p])
	
	final_Eelec = 0
	if float(dic_HADDOCK_Eelec) > dic_HADDOCK_Eelec_Mean[p] + 2 * dic_HADDOCK_Eelec_SD[p]:
		final_Eelec= dic_HADDOCK_Eelec_Max[p]
	
	else:
		final_Eelec = float(dic_HADDOCK_Eelec)
	
	final_norm_Eelec = (final_Eelec - dic_HADDOCK_Eelec_Min[p]) / (dic_HADDOCK_Eelec_Max[p] - dic_HADDOCK_Eelec_Min[p])
	
	final_Edesolv = 0
	if float(dic_HADDOCK_Edesolv) > dic_HADDOCK_Edesolv_Mean[p] + 2 * dic_HADDOCK_Edesolv_SD[p]:
		final_Edesolv = dic_HADDOCK_Edesolv_Max[p]
	
	else:
		final_Edesolv = float(dic_HADDOCK_Edesolv)
	
	final_norm_Edesolv = (final_Edesolv - dic_HADDOCK_Edesolv_Min[p]) / (dic_HADDOCK_Edesolv_Max[p] - dic_HADDOCK_Edesolv_Min[p])
	
	final_BSA = 0
	if float(dic_HADDOCK_BSA) > dic_HADDOCK_BSA_Mean[p] + 2 * dic_HADDOCK_BSA_SD[p]:
		final_BSA = dic_HADDOCK_BSA_Max[p]
	
	else:
		final_BSA = float(dic_HADDOCK_BSA)
	
	final_norm_BSA = (final_BSA - dic_HADDOCK_BSA_Min[p]) / (dic_HADDOCK_BSA_Max[p] - dic_HADDOCK_BSA_Min[p])
	
	print >>fw2, p + ',' + c + ',0,' + str(dic_HADDOCK[p]) + ',' + str(dic_HADDOCK_Evdw) + ',' + str(dic_HADDOCK_Eelec) + ',' + str(dic_HADDOCK_Edesolv) + ',' + str(dic_HADDOCK_BSA) + ',' + str(sum(IC_list_P1 + IC_list_P2) / len(IC_list_P1 + IC_list_P2)) + ',' + str(max(IC_list_P1 + IC_list_P2)) + ',' + str(min(IC_list_P1 + IC_list_P2)) + ',' + str(sum(IC_list_P1) / len(IC_list_P1)) + ',' + str(max(IC_list_P1)) + ',' + str(min(IC_list_P1)) + ',' + str(sum(IC_list_P2) / len(IC_list_P2)) + ',' + str(max(IC_list_P2)) + ',' + str(min(IC_list_P2)) + ',' + str(Inter_Dist_List[9]) + ',' + str(Inter_Dist_List[8]) + ',' + str(Inter_Dist_List[7]) + ',' + str(Inter_Dist_List[6]) + ',' + str(Inter_Dist_List[5]) + ',' + str(Inter_Dist_List[4]) + ',' + str(Inter_Dist_List[3]) + ',' + str(Inter_Dist_List[2]) + ',' + str(Inter_Dist_List[1]) + ',' + str(Inter_Dist_List[0]) + ',' + str(len(interf)) + ',' + str(len(Inter_Res_List)) + ',' + str(float(len(interf)) / float(len(Inter_ResA_List) * len(Inter_ResB_List))) + ','  + str(final_norm_Score) + ',' + str(float(final_norm_Evdw)) + ',' + str(float(final_norm_Eelec)) + ',' + str(float(final_norm_Edesolv)) + ',' + str(float(final_norm_BSA)) + ',' + str(IP_sum) + ',' + str(IP_avg) + ',' + str(IPwPSSM_sum) + ',' + str(IPwPSSM_avg) + ',' + str(dist_IPwPSSM_sum) + ',' + str(dist_IPwPSSM_avg)
	
	
	#	print >>fw2, p + ',' + c + ',0,' + str(dic_HADDOCK_Evdw) + ',' + str(dic_HADDOCK_Eelec) + ',' + str(dic_HADDOCK_Edesolv) + ',' + str(dic_HADDOCK_BSA) + ',' + str(sum(redund_IC_list_P1 + redund_IC_list_P2) / len(redund_IC_list_P1 + redund_IC_list_P2)) + ',' + str(max(redund_IC_list_P1 + redund_IC_list_P2)) + ',' + str(min(redund_IC_list_P1 + redund_IC_list_P2))
	
	dn += 1
pn += 1

fw2.close()


#For Decoys

PSSM_list_U = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/AddResNumPSSM_seqFLs_decoys/*.protein1.ResNumPSSM')
IRMSD_list_U = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/i-rmsd/water/*.irmsd')

PSSM_PDBs_U = []
IRMSD_PDBs_U = []

for pssm_U in PSSM_list_U:
	PSSM_PDBs_U.append(pssm_U.split('/')[-1].replace(".protein1.ResNumPSSM", ""))

PSSM_PDBs_U.sort()

for irmsd_U in IRMSD_list_U: # only Water models are considered
	IRMSD_PDBs_U.append(irmsd_U.split('/')[-1].replace(".irmsd", ""))

IRMSD_PDBs_U.sort()

PDBs_Final_U = list(set(PSSM_PDBs_U).intersection(set(IRMSD_PDBs_U)))
PDBs_Final_U.sort()

outputDIR2_U = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/Decoy_HDScore-HDTerms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormHDScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA/'

os.system('rm -rf ' + outputDIR2_U)

os.system('mkdir ' + outputDIR2_U)

fw2 = open(outputDIR2_U + 'Training.csv', 'w')

pn_U = 1
tn_U = 1
ori_tn_U = 1

for p_U in PDBs_Final_U:
	
	# PSSM for two unbound proteins
	f_PSSM_1_U = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/AddResNumPSSM_seqFLs_decoys/' + p_U + '.protein1.ResNumPSSM', 'r')
	PSSM1_U = f_PSSM_1_U.read().split('\n')
	f_PSSM_1_U.close()

	f_PSSM_2_U = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/AddResNumPSSM_seqFLs_decoys/' + p_U + '.protein2.ResNumPSSM', 'r')
	PSSM2_U = f_PSSM_2_U.read().split('\n')
	f_PSSM_2_U.close()

	dic_PSSM1_IC_U = {}
	dic_PSSM1_AA_U = {}
	dic_PSSM1_Value_U = {}
	
	PSSM_AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	for i in range(4, len(PSSM1_U) - 7):
		if PSSM1_U[i] == '':
			continue
		
		dic_PSSM1_IC_U[PSSM1_U[i].split()[0]] = PSSM1_U[i].split()[-2]
		dic_PSSM1_AA_U[PSSM1_U[i].split()[0]] = PSSM1_U[i].split()[2]
		
		for q in range(3, 23):
			dic_PSSM1_Value_U[PSSM1_U[i].split()[0] + '_' + PSSM_AAs[q - 3]] = round(sigmoid(float(PSSM1_U[i].split()[q])), 3)
	
	dic_PSSM2_IC_U = {}
	dic_PSSM2_AA_U = {}
	dic_PSSM2_Value_U = {}

	for i in range(4, len(PSSM2_U) - 7):
		if PSSM2_U[i] == '':
			continue
		
		dic_PSSM2_IC_U[PSSM2_U[i].split()[0]] = PSSM2_U[i].split()[-2]
		dic_PSSM2_AA_U[PSSM2_U[i].split()[0]] = PSSM2_U[i].split()[2]
		
		for q in range(3, 23):
			dic_PSSM2_Value_U[PSSM2_U[i].split()[0] + '_' + PSSM_AAs[q - 3]] = round(sigmoid(float(PSSM2_U[i].split()[q])), 3)
	
	# I-RMSD for each decoy
	dic_IRMSD_U = {}

	f_IRMSD_U = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/i-rmsd/water/' + p_U + '.irmsd', 'r')
	IRMSD_U = f_IRMSD_U.read().split('\n')
	f_IRMSD_U.close()

	for i in range(0, len(IRMSD_U) - 1):
		dic_IRMSD_U[IRMSD_U[i].split()[0].split('.')[0]] = IRMSD_U[i].split()[1]

	# Each decoy
	decoy_list_U = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/' + p_U + '_*w.pdb')

	decoy_PDBs_U = []

	for decoy_U in decoy_list_U:
		decoy_PDBs_U.append(decoy_U.split('/')[-1])

	decoy_PDBs_U.sort()

	dic_HADDOCK_Evdw_U = {}
	dic_HADDOCK_Eelec_U = {}
	dic_HADDOCK_Edesolv_U = {}
	dic_HADDOCK_BSA_U = {}

	dn_U = 1
	for d in decoy_PDBs_U:
		ori_tn_U += 1

		if d.find('w.pdb') == -1:
			continue

		print p_U + ' - ' + d + '        ' + 'p = ' + str(pn_U) + ', d = ' + str(dn_U) + ', processed_total = ' + str(tn_U) + ', original_total = ' + str(ori_tn_U)

		# Interfacial Residue Extraction for each decoy
		f_interf_U = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/preproc_decoys/Atom_distance_15/CA-CA_dist_8/' + d.replace('.pdb', '.contacts'), 'r')
		interf_U = f_interf_U.read().split('\n')
		f_interf_U.close()

		## Discard decoys which have less than 10 interfacial pairs
		if len(interf_U) - 1 < 10:
			continue

		tn_U += 1

		# Class assignment
		c_U = '0'
		if float(dic_IRMSD_U[d.replace('.pdb', '')]) <= 4:
			c_U = '1'

		# HADDOCK values for each decoy
		HADDOCK_vals_U = os.popen('bash /gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/extract_haddock_terms.sh /gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/' + d).read()
		
		dic_HADDOCK_Evdw_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[1]
		dic_HADDOCK_Eelec_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[2]
		dic_HADDOCK_Edesolv_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[3]
		dic_HADDOCK_BSA_U[d.split('.')[0]] = HADDOCK_vals_U.split('\n')[1].split('\t')[4]

		dic_Inter_P1_U = []
		dic_Inter_P2_U= []
		Inter_Dist_List_U = []
		Inter_ResA_List_U = []
		Inter_ResB_List_U = []
		
		IP_sum_U = 0
		IP_avg_U = 0
		
		IPwPSSM_sum_U = 0
		IPwPSSM_avg_U = 0
		
		dist_IPwPSSM_sum_U = 0
		dist_IPwPSSM_avg_U = 0

		num_error = 0
		
		for i in range(0, len(interf_U) - 1):
			dic_Inter_P1_U.append(interf_U[i].split('\t')[2])
			dic_Inter_P2_U.append(interf_U[i].split('\t')[5])
#			dic_Inter_Dist_U[interf_U[i].split('\t')[2] + '-' + interf_U[i].split('\t')[5]] = interf_U[i].split('\t')[6]
			Inter_Dist_List_U.append(float(interf_U[i].split('\t')[6]))
			Inter_ResA_List_U.append(interf_U[i].split('\t')[1] + '_' + interf_U[i].split('\t')[2])
			Inter_ResB_List_U.append(interf_U[i].split('\t')[4] + '_' + interf_U[i].split('\t')[5])
			
			try:
				IP_sum_U += ppIP_dic[dic_PSSM1_AA_U[interf_U[i].split('\t')[2]] + '_' + dic_PSSM2_AA_U[interf_U[i].split('\t')[5]]]
			
			except KeyError:
				num_error += 1
			
			for aa1 in PSSM_AAs:
				for aa2 in PSSM_AAs:
					
					try:
						IPwPSSM_sum_U += ppIP_dic[aa1 + '_' + aa2] * dic_PSSM1_Value_U[interf_U[i].split('\t')[2] + '_' + aa1] * dic_PSSM2_Value_U[interf_U[i].split('\t')[5] + '_' + aa2]

					except KeyError:
						continue

					try:
						dist_IPwPSSM_sum_U += float(ppIP_dic[aa1 + '_' + aa2] * dic_PSSM1_Value_U[interf_U[i].split('\t')[2] + '_' + aa1] * dic_PSSM2_Value_U[interf_U[i].split('\t')[5] + '_' + aa2]) / float(interf_U[i].split('\t')[6])

					except KeyError:
						continue

		IP_avg_U = IP_sum_U / float(len(interf_U) - num_error)
		IPwPSSM_avg_U = IPwPSSM_sum_U / float(len(interf_U) - num_error)
		dist_IPwPSSM_avg_U = dist_IPwPSSM_sum_U / float(len(interf_U) - num_error)

		Inter_ResA_List_U = list(set(Inter_ResA_List_U))
		Inter_ResB_List_U = list(set(Inter_ResB_List_U))

#		redund_dic_Inter_P1 = dic_Inter_P1
#		redund_dic_Inter_P2 = dic_Inter_P2
		
		Inter_Dist_List_U.sort(reverse = True)
		Inter_Res_List_U = list(set(Inter_ResA_List_U + Inter_ResB_List_U))

		dic_Inter_P1_U = list(set(dic_Inter_P1_U))
		dic_Inter_P2_U = list(set(dic_Inter_P2_U))
		
#		redund_IC_list_P1 = []
#		redund_IC_list_P2 = []
		
		IC_list_P1_U = []
		IC_list_P2_U = []
		
		for int_p1 in dic_Inter_P1_U:
			try:
				IC_list_P1_U.append(float(dic_PSSM1_IC_U[int_p1]))
	
			except KeyError:
				continue

#		for r_int_p1 in redund_dic_Inter_P1:
#			try:
#				redund_IC_list_P1.append(float(dic_PSSM1_IC[r_int_p1]))
#
#			except KeyError:
#				continue

		for int_p2 in dic_Inter_P2_U:
			try:
				IC_list_P2_U.append(float(dic_PSSM2_IC_U[int_p2]))

			except KeyError:
				continue

#		for r_int_p2 in redund_dic_Inter_P2:
#			try:
#				redund_IC_list_P2.append(float(dic_PSSM2_IC[r_int_p2]))
#
#			except KeyError:
#				continue

		IC_list_All_U = IC_list_P1_U + IC_list_P2_U
		IC_list_All_U.sort(reverse = True)

		final_Score = 0
		if float(dic_HADDOCK[d.replace('.pdb', '')]) > dic_HADDOCK_Score_Mean[p_U] + 2 * dic_HADDOCK_Score_SD[p_U]:
			final_Score = dic_HADDOCK_Score_Max[p_U]

		else:
			final_Score = float(dic_HADDOCK[d.replace('.pdb', '')])

		final_norm_Score = (final_Score - dic_HADDOCK_Score_Min[p_U]) / (dic_HADDOCK_Score_Max[p_U] - dic_HADDOCK_Score_Min[p_U])

		final_Evdw = 0
		if float(dic_HADDOCK_Evdw_U[d.split('.')[0]]) > dic_HADDOCK_Evdw_Mean[p_U] + 2 * dic_HADDOCK_Evdw_SD[p_U]:
			final_Evdw= dic_HADDOCK_Evdw_Max[p_U]

		else:
			final_Evdw = float(dic_HADDOCK_Evdw_U[d.split('.')[0]])

		final_norm_Evdw = (final_Evdw - dic_HADDOCK_Evdw_Min[p_U]) / (dic_HADDOCK_Evdw_Max[p_U] - dic_HADDOCK_Evdw_Min[p_U])

		final_Eelec = 0
		if float(dic_HADDOCK_Eelec_U[d.split('.')[0]]) > dic_HADDOCK_Eelec_Mean[p_U] + 2 * dic_HADDOCK_Eelec_SD[p_U]:
			final_Eelec= dic_HADDOCK_Eelec_Max[p_U]

		else:
			final_Eelec = float(dic_HADDOCK_Eelec_U[d.split('.')[0]])

		final_norm_Eelec = (final_Eelec - dic_HADDOCK_Eelec_Min[p_U]) / (dic_HADDOCK_Eelec_Max[p_U] - dic_HADDOCK_Eelec_Min[p_U])

		final_Edesolv = 0
		if float(dic_HADDOCK_Edesolv_U[d.split('.')[0]]) > dic_HADDOCK_Edesolv_Mean[p_U] + 2 * dic_HADDOCK_Edesolv_SD[p_U]:
			final_Edesolv = dic_HADDOCK_Edesolv_Max[p_U]

		else:
			final_Edesolv = float(dic_HADDOCK_Edesolv_U[d.split('.')[0]])

		final_norm_Edesolv = (final_Edesolv - dic_HADDOCK_Edesolv_Min[p_U]) / (dic_HADDOCK_Edesolv_Max[p_U] - dic_HADDOCK_Edesolv_Min[p_U])

		final_BSA = 0
		if float(dic_HADDOCK_BSA_U[d.split('.')[0]]) > dic_HADDOCK_BSA_Mean[p_U] + 2 * dic_HADDOCK_BSA_SD[p_U]:
			final_BSA = dic_HADDOCK_BSA_Max[p_U]

		else:
			final_BSA = float(dic_HADDOCK_BSA_U[d.split('.')[0]])

		final_norm_BSA = (final_BSA - dic_HADDOCK_BSA_Min[p_U]) / (dic_HADDOCK_BSA_Max[p_U] - dic_HADDOCK_BSA_Min[p_U])

		print >>fw2, d + ',' + c_U + ',' + dic_IRMSD_U[d.replace('.pdb', '')] + ',' + str(dic_HADDOCK[d.replace('.pdb', '')]) + ',' + str(dic_HADDOCK_Evdw_U[d.split('.')[0]]) + ',' + str(dic_HADDOCK_Eelec_U[d.split('.')[0]]) + ',' + str(dic_HADDOCK_Edesolv_U[d.split('.')[0]]) + ',' + str(dic_HADDOCK_BSA_U[d.split('.')[0]]) + ',' + str(sum(IC_list_P1_U + IC_list_P2_U) / len(IC_list_P1_U + IC_list_P2_U)) + ',' + str(max(IC_list_P1_U + IC_list_P2_U)) + ',' + str(min(IC_list_P1_U + IC_list_P2_U)) + ',' + str(sum(IC_list_P1_U) / len(IC_list_P1_U)) + ',' + str(max(IC_list_P1_U)) + ',' + str(min(IC_list_P1_U)) + ',' + str(sum(IC_list_P2_U) / len(IC_list_P2_U)) + ',' + str(max(IC_list_P2_U)) + ',' + str(min(IC_list_P2_U)) + ',' + str(Inter_Dist_List_U[9]) + ',' + str(Inter_Dist_List_U[8]) + ',' + str(Inter_Dist_List_U[7]) + ',' + str(Inter_Dist_List_U[6]) + ',' + str(Inter_Dist_List_U[5]) + ',' + str(Inter_Dist_List_U[4]) + ',' + str(Inter_Dist_List_U[3]) + ',' + str(Inter_Dist_List_U[2]) + ',' + str(Inter_Dist_List_U[1]) + ',' + str(Inter_Dist_List_U[0]) + ',' + str(len(interf_U)) + ',' + str(len(Inter_Res_List_U)) + ',' + str(float(len(interf_U)) / float(len(Inter_ResA_List_U) * len(Inter_ResB_List_U))) + ','  + str(final_norm_Score) + ',' + str(float(final_norm_Evdw)) + ',' + str(float(final_norm_Eelec)) + ',' + str(float(final_norm_Edesolv)) + ',' + str(float(final_norm_BSA)) + ',' + str(IP_sum_U) + ',' + str(IP_avg_U) + ',' + str(IPwPSSM_sum_U) + ',' + str(IPwPSSM_avg_U) + ',' + str(dist_IPwPSSM_sum_U) + ',' + str(dist_IPwPSSM_avg_U)
#		print >>fw2, d + ',' + c + ',' + dic_IRMSD[d.replace('.pdb', '')] + ',' + str(dic_HADDOCK_Evdw[d.split('.')[0]]) + ',' + str(dic_HADDOCK_Eelec[d.split('.')[0]]) + ',' + str(dic_HADDOCK_Edesolv[d.split('.')[0]]) + ',' + str(dic_HADDOCK_BSA[d.split('.')[0]]) + ',' + str(sum(redund_IC_list_P1 + redund_IC_list_P2) / len(redund_IC_list_P1 + redund_IC_list_P2)) + ',' + str(max(redund_IC_list_P1 + redund_IC_list_P2)) + ',' + str(min(redund_IC_list_P1 + redund_IC_list_P2))

		dn_U += 1
	pn_U += 1

fw2.close()
