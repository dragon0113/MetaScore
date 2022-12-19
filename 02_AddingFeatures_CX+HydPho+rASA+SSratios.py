#PSSM_IC A & B, HADDOCK score added as features

# Ex of execution: python 01-1_AddingFeatures_CX+HydPho+rASA+SSratios.py 1AK4
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

p = sys.argv[1]

def sigmoid(x):
	return 1 / (1 + math.exp(-x))

#surface_area

aa_sa = {}
fr_sa = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/code_BM4/surface_area.txt', 'r')
sa_lines = fr_sa.read().split('\n')
fr_sa.close()

for i in range(2, len(sa_lines) - 1):
	aa_sa[sa_lines[i].split()[1]] = float(sa_lines[i].split()[3])

# Hydropathy of AAs
hydpho = {}
hydpho['ALA'] = 1.8
hydpho['ARG'] = -4.5
hydpho['ASN'] = -3.5
hydpho['ASP'] = -3.5
hydpho['CYS'] = 2.5
hydpho['GLU'] = -3.5
hydpho['GLN'] = -3.5
hydpho['GLY'] = -0.4
hydpho['HIS'] = -3.2
hydpho['ILE'] = 4.5
hydpho['LEU'] = 3.8
hydpho['LYS'] = -3.9
hydpho['MET'] = 1.9
hydpho['PHE'] = 2.8
hydpho['PRO'] = -1.6
hydpho['SER'] = -0.8
hydpho['THR'] = -0.7
hydpho['TRP'] = -0.9
hydpho['TYR'] = -1.3
hydpho['VAL'] = 4.2

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

ama_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

n = 1
tn = 1
ori_tn = 1

dic_CX_mean = {}
dic_CX_SD = {}
dic_Hydphob = {}
dic_rasa = {}
dic_ss_code_H = {}
dic_ss_code_G = {}
dic_ss_code_I = {}
dic_ss_code_E = {}
dic_ss_code_B = {}
dic_ss_code_T = {}
dic_ss_code_C = {}

print p

#	dn = 1

#	ori_tn += 1

#	print p + '        ' + 'p = ' + str(pn) + ', d = ' + str(dn) + ', processed_total = ' + str(tn) + ', original_total = ' + str(ori_tn)

# Interfacial Residue Extraction for each decoy
f_interf = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/bound/pdbFLs_refined/Atom_distance_15/CA-CA_dist_8/' + p + '.contacts', 'r')
interf = f_interf.read().split('\n')
f_interf.close()
	
f_bound = open("/gpfs/group/vuh14/legacy/YJ_DOCK_PP/BM4_dimers/bound/pdbFLs_refined/" + p + '.pdb', 'r')
bound = f_bound.read().split('\n')
f_bound.close()

# rASA and secondary structure ratio

stride_val = os.popen("/gpfs/group/vuh14/legacy/lixue/bin/stride/stride /gpfs/group/vuh14/legacy/YJ_DOCK_PP/BM4_dimers/bound/pdbFLs_refined/" + p + '.pdb | grep "ASG " | grep -v "SEQ"').read().split('\n')

pro_rasa = {}
pro_ss = {}

for i in range(0, len(stride_val) - 1):
	chain = stride_val[i].split()[2]
	res_idx = stride_val[i].split()[3]
	ss_type = stride_val[i].split()[5]
	rasa = stride_val[i].split()[9]
	aa_three = stride_val[i].split()[1]
	
	print res_idx, ss_type, rasa, aa_three
	
	pro_ss[chain + '_' + res_idx] = ss_type.upper()
	
	try:
		pro_rasa[chain + '_' + res_idx] = float(rasa) / aa_sa[aa_three]
	
	except KeyError:
		continue

# CX for CA

dt = 10.0 # radius 10A
va = 20.1 # atom volume
vs = (4.0 * 3.141592654 * dt * dt * dt / 3.0) # volume of the sphere

dic_CX = {}

for k in range(0, len(bound) - 1):
	if bound[k].find('ATOM') != 0:
		continue
	
	aa1 = bound[k][17:20].strip() # Residue name
	
	if not(aa1 in ama_list):
		continue
	
	at1 = bound[k][12:16].strip() # atom name
	
	if not(at1 == 'CA'):
		continue
	
	r1 = bound[k][22:27].strip() # resSeq number
	x1 = float(bound[k][30:38].strip()) # X coordinate
	y1 = float(bound[k][38:46].strip()) # Y coordinate
	z1 = float(bound[k][46:54].strip()) # Z coordinate
	chain1 = bound[k][71:].strip()
	
	num_atoms = 0
	
	for l in range(0, len(bound) - 1):
		if bound[l].find('ATOM') != 0:
			continue
		
		aa2 = bound[l][17:20].strip()
		
		if not(aa2 in ama_list):
			continue
		
		chain2 = bound[l][71:].strip()
		
		if not(chain1 == chain2):
			continue
		
		at2 = bound[l][12:16].strip() # atom name
		r2 = bound[l][22:27].strip()
		x2 = float(bound[l][30:38].strip())
		y2 = float(bound[l][38:46].strip())
		z2 = float(bound[l][46:54].strip())
		
		d = math.sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2))
#		print d
		if d <= float(dt):
			num_atoms += 1
	
	vi = num_atoms * va
	dic_CX[chain1 + '_' + r1] = (vs - vi) / vi
	
	if dic_CX[chain1 + '_' + r1] < 0:
		dic_CX[chain1 + '_' + r1] = 0

cx_list =[]
hydpho_sum = 0
rasa_sum = 0
ss_code = {}
ss_code['H'] = 0
ss_code['G'] = 0
ss_code['I'] = 0
ss_code['E'] = 0
ss_code['B'] = 0
ss_code['T'] = 0
ss_code['C'] = 0

for i in range(0, len(interf) - 1):
	try:
		cx1 = dic_CX[interf[i].split('\t')[1] + '_' + interf[i].split('\t')[2]] # CX of atoms in chain A
		cx2 = dic_CX[interf[i].split('\t')[4] + '_' + interf[i].split('\t')[5]] # CX of atoms in chain B
		hydpho1 = hydpho[interf[i].split('\t')[0]]
		hydpho2 = hydpho[interf[i].split('\t')[3]]
		rasa1 = pro_rasa[interf[i].split('\t')[1] + '_' + interf[i].split('\t')[2]]
		rasa2 = pro_rasa[interf[i].split('\t')[4] + '_' + interf[i].split('\t')[5]]
		ss_code[pro_ss[interf[i].split('\t')[1] + '_' + interf[i].split('\t')[2]]] += 1
		ss_code[pro_ss[interf[i].split('\t')[4] + '_' + interf[i].split('\t')[5]]] += 1

	except KeyError:
		continue
	
	cxFeat = (pow(max(cx1, cx2), 2) + 1) / (pow(min(cx1, cx2), 2) + 1)
	print cxFeat
	cx_list.append(cxFeat)

	hydpho_sum += (hydpho1 + hydpho2)
	rasa_sum += (rasa1 + rasa2)

cx_list_array = np.array(cx_list)

dic_CX_mean[p] = np.mean(cx_list_array, axis = 0)

dic_CX_SD[p] = np.std(cx_list_array, axis = 0)

dic_Hydphob[p] = hydpho_sum / (len(interf) - 1)

dic_rasa[p] = rasa_sum / (len(interf) - 1) / 2

dic_ss_code_H[p] = float(ss_code['H']) / (len(interf) - 1) / 2
dic_ss_code_G[p] = float(ss_code['G']) / (len(interf) - 1) / 2
dic_ss_code_I[p] = float(ss_code['I']) / (len(interf) - 1) / 2
dic_ss_code_E[p] = float(ss_code['E']) / (len(interf) - 1) / 2
dic_ss_code_B[p] = float(ss_code['B']) / (len(interf) - 1) / 2
dic_ss_code_T[p] = float(ss_code['T']) / (len(interf) - 1) / 2
dic_ss_code_C[p] = float(ss_code['C']) / (len(interf) - 1) / 2


decoy_list_U = glob.glob('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/' + p + '_*w.pdb')

decoy_PDBs_U = []

for decoy_U in decoy_list_U:
	decoy_PDBs_U.append(decoy_U.split('/')[-1])

decoy_PDBs_U.sort()

for d in decoy_PDBs_U:
	print d.replace('.pdb', '')
	
	if d.find('w.pdb') == -1:
		continue
	
	
	# Interfacial Residue Extraction for each decoy
	f_interf_U = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/preproc_decoys/Atom_distance_15/CA-CA_dist_8/' + d.replace('.pdb', '.contacts'), 'r')
	interf_U = f_interf_U.read().split('\n')
	f_interf_U.close()
	
	if len(interf_U) == 1:
		continue
	
	f_unbound = open('/gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/preproc_decoys/' + d, 'r')
	unbound = f_unbound.read().split('\n')
	f_unbound.close()
	
	stride_val_U = os.popen("/gpfs/group/vuh14/legacy/lixue/bin/stride/stride /gpfs/group/vuh14/legacy/YJ_DOCK_PP/GraphKernel_scoringFun/BM4_dimers/restraints_decoys/preproc_decoys/" + d + ' | grep "ASG " | grep -v "SEQ"').read().split('\n')
	
	pro_rasa_U = {}
	pro_ss_U = {}
	
	for i in range(0, len(stride_val_U) - 1):
		chain = stride_val_U[i].split()[2]
		res_idx = stride_val_U[i].split()[3]
		ss_type = stride_val_U[i].split()[5]
		rasa = stride_val_U[i].split()[9]
		aa_three = stride_val_U[i].split()[1]
		
		print res_idx, ss_type, rasa, aa_three
		
		pro_ss_U[chain + '_' + res_idx] = ss_type.upper()
		
		try:
			pro_rasa_U[chain + '_' + res_idx] = float(rasa) / aa_sa[aa_three]
		
		except KeyError:
			continue
	
	dic_CX_U = {}
	
	for k_U in range(0, len(unbound) - 1):
		if unbound[k_U].find('ATOM') != 0:
			continue
		
		aa1_U = unbound[k_U][17:20].strip() # Residue name
		
		if not(aa1_U in ama_list):
			continue
		
		at1_U = unbound[k_U][12:16].strip() # atom name
		
		if not(at1_U == 'CA'):
			continue
		
		r1_U = unbound[k_U][22:27].strip() # resSeq number
		x1_U = float(unbound[k_U][30:38].strip()) # X coordinate
		y1_U = float(unbound[k_U][38:46].strip()) # Y coordinate
		z1_U = float(unbound[k_U][46:54].strip()) # Z coordinate
		chain1_U = unbound[k_U][71:].strip()
		
		num_atoms_U = 0
		
		for l_U in range(0, len(unbound) - 1):
			if unbound[l_U].find('ATOM') != 0:
				continue
			
			aa2_U = unbound[l_U][17:20].strip()
			
			if not(aa2_U in ama_list):
				continue
			
			chain2_U = unbound[l_U][71:].strip()
			
			if not(chain1_U == chain2_U):
				continue
			
			at2_U = unbound[l_U][12:16].strip() # atom name
			r2_U = unbound[l_U][22:27].strip()
			x2_U = float(unbound[l_U][30:38].strip())
			y2_U = float(unbound[l_U][38:46].strip())
			z2_U = float(unbound[l_U][46:54].strip())
			
			d_U = math.sqrt(pow((x1_U - x2_U), 2) + pow((y1_U - y2_U), 2) + pow((z1_U - z2_U), 2))
	#		print d
			if d_U <= float(dt):
				num_atoms_U += 1
		
		vi_U = num_atoms_U * va
		dic_CX_U[chain1_U + '_' + r1_U] = (vs - vi_U) / vi_U
		
		if dic_CX_U[chain1_U + '_' + r1_U] < 0:
			dic_CX_U[chain1_U + '_' + r1_U] = 0
	
	cx_list_U =[]
	hydpho_sum_U = 0
	rasa_sum_U = 0
	ss_code_U = {}
	ss_code_U['H'] = 0
	ss_code_U['G'] = 0
	ss_code_U['I'] = 0
	ss_code_U['E'] = 0
	ss_code_U['B'] = 0
	ss_code_U['T'] = 0
	ss_code_U['C'] = 0
	
	for i in range(0, len(interf_U) - 1):
		try:
			cx1_U = dic_CX_U[interf_U[i].split('\t')[1] + '_' + interf_U[i].split('\t')[2]] # CX of atoms in chain A
			cx2_U = dic_CX_U[interf_U[i].split('\t')[4] + '_' + interf_U[i].split('\t')[5]] # CX of atoms in chain B
			hydpho1_U = hydpho[interf_U[i].split('\t')[0]]
			hydpho2_U = hydpho[interf_U[i].split('\t')[3]]
			rasa1_U = pro_rasa_U[interf_U[i].split('\t')[1] + '_' + interf_U[i].split('\t')[2]]
			rasa2_U = pro_rasa_U[interf_U[i].split('\t')[4] + '_' + interf_U[i].split('\t')[5]]
			ss_code_U[pro_ss_U[interf_U[i].split('\t')[1] + '_' + interf_U[i].split('\t')[2]]] += 1
			ss_code_U[pro_ss_U[interf_U[i].split('\t')[4] + '_' + interf_U[i].split('\t')[5]]] += 1
		
		except KeyError:
			continue
		
		cxFeat_U = (pow(max(cx1_U, cx2_U), 2) + 1) / (pow(min(cx1_U, cx2_U), 2) + 1)
		print cxFeat_U
		cx_list_U.append(cxFeat_U)
		
		hydpho_sum_U += (hydpho1_U + hydpho2_U)
		rasa_sum_U += (rasa1_U + rasa2_U)
	
	cx_list_array_U = np.array(cx_list_U)
	
	dic_CX_mean[d.replace('.pdb', '')] = np.mean(cx_list_array_U, axis = 0)
	
	dic_CX_SD[d.replace('.pdb', '')] = np.std(cx_list_array_U, axis = 0)
	
	dic_Hydphob[d.replace('.pdb', '')] = hydpho_sum_U / (len(interf_U) - 1)
	
	dic_rasa[d.replace('.pdb', '')] = rasa_sum_U / (len(interf_U) - 1) / 2
	
	dic_ss_code_H[d.replace('.pdb', '')] = float(ss_code_U['H']) / (len(interf_U) - 1) / 2
	dic_ss_code_G[d.replace('.pdb', '')] = float(ss_code_U['G']) / (len(interf_U) - 1) / 2
	dic_ss_code_I[d.replace('.pdb', '')] = float(ss_code_U['I']) / (len(interf_U) - 1) / 2
	dic_ss_code_E[d.replace('.pdb', '')] = float(ss_code_U['E']) / (len(interf_U) - 1) / 2
	dic_ss_code_B[d.replace('.pdb', '')] = float(ss_code_U['B']) / (len(interf_U) - 1) / 2
	dic_ss_code_T[d.replace('.pdb', '')] = float(ss_code_U['T']) / (len(interf_U) - 1) / 2
	dic_ss_code_C[d.replace('.pdb', '')] = float(ss_code_U['C']) / (len(interf_U) - 1) / 2

outputDIR1 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormHDScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/'

outputDIR2 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+PSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormHDScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA+CX+HydPho+rASA+SSratios_14.0AforNeg/Undersampled_Neg1x/'

#os.system('rm -rf ' + outputDIR1)

#os.system('mkdir ' + outputDIR1)
#os.system('mkdir ' + outputDIR2)

inputDIR1 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormHDScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA_14.0AforNeg/'

inputDIR2 = '/gpfs/group/vuh14/legacy/YJ_DOCK_PP/RFBased_ScoringFunction/processed_HDScore-HDTerms+InterfacialPSSMIC-Para3-AllP1P2+Dist-Top10+NumLink+NumIntRes+LinkDensity+NormHDScore-NormHDTerms+IP_SA+IPPSSM_SA+IPPSSMdist_SA_14.0AforNeg/Undersampled_Neg1x/'

fr1 = open(inputDIR1 + p + '.csv', 'r')
line1 = fr1.read().split('\n')
fr1.close()

fr2 = open(inputDIR2 + p + '.csv', 'r')
line2 = fr2.read().split('\n')
fr2.close()

fw1 = open(outputDIR1 + p + '.csv', 'w')

for i in range(0, len(line1) - 1):
	if i == 0:
		print >>fw1, line1[i] + ',CX_mean,CX_SD,HydPhoPair_mean,rASA_mean,ssH_ratio,ssG_ratio,ssI_ratio,ssE_ratio,ssB_ratio,ssT_ratio,ssC_ratio'
	
	else:
		print >>fw1, line1[i] + ',' + str(dic_CX_mean[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_CX_SD[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_Hydphob[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_rasa[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_H[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_G[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_I[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_E[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_B[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_T[line1[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_C[line1[i].split(',')[0].split('.')[0]])


fw2 = open(outputDIR2 + p + '.csv', 'w')

for i in range(0, len(line2) - 1):
	if i == 0:
		print >>fw2, line2[i] + ',CX_mean,CX_SD,HydPhoPair_mean,rASA_mean,ssH_ratio,ssG_ratio,ssI_ratio,ssE_ratio,ssB_ratio,ssT_ratio,ssC_ratio'
	
	else:
		print >>fw2, line2[i] + ',' + str(dic_CX_mean[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_CX_SD[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_Hydphob[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_rasa[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_H[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_G[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_I[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_E[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_B[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_T[line2[i].split(',')[0].split('.')[0]]) + ',' + str(dic_ss_code_C[line2[i].split(',')[0].split('.')[0]])
