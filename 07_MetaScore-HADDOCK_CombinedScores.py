# Ex of execution: python 07_ML+HADDOCK_CombinedScores.py /gpfs/home/yuj114/group/YJ_DOCK_PP/RFBased_ScoringFunction/Score_Comparison/HADDOCK-Val/108-fold/10_2_108-foldCV/
# in Biostar server

import glob
import sys
import os
import numpy as np

def calculate_rank(vector):
	a = {}
	rank = 1
	
	for num in sorted(vector):
		if num not in a:
			a[num] = rank
		
		rank += 1
	
	return[a[i] for i in vector]

para = sys.argv[1] #'/gpfs/home/yuj114/group/YJ_DOCK_PP/RFBased_ScoringFunction/Score_Comparison/HADDOCK-Val/108-fold/10_2_108-foldCV/'

outputDIR1 = para + 'ML+HADDOCK_CombinedScores/'

os.system('rm -rf ' + outputDIR1)
os.system('mkdir ' + outputDIR1)

sc_list = glob.glob(para + '*.comp')

for sc in sc_list:
	fw = open(outputDIR1 + sc.split('/')[-1], 'w')

	print >>fw, "Decoy" + '\t' + "Class" + '\t' + "i-RMSD" + '\t' + "ML+HADDOCK_Combined" + '\t' + "HADDOCK-Score" + '\t' + "Rank_i-RMSD" + '\t' + "Rank-ML+HADDOCK" + '\t' + "Rank-HADDOCK"

	if sc.split('/')[-1].split('.')[0] in ['1ZM4', '2OT3', '1GXD', '2G77']:
		print "No decoys. Pass!"
		continue
	
	fr = open(sc, 'r')
	lines = fr.read().split('\n')
	fr.close()
	
	rf_list = []
	hd_list = []
	
	for i in range(1, len(lines) - 1):
		rf_list.append(float(lines[i].split('\t')[3]))
		hd_list.append(float(lines[i].split('\t')[4]))
	
	elements = np.array(hd_list)
	
	mean = np.mean(elements, axis = 0)
	sd = np.std(elements, axis = 0)
	
	s = ''
	num = 0
	
	new_hd_list = []

	for hd in hd_list:
		if hd > mean + 2 * sd:
			s = sc.split('/')[-1]
			num += 1
	
		else:
			new_hd_list.append(hd)

	if num != 0:
		print s.split('.')[0] + '\t' + str(num)

	rf_min = min(rf_list)
	rf_max =  max(rf_list)
	hd_min = min(new_hd_list)
	hd_max = max(new_hd_list)

	rank_RFnHD_list = []

	for i in range(1, len(lines) - 1):
		rf_i = float(lines[i].split('\t')[3])
		hd_i = float(lines[i].split('\t')[4])

		if hd_i > mean + 2 * sd:
			hd_i = hd_max
		
		try:
			norm_rf = (rf_i - rf_min) / (rf_max - rf_min)
	
		except ZeroDivisionError:
			norm_rf = 1

		try:
			norm_hd = (hd_i - hd_min) / (hd_max - hd_min)

		except ZeroDivisionError:
			norm_hd = 1

		rank_RFnHD_list.append((norm_rf + norm_hd) / 2)

	ranked_RFnHD_list = calculate_rank(rank_RFnHD_list)

	for i in range(1, len(lines) - 1):
		rf_i = float(lines[i].split('\t')[3])
		hd_i = float(lines[i].split('\t')[4])

		if hd_i > mean + 2 * sd:
			hd_i = hd_max

		try:
			norm_rf = (rf_i - rf_min) / (rf_max - rf_min)

		except ZeroDivisionError:
			norm_rf = 1

		try:
			norm_hd = (hd_i - hd_min) / (hd_max - hd_min)

		except ZeroDivisionError:
			norm_hd = 1

		print >>fw, lines[i].split('\t')[0] + '\t' + lines[i].split('\t')[1] + '\t' + lines[i].split('\t')[2] + '\t' + str((norm_rf + norm_hd) / 2) + '\t' + lines[i].split('\t')[4] + '\t' + lines[i].split('\t')[5] + '\t' + str(ranked_RFnHD_list[i - 1]) + '\t' + lines[i].split('\t')[7]
