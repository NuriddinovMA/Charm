import os
import numpy as np

def user_coverage_statistic_example(path_to_contact_dump, resolution):
	files = sorted(os.listdir(path_to_contact_dump))
	binCoverageStatistic = {'mean_log_sum':0,'mean_log_mult':0,'mean_log_dif':0}
	H = {}
	M = [[],[]]
	
	for file in files:
		print(file)
		chrm1,chrm2 = file.split('.')[-3:-1]
		with open(path_to_contact_dump + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			bin1,bin2,contact_count = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2])
			try: H[chrm1,bin1] += contact_count
			except KeyError: H[chrm1,bin1] = contact_count
			try: H[chrm2,bin2] += contact_count
			except KeyError: H[chrm2,bin2] = contact_count
	bins = sorted(H.keys())
	k = 0
	for i in range(len(bins)):
		for j in range(i,len(bins)):
			k += 1
			#if k % 1000000 == 0: print('.',end='')
			M[0].append( H[bins[i]] )
			M[1].append( H[bins[j]] )
	M[0] = np.array(M[0])
	M[1] = np.array(M[1])
	binCoverageStatistic['mean_log_sum'] = np.mean(np.log2(M[0]+M[1]))
	binCoverageStatistic['mean_log_mult'] = np.mean(np.log2(M[0]*M[1]))
	binCoverageStatistic['mean_log_dif'] = np.mean(np.abs(np.log2(M[0]/M[1])))
	return binCoverageStatistic

def user_distance_dependent_statistic_example(path_to_contact_dump, resolution):
	files = os.listdir(path_to_contact_dump)
	distanceDependentStatistics = {-1000:[0,0,0],0:[0,0,0],1:[0,0,0]}
	H = {}
	for file in files:
		chrm1,chrm2 = file.split('.')[-3:-1]
		with open(path_to_contact_dump + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			bin1,bin2,contact_count = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2])
			if chrm1 == chrm2 : l = abs(bin2-bin1)
			else: l = -1000
			try: H[l].append(contact_count)
			except KeyError: H[l] = [contact_count]
			#
			#Some operation
			#
	for l in H: distanceDependentStatistics[l] = [np.mean(np.log2(H[l])),np.median(np.log2(H[l]))]
	return distanceDependentStatistics

def user_predict_example(covarege_bin1,covarege_bin2,distance_dependent_statistics, cont_AB, oe_AB):
	predicted_contact_oe = np.log2(covarege_bin1*covarege_bin2)/np.log2(distance_dependent_statistics[2])
	predicted_contact = distance_dependent_statistics[0]*predicted_contact_oe
	return predicted_contact, predicted_contact_oe

def user_pick_contact_example(coverage_bin1_array,coverage_bin2_array,mean_contact_distance_array):
	pick_array = np.log2(coverage_bin1_array*coverage_bin2_array*mean_contact_distance_array)
	return pick_array

def user_randomize_example(contact_array, old_contact_count, new_contact_count):
	loc = contact_array*new_contact_count/old_contact_count
	scale = np.sqrt(contact_array/old_contact_count)
	randomized_contact_array = np.random.logistic(loc*scale)
	return randomized_contact_array