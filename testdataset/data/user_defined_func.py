import os
import numpy as np

#
# This file contains workable examples of custom functions for simulation rearrangements
#


def user_coverage_statistic_example(path_to_contact_dump, resolution):
	#####################################################################
	# There is the function to calculate full genomic statistic of coverage of contacting bins.
	# As input parameters this function MUST USE:
	# (!!!) "path_to_contact_dump" - path to folder contains files with contact dumped from hic-matrices
	# (!!!) "resolution" - resolution of hic-matrices
	# Both these parameters are defined automatically.
	# DONT CHANGE!
	######################################################################
	
	binCoverageStatistic = {'mean_log_sum':0}
	
	##########################################
	# binCoverageStatistic - the hash containing calculated statisitcs
	# if bin_i and bin_j have coverage_i and coverage_j
	# 'mean_log_sum' = mean(log2(coverage_i + coverage_j)) - the average of logarithm of the sum of contacting bin coverage
	##########################################
	
	############################################
	# This block is adapted to read coverage from the intermediate files with contact statistic
	############################################
	coverage_Hash = {} # The hash containing coverage of bins {bin_id:bin_coverage}
	files = sorted(os.listdir(path_to_contact_dump)) # path to folder with intermediate files
	for file in files:
		chrm1,chrm2 = file.split('.')[-3:-1] #defininig chromosome names from the file name
		with open(path_to_contact_dump + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			bin1,bin2,contact_count = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2]) #defining bin numbers from genome coordinate and resolution
			try: coverage_Hash[chrm1,bin1] += contact_count
			except KeyError: coverage_Hash[chrm1,bin1] = contact_count
			try: coverage_Hash[chrm2,bin2] += contact_count
			except KeyError: coverage_Hash[chrm2,bin2] = contact_count
	
	########################
	#Next block calculate statistics, you can change it as wish
	########################
	bins = sorted(coverage_Hash.keys())
	M = []
	for i in range(len(bins)):
		for j in range(i,len(bins)): M.append( coverage_Hash[bins[i]] + coverage_Hash[bins[j]] )
	M = np.array(M)
	binCoverageStatistic['mean_log_sum'] = np.mean(np.log2(M))
	return binCoverageStatistic

def user_distance_dependent_statistic_example(path_to_contact_dump, resolution):
	#####################################################################
	# There is the function to calculate distance dependend statics.
	# As input parameters this function MUST USE:
	# (!!!) "path_to_contact_dump" - path to folder contains files with contact dumped from hic-matrices
	# (!!!) "resolution" - resolution of hic-matrices
	# Both these parameters are defined automatically. DONT CHANGE IT!
	######################################################################
	
	distanceDependentStatistics = {-1000:[0],0:[0],1:[0]}
	##########################################
	# distanceDependentStatistics - the hash containing calculated distance dependent statisitcs {genome_distance_in_bins:[statistics_list]}
	# key "-1000" is reserved for INTERchromosome contact. DONT CHANGE!
	##########################################

	distance_Hash = {} # The hash containing list of contacts for every genomic distance
	files = os.listdir(path_to_contact_dump)
	for file in files:
		chrm1,chrm2 = file.split('.')[-3:-1] #defininig chromosome names from the file name
		with open(path_to_contact_dump + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			bin1,bin2,contact_count = int(parse[0])//resolution,int(parse[1])//resolution,float(parse[2]) #defining bin numbers from genome coordinate and resolution
			if chrm1 == chrm2 : distance_in_bins = abs(bin2-bin1) #defining genomic distance for intrAchromosome contacts
			else: distance_in_bins = -1000 #defining genomic distance for intErchromosome contacts
			try: distance_Hash[distance_in_bins].append(contact_count)
			except KeyError:distance_Hash[distance_in_bins] = [contact_count]
	
	########################
	#Next block calculate statistics, you can change it as wish
	########################
	for distance_in_bins in distance_Hash: distanceDependentStatistics[distance_in_bins] = [np.mean(np.log2(distance_Hash[distance_in_bins]))]
	return distanceDependentStatistics

def user_predict_example(covarege_bin1, covarege_bin2, cont_AB, oe_AB, distance, distance_dependent_statistics, total_contact_statistics):
	##################
	#There is the function to predict contact count instead 0 values ("predict_null_contacts" parametr)
	#covarege_bin1 and covarege_bin1 - the coverages of contacting bins
	#distance_dependent_statistics - the hash containg contact statistics for given distance
	#cont_AB - contact count between pseudo-compartment of this bins
	#oe_AB - observed/expected value for pseudo-compartment of this bins
	# These parameters are defined automatically. DONT CHANGE IT!
	###################
	predicted_contact_oe = np.log2(covarege_bin1*covarege_bin2)/np.log2(distance_dependent_statistics['mean_sum_coverage']) #predicted observed/expected values
	predicted_contact = distance_dependent_statistics['mean']*predicted_contact_oe  #predicted contact count
	return predicted_contact, predicted_contact_oe

def user_pick_contact_example(coverage_bin1_array, coverage_bin2_array, mean_contact_distance_array, **kwargs):
	##################
	#There is the function to distribute contact from low resolution to high. 
	#This function calculate "probability coefficient" to pick given bin pair (in high resolution), in order to put contacts from low resolution
	#covarege_bin1_arrau and covarege_bin1_array - the numpy arrays of coverages of given contacting bins
	#mean_contact_distance_array - the numpy array of mean contact values for distance between  given contacting bins
	# These parameters are defined automatically. DONT CHANGE IT!
	###################
	pick_array = np.log2(coverage_bin1_array*coverage_bin2_array*mean_contact_distance_array)
	# numpy array of "probability" to pick given bin pair
	return pick_array

def user_randomize_example(contact_array, old_contact_count, new_contact_count):
	##################
	#There is the function to randomize given contact values
	#contact_array - numpy array of contacts
	#old_contact_count - total contact count in reference hic-map (see "all_contacts" values in .binCov.stat file )
	#new_contact_count - total contact count in simulated hic-map (see "contact_count" parametr in ini-file)
	# These parameters are defined automatically. DONT CHANGE IT!
	###################
	loc = contact_array*new_contact_count/old_contact_count
	scale = np.sqrt(contact_array/old_contact_count)
	randomized_contact_array = np.random.logistic(loc*scale) # numpy array of randomized contacts
	return randomized_contact_array
