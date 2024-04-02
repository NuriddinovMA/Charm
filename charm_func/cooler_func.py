import cooler
import pandas as pd
import subprocess
import logging
logging.basicConfig(level=logging.INFO)

def create_contacts_from_cool(path_to_hic_map, path_to_hic_dump, resolution, chrom1, chrom2, normalization):
	logging.getLogger(__name__).info(
		f"create contacts from cool file {path_to_mcool_file} for chomosomes {chrom1} and {chrom2}")
	if normalization != True and normalization != False:
		print("Please use True or False for 'normalization' field")
	cool_ref = cooler.Cooler(path_to_mcool_file+"::resolutions/"+str(resolution))
	contacts = cool_ref.matrix(balance=normalization, as_pixels=True).fetch(chrom1, chrom2)
	bins_chrom1 = cool_ref.bins().fetch(chrom1)
	bins_chrom2 = cool_ref.bins().fetch(chrom2)
	BinDict_chrom1 = {index: (index - min(bins_chrom1.index.to_list())) for index in bins_chrom1.index.to_list()}
	BinDict_chrom2 = {index: (index - min(bins_chrom2.index.to_list())) for index in bins_chrom2.index.to_list()}
	contacts["bin1"] = contacts["bin1_id"].apply(lambda x: BinDict_chrom1[x])
	contacts["bin2"] = contacts["bin2_id"].apply(lambda x: BinDict_chrom2[x])
	if normalization == False:
		contacts[["bin1", "bin2", "count"]].to_csv(path_to_hic_dump, index=False, header = False, sep="\t")
	elif normalization == True:
		contacts[["bin1", "bin2", "balanced"]].to_csv(path_to_hic_dump, index=False, header = False, sep="\t")
	logging.getLogger(__name__).info(
			f"succesfully create {path_to_hic_dump}")

def create_cool_from_contacts(path_to_contact_file, path_to_out_mcool, resolution):
	logging.getLogger(__name__).info(
			f"create mcoool file {path_to_out_mcool}")
	contacts_data = pd.read_csv(path_to_contact_file, sep=" ", header=None)
	chrom1, chrom2 = "chr" + str(pd.unique(contacts_data[0])[0]), "chr" + str(pd.unique(contacts_data[2])[0])
	logging.getLogger(__name__).info(
			f"use modeled contacts for {chrom1} and {chrom2}")
	if chrom1 != chrom2:
		bins_indeces_chr1 = list(range(0, max(contacts_data[1])+1))
		bins_chr1 = pd.DataFrame(data = {"chrom":[chrom1]*len(bins_indeces_chr1), "start" : list(range(0, resolution * len(bins_indeces_chr1) , resolution)), "end" : list(range(resolution, resolution*(len(bins_indeces_chr1)+1), resolution))})
		
		bins_indeces_chr2 = list(range(0, max(contacts_data[3])+1))
		bins_chr2 = pd.DataFrame(data = {"chrom":[chrom2]*len(bins_indeces_chr2), "start" : list(range(0, resolution*len(bins_indeces_chr2), resolution)), "end" : list(range(resolution, resolution*(len(bins_indeces_chr2)+1), resolution))})
		bins_indeces_chr2 = [x + max(bins_indeces_chr1) + 1 for x in bins_indeces_chr2]
	
		bins_indeces = bins_indeces_chr1 + bins_indeces_chr2
		bins = pd.concat([bins_chr1, bins_chr2])
		bins["index"] = bins_indeces
		bins.set_index(["index"])
		contacts_data[3] = contacts_data[3].apply(lambda x: x + max(bins_indeces_chr1) + 1)

		
	else:
		bins_indeces = list(range(0, max(max(contacts_data[1]), max(contacts_data[3]))+1))
		bins = pd.DataFrame(data = {"chrom":[chrom1]*len(bins_indeces), "start" : range(0, resolution*len(bins_indeces), resolution), "end" : range(resolution, resolution*(len(bins_indeces)+1), resolution)})
		bins["index"] = bins_indeces
		bins.set_index(["index"])
	
	pixels = pd.DataFrame(data = {"bin1_id" : contacts_data[1], "bin2_id" : contacts_data[3], "count" : contacts_data[4]})
	path_to_out_cool = path_to_out_mcool + "_" + str(res) + "resolution.cool"
	cooler.create_cooler(path_to_out_cool, bins, pixels)
	## cooler zoomify
	cons_command1 = f"cooler zoomify {path_to_out_cool} -r {str(resolution)}N -o {path_to_out_mcool}"
	subprocess.run(cons_command1, shell=True)
	logging.getLogger(__name__).info(
			f"succesfully create {path_to_out_mcool}")

## test these functions
res = 100000
chrom1 = "chr5"
chrom2 = "chr6"
path_ref_cool = "/storage3/polina/P140_1Kb.mcool"
path_contacts = "/storage2/polina/ExoC/in_mut.del.1.1.liftCon"
normalization = False
create_contacts_from_cool(path_ref_cool, "./"+chrom1 + "_" + chrom2 + "_test_contacts.txt", res, chrom1, chrom2, normalization)
create_cool_from_contacts(path_contacts, "./test3.mcool", res)