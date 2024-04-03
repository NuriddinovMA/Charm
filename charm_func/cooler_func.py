import cooler
import pandas as pd
import subprocess
import logging
logging.basicConfig(level=logging.INFO)

def create_contacts_from_cool(path_to_hic_map, path_to_contact_dump, resolution, chrm1, chrm2, normalization):
	logging.getLogger(__name__).info(
		f"create contacts from cool file {path_to_hic_map} for chomosomes {chrm1} and {chrm2}")
	if normalization != True and normalization != False:
		print("Please use True or False for 'normalization' field")
	cool_ref = cooler.Cooler(path_to_hic_map+"::resolutions/"+str(resolution))
	contacts = cool_ref.matrix(balance=normalization, as_pixels=True).fetch(chrm1, chrm2)
	bins_chrm1 = cool_ref.bins().fetch(chrm1)
	bins_chrm2 = cool_ref.bins().fetch(chrm2)
	BinDict_chrm1 = {index: (index - min(bins_chrm1.index.to_list())) for index in bins_chrm1.index.to_list()}
	BinDict_chrm2 = {index: (index - min(bins_chrm2.index.to_list())) for index in bins_chrm2.index.to_list()}
	contacts["bin1"] = contacts["bin1_id"].apply(lambda x: BinDict_chrm1[x])
	contacts["bin2"] = contacts["bin2_id"].apply(lambda x: BinDict_chrm2[x])
	if normalization == False:
		contacts[["bin1", "bin2", "count"]].to_csv(path_to_contact_dump, index=False, header = False, sep="\t")
	elif normalization == True:
		contacts[["bin1", "bin2", "balanced"]].to_csv(path_to_contact_dump, index=False, header = False, sep="\t")
	logging.getLogger(__name__).info(
			f"succesfully create {path_to_contact_dump}")

def create_cool_from_contacts(path_to_contact_file, path_to_out_mcool, chrom_sizes, hic_resolutions):
	logging.getLogger(__name__).info(
			f"create mcoool file {path_to_out_mcool}")
	contacts_data = pd.read_csv(path_to_contact_file, sep="\t", names=['chrm1','bin1','chrm2','bin2','counts'])
	chroms = sorted(chrom_sizes.keys())
	for resolution in hic_resolutions:
		for i in range(len(chroms)):
			for j in range(i,len(chroms)):
				chrm1,chrm2 = chroms[i],chroms[j]
				logging.getLogger(__name__).info(f"use modeled contacts for {chrm1} and {chrm2}")
				lnc1,lnc2 = chrom_sizes[chrm1]//resolution+1,chrom_sizes[chrm2]//resolution+1
				chrom_pair = contacts_data[(contacts_data['chrm1'] == chrm1) & (contacts_data['chrm2'] == chrm2)].copy()
				chrom_pair['bin1']//= resolution
				chrom_pair['bin2']//= resolution
				chrom_pair.to_csv('chrom_pair1.csv', index=False,sep='\t')
				chrom_pair = chrom_pair.groupby(['chrm1','bin1','chrm2','bin2']).agg({'counts':'sum'})
				chrom_pair = chrom_pair.reset_index()
				chrom_pair.to_csv('chrom_pair2.csv', index=False,sep='\t')
				#
				if chrm1 != chrm2:
					bins_indeces_chr1 = list(range(lnc1))
					bins_chr1 = pd.DataFrame(data = {"chrom":[chrm1]*lnc1, "start" : range(0, lnc1*resolution, resolution ), "end" : range(resolution,(lnc1+1)*resolution, resolution)})
					
					bins_indeces_chr2 = list(range(lnc2))
					bins_chr2 = pd.DataFrame(data = {"chrom":[chrm2]*lnc2, "start" : range(0, lnc2*resolution, resolution ), "end" : range(resolution,(lnc2+1)*resolution, resolution)})
					bins_indeces_chr2 = [x + lnc1 + 1 for x in bins_indeces_chr2]
					bins_indeces = bins_indeces_chr1 + bins_indeces_chr2
					bins = pd.concat([bins_chr1, bins_chr2])
					bins["index"] = bins_indeces
					bins.set_index(["index"])
					chrom_pair['bin2'] = chrom_pair['bin2'].apply(lambda x: x + lnc1 + 1)
				else:
					bins_indeces = list(range(lnc1))
					print(len(bins_indeces))
					bins = pd.DataFrame(data = {"chrom":[chrm1]*lnc1, "start" : range(0, lnc1*resolution, resolution), "end" : range(resolution, (lnc1+1)*resolution, resolution)})
					bins["index"] = bins_indeces
					bins.set_index(["index"])
				
				pixels = pd.DataFrame(data = {"bin1_id" : chrom_pair['bin1'], "bin2_id" : chrom_pair['bin2'], "count" : chrom_pair['counts']})
				path_to_out_cool = path_to_out_mcool + "_" + str(resolution) + "resolution.cool"
				bins.to_csv('bins.csv', index=True,sep='\t')
				pixels.to_csv('pixels.csv', index=True,sep='\t')
				cooler.create_cooler(path_to_out_cool, bins, pixels)
	## cooler zoomify
	cons_command1 = f"cooler zoomify {path_to_out_cool} -r {str(resolution)}N -o {path_to_out_mcool}"
	subprocess.run(cons_command1, shell=True)
	logging.getLogger(__name__).info(
			f"succesfully create {path_to_out_mcool}")

## test these functions
# res = 100000
# chrm1 = "chr5"
# chrm2 = "chr6"
# path_ref_cool = "/storage3/polina/P140_1Kb.mcool"
# path_contacts = "/storage2/polina/ExoC/in_mut.del.1.1.liftCon"
# normalization = False
# create_contacts_from_cool(path_ref_cool, "./"+chrm1 + "_" + chrm2 + "_test_contacts.txt", res, chrm1, chrm2, normalization)
# create_cool_from_contacts(path_contacts, "./test3.mcool", res)