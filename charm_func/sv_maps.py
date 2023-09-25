import numpy as np
import os
from charm_func import sv_func as svf
from charm_func import global_func as gf
from copy import deepcopy

def generate_SV_map(chrom_sizes, resolution, rearrangement_list, work_dir, rname, stand_alone,log_file):
	
	outdir = work_dir + '/rear'
	try: os.makedirs(outdir)
	except OSError: pass
	resolution = int(resolution)
	shift = resolution//2
	ChrmSzs = gf.ChromSizes(chrom_sizes,resolution)
	ChrmIdx = gf.ChromIndexing(chrom_sizes)
	
	with open(rearrangement_list, 'r') as f: lines = f.readlines()
	for line in lines:
		cnt,mut,c1,p11,p12,a,c2,cnv2 = line.split()[:8]
		if (rname == False or rname == mut or stand_alone) and (line[0] != '#'):
			gf.printlog('\t\tstructural variation %s' % line.strip(),log_file)
			
			if a in ['->','!>']: 
				pointviews = ''
				chrm_from, chrm_to, cnv, add_from, add_to = set([]),set([]),set([]),set([]),set([])
				cc_from,cc_to,MT,mutChrmSzs = {},{},{},{}
				mutChrmIdx = {i:ChrmIdx[i] for i in ChrmIdx}
				mutChrmSzs = {c2:ChrmSzs[c2] for c2 in ChrmSzs}
				
				
			idxStart = len(mutChrmIdx) + 1
			try: mutChrmIdx[c2]
			except KeyError:
				mutChrmIdx[c2] = idxStart
				mutChrmIdx[idxStart] = c2

			pointviews += '%s %s %s\n' % (c1,p11,p12)
			cnv2 = int(cnv2)
			
			try: cc_from[c1].add(c2)
			except KeyError: cc_from[c1] = set([c2,])
			try: cc_to[c2].add(c1)
			except KeyError: cc_to[c2] = set([c1,])
			
			p11 = (int(p11)+shift)//resolution
			try: p12 = (int(p12)+shift)//resolution
			except ValueError: p12 = ChrmSzs[c1]
			
			if cnv2 > 0: 
				try: MT[c2] += [(c1,i) for i in range(p11,p12)]*cnv2
				except KeyError: MT[c2] = [(c1,i) for i in range(p11,p12)]*cnv2
			elif cnv2 < 0: 
				try: MT[c2] += [(c1,i) for i in range(p12-1,p11-1,-1)]*abs(cnv2)
				except KeyError: MT[c2] = [(c1,i) for i in range(p12-1,p11-1,-1)]*abs(cnv2)
			else: pass
			
			try: del mutChrmSzs[c1]
			except KeyError: pass
			try: mutChrmSzs[c2] += (p12-p11)*abs(cnv2)
			except KeyError: mutChrmSzs[c2] = (p12-p11)*abs(cnv2)
			
			if a in ['->','>!']:
				svf.markgenerate(MT,resolution,(cnt,mut),outdir)
				
				chrm_from, chrm_to = set([]),set([])
				for c1 in cc_from:
					lenc = len(cc_from[c1])
					cc_from[c1] = list(cc_from[c1])
					for i in range(lenc):
						for j in range(i,lenc):
							c2i,c2j = cc_from[c1][i],cc_from[c1][j]
							if mutChrmIdx[c2i] <= mutChrmIdx[c2j]: chrm_to.add((c2i,c2j))
							else: chrm_to.add((c2j,c2i))
				print('to',chrm_to)
				for c2 in cc_to:
					lenc = len(cc_to[c2])
					cc_to[c2] = list(cc_to[c2])
					for i in range(lenc):
						for j in range(i,lenc): 
							c1i,c1j = cc_to[c2][i],cc_to[c2][j]
							if ChrmIdx[c1i] <= ChrmIdx[c1j]: chrm_from.add((c1i,c1j))
							else: chrm_from.add((c1j,c1i))
				print('from',chrm_from)
				
				full_bins_from,full_bins_to = [],[]
				for c1 in cc_from: full_bins_from = [(c1,i) for i in range( ChrmSzs[c1] )]
				for c2 in cc_to: full_bins_to.extend(MT[c2])
				unique = np.unique(full_bins_to,return_counts=True,axis=0)
				cnvBins = unique[0][unique[1]>1].tolist()+list(set(full_bins_from)-set(full_bins_to))
				
				if len(cnvBins) > 0:
					print('!!!CNV!!!')
					for i in cnvBins: cnv.add(i[0])
					for c1 in cnv:
						for chrm in ChrmSzs: 
							if ChrmIdx[c1] <= ChrmIdx[chrm]: add_from.add((c1,chrm))
							else: add_from.add((chrm,c1))
						for c2 in cc_from[c1]: 
							for chrm in mutChrmSzs: 
								if mutChrmIdx[c2] <= mutChrmIdx[chrm]: add_to.add((c2,chrm))
								else: add_to.add((chrm,c2))
				else: print('!!!NO-CNV!!!')
				print('from',add_from)
				print('to',add_to)
				
				with open('%s/%s.%s.chr.sizes' % (outdir, mut, cnt), 'w') as f: 
					Keys = sorted(mutChrmSzs)
					for key in Keys: f.write( '%s %i\n' % (key, (mutChrmSzs[key]+1)*resolution) )
					map_SV_from_ref ='%s/%s.%s.%i.mark' % (outdir, cnt, mut, resolution)
					map_SV_to_ref ='%s/%s.%s.%i.mark' % (outdir, mut, cnt, resolution)
					chrom_sizes_SV = '%s/%s.%s.chr.sizes' % (outdir, mut,cnt)
				
				chrm_from,add_from,chrm_to,add_to = list(chrm_from),list(add_from),list(chrm_to),list(add_to)
				chrm_from.sort(),add_from.sort(),chrm_to.sort(),add_to.sort()
				gf.printlog('\t\tchrm_from %s' % chrm_from,log_file)
				gf.printlog('\t\tadd_from %s' % add_from,log_file)
				gf.printlog('\t\tchrm_to %s' % chrm_to,log_file)
				gf.printlog('\t\tadd_to %s' % add_to,log_file)
				
				cf,af,ct,at = '','','',''
				
				for i in chrm_from: cf = '%s,%s;%s' % (i[0],i[1],cf)
				for i in add_from: af = '%s,%s;%s' % (i[0],i[1],af)
				for i in chrm_to: ct = '%s,%s;%s' % (i[0],i[1],ct)
				for i in add_to: at = '%s,%s;%s' % (i[0],i[1],at)
				if stand_alone == False: return (cf[:-1],ct[:-1]),(af[:-1],at[:-1]),map_SV_from_ref,pointviews,map_SV_to_ref,chrom_sizes_SV


