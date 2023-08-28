import numpy as np
import os
from charm_func import sv2_func as svf
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
		cnt,mut,c1,p11,p12,a,c2,cnv1,cnv2 = line.split()[:9]
		if rname == False or rname == mut or stand_alone:
			gf.printlog('\t\tstructural variation %s' % line.strip(),log_file)
			
			if a in ['->','!>']: 
				pointviews = ''
				chrm_from, chrm_to, add_from, add_to = set([]),set([]),set([]),set([])
				cc_from,cc_to,MT,mutChrmSzs = {},{},{},{}
			pointviews += '%s %s %s\n' % (c1,p11,p12)
			
			cnv1,cnv2 = int(cnv1),int(cnv2)
			
			try: cc_from[c2].add(c1)
			except KeyError: cc_from[c2] = set([c1,])
			
			try: cc_to[c1].add(c2)
			except KeyError: cc_to[c1] = set([c2,])

			if (c2 not in ChrmSzs) or ((cnv1+abs(cnv2)) > 1):
				for chrm in ChrmSzs: 
					add_from.add((c1,chrm))
					add_to.add((c2,chrm))
			
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
			
			try: mutChrmSzs[c2] += (p12-p11)*abs(cnv2)
			except KeyError: mutChrmSzs[c2] = (p12-p11)*abs(cnv2)
			
			if a in ['->','>!']:
				svf.markgenerate(MT,resolution,(cnt,mut),outdir)
				
				print(cc_from)
				for c2 in cc_from:
					print('1',cc_from[c2])
					cc_from[c2] = list(cc_from[c2])
					cc_from[c2].sort()
					print('1',cc_from[c2])
					for i in range(len(cc_from[c2])):
						for j in range(i,len(cc_from[c2])): 
							chrm_from.add((cc_from[c2][i],cc_from[c2][j]))
							if i != j: 
								print('2',(cc_from[c2][j],cc_from[c2][i]))
								chrm_from -= set([(cc_from[c2][j],cc_from[c2][i])])
							add_from -= set([(cc_from[c2][i],cc_from[c2][j]),(cc_from[c2][j],cc_from[c2][i])])
				
				print(cc_to)
				for c1 in cc_to:
					cc_to[c1] = list(cc_to[c1])
					print('3',cc_to[c1])
					for i in range(len(cc_to[c1])):
						for j in range(i,len(cc_to[c1])):
							chrm_to.add((cc_to[c1][i],cc_to[c1][j]))
							if i != j: chrm_to -= set([(cc_to[c1][j],cc_to[c1][i])])
							add_to -= set([(cc_to[c1][i],cc_to[c1][j]),(cc_to[c1][j],cc_to[c1][i])])
				
				for c2 in ChrmSzs:
					try: mutChrmSzs[c2]
					except KeyError: mutChrmSzs[c2] = ChrmSzs[c2]
				
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


