import numpy as np
import os
import sv_func as svf
import global_func as gf
from copy import deepcopy

def generate_SV_map(chrom_sizes, resolution, rearrangement_list, work_dir, stand_alone):
	
	outdir = work_dir + '/rear'
	try: os.makedirs(outdir)
	except OSError: pass
	resolution = int(resolution)
	ChrmSzs = gf.ChromSizes(chrom_sizes,resolution)
	ChrmIdx = gf.ChromIndexing(chrom_sizes)
	mutChrmSzs = {}
	WT = { ChrmIdx[i+1]:[(ChrmIdx[i+1],j) for j in range(ChrmSzs[ChrmIdx[i+1]])] for i in range(len(ChrmSzs)) }
	cc_f, cc_t, pointviews = set([]),set([]),''
	with open(rearrangement_list, 'r') as f: lines = f.readlines()
	for line in lines:
		print( line )
		cnt,mut,c1,p11,p12,a,c2,p2,cnv1,cnv2 = line.split()[:10]
		cnv1,cnv2 = int(cnv1),int(cnv2)
		cc_f.add(c1)
		cc_t.add(c2)
		pointviews += '%s %s %s\n' % (c1,p11,p12)
		if cnv2 < 0: invert = True
		else: invert = False
		
		p11 = int(p11)//resolution
		if a in ['->','!>']: 
			print('start')
			chrms = [c1,c2]
			try: p12 = int(p12)//resolution
			except ValueError:
				try: p12 = ChrmSzs[c1]
				except KeyError: p12 = 0
			try: p2 = int(p2)//resolution
			except ValueError:
				try: p2 = ChrmSzs[c2]
				except KeyError: p2 = 0
		if a in ['>>','>!']:
			print('+++')
			chrms.extend([c1,c2])
			try: p12 = int(p12)//resolution
			except ValueError: 
				try: p12 = ChrmSzs[c1]
				except KeyError: p12 = mutChrmSzs[c1]
			try: p2 = int(p2)//resolution
			except ValueError:
				try: p2 = ChrmSzs[c2]
				except KeyError:
					try: p2 = mutChrmSzs[c2]
					except KeyError: p2 = 0
		else:
			print( a,'generate MT' )
			MT = { ChrmIdx[i+1]:[(ChrmIdx[i+1],j) for j in range(ChrmSzs[ChrmIdx[i+1]])] for i in range(len(ChrmSzs)) }
		if a in ['>>','>!']:
			try: 
				if (c2,p2) == markHash[c2,p2]: c2t,p2t = c2,p2
				else:
					c2t,p2t = markHash[c2,p2-1]
					p2t+=1
			except KeyError: c2t,p2t = c2,p2
			print( c2t,p2t )
			ln = p12-p11
			print( c1,p11,p12,ln )
			if cnv1 and cnv2: c1t,p11t = c1,p11
			else: c1t,p11t = markHash[c1,p11]
			p12t = p11t+ln
			print( c1t,p11t,p12t )
			print( cnv1, cnv2 )
			MT = svf.mutation(MT,(c1t,p11t,p12t),(c2t,p2t),cnv1=cnv1,cnv2=cnv2)
		else: MT = svf.mutation(MT,(c1,p11,p12),(c2,p2),cnv1=cnv1,cnv2=cnv2)
		for key in MT: mutChrmSzs[key] = len(MT[key])
		if a in ['>>','!>']: markHash = svf.hashgenerate(MT)
		if a in ['!>','->']:
			try: os.remove('%s/%s.%s.%i.mark' % (outdir, cnt, mut, resolution))
			except OSError: pass 
			try: os.remove('%s/%s.%s.%i.mark' % (outdir, mut, cnt, resolution))
			except OSError: pass 
		if a in ['->','>!']:
			chrms = list(set(chrms))
			chrms.sort()
			for c1 in range(len(chrms)): svf.markgenerate(MT,chrms[c1],outdir,resolution,(cnt,mut),c1)
		if a in ['->','>!']:
			with open('%s/%s.%s.chr.sizes' % (outdir, mut, cnt), 'w') as f: 
				Keys = sorted(MT)
				for key in Keys: print( key, (mutChrmSzs[key]+1)*resolution, file=f)
				map_SV_from_ref ='%s/%s.%s.%i.mark' % (outdir, cnt, mut, resolution)
				map_SV_to_ref ='%s/%s.%s.%i.mark' % (outdir, mut, cnt, resolution)
				chrom_sizes_SV = '%s/%s.%s.chr.sizes' % (outdir, cnt,mut)
			print('end')
			
			if stand_alone == False:
				chosen_chroms_from, chosen_chroms_to = '',''
				for c in cc_f: chosen_chroms_from += '%s,' % c
				for c in cc_t: chosen_chroms_to += '%s,' % c
				return map_SV_from_ref,chosen_chroms_from[:-1],pointviews,map_SV_to_ref,chosen_chroms_to[:-1],chrom_sizes_SV
				exit()


