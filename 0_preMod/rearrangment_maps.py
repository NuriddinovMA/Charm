import numpy as np
import os
import argparse
import lift_func as lf
from copy import deepcopy

parser = argparse.ArgumentParser(description='Parameters for generation of rearrangement maps')
parser.add_argument('vars', metavar='rearrangements', help='path to file with the descriptions of rearrangements')
parser.add_argument('-c', dest='chrmsizes', metavar='chromosomes', help='path to the chromosome size file')
parser.add_argument('-o', dest='outdir',metavar='output', help='path to the output directory')
parser.add_argument('-r', dest='resolution', metavar='resolution', type=int, help='resolution of generated maps')
args = parser.parse_args()
ChrmSzs = lf.ChromSizes(args.chrmsizes,args.resolution)
mutChrmSzs = {}
ChrmIdx = lf.ChromIndexing(args.chrmsizes)

WT = { ChrmIdx[i+1]:[(ChrmIdx[i+1],j) for j in range(ChrmSzs[ChrmIdx[i+1]])] for i in range(len(ChrmSzs)) }

with open(args.vars, 'r') as f: lines = f.readlines()
for line in lines[100:]:
	print( line )
	cnt,mut,c1,p11,p12,a,c2,p2,cnv1,cnv2 = line.split()[:10]
	cnv1,cnv2 = int(cnv1),int(cnv2)
	
	if cnv2 < 0: invert = True
	else: invert = False
	
	p11 = int(p11)//args.resolution
	if a in ['->','!>']: 
		print('start')
		chrms = [c1,c2]
		try: p12 = int(p12)//args.resolution
		except ValueError:
			try: p12 = ChrmSzs[c1]
			except KeyError: p12 = 0
		try: p2 = int(p2)//args.resolution
		except ValueError:
			try: p2 = ChrmSzs[c2]
			except KeyError: p2 = 0
	if a in ['>>','>!']:
		print('+++')
		chrms.extend([c1,c2])
		try: p12 = int(p12)//args.resolution
		except ValueError: 
			try: p12 = ChrmSzs[c1]
			except KeyError: p12 = mutChrmSzs[c1]
		try: p2 = int(p2)//args.resolution
		except ValueError:
			try: p2 = ChrmSzs[c2]
			except KeyError: p2 =  mutChrmSzs[c2]
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
		MT = lf.mutation(MT,(c1t,p11t,p12t),(c2t,p2t),cnv1=cnv1,cnv2=cnv2)
	else: MT = lf.mutation(MT,(c1,p11,p12),(c2,p2),cnv1=cnv1,cnv2=cnv2)
	for key in MT: mutChrmSzs[key] = len(MT[key])
	if a in ['>>','!>']: markHash = lf.hashgenerate(MT)
	if a in ['!>','->']:
		try: os.remove('%s/%s.%s.%i.mark' % (args.outdir, cnt, mut, args.resolution))
		except OSError: pass 
		try: os.remove('%s/%s.%s.%i.mark' % (args.outdir, mut, cnt, args.resolution))
		except OSError: pass 
	if a in ['->','>!']:
		chrms = list(set(chrms))
		chrms.sort()
		for c1 in range(len(chrms)): lf.markgenerate(MT,chrms[c1],args.outdir,args.resolution,(cnt,mut),c1)
		
	if a in ['->','>!']:
		with open('%s/%s.%s.chr.sizes' % (args.outdir, cnt,mut), 'w') as f: 
			Keys = sorted(MT)
			for key in Keys: print( key, (mutChrmSzs[key]+1)*args.resolution, file=f)
		print('end')


