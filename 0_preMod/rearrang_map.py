import numpy as np
import lift_func as lf
from copy import deepcopy

path_c = '/mnt/storage/home/manuriddinov/MyData/Vertebrates/genome/mdl'
chrm ='hg19.chr.sizes'
path_m = '/mnt/scratch/ws/manuriddinov/202207111350EVL29/mlt/Mdl'
resolution = 10000
ChrmSzs = lf.ChromSizes(path_c+'/'+chrm,resolution)
ChrmIdx = lf.ChromIndexing(path_c+'/'+chrm)
# print ChrmSzs
WT = { ChrmIdx[i+1]:[(ChrmIdx[i+1],j) for j in range(ChrmSzs[i+1])] for i in range(len(ChrmSzs)/2) }
var = open('patient.txt','r')
lines = var.readlines()
var.close()

m1 = ''
for i in range(len(lines)):
	print lines[i]
	
	parse = lines[i].split()
	m,c1,p11,p12,a,c2,p2,invert,type = parse[1:]
	if invert == 'inverted': invert = True
	else: invert = False
	if type == 'duplication': insert,delete = True,False
	elif type == 'deletion': insert,delete = False,True
	else: insert,delete = True,True
	
	p11 = int(p11)/resolution
	try: p12 = int(p12)/resolution
	except ValueError: p12 = ChrmSzs[c1]
	try: p2 = int(p2)/resolution
	except ValueError: p2 = ChrmSzs[c2]
	
	print m,c1,p11,p12,c2,p2
	if m[-1] in ['+','*','^']: m0 = m[:-1]
	else:  m0 = m
	if m[-1] in ['+','*']: pass
	else:
		MT = { ChrmIdx[i+1]:[(ChrmIdx[i+1],j) for j in range(ChrmSzs[i+1])] for i in range(len(ChrmSzs)/2) }
	if m[-1] in ['+','*']:
		print c2,p2
		c2t,p2t = markHash[c2,p2-1]
		p2t+=1
		print c2t,p2t
		ln = p12-p11
		print c1,p11,p12,ln
		c1t,p11t = markHash[c1,p11]
		p12t = p11t+ln
		print c1t,p11t,p12t
		MT = lf.mutation(MT,(c1t,p11t,p12t),(c2t,p2t),invert=invert,insert=insert,delete=delete)
	else: MT = lf.mutation(MT,(c1,p11,p12),(c2,p2),invert=invert,insert=insert,delete=delete)
	if m[-1] in ['+','^']:
		markHash = lf.hashgenerate(MT)
	else:
		print m
		lf.markgenerate(MT,c1,path_m,path_c,resolution,m0,'w')
		if c2 != c1: lf.markgenerate(MT,c2,path_m,path_c,resolution,m0,'a')
	print 'end'


