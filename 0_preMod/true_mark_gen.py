import numpy as np
import lift_func as lf
from copy import deepcopy

path_c = '/mnt/storage/home/manuriddinov/MyData/Vertebrates/genome/mm10'
chrm ='mm10.chr.sizes'
path_m = '/mnt/scratch/ws/manuriddinov/202203131447EVL27/mlt/mm10'
resolution = 5000
ChrmSzs = lf.ChromSizes(path_c+'/'+chrm,resolution)
ChrmIdx = lf.ChromIndexing(path_c+'/'+chrm)
print ChrmSzs
WT = { ChrmIdx[i+1]:[(ChrmIdx[i+1],j) for j in range(ChrmSzs[i+1])] for i in range(len(ChrmSzs)/2) }
var = open('mm10cap.txt','r')
lines = var.readlines()
var.close()

for line in lines:
	parse = line.split()
	m,c1,p11,p12,a,c2,p2,inv,dups = parse[1:]
	p11,p12 = int(p11)/resolution,int(p12)/resolution
	l = p12 - p11
	try: p2 = int(p2)/resolution
	except ValueError: p2,delete = (p12+l),True
	if inv == 'inverted': inv = True
	else: inv = False
	if dups == 'dups': dups = True
	else: dups = False
	MT = deepcopy(WT)
	print c1,p11,p12,c2,p2
	MT = lf.mutation(MT,(c1,p11),(c2,p2),l,delete=delete,inv=inv,dups=dups)
	lf.markgenerate(MT,c1,path_m,path_c,m,'w')
	if  c2 != c1: lf.markgenerate(MT,c2,path_m,path_c,m,'a')


