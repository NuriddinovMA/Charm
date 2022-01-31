import numpy as np
import lift_func as lf
from copy import deepcopy

def mutation(G,b1,b2,ln,**kwargs):
	try: inv = kwargs['inv']
	except KeyError: inv = False
	try: delete = kwargs['delete']
	except KeyError: delete = False
	try: insert = kwargs['insert']
	except KeyError: insert = False
	reg = G[b1[0]][b1[1]:(b1[1]+ln)]
	if inv: reg.reverse()
	else: pass
	if insert: G[b2[0]] = G[b2[0]][:b2[1]]+reg+G[b2[0]][b2[1]:]
	else: pass
	if delete: 
		if b1[0] == b2[0] and b2[1] < b1[1]: del G[b1[0]][(b1[1]+ln):(b1[1]+2*ln)]
		else: del G[b1[0]][b1[1]:(b1[1]+ln)]
	else: pass
	del reg
	return G

def markgenerate(M,c1,path_mark,path_chrm,name,of):
	f1 = open('%s/cnt.mut.%s.mark' % (path_mark, name), of)
	f2 = open('%s/mut.%s.cnt.mark' % (path_mark, name), of)
	k = 0
	for i in range(len(M[c1])):
		k += 1
		print >> f1, M[c1][i][0], M[c1][i][1]*5000+2400, M[c1][i][1]*5000+2600, c1, i*5000+2400, i*5000+2600, 1, k
		print >> f2, c1, i*5000+2400, i*5000+2600, M[c1][i][0], M[c1][i][1]*5000+2400, M[c1][i][1]*5000+2600, 1, k
	f1.close()
	f2.close()
	with open('%s/mut.%s.chr.sizes' % (path_chrm, name), of) as f: print >> f, c1,len(M[c1])*5000

path_c = '/mnt/storage/home/manuriddinov/MyData/Vertebrates/genome/mdl'
chrm ='hg19.chr.sizes'
path_m = '/mnt/scratch/ws/manuriddinov/202203131447EVL27/mlt/mdl'
ChrmSzs = lf.ChromSizes(path_c+'/'+chrm,5000)
ChrmIdx = lf.ChromIndexing(path_c+'/'+chrm)
print ChrmSzs
WT = { ChrmIdx[i+1]:[(ChrmIdx[i+1],j) for j in range(ChrmSzs[i+1])] for i in range(24) }
length = [5000,10000,20000,50000,200000,500000,2000000]
tlength = [5000,100000,2000000,10000000,40000000]
dir = ['intact','inverted']
m =10570
var = open('mutvar2.txt','a')
var.close()
for s in range(29):
	for lng in length:
		l = lng/5000
		d = s % 4
		c1,ci = np.random.choice(24,2,replace=False)
		c1,ci = ChrmIdx[c1+1],ChrmIdx[ci+1]
		inv = np.random.randint(2,size=8)
		tdir = np.random.randint(2,size=5)
		sz1,sz2 = ChrmSzs[c1], ChrmSzs[ci]
		p3,p4,p5,p6,p7 = np.random.randint(sz1,size=5)
		
		tlength = [1,20,400,2000,8000,sz1]
		for tlow in range(5):
			tl = np.random.randint(tlength[tlow],high=tlength[tlow+1])
			if tdir[tlow]: 
				p1 = np.random.randint(sz1-tl)
				p2 = p1+tl
			else:
				p1 = np.random.randint(tl,high=sz1)
				p2 = p1-tl
			MT = deepcopy(WT)
			MT = mutation(MT,(c1,p1),(c1,p2),l,inv=inv[tlow],insert=True,delete=True)
			m+=1
			markgenerate(MT,c1,path_m,path_c,m,'w')
			with open('mutvar2.txt','a') as var: 
				print >> var, 'variant', m, c1,p1*5000,l*5000,'->',c1,p2*5000,dir[inv[tlow]],'trans'
		#interchromosome
		pi=np.random.randint(sz2)
		MT = deepcopy(WT)
		MT = mutation(MT,(c1,p3),(ci,pi),l,inv=inv[5],insert=True,delete=True)
		m+=1
		markgenerate(MT,c1,path_m,path_c,m,'w')
		markgenerate(MT,ci,path_m,path_c,m,'a')
		with open('mutvar2.txt','a') as var: 
			print >> var, 'variant', m, c1,p3*5000,l*5000,'->',ci,pi*5000,dir[inv[5]],'trans'
		#duplication
		MT = deepcopy(WT)
		MT = mutation(MT,(c1,p4),(c1,p4+l),l,inv=inv[6],insert=True,delete=False)
		m+=1
		markgenerate(MT,c1,path_m,path_c,m,'w')
		with open('mutvar2.txt','a') as var: 
			print >> var, 'variant', m, c1,p4*5000,l*5000,'->',c1,(p4+l)*5000,dir[inv[6]],'dups'
		#invers
		MT = deepcopy(WT)
		MT = mutation(MT,(c1,p5),(c1,p5+l),l,inv=True,insert=True,delete=True)
		m+=1
		markgenerate(MT,c1,path_m,path_c,m,'w')
		with open('mutvar2.txt','a') as var: 
			print >> var, 'variant', m, c1,p5*5000,l*5000,'->',c1,(p5+l)*5000,'inverted','trans'
		#delete
		MT = deepcopy(WT)
		MT = mutation(MT,(c1,p6),(c1,p6),l,inv=inv[7],insert=False,delete=True)
		m+=1
		markgenerate(MT,c1,path_m,path_c,m,'w')
		with open('mutvar2.txt','a') as var: 
			print >> var, 'variant', m, c1,p6*5000,l*5000,'->',c1,'delete',dir[inv[7]],'trans'
		#control
		MT = deepcopy(WT)
		MT = mutation(MT,(c1,p7),(c1,p7+l),l,inv=False,insert=True,delete=True)
		m+=1
		markgenerate(MT,c1,path_m,path_c,m,'w')
		with open('mutvar2.txt','a') as var: 
			print >> var, 'variant', m, c1,p7*5000,l*5000,'->',c1,(p7+l)*5000,'intact','trans'

	print s,m
