import numpy as np
from charm_func import global_func as lf

path_c = '/mnt/storage/home/manuriddinov/MyData/Vertebrates/genome/mdl/hg19.somchr.sizes'
path_m = '/mnt/storage/home/manuriddinov/scratch/mlt/mdl'
path_cov = '/mnt/storage/home/manuriddinov/scratch/ArChER/ExoC/pre/ecSum/ecSum.5000.binCov'
R = 5000

ChrmSzs = lf.ChromSizes(path_c,R)
ChrmIdx = lf.ChromIndexing(path_c)
dir = ['intact','inverted']

covH,covHT = {},{}
f = open(path_cov,'r')
lines = f.readlines()
f.close()
for line in lines[4:]:
	parse = line.split()
	try:
		chrm,b1,cov = parse[0],int(parse[1]),int(float(parse[2]))
		covHT[chrm,b1] = cov
		if cov != 0: 
			try: covH[chrm].append((b1,cov))
			except KeyError: covH[chrm] = [(b1,cov),]
	except IndexError: pass
	except ValueError: pass
del lines
for chrm in covH: 
	covH[chrm] = np.array(covH[chrm])
	covH[chrm] = covH[chrm][covH[chrm][:,1]>=1000]
inv = (2*np.random.randint(2,size=200)-1).tolist()
c = np.random.choice(22,200,replace=True).tolist()
brk = []
for i in c:
	k = str(i+1)
	brk.append(covH[k][np.random.choice(len(covH[k])),1])
with open('break.txt','w') as f:
	for i in range(200): f.write('%i %i %i\n' % (c[i]+1,brk[i]*5000,inv[i]))

