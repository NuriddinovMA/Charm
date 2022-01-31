import lift_func2 as lf

res=5000
path_c = '/mnt/storage/home/manuriddinov/MyData/Vertebrates/genome/mdl/hg19.chr.sizes'
path_m = '/mnt/scratch/ws/manuriddinov/202111211736EVL25/mlt/mdl'
ChrmSzs = lf.ChromSizes(path_c,res)
ChrmIdx = lf.ChromIndexing(path_c)
print ChrmSzs
t1 = 'chr5', 63400000/res
t2 = 'chr10', 56340000/res
print t1,t2
H = []
for chrm in ChrmSzs:
	if ChrmIdx[chrm] in ['chr5','chr10']:
		H += [(t1[0],b,t1[0],b) for b in range(t1[1])]
		H += [( t2[0], b, t1[0], (t1[1]+b-t2[1]) ) for b in range(t2[1],(ChrmSzs[chrm]))]
		H += [(t2[0],b,t2[0],b) for b in range(t2[1])]
		H += [( t1[0], b, t2[0], (t2[1]+b-t1[1]) ) for b in range(t1[1],(ChrmSzs[chrm]))]
	else: H += [(ChrmIdx[chrm],b,ChrmIdx[chrm],b) for b in range(ChrmSzs[chrm])]
f1 = open(path_m+'/p10.m10.full.mark','w')
f2 = open(path_m+'/m10.p10.full.mark','w')
for h in range(len(H)): print >> f1, H[h][0],H[h][1]*5000 + 2400,H[h][1]*5000 + 2600,H[h][2],H[h][3]*5000 + 2400,H[h][3]*5000 + 2600,1,h
for h in range(len(H)): print >> f2, H[h][2],H[h][3]*5000 + 2400,H[h][3]*5000 + 2600,H[h][0],H[h][1]*5000 + 2400,H[h][1]*5000 + 2600,1,h
f1.close()
f2.close()
# f = open('m10.full.chr.sizes','w')
# print >> f, t1[0],(t1[1]+ChrmSzs[2]-t2[1])*5000
# print >> f, t2[0],(t2[1]+ChrmSzs[1]-t1[1])*5000
# f.close()