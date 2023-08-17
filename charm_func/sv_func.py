def mutation(G,b1,b2,**kwargs):
	try: cnv1 = kwargs['cnv1']
	except KeyError: cnv1 = False
	try: cnv2 = kwargs['cnv2']
	except KeyError: cnv2 = False
	ln = b1[2]-b1[1]
	reg = G[b1[0]][b1[1]:(b1[1]+ln)]
	#print(b1, len(reg))
	if cnv2 < 0: reg.reverse()
	else: pass
	if cnv2:
		reg = reg*abs(cnv2)
		try: G[b2[0]] = G[b2[0]][:b2[1]]+reg+G[b2[0]][b2[1]:]
		except KeyError: G[b2[0]] = reg
	else: pass
	if cnv1: pass
	else: 
		#print(print b1,b2)
		if b1[0] == b2[0] and b2[1] <= b1[1]: del G[b1[0]][(b1[1]+ln):(b1[1]+2*ln)]
		else: del G[b1[0]][b1[1]:(b1[1]+ln)]
	return G

def markgenerate(M,c1,output_dir,resolution,name, header):
	#print('resolution',resolution)
	f1 = open('%s/%s.%s.%i.mark' % (output_dir, name[0],name[1],resolution), 'a')
	f2 = open('%s/%s.%s.%i.mark' % (output_dir, name[1],name[0],resolution), 'a')
	if header == 0:
		f1.write('chrm_ref\tpos_ref_1\tpos_ref_2\tchrm_mut\tpos_mut_1\tpos_mut_2\tpoint number\n')
		f2.write('chrm_mut\tpos_mut_1\tpos_mut_2\tchrm_ref\tpos_ref_1\tpos_ref_2\tpoint number\n')
	k = 0
	r = resolution//2
	for i in range(len(M[c1])):
		k += 1
		f1.write('%s\t%i\t%i\t%s\t%i\t%i\t%i\n' % (M[c1][i][0], M[c1][i][1]*resolution+r-50,M[c1][i][1]*resolution+r+50, c1, i*resolution+r-50,i*resolution+r+50, k) )
		f2.write('%s\t%i\t%i\t%s\t%i\t%i\t%i\n' % (c1, i*resolution+r-50,i*resolution+r+50, M[c1][i][0], M[c1][i][1]*resolution+r-50,M[c1][i][1]*resolution+r+50, k) )
	f1.close()
	f2.close()
	#with open('%s/%s.%s.chr.sizes' % (output_dir, name[0],name[1]), of) as f: print( c1,len(M[c1])*resolution, file=f)

def hashgenerate(M):
	markHash = {}
	for c in M:
		for i in range(len(M[c])): markHash[M[c][i][0],M[c][i][1]] = c, i
		try: markHash[M[c][i][0],M[c][i][1]+1] = c, i
		except IndexError: print('!!!',c,len(M[c]), i)
	return markHash
