import math
import numpy as np
import timeit
import sys
import os

def printlog(str, logname):
	print str
	if logname:
		with open(logname, 'a') as log: print >> log, str

def colorList():
	color = []
	for i in range(0,256,85):
		for j in range(0,256,85):
			for k in range(0,256,85): color.append('%i,%i,%i' % (i,j,k))
	return color

def HashTry(Hash, key):
	t = True
	try: Hash[key]
	except KeyError: t = False
	return t

def FindKey(Hash, key):
	if key != 'interchromosome':
		for k in range(key,-1,-1):
			try: 
				Hash[k]
				return k
			except KeyError: pass
 	else: return key

def boolean(x):
	if x in ['False', 'false', 'no','No','NO']: return False
	elif x in [ 'True', 'true','yes','Yes','YES' ]: return True
	else: return x

def alog2(x):
	try: return abs(math.log(x,2))
	except ValueError: return -10

def readBedGraph(fname, resolution,ChrIdxs):
	bG = {}
	f = open(fname,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)-1,-1,-1):
		parse = lines[i].split()
		key = ChrIdxs[parse[0]], int ((int(parse[1])+int(parse[2])) / (2*resolution))
		bG[key] = float(parse[3])
		del lines[i]
	return bG

def ChromIndexing(path):
	ChrInd = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)):
		parse = lines[i].split()
		try: 
			ChrInd[parse[0]] = i+1
			ChrInd[i+1] = parse[0]
		except IndexError: break
	del lines
	return ChrInd

def ChromSizes(path,resolution):
	ChrSzs = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)):
		parse = lines[i].split()
		try:
			ChrSzs[i+1] = int(parse[1])/resolution
			ChrSzs[parse[0]] = int(parse[1])/resolution
		except IndexError: break
	del lines
	return ChrSzs

def iReadInitialContact(path,ChrIdxs,**kwargs): #Reading Contact from database file
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: chrms = kwargs['chrms']
	except KeyError: chrms = ChrIdxs.keys()
	start_time = timeit.default_timer()
	printlog('\tstart reading contacts', logname)
	contactHash = {}
	files = os.listdir(path)
	lf,lc = len(files),0
	for i in range(lf):
		stf = timeit.default_timer()
		parse = files[i].split('.')
		if parse[-2] in chrms and parse[-3] in chrms:
			printlog('\t\topen %s, %s/%i' % (files[i],i,lf), logname)
			start_time2 = timeit.default_timer()
			f = open(path + '/' +files[i],'r')
			lines = f.readlines()
			f.close()
			elpf = timeit.default_timer() - stf
			printlog('\t\t...close, %.2f' % (elpf) , logname)
			ln = len(lines)
			for j in range(ln-1,-1,-1):
				parse = lines[j].split()
				del lines[j]
				try:
					c1,b1,c2,b2 = ChrIdxs[parse[0]], int(parse[1]), ChrIdxs[parse[2]], int(parse[3])
					p,oe,d = float(parse[4]),float(parse[-2]),float(parse[-1])
					contactHash[c1,b1,c2,b2] = p,oe,d
					lc +=1
				except IndexError: print j, parse
				except ValueError: pass
				if (ln-j) % 1000000 == 0:
					elpf = timeit.default_timer() - stf
					printlog('\t\t\tcontact reading: %i, time elapsed: %.2fs' % (ln-j,elpf), logname)
			elpf = timeit.default_timer() - stf
			elp = timeit.default_timer() - start_time
			printlog('\t\tfile is readed, %.2fs (%.2fs)' % (elpf,elp), logname)
	return contactHash

def readMeanHash(fname,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	psHash = {}
	printlog('\tstart reading distance statistics '+fname, logname)
	start_time = timeit.default_timer()
	f = open(fname,'r')
	lines = f.readlines()
	f.close()
	for i in range(1,len(lines)):
		parse = lines[i].split()
		if len(parse[2:]) == 0: 
			psHash = False
			break
		else: 
			try: psHash[int(parse[0])] = float(parse[4])
			except ValueError: psHash[-1000] = float(parse[4])
	elp = timeit.default_timer() - start_time
	printlog('\t...end reading distance statistics %.2fs' % elp, logname)
	return psHash

def readCovHash(fname,ChrIdxs,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	covHash = {}
	printlog('\tstart reading coverage statistics '+fname, logname)
	start_time = timeit.default_timer()
	f = open(fname,'r')
	lines = f.readlines()
	f.close()
	covHash['all'] = float(lines[0].split()[1])
	covHash['mean'] = float(lines[1].split()[1])
	for i in range(3,len(lines)):
		parse = lines[i].split()
		covHash[ChrIdxs[parse[0]],int(parse[1])] = float(parse[2])
	elp = timeit.default_timer() - start_time
	printlog('\t...end reading coverage statistics %.2fs' % elp, logname)
	return covHash


def iReadingMarkPoints(fname, resolution, ChrIdxs,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: chrms = kwargs['chrms']
	except KeyError: chrms = ChrIdxs.keys()
	ObjCoorMPH = {},{},{}
	ObjCoorMP = {}
	printlog('\tstart reading markpoints '+fname, logname)
	start_time = timeit.default_timer()
	f = open(fname, 'r')
	lines = f.readlines()
	f.close()
	ln = len(lines)
	step = int(ln/10)
	for i in range(ln-1,-1,-1):
		parse = lines[i].split()
		c1,p11,p12,c2,p21,p22 = parse[:6]
		if c1 in chrms:
			c1,c2 = ChrIdxs[c1], ChrIdxs[c2]
			p11,p12,p21,p22 = int(p11),int(p12),int(p21),int(p22)
			p1,p2 = (p11+p12)/2,(p21+p22)/2 
			b1 = c1,p1/resolution
			b2 = c2,p2/resolution
			try: ObjCoorMPH[0][b1].add(p1)
			except KeyError:ObjCoorMPH[0][b1] = set([p1])
			try: ObjCoorMPH[1][b2].add(p2)
			except KeyError:ObjCoorMPH[1][b2] = set([p2])
			try: ObjCoorMPH[2][b1+b2] += 1
			except KeyError:ObjCoorMPH[2][b1+b2] = 1
			if (ln - i) % step == 0:
				elp = timeit.default_timer() - start_time
				printlog('\t\treaded %i/%i markpoins, %.2fs' % ((ln-i),ln,elp), logname)
		del lines[i]
	elp = timeit.default_timer() - start_time
	printlog('\t...end reading markpoints %.2fs' % elp, logname)
	Keys = ObjCoorMPH[2].keys()
	printlog('\tstart remapping coefficient calculating from %i markpoints' % len(Keys), logname)
	for key in Keys:
		c0,c1,cx = len(ObjCoorMPH[0][key[:2]]),len(ObjCoorMPH[1][key[2:]]),ObjCoorMPH[2][key]
		try: ObjCoorMP[key[2:]][key[:2]] = 1.*cx/c1, 1.*cx/c0
		except KeyError: ObjCoorMP[key[2:]] = { key[:2]:(1.*cx/c1, 1.*cx/c0) }
		del ObjCoorMPH[2][key]
	del ObjCoorMPH
	elp = timeit.default_timer() - start_time
	printlog('\t...end remapping coefficient calculating %.2fs' % elp, logname)
	return ObjCoorMP

def iDuplicateContact(Contact_disp_L, ObjCoorMP1, ObjCoorMP2):
	c = [0,0,0,0]
	c[-1] = 1.*ObjCoorMP1[0]*ObjCoorMP2[0]*ObjCoorMP1[1]*ObjCoorMP2[1]
	c[0] = c[-1]*Contact_disp_L[0]
	c[1] = c[-1]*Contact_disp_L[1]
	c[2] = c[-1]*Contact_disp_L[2]
	return c
	

def iDifferContact(Contact_disp, covHash, ObjCoorMP, ChrIdxs, out_name, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: scoring = kwargs['scoring']
	except KeyError: scoring = False
	try: coef = kwargs['contact_coef']
	except KeyError: 
		coef = False
		rescale = 1,1
	try: rescale = kwargs['rescale']
	except KeyError: 
		rescale = 1,1
		coef = False
	try: model = kwargs['model']
	except KeyError: model = 'easy'
	try: random = kwargs['random']
	except KeyError: random = False
	try: predict = kwargs['predict']
	except KeyError: predict = False
	DifferContact = {}
	start_time = timeit.default_timer()
	printlog('\tstart liftovering', logname)
	if model == 'oe': data = 1
	else: data = 0
	elp = timeit.default_timer() - start_time
	printlog('\taccount number of interactions and contacts, %.2f' % elp, logname)
	Keys = ObjCoorMP.keys()
	Keys.sort()
	lnk = len(Keys)
	total,num,S,cl = (lnk+1)*lnk/2,0,covHash['all'],[]
	x = int(total/2000000000+1)
	step = int(total/(100*x))
	if total <= 20000000: step = total
	try: regression = kwargs['regression']
	except KeyError: regression = S
	if regression == 0: regression = S
	elp = timeit.default_timer() - start_time
	printlog('\tmodel %s' % model, logname)
	printlog('\trandomize %s' % random, logname)
	if coef: printlog('\tadditional coef reading with rescale %i %i' % (rescale[0],rescale[1]), logname)
	else: printlog('\tno additional coef reading', logname)
	if regression != S: printlog('\t%i contacts are regressed to %i, %.4f, %.2f' % (S,regression, 1.*regression/S,elp), logname)
	else: printlog('\tno regression, %.2f' % (elp), logname)
	printlog('\tstart processing %i real and %i potential interaction' % (len(Contact_disp),total), logname)
	out = open(out_name,'w')
	for i in range(lnk):
		for j in range(i,lnk):
			num += 1
			key1,key2 = Keys[i],Keys[j]
			if (key1[0] != key2[0]): dk = -1000
			else: dk = abs(key1[1]-key2[1])
			c = np.array([0.,0.,0.,0.])
			for k1 in ObjCoorMP[key1]:
				for k2 in ObjCoorMP[key2]:
					if k1[0] == k2[0] or predict == False: real = True
					else: real = False
					if real:
						r = 1
						if (k1+k2) in Contact_disp: c += iDuplicateContact(Contact_disp[k1+k2], ObjCoorMP[key1][k1], ObjCoorMP[key2][k2])
						elif (k2+k1) in Contact_disp: c += iDuplicateContact(Contact_disp[k2+k1], ObjCoorMP[key1][k1], ObjCoorMP[key2][k2])
						else: real = False
					if real == False and predict == True:
						r = 0
						if coef:
							c1,b1,c2,b2 = k1+k2
							rb1 = (b1*rescale[0]+rescale[1]/2)/rescale[1]
							rb2 = (b2*rescale[0]+rescale[1]/2)/rescale[1]
							if (c1,rb1,c2,rb2) in coef: mlt = coef[c1,rb1,c2,rb2][1]
							elif (c2,rb2,c1,rb1) in coef: mlt = coef[c2,rb2,c1,rb1][1]
							else: mlt = 0
						else: mlt = 1
						try:
							if (k1[0] != k2[0]): mdk = -1000
							else: mdk = abs(k1[1]-k2[1])
							poe = mlt*covHash[k1]*covHash[k2]/covHash['mean']
							pc = scoring[mdk]*poe
							c += iDuplicateContact([pc,poe,mdk], ObjCoorMP[key1][k1], ObjCoorMP[key2][k2])
						except KeyError: pass
					else: pass 
			if c[-1] > 0:
				if model == 'balanced': norm, balance = 1.,c[-1]
				elif model == 'oe': norm, balance = scoring[dk], c[-1]
				elif model == 'strong': 
					try: norm, balance = np.random.choice(2,1,p=[1-c[-1],c[-1]]), c[-1]
					except ValueError: print key1[0], key1[1], key2[0], key2[1], ObjCoorMP[key1],ObjCoorMP[key2],c[-1]
				elif model == 'easy':norm,balance = 1.,1.
				else: norm,balance = 1.,1.
				cc = round(c[data]*norm/balance,4)
				cl.append([key1[0], key1[1], key2[0], key2[1], cc, r, c[0], c[1], c[2]])
			else: pass
			if num % step == 0:
				try:
					np.random.seed()
					cl = np.array(cl)
					if random == 'binomial': cl[:,4] = np.random.binomial(regression,cl[:,4]/S)
					else: cl[:,4] = np.round(regression*cl[:,4]/S)
					cl = cl[cl[:,4] != 0]
					cl = cl.tolist()
					for c in cl: print >> out, '%s %i %s %i %.4f %i %i %.4f %i'% (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8])
				except IndexError: pass
				cl = []
				elp = timeit.default_timer() - start_time
				printlog('\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
	if len(cl) != 0:
		cl = np.array(cl)
		if random == 'binomial': cl[:,4] = np.random.binomial(regression,cl[:,4]/S)
		else: cl[:,4] = np.round(regression*cl[:,4]/S)
		cl = cl[cl[:,4] != 0]
		cl = cl.tolist()
		for c in cl: print >> out, '%s %i %s %i %.4f %i %i %.4f %i'% (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8])
		elp = timeit.default_timer() - start_time
		printlog('\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)
	out.close()
	elp = timeit.default_timer() - start_time
	printlog('\t%i interactions are lifovered, %.2fs.' % (total,elp), logname)

def _random(h,r,S):
	np.random.seed()
	h = np.array(h)
	h[:,4] = np.random.binomial(r,h[:,4]/S)
	h = h[h[:,4] != 0]
	h = h.tolist()
	return h

def iContactRegression(fname,covHash,ChrIdxs,ChrSizes,resolution,out_name,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: scoring = kwargs['scoring']
	except KeyError: scoring = False
	try: 
		coef = kwargs['coef']
		rescale = kwargs['rescale']
	except KeyError: 
		coef = False
		rescale = 1,1
	try: scoring = kwargs['scoring']
	except KeyError: scoring = False
	S,M = covHash['all'],covHash['mean']
	try: regression = kwargs['regression']
	except KeyError: regression = S
	if regression == 0: regression = S
	start_time = timeit.default_timer()
	printlog('\tstart randomize', logname)
	parse = fname.split('.')
	chrm1,chrm2 = parse[-3],parse[-2]
	num,H,Keys = 0, [], set([])
	if chrm1 != chrm2: total,type = ChrSizes[chrm1]*ChrSizes[chrm2],0
	else: total,type = (ChrSizes[chrm1]+1)*ChrSizes[chrm1]/2,1
	out = open(out_name,'w')
	if type==1:
		with open(fname,'r') as f: lines = f.readlines()
		for i in range(len(lines)-1,0,-1):
			num += 1
			parse = lines[i].split()
			del lines[i]
			c1,b1,c2,b2,p = ChrIdxs[parse[0]],int(parse[1]),ChrIdxs[parse[2]],int(parse[3]),float(parse[4])
			Keys.add((c1,b1,c2,b2))
			H.append([c1,b1,c2,b2,p])
			if num % 10000000 == 0:
				try:
					H = _random(H,regression,S)
					for c in H: print >> out, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%i' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4])
				except IndexError: pass
				H = []
				elp = timeit.default_timer() - start_time
				printlog('\t\t %i interactions from %i are processed, %.2fs.' % (num,total,elp), logname)
		elp = timeit.default_timer() - start_time
		printlog('\t%i real interactions from are processed, %.2fs.' % (num,elp), logname)
	printlog('\tstart processing potential contacts', logname)
	for i in range(ChrSizes[chrm1]):
		k = i*type
		for j in range(k,ChrSizes[chrm2]):
			c1,b1,c2,b2 = ChrIdxs[chrm1],i,ChrIdxs[chrm2],j
			if (c1,b2,c2,b2) in Keys: pass
			else:
				num += 1
				if scoring:
					dk = abs(b2-b1)*type + 1000*(type-1)
					score = scoring[dk]
				else: score = 1.0
				if coef:
					rb1 = (b1*rescale[0]+rescale[1]/2)/rescale[1]
					rb2 = (b2*rescale[0]+rescale[1]/2)/rescale[1]
					if (c1,rb1,c2,rb2) in coef: mlt = coef[c1,rb1,c2,rb2][1]
					elif (c2,rb2,c1,rb1) in coef: mlt = coef[c2,rb2,c1,rb1][1]
					else: mlt = 0
				else: mlt=1
				try:
					p = score*mlt*covHash[c1,b1]*covHash[c2,b2]/M
					H.append([c1,b1,c2,b2, p])
				except KeyError: pass
				if num % 10000000 == 0:
					try:
						H = _random(H,regression,S)
						for c in H: print >> out, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%i' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4])
					except IndexError: pass
					H = []
					elp = timeit.default_timer() - start_time
					printlog('\t\t %i interactions from %i are processed, %.2fs.' % (num,total,elp), logname)
	if len(H) != 0:
		try:
			H = _random(H,regression,S)
			for c in H: print >> out, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%i' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4])
		except IndexError: pass
		H = []
		elp = timeit.default_timer() - start_time
		printlog('\t\t %i interactions are processed, %.2fs.' % (num,elp), logname)
	out.close()
	elp = timeit.default_timer() - start_time
	os.system('sort -k 3n -k 7n %s -o %s' % (out_name,out_name))
	printlog('\t%i interactions are randomized, %.2fs.' % (total,elp), logname)

def JuiceboxPre(dir,res,out_path,out_name,**kwargs):
	try: os.mkdir('%s/%s' % (out_path,out_name))
	except OSError: pass
	files = os.listdir(dir)
	for file in files:
		parse = file.split('.')
		chr1,chr2 = parse[-3],parse[-2]
		with open(dir+'/'+file, 'r') as f: lines = f.readlines()
		f1 = open('%s/%s/%s.%s.%s.pre' % (out_path,out_name,out_name,chr1,chr2),'w')
		for line in lines:
			try:
				c1,b1,c2,b2,p = line.split()[:5]
				b1,b2,p = int(b1)*res,int(b2)*res, float(p)
				print >> f1, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%i' % (c1,b1,c2,b2,p)
			except ValueError: print line
		f1.close()

def SummingPre(dir_mut,dir_wt1,dir_wt2,out_path,out_name,**kwargs):
	try: order = kwargs['order']
	except KeyError: order = False
	try: format = kwargs['format']
	except KeyError: format = 'pre'
	print 'order',order
	H,F = {},{}
	names = set([])
	files = os.listdir(dir_mut)
	for file in files: 
		parse = file.split('.')
		F[parse[-3],parse[-2]] = [dir_mut+'/'+file,]
	files = os.listdir(dir_wt1)
	for file in files: 
		parse = file.split('.')
		if (parse[-3],parse[-2]) in F: pass
		else: F[parse[-3],parse[-2]] = [dir_wt1+'/'+file,]
	files = os.listdir(dir_wt2)
	for file in files: 
		parse = file.split('.')
		F[parse[-3],parse[-2]].append(dir_wt2+'/'+file)
	try: os.mkdir('%s/%s.summ' % (out_path,out_name))
	except OSError: pass
	for name in F.keys():
		print F[name]
		H = {}
		with open(F[name][0], 'r') as f: lines = f.readlines()
		for i in range(len(lines)-1,-1,-1):
			x1,c1,b1,y1,x2,c2,b2,y2,p = lines[i].split()
			b1,b2,p = int(b1),int(b2), float(p)
			key = c1,b1,c2,b2 
			H[key] = p
			del lines[i]
		with open(F[name][1], 'r') as f: lines = f.readlines()
		for i in range(len(lines)-1,-1,-1):
			x1,c1,b1,y1,x2,c2,b2,y2,p = lines[i].split()
			b1,b2,p = int(b1),int(b2), float(p)
			key = c1,b1,c2,b2 
			try: H[key] += p
			except KeyError: H[key] = p
			del lines[i]
		del F[name]
		Keys = H.keys()
		if order: 
			print 'order',name
			Keys.sort(key=lambda k:(order[k[0]],order[k[2]],k[1],k[3]))
		else: Keys.sort(key=lambda k:(k[0],k[2],k[1],k[3]))
		fname='%s/%s.summ/%s.summ.%s.%s.pre' % (out_path,out_name,out_name,name[0],name[1])
		f = open(fname,'w')
		if format == 'short':
			for key in Keys: print >> f, '%s\t%i\t%s\t%i\t%i' % (key[0],key[1],key[2],key[3],H[key])
		else:
			for key in Keys: print >> f, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%i' % (key[0],key[1],key[2],key[3],H[key])
		f.close()
	files = os.listdir('%s/%s.summ' % (out_path,out_name))
	if order: files.sort(key=lambda k: (order[k.split('.')[-3]],order[k.split('.')[-2]]))
	else: files.sort(key=lambda k: k.split('.')[-3:-1])
	if format == 'short': end = 'pre.short'
	else:  end = 'pre'
	f = open('%s/%s.summ.%s' % (out_path,out_name,end),'w')
	f.close()
	for file in files:
		os.system( 'cat %s/%s.summ/%s >> %s/%s.summ.%s' % (out_path,out_name,file,out_path,out_name,end) )