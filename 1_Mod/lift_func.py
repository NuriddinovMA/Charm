import math
import numpy as np
import timeit
import sys
import os

def printlog(str, logname):
	
	if logname:
		print str
		if logname == 'stdout': pass
		else:
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
	
	if x in ['False', 'false', 'no','No','NO']: var = False
	elif x in [ 'True', 'true','yes','Yes','YES' ]: var = True
	else: var = x
	try:
		if var[0] in ['False', 'false', 'no','No','NO']: return False
		elif var[0] in [ 'True', 'true','yes','Yes','YES' ]: return True
		else:var = x
	except: pass
	return var

def alog2(x):
	try: return abs(math.log(x,2))
	except ValueError: return -10

def mutation(G,b1,b2,**kwargs):
	try: cnv1 = kwargs['cnv1']
	except KeyError: cnv1 = False
	try: cnv2 = kwargs['cnv2']
	except KeyError: cnv2 = False
	ln = b1[2]-b1[1]
	reg = G[b1[0]][b1[1]:(b1[1]+ln)]
	print b1, len(reg)
	if cnv2 < 0: reg.reverse()
	else: pass
	if cnv2:
		reg = reg*abs(cnv2)
		try: G[b2[0]] = G[b2[0]][:b2[1]]+reg+G[b2[0]][b2[1]:]
		except KeyError: G[b2[0]] = reg
	else: pass
	if cnv1: pass
	else: 
		print b1,b2
		if b1[0] == b2[0] and b2[1] < b1[1]: del G[b1[0]][(b1[1]+ln):(b1[1]+2*ln)]
		else: del G[b1[0]][b1[1]:(b1[1]+ln)]
	return G

def markgenerate(M,c1,output_dir,resolution,name, header):
	print 'resolution',resolution
	f1 = open('%s/%s.%s.%i.mark' % (output_dir, name[0],name[1],resolution), 'a')
	f2 = open('%s/%s.%s.%i.mark' % (output_dir, name[1],name[0],resolution), 'a')
	if header == 0:
		print >> f1, 'chrm_ref\tpos_ref_1\tpos_ref_2\tchrm_mut\tpos_mut_1\tpos_mut_2\tpoint number'
		print >> f2, 'chrm_mut\tpos_mut_1\tpos_mut_2\tchrm_ref\tpos_ref_1\tpos_ref_2\tpoint number'
	k = 0
	r = resolution/2
	for i in range(len(M[c1])):
		k += 1
		print >> f1, '%s\t%i\t%i\t%s\t%i\t%i\t%i' % (M[c1][i][0], M[c1][i][1]*resolution+r-50,M[c1][i][1]*resolution+r+50, c1, i*resolution+r-50,i*resolution+r+50, k)
		print >> f2, '%s\t%i\t%i\t%s\t%i\t%i\t%i' % (c1, i*resolution+r-50,i*resolution+r+50, M[c1][i][0], M[c1][i][1]*resolution+r-50,M[c1][i][1]*resolution+r+50, k)
	f1.close()
	f2.close()
	#with open('%s/%s.%s.chr.sizes' % (output_dir, name[0],name[1]), of) as f: print >> f, c1,len(M[c1])*resolution

def hashgenerate(M):
	markHash = {}
	f = open('test.txt','w')
	for c in M:
		#print(M[c])
		for i in range(len(M[c])):
			markHash[M[c][i][0],M[c][i][1]] = c, i
			print >> f, M[c][i][0], M[c][i][1], c, i
		try: markHash[M[c][i][0],M[c][i][1]+1] = c, i
		except IndexError: print('!!!',c,len(M[c]), i)
	f.close()
	return markHash

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
		try: ChrSzs[parse[0]] = int(parse[1])/resolution
		except IndexError: break
	del lines
	return ChrSzs

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
			try: psHash[int(parse[0])] = float(parse[4]),float(parse[6]),float(parse[7])
			except ValueError: psHash[-1000] = float(parse[4]),float(parse[6]),float(parse[7])
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
	covHash['mult_mean'] = float(lines[2].split()[1])
	for i in range(4,len(lines)):
		parse = lines[i].split()
		try: covHash[ChrIdxs[parse[0]],int(parse[1])] = float(parse[2])
		except KeyError: pass
	elp = timeit.default_timer() - start_time
	printlog('\t...end reading coverage statistics %.2fs' % elp, logname)
	return covHash

def covAlignment(covHash1,covHash2):
	covAli = {}
	for key in covHash1:
		try: covAli[key] = 1.*(covHash2[key]/covHash2['mean'])/(covHash1[key]/covHash1['mean'])
		except KeyError: covAli[key] = 0
		except ZeroDivisionError: covAli[key] = 0
	return covAli

def iReadInitialContact(path,ChrIdxs,**kwargs): #Reading Contact from database file
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: norm = kwargs['norm']
	except KeyError: norm = {}
	try: chrms = kwargs['chrms']
	except KeyError: chrms = ChrIdxs.keys()
	try: scale = kwargs['scale']
	except KeyError: scale = 1
	try: totality = kwargs['totality']
	except KeyError: totality = 'all'
	start_time = timeit.default_timer()
	printlog('\tstart reading contacts', logname)
	contactHash = {}
	files = os.listdir(path)
	lf,lc = len(files),0
	for i in range(lf):
		stf = timeit.default_timer()
		parse = files[i].split('.')
		T = False
		if totality == 'intra': T = ((parse[-2] in chrms) and (parse[-2] == parse[-3]))
		elif totality == 'inter': T = ((parse[-2] in chrms) and (parse[-3] in chrms) and (parse[-2] != parse[-3]))
		else: T = ((parse[-2] in chrms) and (parse[-3] in chrms))
		if T:
			printlog('\t\topen %s, %s/%i' % (files[i],i,lf), logname)
			start_time2 = timeit.default_timer()
			f = open(path + '/' +files[i],'r')
			lines = f.readlines()
			f.close()
			elpf = timeit.default_timer() - stf
			printlog('\t\t...close, %.2f' % (elpf) , logname)
			ln = len(lines)
			Keys = []
			for j in range(ln-1,0,-1):
				parse = lines[j].split()
				del lines[j]
				try:
					c1,b1,c2,b2 = ChrIdxs[parse[0]], int(parse[1])//scale, ChrIdxs[parse[2]], int(parse[3])//scale
					p,oe,mult_cov = float(parse[4]),float(parse[5]),float(parse[6])
					Keys.append((c1,b1,c2,b2))
					try: k = norm[c1,b1]*norm[c2,b2]
					except KeyError: k = 1
					try: 
						contactHash[c1,b1,c2,b2][0] += p*k 
						contactHash[c1,b1,c2,b2][1].append(oe)
						contactHash[c1,b1,c2,b2][2] += mult_cov
					except KeyError: contactHash[c1,b1,c2,b2] = [p*k,[oe*k,],mult_cov]
					lc +=1
				except IndexError: print 'iReadInitialContact IndexError', j, parse
				except ValueError: print 'iReadInitialContact ValueError', j, parse
				if (ln-j) % 1000000 == 0:
					elpf = timeit.default_timer() - stf
					printlog('\t\t\tcontact reading: %i, time elapsed: %.2fs' % (ln-j,elpf), logname)
			for key in Keys: contactHash[key][1] = np.mean(contactHash[key][1])
			del Keys
			elpf = timeit.default_timer() - stf
			elp = timeit.default_timer() - start_time
			printlog('\t\tfile is readed, %.2fs (%.2fs)' % (elpf,elp), logname)
	return contactHash

def iReadingMarkPoints(fname, resolution, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: ChrIdxs_from = kwargs['chrm_index_from']
	except KeyError:
		print 'ERROR!!! Dont exist control chromosome order (chrom_order_from)'
		lf.printlog('ERROR!!! Dont exist control chromosome order (chrom_order_from)', logname)
		return False
	try: ChrIdxs_to = kwargs['chrm_index_to']
	except KeyError: ChrIdxs_to = ChrIdxs_from
	try: chrms_from = kwargs['chrms_from']
	except KeyError: chrms = ChrIdxs_from.keys()
	try: chrms_to = kwargs['chrms_to']
	except KeyError: chrms = ChrIdxs_to.keys()
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
		if c1 in chrms_from and c2 in chrms_to:
			c1,c2 = ChrIdxs_from[c1], ChrIdxs_to[c2]
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

def _recalculateContact(ContactsHash_L, relation1, relation2):
	
	rc = [0,0,0,0]
	rc[-1] = 1.*relation1[0]*relation2[0]*relation1[1]*relation2[1]
	rc[0] = rc[-1]*ContactsHash_L[0]
	rc[1] = rc[-1]*ContactsHash_L[1]
	rc[2] = rc[-1]*ContactsHash_L[2]
	return rc

def _randomizing(h, new_counts, old_counts, random):
	np.random.seed()
	h = np.array(h)
	if random == 'binomial': h[:,4] = np.random.binomial(new_counts,h[:,4]/old_counts)
	elif random == 'normal': 
		loc = np.sqrt(h[:,4]/S)
		loc[loc == 0] = 1
		h[:,4] = np.random.normal(h[:,4],loc)
	elif random == 'hypergeometric':
		ngood = h[:,5]
		ngood[ngood < 1] = 1
		nsample = h[:,6]
		nsample[nsample < 1] = 1
		nbad = old_counts - ngood
		nsample,ngood,nbad = np.int32(np.round(nsample)),np.int32(np.round(ngood)),np.int32(np.round(nbad))
		h[:,4] = np.random.hypergeometric(ngood,nbad,nsample)
		h[:,4] = np.round(new_counts*h[:,4]/old_counts)
	elif random == 'round': h[:,4] = np.round(1.*new_counts*h[:,4]/old_counts,decimals=8)
	elif random == 'choice':
		h[:,4] = 1.
		idx = np.random.choice(len(h),size=int(new_counts),p=old_counts)
		h,h_c = np.unique(h[idx],return_counts=True,axis=0)
		h[:,4] = h[:,4]*h_c
	else: pass
	h = h[h[:,4] > 0] 
	h = h.tolist()
	return h

def _enumerateContacts(coor_to, contH, covH, meanH, mapH, res, params_ec): #model, nullc, ab_cont, ab_cov, ab_res
	c1,b1,c2,b2 = coor_to
	model, nullc, ab_cont, ab_cov, ab_res = params_ec
	pc,poe,mcov = 0.,0.,0.
	lifted_c = np.array([0.,0.,0.,0.])
	if (c1 != c2): real_dist = -1000
	else: real_dist = abs(b2-b1)
	y,r = -1,0
	for k1 in mapH[c1,b1]:
		for k2 in mapH[c2,b2]:
			if k1[0] == k2[0] or nullc == False: real = True
			else: real = False
			if real:
				r = 1
				if (k1+k2) in contH:
					lifted_c += _recalculateContact(contH[k1+k2], mapH[c1,b1][k1], mapH[c2,b2][k2])
				elif (k2+k1) in contH:
					lifted_c += _recalculateContact(contH[k2+k1], mapH[c1,b1][k1], mapH[c2,b2][k2])
				else: real = False
			if real == False and nullc:
				r = 0
				pc,poe,mcov = _predictContacts(k1+k2, covH, meanH, res, params_ec[1:])# nullc, ab_cont, ab_cov, ab_res
				lifted_c += _recalculateContact([pc,poe,mcov], mapH[c1,b1][k1],  mapH[c2,b2][k2])
			else: pass
	data, norm, balance = 0, -1.,-1.
	try: meanH[real_dist]
	except KeyError: real_dist = max(meanH.keys())
	if lifted_c[-1] > 0:
		if model == 'balanced': data,norm, balance = 1, meanH[real_dist][0], lifted_c[-1]
		elif model == 'align_sensitive': data,norm, balance = 1, meanH[real_dist][0], 1.
		elif model == 'distance_sensitive': data, norm, balance = 0, 1.,lifted_c[-1]
		elif model == 'easy': data, norm, balance = 0, 1.,1.
		else: pass
		cc = round(lifted_c[data]*norm/balance,8)
	else: cc = 0
	x = c1, b1, c2, b2, cc, cc/meanH[real_dist][0], round(lifted_c[2]*norm/balance,8), lifted_c[0], lifted_c[1], lifted_c[2], r, meanH[real_dist][0], norm, balance
	return x 

def _predictContacts(coor_from, covH, meanH, res, params_pc):
	
	meanMCov = covH['mult_mean']
	c1,b1,c2,b2 = coor_from
	nullc, ab_cont, ab_cov, ab_res = params_pc
	
	if (c1 != c2): modeled_dist = -1000
	else: modeled_dist = abs(b1-b2)
	if meanH: score = meanH[modeled_dist]
	else: score = 1.0,1.0,1.0
	if ab_cont:
		rb1 = (b1*res+ab_res/2)/ab_res
		rb2 = (b2*res+ab_res/2)/ab_res
		if (c1,rb1,c2,rb2) in ab_cont: cont_AB,mult_AB = ab_cont[c1,rb1,c2,rb2][:2]
		elif (c2,rb2,c1,rb1) in ab_cont: cont_AB,mult_AB = ab_cont[c2,rb2,c1,rb1][:2]
		else: cont_AB,mult_AB = 1,0
	if nullc == 'mixed':
		if modeled_dist == -1000:
			poe = (covH[c1,b1]*covH[c2,b2])/meanMCov
			pc = mult_AB*poe*score[0]
		else:
			try: pc = cont_AB*(covH[c1,b1]*covH[c2,b2])/(ab_cov[c1,rb1]*ab_cov[c2,rb2])**2
			except ZeroDivisionError: pc = 0
			except KeyError: pc = 0
			try: poe = mult_AB*pc/score[0]
			except ZeroDivisionError: poe = 0
			except KeyError: poe = 0
	elif nullc == 'ab':
		try: pc = cont_AB*(covH[c1,b1]*covH[c2,b2])/(ab_cov[c1,rb1]*ab_cov[c2,rb2])**2
		except ZeroDivisionError: pc = 0
		try: poe = mult_AB*pc/score[0]
		except ZeroDivisionError: poe = 0
	elif nullc == 'sq_ab':
		poe = mult_AB
		try: pc = cont_AB*np.sqrt((covH[c1,b1]*covH[c2,b2])/(ab_cov[c1,rb1]*ab_cov[c2,rb2]))
		except ZeroDivisionError: pc = 0
	elif nullc == 'mc':
		poe = (covH[c1,b1]*covH[c2,b2])/meanMCov
		pc = score[0]*poe
	elif nullc == 'sq_mc':
		poe = np.sqrt((covH[c1,b1]*covH[c2,b2])/meanMCov)
		pc = score[0]*poe
	else:
		poe = (covH[c1,b1]*covH[c2,b2])/meanMCov
		pc = score[0]*poe
	modelled_cont = pc, poe, covH[c1,b1]*covH[c2,b2]
	return modelled_cont

def _rescaleContacts(high_cl, high_res, contH, covH, meanH, mapH, res, params_rand, params_rc): # params_rC: model, nullc, ab_cont, ab_cov, ab_res
	_cl = []
	regression,allCon,random = params_rand
	for hi in high_cl:
		c1,b1,c2,b2,count = hi[:5]
		ri1,ri2,rj1,rj2 = int(b1*high_res/res),int((b1+1)*high_res/res),int(b2*high_res/res),int((b2+1)*high_res/res)
		_s,_clh = [],[]
		for ri in range(ri1,ri2):
			if (c1,b1) == (c2,b2): rj1 = ri
			for rj in range(rj1,rj2):
				try:
					mapH[c1,ri],mapH[c2,rj],
					y = _enumerateContacts((c1,ri,c2,rj), contH, covH, meanH, mapH, res, params_rc)
					_clh.append(y)
				except KeyError: pass
		_clh = np.array(_clh)
		try:
			if random: 
				p = 1.*_clh[:,9]*_clh[:,11]
				p = p/np.sum(p)
				_clh = _randomizing(_clh, count, p, 'choice')
			else: _clh = _randomizing(_clh, regression,allCon,random)
			_clh.sort()
			_cl.extend(_clh)
		except IndexError: pass#print 'PRCC IndexError:' #_cl, _counts, _sum_count, high_cl[:10], high_cl[-10:]
		except ValueError: pass
	return _cl

def iLiftOverContact(ContactsHash, covHash, ObjCoorMP, resolution, ChrIdxs, out_name, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: scoring = kwargs['scoring']
	except KeyError: scoring = False
	try: CoefObjCoorMP = kwargs['markpoints_coef']
	except KeyError: CoefObjCoorMP = False
	try: cov_coef = kwargs['coverage_coef']
	except KeyError: cov_coef = False
	try: scoring_coef = kwargs['scoring_coef']
	except KeyError: scoring_coef = False
	try: contact_coef = kwargs['contact_coef']
	except KeyError:
		contact_coef = False
		cov_coef = False
		scoring_coef = False
		CoefObjCoorMP = False
	try: coef_res = kwargs['coef_resolution']
	except KeyError: c_res = 1
	try: ab_cov = kwargs['coverage_ab']
	except KeyError: ab_cov = False
	try: ab_cont = kwargs['contact_ab']
	except KeyError: ab_cont,ab_cov = False,False
	try: ab_res = kwargs['resolution_ab']
	except KeyError: ab_res = 1
	try: model = kwargs['model']
	except KeyError: model = 'easy'
	try: random = kwargs['random']
	except KeyError: random = False
	try: nullc = kwargs['null_contacts']
	except KeyError: nullc = False
	
	DifferContact = {}
	start_time = timeit.default_timer()
	printlog('\tstart liftovering', logname)
	elp = timeit.default_timer() - start_time
	printlog('\taccount number of interactions and contacts, %.2f' % elp, logname)
	Keys = ObjCoorMP.keys()
	Keys.sort()
	lnk = len(Keys)
	total,num = (lnk+1)*lnk/2,0,
	if CoefObjCoorMP:
		CoefKeys = CoefObjCoorMP.keys()
		CoefKeys.sort()
		coef_lnk = len(CoefKeys)
		total = (coef_lnk+1)*coef_lnk/2
	allCon,meanMCov = covHash['all'],covHash['mult_mean']
	try: regression = kwargs['regression']
	except KeyError: regression = allCon
	if regression == 0: regression = allCon
	elp = timeit.default_timer() - start_time
	printlog('\tmodel %s' % model, logname)
	printlog('\trandomize %s' % random, logname)
	printlog('\tpredict null contact %s' % nullc, logname)
	if contact_coef and coef_res: printlog('\tadditional coef reading with rescale %i %i' % (resolution,coef_res), logname)
	else: printlog('\tno additional coef reading', logname)
	if regression != allCon: printlog('\t%i contacts are regressed to %i, %.4f, %.2f' % (allCon,regression, 1.*regression/allCon,elp), logname)
	else: printlog('\tno regression, %.2f' % (elp), logname)
	out = open(out_name,'w')
	
	params_model = model, nullc, ab_cont, ab_cov, ab_res
	params_random = regression, allCon, random
	
	if CoefObjCoorMP:
		coef_cl = []
		printlog('\tstart processing %i real and %i potential interaction by Scaled Simulation' % (len(contact_coef),total), logname)
		for i in range(coef_lnk):
			for j in range(i,coef_lnk):
				num += 1
				key1,key2 = CoefKeys[i],CoefKeys[j]
				coef_cl.append( _enumerateContacts(key1+key2, contact_coef, cov_coef, scoring_coef, CoefObjCoorMP, coef_res, params_model ) )#, model, nullc, ab_cont, ab_cov, ab_res) )
				if num % 2000000 == 0:
					try: coef_cl = _randomizing(coef_cl, regression, allCon, random)
					except IndexError: print 'iLOF-1, randomize, IndexError:', coef_cl
					cl = _rescaleContacts(coef_cl, coef_res, ContactsHash, covHash, scoring, ObjCoorMP, resolution, params_random, params_model)#, model, nullc, ab_cont, ab_cov, ab_res, regression, random)
					# hash = {}
					# for c in range(len(cl)-1,-1,-1):
						# try: hash[cl[c][0],cl[c][1],cl[c][2],cl[c][3]][0] += cl[c][4]
						# except KeyError: hash[cl[c][0],cl[c][1],cl[c][2],cl[c][3]] = cl[c][4:]
						# del cl[c]
					# del cl
					# for c in hash.keys():
						# print >> out, '%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], hash[c][0], hash[c][1], hash[c][2], hash[c][3], hash[c][4], hash[c][5], hash[c][6], hash[c][7], hash[c][8], hash[c][9])
					for c in cl:
						print >> out, '%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13])
					coef_cl,cl,hash = [],[],{}
					elp = timeit.default_timer() - start_time
					printlog('\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
		if len(coef_cl) != 0:
			try: coef_cl = _randomizing(coef_cl, regression, allCon, random)
			except IndexError: print 'iLOF-1, randomize, IndexError:', coef_cl
			cl = _rescaleContacts(coef_cl, coef_res, ContactsHash, covHash, scoring, ObjCoorMP, resolution, params_random, params_model)#, model, nullc, ab_cont, ab_cov, ab_res, regression, random)
			# hash = {}
			# for c in range(len(cl)-1,-1,-1):
				# try: hash[cl[c][0],cl[c][1],cl[c][2],cl[c][3]][0] += cl[c][4]
				# except KeyError: hash[cl[c][0],cl[c][1],cl[c][2],cl[c][3]] = cl[c][4:]
				# del cl[c]
			# del cl
			# for c in hash.keys():
				# print >> out, '%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], hash[c][0], hash[c][1], hash[c][2], hash[c][3], hash[c][4], hash[c][5], hash[c][6], hash[c][7], hash[c][8], hash[c][9])
			for c in cl:
				print >> out, '%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13])
			coef_cl,cl,hash = [],[],{}
			elp = timeit.default_timer() - start_time
			printlog('\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)
	else:
		cl = []
		printlog('\tstart processing %i real and %i potential interaction' % (len(ContactsHash),total), logname)
		for i in range(lnk):
			for j in range(i,lnk):
				num += 1
				key1,key2 = Keys[i],Keys[j]
				cl.append( _enumerateContacts(key1+key2, ContactsHash, covHash, scoring, ObjCoorMP, resolution, params_model) )#, model, nullc, ab_cont, ab_cov, ab_res))
				if num % 2000000 == 0:
					try: cl = _randomizing(cl, regression, allCon, random)
					except IndexError: print 'iLOF-3, randomize, IndexError:',cl
					for c in cl:
						print >> out, '%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13])
					cl = []
					elp = timeit.default_timer() - start_time
					printlog('\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
		if len(cl) != 0:
			cl = _randomizing(cl, regression, allCon, random)
			for c in cl: 
				print >> out, '%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13])
			cl = []
			elp = timeit.default_timer() - start_time
			printlog('\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)
	out.close()
	elp = timeit.default_timer() - start_time
	printlog('\t%i interactions are lifovered, %.2fs.' % (total,elp), logname)

def _listing(coor, contactHash, covH, meanH, res, params_list):
	
	c1,b1,c2,b2 = coor
	nullc, ab_cont, ab_cov, ab_res = params_pc = params_list
	pc,poe,mcov = 0,0,0
	if (c1 != c2): real_dist = -1000
	else: real_dist = abs(b1-b2)
	try: mean_pc = meanH[real_dist][0]
	except KeyError: mean_pc = meanH[max(meanH.keys())][0]
	
	if (c1,b1,c2,b2) in contactHash: pc,poe,mcov = contactHash[c1,b1,c2,b2]
	elif (c2,b2,c1,b1) in contactHash: pc,poe,mcov = contactHash[c2,b2,c1,b1]
	elif nullc:
		try: pc,poe,mcov = _predictContacts((c1,b1,c2,b2), covH, meanH, res, params_list)
		except KeyError: pass
	else: pass
	_h = (c1,b1,c2,b2,pc,poe,mcov,mean_pc)
	return _h

def _easyRescale(high_cl, high_res, covH, meanH, res, params_rand):
	regression, allCon, random = params_rand
	_cl,_clh = [],[]
	for hi in range(len(high_cl)-1,-1,-1):
		c1,b1,c2,b2,count = high_cl[hi][:5]
		del high_cl[hi]
		ri1,ri2,rj1,rj2 = int(b1*high_res/res),int((b1+1)*high_res/res),int(b2*high_res/res),int((b2+1)*high_res/res)
		_clh = []
		for ri in range(ri1,ri2):
			if (c1,b1) == (c2,b2): rj1 = ri
			for rj in range(rj1,rj2):
				try: 
					covH[c1,ri],covH[c2,rj]
					if (c1 != c2): dist = -1000
					else: dist = abs(ri-rj)
					try: mean_con = meanH[dist][0]
					except KeyError: mean_con =  meanH[max(meanH.keys())][0]
					_clh.append((c1,ri,c2,rj,1.,1.,covH[c1,ri]*covH[c2,rj],mean_con))
				except KeyError: pass
		try:
			if random:
				_clh = np.array(_clh)
				p = 1.*_clh[:,-2]*_clh[:,-1]
				p = p/np.sum(p)
				_clh = _randomizing(_clh, count, p, 'choice')
			else: _clh = _randomizing(_clh, regression, allCon, random)
		except IndexError: pass#
		except ValueError: pass
		_cl.extend(_clh)
		_clh = []
	del _clh
	# hash = {}
	# for c in range(len(_cl)-1,-1,-1):
		# try: hash[_cl[c][0],_cl[c][1],_cl[c][2],_cl[c][3]] += _cl[c][4]
		# except KeyError: hash[_cl[c][0],_cl[c][1],_cl[c][2],_cl[c][3]] = _cl[c][4]
		# del _cl[c]
	# del _cl
	return _cl

def iContactRegression(ContactHash,covHash,resolution,chroms,ChrIdxs,coef_ChrSizes,out_name,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: scoring = kwargs['scoring']
	except KeyError: scoring = False
	try: CoefObjCoorMP = kwargs['markpoints_coef']
	except KeyError: CoefObjCoorMP = False
	try: cov_coef = kwargs['coverage_coef']
	except KeyError:
		cov_coef = False
		CoefObjCoorMP = False
	try: scoring_coef = kwargs['scoring_coef']
	except KeyError:
		scoring_coef = False
		CoefObjCoorMP = False
	try: contact_coef = kwargs['contact_coef']
	except KeyError:
		contact_coef = False
		cov_coef = False
		scoring_coef = False
		CoefObjCoorMP = False
	try: coef_res = kwargs['coef_resolution']
	except KeyError: c_res = 1
	try: ab_cov = kwargs['coverage_ab']
	except KeyError: ab_cov = False
	try: ab_cont = kwargs['contact_ab']
	except KeyError: ab,ab_cov = False,False
	try: ab_res = kwargs['resolution_ab']
	except KeyError: ab_res = 1
	try: model = kwargs['model']
	except KeyError: model = 'easy'
	try: random = kwargs['random']
	except KeyError: random = False
	try: predict = kwargs['predict']
	except KeyError: predict = False
	try: nullc = kwargs['null_contacts']
	except KeyError: nullc = False
	try: regression = kwargs['regression']
	except KeyError: regression = S
	if regression == 0: regression = S
	total = 0
	allCon,meanMCov = covHash['all'],covHash['mult_mean']
	params_model = nullc, ab_cont, ab_cov, ab_res
	params_random = regression, allCon, random
	for ci in range(len(chroms)):
		for cj in range((ci+1),len(chroms)):
			chrm1,chrm2 = chroms[ci],chroms[cj]
			if chrm1 == chrm2: total += (coef_ChrSizes[chrm1]+1)*coef_ChrSizes[chrm1]/2
			else: total += coef_ChrSizes[chrm1]*coef_ChrSizes[chrm2]
	num,H,cl = 0, [], []
	start_time = timeit.default_timer()
	elp = timeit.default_timer() - start_time
	printlog('\t%i contacts are regressed to %i, %.2fs.' % (allCon,regression,elp), logname)
	print chroms,coef_ChrSizes
	for ci in range(len(chroms)):
		for cj in range((ci+1),len(chroms)):
			chrm1,chrm2 = chroms[ci],chroms[cj]
			full_name = '%s.%s.%s.allCon'% (out_name,chrm1,chrm2)
			out = open( full_name,'w')
			elp = timeit.default_timer() - start_time
			printlog('\t\tregression %s %s chromosome pair, %.2fs.' % (chrm1,chrm2,elp), logname)
			if chrm1 == chrm2: k = 1
			else: k = 0
			for i in range(coef_ChrSizes[chrm1]):
				j1 = i*k
				for j in range(j1,coef_ChrSizes[chrm2]):
					num += 1.
					c1,b1,c2,b2 = ChrIdxs[chrm1],i,ChrIdxs[chrm2],j
					H.append(_listing((c1,b1,c2,b2), contact_coef, cov_coef, scoring_coef, coef_res, params_model))
					if num % 200000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t\t%s randomization, %.2fs.' % (random,elp), logname)
						H = _randomizing(H, regression, allCon, random)
						H = _easyRescale(H, coef_res, covHash, scoring, resolution, params_random)
						elp = timeit.default_timer() - start_time
						printlog('\t\twriting, %.2fs.' % (elp), logname)
						for c in H:
							try: print >> out, '%s\t%i\t%s\t%i\t%.4f' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4])
							except IndexError: print 'IndexError', c
						elp = timeit.default_timer() - start_time
						printlog('\t\t%i interactions from %i are processed, %.2fs.' % (num,total,elp), logname)
						H = []
			if len(H) > 0:
				H = _randomizing(H, regression, allCon, random)
				H = _easyRescale(H, coef_res, covHash, scoring, resolution, params_random)
				elp = timeit.default_timer() - start_time
				printlog('\t\twriting, %.2fs.' % (elp), logname)
				for c in H:
					try: print >> out, '%s\t%i\t%s\t%i\t%.4f' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4])
					except IndexError:print 'IndexError', c
			del H
			elp = timeit.default_timer() - start_time
			printlog('\t\tend regression %s %s chromosome pair, %.2fs.' % (chrm1,chrm2,elp), logname)
			out.close()
			#os.system('sort -k 2n -k 4n %s -o %s' % (full_name, full_name))
	elp = timeit.default_timer() - start_time
	printlog('\t%i interactions are randomized, %.2fs.' % (total,elp), logname)

def SummingPre(dir_mut,dir_wt1,dir_wt2,out_path,out_name,**kwargs):
	try: order = kwargs['order']
	except KeyError: order = False
	try: format = kwargs['format']
	except KeyError: format = 'pre'
	try: res = kwargs['out_res']
	except KeyError: res = 1
	print 'order',order
	H,F = {},{}
	names = set([])
	if dir_mut:
		files = os.listdir(dir_mut)
		for file in files: 
			parse = file.split('.')
			F[parse[-3],parse[-2]] = [dir_mut+'/'+file,0,1]
	if dir_wt1:
		files = os.listdir(dir_wt1)
		for file in files: 
			parse = file.split('.')
			if (parse[-3],parse[-2]) in F: pass
			else: F[parse[-3],parse[-2]] = [dir_wt1+'/'+file,0,res]
	if dir_wt2:
		files = os.listdir(dir_wt2)
		for file in files: 
			parse = file.split('.')
			try: F[parse[-3],parse[-2]][1] = dir_wt2+'/'+file
			except KeyError: pass
	os.system('mkdir %s/%s.summ' % (out_path,out_name))
	for name in F.keys():
		print F[name]
		H = {}
		with open(F[name][0], 'r') as f: lines = f.readlines()
		for i in range(len(lines)-1,-1,-1):
			div = F[name][2]
			try:
				c1,b1,c2,b2,p = lines[i].split()[:5]
				b1,b2,p = int(b1),int(b2), float(p)
				key = c1,b1//div*res,c2,b2//div*res
				if np.isnan(p): p = 0
				try: H[key] += p
				except KeyError: H[key] = p
				del lines[i]
			except ValueError:pass
		if dir_wt2:
			with open(F[name][1], 'r') as f: lines = f.readlines()
			print F[name][1]
			for i in range(len(lines)-1,-1,-1):
				try:
					c1,b1,c2,b2,p = lines[i].split()[:5]
					b1,b2,p = int(b1),int(b2), float(p)
					key = c1,b1//res*res,c2,b2//res*res
					if np.isnan(p): p = 0
					try: H[key] += p
					except KeyError: H[key] = p
					del lines[i]
				except ValueError:pass
		del F[name]
		Keys = H.keys()
		if order: Keys.sort(key=lambda k:(order[k[0]],order[k[2]],k[1],k[3]))
		else: Keys.sort(key=lambda k:(k[0],k[2],k[1],k[3]))
		fname='%s/%s.summ/%s.summ.%s.%s.pre' % (out_path,out_name,out_name,name[0],name[1])
		f = open(fname,'w')
		if format == 'short':
			for key in Keys: print >> f, '%s\t%i\t%s\t%i\t%.8f' % (key[0],key[1],key[2],key[3],H[key])
		else:
			for key in Keys: print >> f, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%.8f' % (key[0],key[1],key[2],key[3],H[key])
		f.close()
	files = os.listdir('%s/%s.summ' % (out_path,out_name))
	if order: files.sort(key=lambda k: (order[k.split('.')[-3]],order[k.split('.')[-2]]))
	else: files.sort(key=lambda k: k.split('.')[-3:-1])
	if format == 'short': end = 'pre.short'
	else: end = 'pre'
	f = open('%s/%s.summ.%s' % (out_path,out_name,end),'w')
	f.close()
	for file in files:
		os.system( 'cat %s/%s.summ/%s >> %s/%s.summ.%s' % (out_path,out_name,file,out_path,out_name,end) )

# def SummingPointview(dir_mut,dir_wt1,dir_wt2,pointview,out_path,out_name,**kwargs):
	# try: order = kwargs['order']
	# except KeyError: order = False
	# try: format = kwargs['format']
	# except KeyError: format = 'pre'
	# try: res = kwargs['out_res']
	# except KeyError: res = 1
	# print pointview
	# H,F = {},{}
	# names = set([])
	# files = os.listdir(dir_mut)
	# for file in files: 
		# parse = file.split('.')
		# F[parse[-3],parse[-2]] = [dir_mut+'/'+file,]
	# if dir_wt1:
		# files = os.listdir(dir_wt1)
		# for file in files: 
			# parse = file.split('.')
			# if (parse[-3],parse[-2]) in F: pass
			# elif pointview[0] in (parse[-3],parse[-2]): F[parse[-3],parse[-2]] = [dir_wt1+'/'+file,]
			# else: pass
	# if dir_wt2:
		# files = os.listdir(dir_wt2)
		# for file in files: 
			# parse = file.split('.')
			# if (parse[-3],parse[-2]) in F: F[parse[-3],parse[-2]].append(dir_wt2+'/'+file)
			# else: pass
	# for key in F: print F[key]
	# os.system('mkdir %s/%s.summ' % (out_path,out_name))
	# print F
	# for name in F.keys():
		# print F[name]
		# H = {}
		# with open(F[name][0], 'r') as f: lines = f.readlines()
		# for i in range(len(lines)-1,-1,-1):
			# x1,c1,b1,y1,x2,c2,b2,y2,p = lines[i].split()
			# b1,b2,p = int(b1),int(b2), float(p)
			# if np.isnan(p): p = 0
			# if ( c1 == pointview[0] and pointview[1] <= b1 < pointview[2] ) or ( c2 == pointview[0] and pointview[1] <= b2 < pointview[2] ):
				# key = c1,b1//res*res,c2,b2//res*res
				# try: H[key] += p
				# except KeyError: H[key] = p
			# del lines[i]
		# if dir_wt2:
			# with open(F[name][1], 'r') as f: lines = f.readlines()
			# for i in range(len(lines)-1,-1,-1):
				# x1,c1,b1,y1,x2,c2,b2,y2,p = lines[i].split()
				# b1,b2,p = int(b1),int(b2), float(p)
				# if np.isnan(p): p = 0
				# if ( c1 == pointview[0] and pointview[1] <= b1 < pointview[2] ) or ( c2 == pointview[0] and pointview[1] <= b2 < pointview[2] ):
					# key = c1,b1//res*res,c2,b2//res*res
					# try: H[key] += p
					# except KeyError: H[key] = p
				# del lines[i]
		# del F[name]
		# Keys = H.keys()
		# if order: 
			# print 'order',name
			# Keys.sort(key=lambda k:(order[k[0]],order[k[2]],k[1],k[3]))
		# else: Keys.sort(key=lambda k:(k[0],k[2],k[1],k[3]))
		# fname='%s/%s.summ/%s.summ.%s.%s.pre' % (out_path,out_name,out_name,name[0],name[1])
		# f = open(fname,'w')
		# if format == 'short':
			# for key in Keys: print >> f, '%s\t%i\t%s\t%i\t%.4f' % (key[0],key[1],key[2],key[3],H[key])
		# else:
			# for key in Keys: print >> f, '0\t%s\t%i\t0\t1\t%s\t%i\t1\t%.4f' % (key[0],key[1],key[2],key[3],H[key])
		# f.close()
	# files = os.listdir('%s/%s.summ' % (out_path,out_name))
	# if order: files.sort(key=lambda k: (order[k.split('.')[-3]],order[k.split('.')[-2]]))
	# else: files.sort(key=lambda k: k.split('.')[-3:-1])
	# if format == 'short': end = 'pre.short'
	# else: end = 'pre'
	# f = open('%s/%s.summ.%s' % (out_path,out_name,end),'w')
	# f.close()
	# for file in files:
		# os.system( 'cat %s/%s.summ/%s >> %s/%s.summ.%s' % (out_path,out_name,file,out_path,out_name,end) )
	
def AddNormVector(path,Rlist,hic,norm):
	parse = Rlist.split(',')
	f = open(hic+'.norm','w')
	for r in parse:
		resolution = int(r)
		chrSizes = ChromSizes(path,resolution)
		Keys = chrSizes.keys()
		Keys.sort()
		for key in Keys:
			try:
				print >> f, 'vector ' + norm + ' ' + key + ' ' + r + ' BP'
				for i in range(chrSizes[key]): print >> f, '1.0'
			except TypeError: pass
			

def netParser(name,**kwargs):
	try: min_length = kwargs['min_length']
	except KeyError: min_length = 0
	try: gap_length = kwargs['gap_length']
	except KeyError: gap_length = 100000
	try: chr_from = kwargs['chr_from']
	except KeyError: chr_from = False
	try: chr_to = kwargs['chr_to']
	except KeyError: chr_to = False
	f = open(name, 'r')
	lines = f.readlines()
	f.close()
	parsedNet,parsedId = {}, {}
	id = 0
	for line in lines:
		parse = line.split()
		if parse[0] == 'net': 
			chrm = parse[1]
			if chr_from == chrm or chr_from == False: parsedNet[chrm] = []
			else: chrm = False
			i = 0
		elif parse[0] == 'fill':
			start1,ln1,start2,ln2 = int(parse[1]),int(parse[2]),int(parse[5]),int(parse[6])
			dir = int(parse[4]+'1')
			fi1,fi2 = start1+ln1,start2+ln2
			if (chr_to == parse[3] or chr_to == False) and chrm and ln1 > min_length and ln2 > min_length:
				id += 1
				i += 1
				data = chrm,start1,fi1,parse[3],start2,fi2,dir,id
				parsedNet[chrm].append( data )
	del lines
	print 'lines readed'
	chrms = parsedNet.keys()
	f = open('test.csv','w')
	for chrm in chrms:
		list(set(parsedNet[chrm]))
		parsedNet[chrm].sort()
		ln = len(parsedNet[chrm])
		for i in range(ln-1):
			chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0 = parsedNet[chrm][i]
			nid = id0
			try:
				parsedId[nid]
				print >> f, 'repeat', chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0
			except KeyError:
				print >> f, 'start', chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0
				parsedId[nid] = set([parsedNet[chrm][i],])
				for j in range(i+1,ln):
					chrm_from1,start_from1,fi_from1,chrm_to1,start_to1,fi_to1,dir1,id1 = parsedNet[chrm][j]
					if start_from1 - fi_from0 > gap_length or dir1 != dir0:
						break
						print >> f, 'break_0', chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0
						print >> f, 'break', chrm_from1,start_from1,fi_from1,chrm_to1,start_to1,fi_to1,dir1,id1
					elif chrm_to1 == chrm_to0:
						gap = start_from1 - fi_from0
						if dir1 == dir0 == 1 and start_to1 >= fi_to0 and 1000 < (start_to1 - fi_to0) < gap_length:
							data = chrm_from1,(fi_from0+250),(start_from1-250),chrm_to1,(fi_to0+250),(start_to1-250),dir1,-nid
							parsedId[nid].add(data)
							parsedNet[chrm][j] = parsedNet[chrm][j][:-1] + (nid,)
							parsedId[nid].add(parsedNet[chrm][j])
							print >> f, 'grow', chrm_from1,start_from1,fi_from1,chrm_to1,start_to1,fi_to1,dir1,id1
							chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0 = parsedNet[chrm][j]
						elif dir1 == dir0 == 1 and start_to1 >= fi_to0 and (start_to1 - fi_to0) < gap_length:
							parsedNet[chrm][j] = parsedNet[chrm][j][:-1] + (nid,)
							parsedId[nid].add(parsedNet[chrm][j])
							print >> f, 'grow', chrm_from1,start_from1,fi_from1,chrm_to1,start_to1,fi_to1,dir1,id1
							chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0 = parsedNet[chrm][j]
						elif dir1 == dir0 == -1 and start_to1 <= fi_to0 and 1000 < (start_to0 - fi_to1) < gap_length:
							data = chrm_from1,(fi_from0+250),(start_from1-250),chrm_to1,(fi_to1+250),(start_to0-250),dir1,-nid
							parsedId[nid].add(data)
							parsedNet[chrm][j] = parsedNet[chrm][j][:-1] + (nid,)
							parsedId[nid].add(parsedNet[chrm][j])
							print >> f, 'grow', chrm_from1,start_from1,fi_from1,chrm_to1,start_to1,fi_to1,dir1,id1
							chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0 = parsedNet[chrm][j]
						elif dir1 == dir0 == -1 and start_to1 <= fi_to0 and (start_to0 - fi_to1) < gap_length:
							parsedNet[chrm][j] = parsedNet[chrm][j][:-1] + (nid,)
							parsedId[nid].add(parsedNet[chrm][j])
							print >> f, 'grow', chrm_from1,start_from1,fi_from1,chrm_to1,start_to1,fi_to1,dir1,id1
							chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0 = parsedNet[chrm][j]
						else:
							print >> f, 'pass_0', chrm_from0,start_from0,fi_from0,chrm_to0,start_to0,fi_to0,dir0,id0
							print >> f, 'pass', chrm_from1,start_from1,fi_from1,chrm_to1,start_to1,fi_to1,dir1,id1
					else: pass
	f.close()
	
	pN = []
	Keys = parsedId.keys()
	for key in Keys:
		if len(parsedId[key])> 0:
			parsedId[key] = list(parsedId[key])
			parsedId[key].sort()
			pN.append( parsedId[key] )
		del parsedId[key]
	del Keys
	del parsedNet
	pN.sort()
	return pN

def net2pre(parsedNet,out,**kwargs):
	markPoints = []
	try: bin_length = kwargs['bin_length']
	except KeyError: bin_length = 1000
	f1 = open(out+'.pre.mark', 'w')
	f2 = open(out+'.2D.ann', 'w')
	
	print >> f1, 'chr1\tstart1\tend1\tchr2\tstart2\tend2\tid'
	print >> f2, 'chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tcomment'
	for i in parsedNet:
		a = 0
		if i[0][6] > 1: st,fi = i[0][4],i[-1][5]
		else: st,fi = i[-1][4],i[0][5]
		print >> f2, '%s\t%i\t%i\t%s\t%i\t%i\t%s\t%s:%i-%i:%i' % (i[0][0],i[0][1],i[-1][2],i[0][0],i[0][1],i[-1][2],'0,255,0',i[0][3],st,fi,i[0][6])
		for j in i:
			try:
				if j[6] > 0: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i\t%i\t%i' % (j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7])
				else: print >> f1, '%s\t%i\t%i\t%s\t%i\t%i\t%i\t%i' % (j[0],j[1],j[2],j[3],j[5],j[4],j[6],j[7])
			except TypeError: a = 1
		if a == 1: print i
	f1.close()
	f2.close()
	
	f3 = open(out+'.mark', 'w')
	for i in parsedNet:
		dir = i[0][6]
		if dir > 0:
			for j in range(len(i)):
				ln1,ln2 = i[j][2]-i[j][1],i[j][5]-i[j][4]
				if ln1 > 300: 
					k = ln1/150.0
					sp1,sp2 = ln1/k,ln2/k
					h = []
					for n in range(int(k)): h.append((i[j][1]+sp1*n,i[j][4]+sp2*n))
					h.append( (i[j][2],i[j][5] ) )
					for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]),int(h[n+1][0]),i[j][3],int(h[n][1]),int(h[n+1][1]),i[j][6],i[j][7]) )
				else: markPoints.append(i[j])
				try:
					gap1,gap2 = i[j+1][1]-i[j][2],i[j+1][4]-i[j][5]
					if gap1 > 1500:
						k = gap1/1000.0
						sp1,sp2 = gap1/k,gap2/k
						h = []
						for n in range(int(k)): h.append( (i[j][2]+sp1*(n+.5),i[j][5]+sp2*(n+.5)) )
						for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]-100),int(h[n][0]+100),i[j][3],int(h[n][1]-100),int(h[n][1]+100),i[j][6],str(i[j][7])+'_gap') )
					elif gap1 < 500: pass
					else: markPoints.append( (i[j][0],i[j][2]+gap1/2-100,i[j][2]+gap1/2+100,i[j][3],i[j][5]+gap2/2-100,i[j][5]+gap2/2+100,i[j][6],str(i[j][7])+'_gap') )
				except IndexError: pass
		else:
			for j in range(len(i)):
				ln1,ln2 = i[j][2]-i[j][1],i[j][5]-i[j][4]
				if ln1 > 300: 
					k = ln1/150.0
					sp1,sp2 = ln1/k,ln2/k
					h = []
					for n in range(int(k)): h.append((i[j][1]+sp1*n,i[j][5]-sp2*n))
					h.append( (i[j][2],i[j][4] ) )
					for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]),int(h[n+1][0]),i[j][3],int(h[n+1][1]),int(h[n][1]),i[j][6],i[j][7]) )
				else: markPoints.append(i[j])
				try:
					gap1,gap2 = i[j+1][1]-i[j][2],i[j][4]-i[j+1][5]
					if gap1 > 1500: 
						k = ln1/1000.0
						sp1,sp2 = gap1/k,gap2/k
						h = []
						for n in range(int(k)): h.append((i[j][2]+sp1*n,i[j][4]-sp2*n))
						h.append( (i[j+1][1],i[j+1][5] ) )
						for n in range(int(k)): markPoints.append( (i[j][0],int(h[n][0]-100),int(h[n][0]+100),i[j][3],int(h[n][1]-100),int(h[n][1]+100),i[j][6],str(i[j][7])+'_gap') )
					elif gap1 < 500: pass
					else: markPoints.append( (i[j][0],i[j][2]+gap1/2-100,i[j][1]+gap1/2+100,i[j][3],i[j][4]-gap1/2-100,i[j][4]-gap1/2+100,i[j][6],str(i[j][7])+'_gap') )
				except IndexError: pass

	for i in markPoints: print >> f3,'%s\t%i\t%i\t%s\t%i\t%i\t%i\t%s' % (i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7])
	f3.close()
	
'''
def iPointviewContact(ContactsHash, covHash, ObjCoorMP, ChrIdxs, pointviews, out_name, **kwargs):
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
	printlog('\tstart processing %i real and %i potential interaction' % (len(ContactsHash),total), logname)
	out = open(out_name,'w')
	
	keyj,skip = [],[]
	for j in range(lnk):
		key1 = Keys[j]
		for pointview in pointviews:
			if pointviews[0] == 'from':
				for k1 in ObjCoorMP[key1]:
					if k1[0] == pointview[0] and pointview[1]<=k1[1]<pointview[1]:
						keyj.append(j)
						# print 'from',j, k1, Keys[j]
						break
			else:
				if key1[0] == pointview[0] and pointview[1]<=key1[1]<pointview[2]:
					keyj.append(j)
					# print 'to',j, key1, Keys[j]

	keyj.sort()
	for i in range(lnk):
		for j in keyj:
			if (i in keyj) and (j<i): pass
			else:
				num += 1
				key1,key2 = Keys[i],Keys[j]
				# print key1,key2
				if (key1[0] != key2[0]): dk = -1000
				else: dk = abs(key1[1]-key2[1])
				c = np.array([0.,0.,0.,0.])
				tst = 0
				for k1 in ObjCoorMP[key1]:
					for k2 in ObjCoorMP[key2]:
						if k1[0] == k2[0] or predict == False: real = True
						else: real = False
						if real:
							r = 1
							if (k1+k2) in ContactsHash: 
								tst += 1
								c += iDuplicateContact(ContactsHash[k1+k2], ObjCoorMP[key1][k1], ObjCoorMP[key2][k2])
							elif (k2+k1) in ContactsHash: 
								tst += 1
								c += iDuplicateContact(ContactsHash[k2+k1], ObjCoorMP[key1][k1], ObjCoorMP[key2][k2])
							else: real = False
							tst += 1
						if real == False and predict == True:
							r = 0
							try:
								pc,poe,mdk = ModelContacts(k1+k2, scoring, covHash, nullc, coef, cov_coef, rescale)
								c += iDuplicateContact([pc,poe,mdk], ObjCoorMP[key1][k1], ObjCoorMP[key2][k2])
							except KeyError: pass
						else: pass 
				if c[-1] > 0:
					try: scoring[dk]
					except KeyError: dk = max(scoring.keys())
					if model == 'balanced': norm, balance = 1.,c[-1]
					elif model == 'oe': norm, balance = scoring[dk], c[-1]
					elif model == 'shifted': norm,balance,mult = 1.,1.,mult*c[-1]
					elif model == 'easy':norm,balance = 1.,1.
					else: norm,balance = 1.,1.
					cc = round(c[data]*norm/balance,8)
					cl.append([key1[0], key1[1], key2[0], key2[1], cc, cc/scoring[dk], dk, r, c[0], c[1], c[2], norm, balance])
				else: pass
				if num % step == 0:
					try:
						cl = _randomizing(cl,regression,S, random)
						for c in cl: print >> out, '%s %i %s %i %.8f %.8f %i %i %.8f %.8f %.8f %.8f %.8f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12])
					except IndexError: pass
					cl = []
					elp = timeit.default_timer() - start_time
					printlog('\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
	if len(cl) != 0:
		cl = _randomizing(cl,regression,S, random)
		for c in cl: print >> out, '%s %i %s %i %.8f %.8f %i %i %.8f %.8f %.8f %.8f %.8f' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12])
		elp = timeit.default_timer() - start_time
		printlog('\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)
	out.close()
	elp = timeit.default_timer() - start_time
	printlog('\t%i interactions are lifovered, %.2fs.' % (total,elp), logname)
'''