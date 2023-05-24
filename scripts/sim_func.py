import math
import numpy as np
import timeit
import sys
import os
import global_func as gf

def readBedGraph(fname, resolution,ChrIdxs):
	bG = {}
	f = open(fname,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)-1,-1,-1):
		parse = lines[i].split()
		key = ChrIdxs[parse[0]], (int(parse[1])+int(parse[2]))//(2*resolution)
		bG[key] = float(parse[3])
		del lines[i]
	return bG


def readMeanHash(fname,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	psHash = {}
	gf.printlog('\tstart reading distance statistics '+fname, logname)
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
	gf.printlog('\t...end reading distance statistics %.2fs' % elp, logname)
	return psHash

def readCovHash(fname,ChrIdxs,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	covHash = {}
	gf.printlog('\tstart reading coverage statistics '+fname, logname)
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
	gf.printlog('\t...end reading coverage statistics %.2fs' % elp, logname)
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
	try: chrms = kwargs['chrms']
	except KeyError: chrms = sorted(ChrIdxs)
	try: scale = kwargs['scale']
	except KeyError: scale = 1
	start_time = timeit.default_timer()
	gf.printlog('\tstart reading contacts', logname)
	contactHash = {}
	files = os.listdir(path)
	lf,lc = len(files),0
	for i in range(lf):
		stf = timeit.default_timer()
		parse = files[i].split('.')
		gf.printlog('\t\topen %s, %s/%i' % (files[i],i,lf), logname)
		start_time2 = timeit.default_timer()
		f = open(path + '/' +files[i],'r')
		lines = f.readlines()
		f.close()
		elpf = timeit.default_timer() - stf
		gf.printlog('\t\t...close, %.2f' % (elpf) , logname)
		ln = len(lines)
		Keys = []
		for j in range(ln-1,0,-1):
			parse = lines[j].split()
			del lines[j]
			try:
				c1,b1,c2,b2 = ChrIdxs[parse[0]], int(parse[1])//scale, ChrIdxs[parse[2]], int(parse[3])//scale
				p,oe,mult_cov = float(parse[4]),float(parse[5]),float(parse[6])
				Keys.append((c1,b1,c2,b2))
				try: 
					contactHash[c1,b1,c2,b2][0] += p 
					contactHash[c1,b1,c2,b2][1].append(oe)
					contactHash[c1,b1,c2,b2][2] += mult_cov
				except KeyError: contactHash[c1,b1,c2,b2] = [p,[oe,],mult_cov]
				lc +=1
			except IndexError: print( 'iReadInitialContact IndexError', j, parse )
			except ValueError: print( 'iReadInitialContact ValueError', j, parse )
			if (ln-j) % 1000000 == 0:
				elpf = timeit.default_timer() - stf
				gf.printlog('\t\t\tcontact reading: %i, time elapsed: %.2fs' % (ln-j,elpf), logname)
		for key in Keys: contactHash[key][1] = np.mean(contactHash[key][1])
		del Keys
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\tfile is readed, %.2fs (%.2fs)' % (elpf,elp), logname)
	return contactHash

def iReadingMarkPoints(fname, resolution, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: ChrIdxs_from = kwargs['chrm_index_from']
	except KeyError:
		gf.printlog('iReadingMarkPoints Error! Dont exist control chromosome order (chrom_order_from)', logname)
		exit()
	try: ChrIdxs_to = kwargs['chrm_index_to']
	except KeyError: ChrIdxs_to = ChrIdxs_from
	try: chrms_from = kwargs['chrms_from']
	except KeyError: chrms = sorted(ChrIdxs_from)
	try: chrms_to = kwargs['chrms_to']
	except KeyError: chrms = sorted(ChrIdxs_to)
	ObjCoorMPH = {},{},{}
	ObjCoorMP = {}
	gf.printlog('\tstart reading markpoints '+fname, logname)
	start_time = timeit.default_timer()
	f = open(fname, 'r')
	lines = f.readlines()
	f.close()
	ln = len(lines)
	step = ln//10
	for i in range(ln-1,-1,-1):
		parse = lines[i].split()
		c1,p11,p12,c2,p21,p22 = parse[:6]
		if c1 in chrms_from and c2 in chrms_to:
			c1,c2 = ChrIdxs_from[c1], ChrIdxs_to[c2]
			p11,p12,p21,p22 = int(p11),int(p12),int(p21),int(p22)
			p1,p2 = (p11+p12)//2,(p21+p22)//2 
			b1 = c1,p1//resolution
			b2 = c2,p2//resolution
			try: ObjCoorMPH[0][b1].add(p1)
			except KeyError:ObjCoorMPH[0][b1] = set([p1])
			try: ObjCoorMPH[1][b2].add(p2)
			except KeyError:ObjCoorMPH[1][b2] = set([p2])
			try: ObjCoorMPH[2][b1+b2] += 1
			except KeyError:ObjCoorMPH[2][b1+b2] = 1
			if (ln - i) % step == 0:
				elp = timeit.default_timer() - start_time
				gf.printlog('\t\treaded %i/%i markpoins, %.2fs' % ((ln-i),ln,elp), logname)
		del lines[i]
	elp = timeit.default_timer() - start_time
	gf.printlog('\t...end reading markpoints %.2fs' % elp, logname)
	Keys = sorted(ObjCoorMPH[2])
	gf.printlog('\tstart remapping lowficient calculating from %i markpoints' % len(Keys), logname)
	for key in Keys:
		c0,c1,cx = len(ObjCoorMPH[0][key[:2]]),len(ObjCoorMPH[1][key[2:]]),ObjCoorMPH[2][key]
		try: ObjCoorMP[key[2:]][key[:2]] = 1.*cx/c1, 1.*cx/c0
		except KeyError: ObjCoorMP[key[2:]] = { key[:2]:(1.*cx/c1, 1.*cx/c0) }
		del ObjCoorMPH[2][key]
	del ObjCoorMPH
	elp = timeit.default_timer() - start_time
	gf.printlog('\t...end remapping lowficient calculating %.2fs' % elp, logname)
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

def _enumerateContacts(coor_to, contH, covH, meanH, mapH, res, params_ec): #model, nullc, pab_cont, pab_cov, pab_res
	c1,b1,c2,b2 = coor_to
	model, nullc, pab_cont, pab_cov, pab_res = params_ec
	pc,poe,mcov = 0.,0.,0.
	lifted_c = np.array([0.,0.,0.,0.])
	if (c1 != c2): real_dist = -1000
	else: real_dist = abs(b2-b1)
	real = 0
	try:
		for k1 in mapH[c1,b1]:
			for k2 in mapH[c2,b2]:
				if k1[0] == k2[0] or nullc == False: real = 1
				else: real = 0
				if real == 1:
					if (k1+k2) in contH:
						lifted_c += _recalculateContact(contH[k1+k2], mapH[c1,b1][k1], mapH[c2,b2][k2])
					elif (k2+k1) in contH:
						lifted_c += _recalculateContact(contH[k2+k1], mapH[c1,b1][k1], mapH[c2,b2][k2])
					else: real = 0
				if real == 0 and nullc:
					pc,poe,mcov = _predictContacts(k1+k2, covH, meanH, res, params_ec[1:])# nullc, pab_cont, pab_cov, pab_res
					lifted_c += _recalculateContact([pc,poe,mcov], mapH[c1,b1][k1], mapH[c2,b2][k2])
				else: pass
	except KeyError: pass
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
	x = c1, b1, c2, b2, cc, cc/meanH[real_dist][0], round(lifted_c[2]*norm/balance,8), lifted_c[0], lifted_c[1], lifted_c[2], real, meanH[real_dist][0], norm, balance
	return x 

def _predictContacts(coor_from, covH, meanH, res, params_pc):
	
	meanMCov = covH['mult_mean']
	c1,b1,c2,b2 = coor_from
	nullc, pab_cont, pab_cov, pab_res = params_pc
	
	if (c1 != c2): modeled_dist = -1000
	else: modeled_dist = abs(b1-b2)
	if meanH: score = meanH[modeled_dist]
	else: score = 1.0,1.0,1.0
	if pab_cont:
		rb1 = (b1*res+pab_res//2)//pab_res
		rb2 = (b2*res+pab_res//2)//pab_res
		if (c1,rb1,c2,rb2) in pab_cont: cont_AB,mult_AB = pab_cont[c1,rb1,c2,rb2][:2]
		elif (c2,rb2,c1,rb1) in pab_cont: cont_AB,mult_AB = pab_cont[c2,rb2,c1,rb1][:2]
		else: cont_AB,mult_AB = 1,0
	if nullc == 'mixed':
		if modeled_dist == -1000:
			poe = (covH[c1,b1]*covH[c2,b2])/meanMCov
			pc = mult_AB*poe*score[0]
		else:
			try: pc = cont_AB*(covH[c1,b1]*covH[c2,b2])/(pab_cov[c1,rb1]*pab_cov[c2,rb2])**2
			except ZeroDivisionError: pc = 0
			except KeyError: pc = 0
			try: poe = mult_AB*pc/score[0]
			except ZeroDivisionError: poe = 0
			except KeyError: poe = 0
	elif nullc == 'ab':
		try: pc = cont_AB*(covH[c1,b1]*covH[c2,b2])/(pab_cov[c1,rb1]*pab_cov[c2,rb2])**2
		except ZeroDivisionError: pc = 0
		try: poe = mult_AB*pc/score[0]
		except ZeroDivisionError: poe = 0
	elif nullc == 'sq_ab':
		poe = mult_AB
		try: pc = cont_AB*np.sqrt((covH[c1,b1]*covH[c2,b2])/(pab_cov[c1,rb1]*pab_cov[c2,rb2]))
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

def _rescaleContacts(high_cl, high_res, contH, covH, meanH, mapH, res, params_rand, params_rc): # params_rC: model, nullc, pab_cont, pab_cov, pab_res
	_cl = []
	contact_count,allCon,random = params_rand
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
			else: _clh = _randomizing(_clh, contact_count,allCon,random)
			_clh.sort()
			_cl.extend(_clh)
		except IndexError: pass
		except ValueError: pass
	return _cl

def iLiftOverContact(ContactsHash, covHash, ObjCoorMP, resolution, ChrIdxs, out_name, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: scoring = kwargs['scoring']
	except KeyError: scoring = False
	try: LowObjCoorMP = kwargs['markpoints_low']
	except KeyError: LowObjCoorMP = False
	try: cov_low = kwargs['coverage_low']
	except KeyError: cov_low = False
	try: scoring_low = kwargs['scoring_low']
	except KeyError: scoring_low = False
	try: contact_low = kwargs['contact_low']
	except KeyError:
		contact_low = False
		cov_low = False
		scoring_low = False
		LowObjCoorMP = False
	try: low_res = kwargs['resolution_low']
	except KeyError: low_res = 1
	try: pab_cov = kwargs['coverage_pab']
	except KeyError: pab_cov = False
	try: pab_cont = kwargs['contact_pab']
	except KeyError: pab_cont,pab_cov = False,False
	try: pab_res = kwargs['resolution_pab']
	except KeyError: pab_res = 1
	try: model = kwargs['model']
	except KeyError: model = 'easy'
	try: random = kwargs['random']
	except KeyError: random = False
	try: nullc = kwargs['null_contacts']
	except KeyError: nullc = False
	try: pointviews = kwargs['pointviews']
	except KeyError: pointviews = False
	
	DifferContact = {}
	start_time = timeit.default_timer()
	gf.printlog('\tstart liftovering', logname)
	elp = timeit.default_timer() - start_time
	gf.printlog('\taccount number of interactions and contacts, %.2f' % elp, logname)
	aKeys = sorted(ObjCoorMP)
	alnk = len(aKeys)
	total,num,bKeys,blnk,type = (alnk+1)*alnk//2,0,aKeys,alnk,1
	if LowObjCoorMP:
		LowKeys = sorted(LowObjCoorMP)
		low_lnk = len(LowKeys)
		total = (low_lnk+1)*low_lnk//2
	allCon,meanMCov = covHash['all'],covHash['mult_mean']
	try: contact_count = kwargs['contact_count']
	except KeyError: contact_count = allCon
	if contact_count == 0: contact_count = allCon
	elp = timeit.default_timer() - start_time
	gf.printlog('\tmodel %s' % model, logname)
	gf.printlog('\trandomize %s' % random, logname)
	gf.printlog('\tpredict null contact %s' % nullc, logname)
	if contact_low and low_res: gf.printlog('\tadditional low reading with rescale %i %i' % (resolution,low_res), logname)
	else: gf.printlog('\tno additional low reading', logname)
	if contact_count != allCon: gf.printlog('\t%i contacts are regressed to %i, %.4f, %.2f' % (allCon,contact_count, 1.*contact_count/allCon,elp), logname)
	else: gf.printlog('\tno contact_count, %.2f' % (elp), logname)
	out = open(out_name,'w')
	
	params_model = model, nullc, pab_cont, pab_cov, pab_res
	params_random = contact_count, allCon, random
	pKeys = set([])
	if pointviews: 
		s = ''
		for p in pointviews: s += ('%s %i %i ' % (p[0],p[1]*low_res,(p[1]+1)*low_res))
		gf.printlog('\tPointviews: %s'% s, logname) 
	else: gf.printlog('\tNo pointviews', logname) 
	
	if LowObjCoorMP:
		low_cl = []
		gf.printlog('\tstart processing %i real and %i potential interaction by model rescaling' % (len(contact_low),total), logname)
		for i in range(low_lnk):
			for j in range(i,low_lnk):
				num += 1
				key1,key2 = LowKeys[i],LowKeys[j]
				if set(LowObjCoorMP[key1].keys()) & set(pointviews): pKeys.add(key1)
				elif set(LowObjCoorMP[key2].keys()) & set(pointviews): pKeys.add(key2)
				else:
					low_cl.append( _enumerateContacts(key1+key2, contact_low, cov_low, scoring_low, LowObjCoorMP, low_res, params_model ) )#, model, nullc, pab_cont, pab_cov, pab_res) )
					if num % 2000000 == 0:
						try: low_cl = _randomizing(low_cl, contact_count, allCon, random)
						except IndexError: pass
						cl = _rescaleContacts(low_cl, low_res, ContactsHash, covHash, scoring, ObjCoorMP, resolution, params_random, params_model)#, model, nullc, pab_cont, pab_cov, pab_res, contact_count, random)
						for c in cl:
							out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13]) )
						low_cl,cl,hash = [],[],{}
						elp = timeit.default_timer() - start_time
						gf.printlog('\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
		if len(low_cl) != 0:
			try: low_cl = _randomizing(low_cl, contact_count, allCon, random)
			except IndexError: pass
			cl = _rescaleContacts(low_cl, low_res, ContactsHash, covHash, scoring, ObjCoorMP, resolution, params_random, params_model)#, model, nullc, pab_cont, pab_cov, pab_res, contact_count, random)
			for c in cl:
				out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13]) )
			low_cl,cl,hash = [],[],{}
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)
		bKeys = []
		for key in pKeys:
			for i in range(key[1]*low_res//resolution,(key[1]+1)*low_res//resolution): bKeys.append((key[0],i))
		blnk,type = len(bKeys),0
		bKeys.sort()
		aKeys = list(set(aKeys) - set(bKeys))
		aKeys.sort()
		alnk = len(aKeys)
		aKeys += bKeys

	if LowObjCoorMP == False or pKeys:
		if pKeys: gf.printlog('\tPointview edges modelling', logname)
		else: gf.printlog('\tModelling without rescaling', logname)
		cl = []
		gf.printlog('\tstart processing %i real and %i potential interaction' % (len(ContactsHash),total), logname)
		for i in range(alnk+blnk*(1-type)):
			j0 = max((type*i),(i-alnk))
			for j in range(j0,blnk):
				num += 1
				key1,key2 = aKeys[i],bKeys[j]
				cl.append( _enumerateContacts(key1+key2, ContactsHash, covHash, scoring, ObjCoorMP, resolution, params_model) )#, model, nullc, pab_cont, pab_cov, pab_res))
				if num % 2000000 == 0:
					try: cl = _randomizing(cl, contact_count, allCon, random)
					except IndexError: pass
					for c in cl:
						out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13]) )
					cl = []
					elp = timeit.default_timer() - start_time
					gf.printlog('\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
		if len(cl) != 0:
			cl = _randomizing(cl, contact_count, allCon, random)
			for c in cl: 
				out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13]) )
			cl = []
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)
	out.close()
	elp = timeit.default_timer() - start_time
	gf.printlog('\t%i interactions are lifovered, %.2fs.' % (total,elp), logname)

def _listing(coor, contactHash, covH, meanH, res, params_list):
	
	c1,b1,c2,b2 = coor
	nullc, pab_cont, pab_cov, pab_res = params_pc = params_list
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
	contact_count, allCon, random = params_rand
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
					_clh.append((c1,ri,c2,rj,1.,1.,( covH[c1,ri]*covH[c2,rj] ),mean_con))
				except KeyError: pass
				except ZeroDivisionError: pass
		try:
			if random:
				_clh = np.array(_clh)
				p = 1.*_clh[:,-2]*_clh[:,-1]
				p = p/np.sum(p)
				_clh = _randomizing(_clh, count, p, 'choice')
			else: _clh = _randomizing(_clh, contact_count, allCon, random)
		except IndexError: pass#
		except ValueError: pass
		_cl.extend(_clh)
		_clh = []
	del _clh
	return _cl

def iContactRegression(covHash,resolution,chosen_chroms,ChrIdxs,low_ChrSizes,out_name,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: scoring = kwargs['scoring']
	except KeyError: scoring = False
	try: LowObjCoorMP = kwargs['markpoints_low']
	except KeyError: LowObjCoorMP = False
	try: cov_low = kwargs['coverage_low']
	except KeyError:
		cov_low = False
		LowObjCoorMP = False
	try: scoring_low = kwargs['scoring_low']
	except KeyError:
		scoring_low = False
		LowObjCoorMP = False
	try: contact_low = kwargs['contact_low']
	except KeyError:
		contact_low = False
		cov_low = False
		scoring_low = False
		LowObjCoorMP = False
	try: low_res = kwargs['resolution_low']
	except KeyError: low_res = 1
	try: pab_cov = kwargs['coverage_pab']
	except KeyError: pab_cov = False
	try: pab_cont = kwargs['contact_pab']
	except KeyError: pab_cont,pab_cov = False,False
	try: pab_res = kwargs['resolution_pab']
	except KeyError: pab_res = 1
	try: model = kwargs['model']
	except KeyError: model = 'easy'
	try: random = kwargs['random']
	except KeyError: random = False
	try: nullc = kwargs['null_contacts']
	except KeyError: nullc = False
	allCon,meanMCov = covHash['all'],covHash['mult_mean']
	try: contact_count = kwargs['contact_count']
	except KeyError: contact_count = allCon
	if contact_count == 0: contact_count = allCon
	total = 0
	params_model = nullc, pab_cont, pab_cov, pab_res
	params_random = contact_count, allCon, random
	chroms = chosen_chroms.split(',')
	if len(chroms) == 1: chroms *= 2
	for ci in range(len(chroms)):
		for cj in range((ci+1),len(chroms)):
			chrm1,chrm2 = chroms[ci],chroms[cj]
			if chrm1 == chrm2: total += (low_ChrSizes[chrm1]+1)*low_ChrSizes[chrm1]//2
			else: total += low_ChrSizes[chrm1]*low_ChrSizes[chrm2]
	num,H,cl = 0, [], []
	start_time = timeit.default_timer()
	elp = timeit.default_timer() - start_time
	gf.printlog('\t%i contacts are regressed to %i, %.2fs.' % (allCon,contact_count,elp), logname)
	for ci in range(len(chroms)):
		for cj in range((ci+1),len(chroms)):
			chrm1,chrm2 = chroms[ci],chroms[cj]
			full_name = '%s.%s.%s.allCon'% (out_name,chrm1,chrm2)
			out = open( full_name,'w')
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\tcontact_count %s %s chromosome pair, %.2fs.' % (chrm1,chrm2,elp), logname)
			if chrm1 == chrm2: k = 1
			else: k = 0
			for i in range(low_ChrSizes[chrm1]):
				j1 = i*k
				for j in range(j1,low_ChrSizes[chrm2]):
					num += 1.
					c1,b1,c2,b2 = ChrIdxs[chrm1],i,ChrIdxs[chrm2],j
					H.append(_listing((c1,b1,c2,b2), contact_low, cov_low, scoring_low, low_res, params_model))
					if num % 200000 == 0:
						elp = timeit.default_timer() - start_time
						gf.printlog('\t\t%s randomization, %.2fs.' % (random,elp), logname)
						H = _randomizing(H, contact_count, allCon, random)
						H = _easyRescale(H, low_res, covHash, scoring, resolution, params_random)
						elp = timeit.default_timer() - start_time
						gf.printlog('\t\twriting, %.2fs.' % (elp), logname)
						for c in H:
							try: out.write('%s\t%i\t%s\t%i\t%.4f\n' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4]))
							except IndexError: pass
						elp = timeit.default_timer() - start_time
						gf.printlog('\t\t%i interactions from %i are processed, %.2fs.' % (num,total,elp), logname)
						H = []
			if len(H) > 0:
				H = _randomizing(H, contact_count, allCon, random)
				H = _easyRescale(H, low_res, covHash, scoring, resolution, params_random)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t\twriting, %.2fs.' % (elp), logname)
				for c in H:
					try: out.write('%s\t%i\t%s\t%i\t%.4f\n' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4]))
					except IndexError: pass
			del H
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\tend contact_count %s %s chromosome pair, %.2fs.' % (chrm1,chrm2,elp), logname)
			out.close()
	elp = timeit.default_timer() - start_time
	gf.printlog('\t%i interactions are randomized, %.2fs.' % (total,elp), logname)

def AddNormVector(path,Rlist,hic,norm):
	parse = Rlist.split(',')
	f = open(hic+'.norm','w')
	for r in parse:
		resolution = int(r)
		chrSizes = gf.ChromSizes(path,resolution)
		Keys = sorted(chrSizes)
		for key in Keys:
			try:
				f.write( 'vector ' + norm + ' ' + key + ' ' + r + ' BP\n' )
				for i in range(chrSizes[key]): f.write( '1.0\n' )
			except TypeError: pass
			