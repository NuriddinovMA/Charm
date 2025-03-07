import os
import sys
import timeit
from charm_func import global_func as gf

try: import numpy as np
except ModuleNotFoundError:
	print('Lethal Error! NumPy not found!')
	exit()

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
	gf.printlog('\t\t\tstart reading distance statistics '+fname, logname)
	start_time = timeit.default_timer()
	f = open(fname,'r')
	lines = f.readlines()
	f.close()
	Keys = lines[0].split()
	lnk = len(Keys)
	print(Keys)
	for i in range(1,len(lines)):
		parse = lines[i].split()
		if len(parse[2:]) == 0: 
			psHash = False
			break
		else:
			psHash[int(parse[0])] = { Keys[i]:float(parse[i]) for i in range(2,lnk) }
			# try: 
			# except ValueError: psHash[-1000] = k
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t...end reading distance statistics %.2fs' % elp, logname)
	return psHash

def readCovHash(fname,ChrIdxs,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: treshold = kwargs['treshold']
	except KeyError: treshold = 0
	print(treshold)
	covHash = {}
	gf.printlog('\t\t\tstart reading coverage statistics '+fname, logname)
	start_time = timeit.default_timer()
	with open(fname,'r') as f: lines = f.readlines()
	for i in range(1,len(lines)):
		parse = lines[i].split()
		coverage = float(parse[2])
		if coverage <= treshold: coverage = 0
		try: covHash[ChrIdxs[parse[0]],int(parse[1])] = coverage
		except KeyError: pass
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t...end reading coverage statistics %.2fs' % elp, logname)
	return covHash

def readTotalStatHash(fname,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	totalHash = {}
	gf.printlog('\t\t\tstart reading distance independing statistics '+fname, logname)
	start_time = timeit.default_timer()
	with open(fname,'r') as f: lines = f.readlines()
	for line in lines:
		parse = line.split()
		totalHash[parse[0]] = float(parse[1])
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t...end reading coverage statistics %.2fs' % elp, logname)
	return totalHash

def iReadInitialContact(path,ChrIdxs,**kwargs): #Reading Contact from database file
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: chrms = kwargs['chrms']
	except KeyError: chrms = sorted(ChrIdxs)
	try: scale = kwargs['scale']
	except KeyError: scale = 1
	try: index = kwargs['index']
	except KeyError: index = 4
	start_time = timeit.default_timer()
	gf.printlog('\t\t\tstart reading contacts', logname)
	contactHash = {}
	files = os.listdir(path)
	lf,lc = len(files),0

	for i in range(lf):
		parse = files[i].split('.')
		Keys = []
		if ([parse[-2],parse[-3]] in chrms) or ([parse[-3],parse[-2]] in chrms):
			stf = timeit.default_timer()
			gf.printlog('\t\t\t\topen %s, %s/%i' % (files[i],i,lf), logname)
			f = open(path + '/' +files[i],'r')
			lines = f.readlines()
			f.close()
			elpf = timeit.default_timer() - stf
			gf.printlog('\t\t\t\t...close, %.2f' % (elpf) , logname)
			ln = len(lines)
			
			for j in range(ln-1,0,-1):
				parse = lines[j].split()
				del lines[j]
				try:
					c1,b1,c2,b2 = ChrIdxs[parse[0]], int(parse[1])//scale, ChrIdxs[parse[2]], int(parse[3])//scale
					p,oe = float(parse[index]),float(parse[5])
					Keys.append((c1,b1,c2,b2))
					try: 
						contactHash[c1,b1,c2,b2][0] += p 
						contactHash[c1,b1,c2,b2][1].append(oe)
					except KeyError: contactHash[c1,b1,c2,b2] = [p,[oe,]]
					except AttributeError: print('iReadInitialContact AttributeError', parse, c1,b1,c2,b2, contactHash[c1,b1,c2,b2],p,oe)
					lc +=1
				except IndexError: print( 'iReadInitialContact IndexError', j, parse )
				except ValueError: print( 'iReadInitialContact ValueError', j, parse )
				if (ln-j) % 1000000 == 0:
					elpf = timeit.default_timer() - stf
					gf.printlog('\t\t\t\t\tcontact reading: %i, time elapsed: %.2fs' % (ln-j,elpf), logname)
			elpf = timeit.default_timer() - stf
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t\t\tfile is readed, %.2fs (%.2fs)' % (elpf,elp), logname)
		for key in Keys: contactHash[key][1] = np.mean(contactHash[key][1])
		del Keys
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
	except KeyError: chrms_from = ChrIdxs_from
	try: chrms_to = kwargs['chrms_to']
	except KeyError: chrms_to = ChrIdxs_to
	ObjCoorMPH = {},{},{}
	ObjCoorMP = {}
	chosen_from = set([])
	gf.printlog('\t\t\tstart reading markpoints '+fname, logname)
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
			chosen_from.add(c1)
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
				gf.printlog('\t\t\t\treaded %i/%i markpoins, %.2fs' % ((ln-i),ln,elp), logname)
		del lines[i]
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t...end reading markpoints %.2fs' % elp, logname)
	Keys = sorted(ObjCoorMPH[2])
	gf.printlog('\t\t\tstart remapping lowficient calculating from %i markpoints' % len(Keys), logname)
	for key in Keys:
		c0,c1,cx = len(ObjCoorMPH[0][key[:2]]),len(ObjCoorMPH[1][key[2:]]),ObjCoorMPH[2][key]
		try: ObjCoorMP[key[2:]][key[:2]] = 1.*cx/c1, 1.*cx/c0
		except KeyError: ObjCoorMP[key[2:]] = { key[:2]:(1.*cx/c1, 1.*cx/c0) }
		del ObjCoorMPH[2][key]
	del ObjCoorMPH
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t...end remapping lowficient calculating %.2fs' % elp, logname)
	chosen_from = list(chosen_from)
	chosen_from.sort()
	return ObjCoorMP,chosen_from

def createUntouchedMarkPoints(chroms,chrSizes,ChrIdxs_to,ChrIdxs_from):
	UnCoorMP = {}
	for chrom in chroms:
		for i in range(chrSizes[chrom]): UnCoorMP[(ChrIdxs_to[chrom],i)] = {(ChrIdxs_from[chrom],i):(1.,1.)}
	return UnCoorMP

def _recalculateContact(ContactsHash_L, relation1, relation2):
	
	rc = [0,0,0,0,0]
	rc[-1] = 1.*relation1[0]*relation2[0]*relation1[1]*relation2[1]
	rc[0] = rc[-1]*ContactsHash_L[0]
	rc[1] = rc[-1]*ContactsHash_L[1]
	rc[2] = rc[-1]*ContactsHash_L[2]
	rc[3] = rc[-1]*ContactsHash_L[3]
	return rc

def _randomizing(h, old_counts, new_counts, random_func,**kwargs):
	np.random.seed()
	h = np.array(h)
	h = h[np.isfinite(h[:,4])]
	if random_func: h[:,4] = random_func(h[:,4], old_counts, new_counts)
	h = h[h[:,4] > 0]
	h = h.tolist() 
	return h

def _enumerateContacts(coor_to, contH, covHash, totalH, distanceH, mapH1, mapH2, res, params_ec,**kwargs): #model, pred_null_model, pab_cont, pab_cov, pab_res
	try: out = kwargs['out']
	except KeyError: out = False
	
	c1,b1,c2,b2 = coor_to
	model, noised, pred_null_model, pab_cont, pab_cov, pab_res = params_ec
	
	pc,poe,cov1,cov2,real = 0.,0.,0.,0.,0
	lifted_c = np.array([0.,0.,0.,0.,0.,])
	if (c1 != c2): modeled_dist = -1000
	else: modeled_dist = abs(b2-b1)
	try: distanceH[modeled_dist]
	except KeyError: modeled_dist = max(list(distanceH.keys()))
	if noised: noise = 1
	else: noise = 0
	try:
		for k1 in mapH1[c1,b1]:
			for k2 in mapH2[c2,b2]:
				try: mcov = covHash[k1]*covHash[k2]
				except TypeError: mcov = 1
				cont = 0,0,-1
				if k1[0] == k2[0] or pred_null_model == False: real = 1
				else: real = 0
				if (k1[0] != k2[0]): from_dist = -1000
				else: from_dist = abs(k1[1]-k2[1])
				if real == 1:
					if (k1+k2) in contH and contH[k1+k2][0] > noise:
						try: cont = contH[k1+k2][0],contH[k1+k2][1],covHash[k1],covHash[k2]
						except TypeError: cont = contH[k1+k2][0],contH[k1+k2][1],0,0
						lifted_c += _recalculateContact(cont, mapH1[c1,b1][k1], mapH2[c2,b2][k2])
					elif (k2+k1) in contH and contH[k2+k1][0] > noise:
						try: cont = contH[k2+k1][0],contH[k2+k1][1],covHash[k2],covHash[k1]
						except TypeError: cont = contH[k2+k1][0],contH[k2+k1][1],0,0
						lifted_c += _recalculateContact(cont, mapH1[c1,b1][k1], mapH2[c2,b2][k2])
					else: real = 0
				if real == 0 and mcov > 0 and pred_null_model:
					cont = _predictContacts(k1+k2, modeled_dist, covHash, totalH, distanceH, res, params_ec[2:], out=out)# pred_null_model, pab_cont, pab_cov, pab_res
					lifted_c += _recalculateContact(cont, mapH1[c1,b1][k1], mapH2[c2,b2][k2])
					real == 2
				else: pass
				
	except KeyError: pass
	data, norm, balance = 0, -1.,-1.
	try: distanceH[modeled_dist]
	except KeyError: modeled_dist = max(distanceH.keys())
	if lifted_c[-1] != 0:
		if model == 'balanced': data,norm, balance = 1, distanceH[modeled_dist]['mean_contact'], lifted_c[-1]
		elif model == 'align_sensitive': data,norm, balance = 1, distanceH[modeled_dist]['mean_contact'], 1.
		elif model == 'distance_sensitive': data, norm, balance = 0, 1.,lifted_c[-1]
		elif model == 'easy': data, norm, balance = 0, 1.,1.
		else: pass
		cc = round(lifted_c[data]*norm/balance,8)
	else: cc = 0
	x = c1, b1, c2, b2, cc, cc/distanceH[modeled_dist]['mean_contact'], round(lifted_c[2]*norm/balance,8), lifted_c[0], lifted_c[1], lifted_c[2], lifted_c[3], real, distanceH[modeled_dist]['mean_contact'], modeled_dist, norm, balance
	return x

def _predictContacts(coor_from, modeled_dist, covHash, totalH, distanceH, res, params_pc,**kwargs):
	try: out = kwargs['out']
	except KeyError: out = False

	c1,b1,c2,b2 = coor_from
	pred_null_model, pab_cont, pab_cov, pab_res_list = params_pc
	cont_AB,oe_AB = 1.,1.
	if (c1 != c2): from_dist = -1000
	else: from_dist = abs(b2-b1)
	if distanceH: distance_stats = distanceH[modeled_dist]
	else: distance_stats = 1.0,1.0,1.0,1.0,1.0,1.0,1.0

	if pab_cont:
		pab_res = pab_res_list[-1]
		for i in range(len(pab_res_list)):
			k = (1-distanceH[from_dist]['contacts=0']/distanceH[from_dist]['all_contacts'])*pab_res_list[i]/res
			if k > 0.5 and pab_res_list[i] >= res:
				pab_res = pab_res_list[i]
				break
		pab_cont = pab_cont[pab_res]
		pab_cov = pab_cov[pab_res]
		rb1 = (b1*res-pab_res//2)//pab_res
		rb2 = (b2*res-pab_res//2)//pab_res
		if rb1 < 0: rb1 = 0
		if rb2 < 0: rb2 = 0
		try: cont_AB,oe_AB = pab_cont[c1,rb1,c2,rb2][:2]
		except KeyError:
			try: cont_AB,oe_AB = pab_cont[c2,rb2,c1,rb1][:2]
			except KeyError: cont_AB,oe_AB = 1,0
		k = (2*pab_res/res)
		cov1 = pab_cov[c1,rb1]
		try: cov1 += pab_cov[c1,rb1+1]
		except KeyError: pass
		cov2 = pab_cov[c2,rb2]
		try: cov2 += pab_cov[c2,rb2+1]
		except KeyError: pass
	pc,poe = pred_null_model(covHash[c1,b1], covHash[c2,b2], cont_AB, oe_AB, modeled_dist, distanceH, totalH)
	modelled_cont = pc, poe, covHash[c1,b1], covHash[c2,b2]
	return modelled_cont

def _rescaleContacts(high_cl, high_res, contH, covHash, totalH, distanceH, mapH1, mapH2,res, params_rand, params_rc,**kwargs): # params_rC: model, pred_null_model, pick_extra, pick_intra, pab_cont, pab_cov, pab_res
	try: out = kwargs['out']
	except KeyError: out = False

	_cl = []
	contact_count, allCon, random_func = params_rand
	for hi in high_cl:
		c1,b1,c2,b2,count,oe = hi[:6]
		ri1,ri2,rj1,rj2 = int(b1*high_res/res),int((b1+1)*high_res/res),int(b2*high_res/res),int((b2+1)*high_res/res)
		_s,_clh = [],[]
		for ri in range(ri1,ri2):
			if (c1,b1) == (c2,b2): rj1 = ri
			for rj in range(rj1,rj2):
				try:
					mapH1[c1,ri],mapH2[c2,rj]
					y = _enumerateContacts((c1,ri,c2,rj), contH, covHash, totalH, distanceH, mapH1, mapH2, res, params_rc[:6], out=out)
					_clh.append(y)
				except KeyError: pass
		_clh = np.array(_clh)
		try:
			if random_func:
				k = _clh[:,9]*_clh[:,10]
				k[k>0] = 1.
				pick_probality_func = params_rc[6]
				#print(k,pick_probality_func.__name__)
				#print('!',k*_clh[:,9],k*_clh[:,10],_clh[:,12],'!')
				p = pick_probality_func(k*_clh[:,9],k*_clh[:,10],_clh[:,12],contact_distance=_clh[:,13])
				p = p/np.sum(p)
				#print('1',_clh[:,4],p,count)
				_clh = _randomizing(_clh, p, count, gf._choice)
				#print('4',np.array(_clh)[:,4])
				#raise NameError()
			else: _clh = _randomizing(_clh, allCon, contact_count, False)
			_clh.sort()
			_cl.extend(_clh)
		except IndexError: pass
		except ValueError: pass
	return _cl

def iLiftOverContact(ContactsHash, covHash, ObjCoorMP, resolution, ChrIdxs, out_name, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: distance_dependence = kwargs['distance_dependence']
	except KeyError: distance_dependence = False
	try: total_statistics = kwargs['total_statistics']
	except KeyError: total_statistics = False
	try: LowObjCoorMP = kwargs['markpoints_low']
	except KeyError: LowObjCoorMP = False
	try: cov_low = kwargs['coverage_low']
	except KeyError: cov_low = False
	try: distance_dependence_low = kwargs['distance_dependence_low']
	except KeyError: distance_dependence_low = False
	try: total_statistics_low = kwargs['total_statistics_low']
	except KeyError: total_statistics_low = False
	try: contact_low = kwargs['contact_low']
	except KeyError:
		contact_low = False
		cov_low = False
		distance_dependence_low = False
		LowObjCoorMP = False
	try: low_res = int(kwargs['resolution_low'])
	except KeyError: low_res = False
	try: pab_cov = kwargs['coverage_pab']
	except KeyError: pab_cov = False
	try: pab_cont = kwargs['contact_pab']
	except KeyError: pab_cont,pab_cov = False,False
	try: pab_res = kwargs['resolution_pab']
	except KeyError: pab_res = False
	try: model = kwargs['model']
	except KeyError: model = 'easy'
	try: random_func = kwargs['random_func']
	except KeyError: random_func = False
	
	try: pred_null_model = kwargs['predict_null_contacts']
	except KeyError: pred_null_model = False
	try: pick_contacts = kwargs['pick_contacts']
	except KeyError: pick_contacts = False
	
	try: noised = kwargs['noised']
	except KeyError: noised = False
	try: pointviews = kwargs['pointviews']
	except KeyError: pointviews = False
	try: untouched_low = kwargs['untouched_low']
	except KeyError: untouched_low = False
	try: untouched = kwargs['untouched']
	except KeyError:
		untouched = False
		untouched_low = False
	if LowObjCoorMP and untouched_low == False: untouched = False
	
	DifferContact = {}
	start_time = timeit.default_timer()
	gf.printlog('\t\t\tstart liftovering', logname)
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\taccount number of interactions and contacts, %.2f' % elp, logname)
	aKeys = sorted(ObjCoorMP)
	alnk = len(aKeys)
	total,num = (alnk+1)*alnk//2,0
	
	if LowObjCoorMP:
		LowKeys = sorted(LowObjCoorMP)
		low_lnk = len(LowKeys)
		LowKeys2 = LowKeys
		low_lnk2 = low_lnk
		LowObjCoorMP2 = LowObjCoorMP
		total = (low_lnk+1)*low_lnk//2
		ii = 1
		if untouched_low:
			LowKeys2 = sorted(untouched_low)
			low_lnk2 = len(LowKeys2)
			total = low_lnk2*low_lnk
			LowObjCoorMP2 = untouched_low
			ii = 0
	
	if untouched: ObjCoorMP2 = untouched
	else: ObjCoorMP2 = ObjCoorMP
	
	try: allCon,meanCov,meanMCov = total_statistics['all_contacts'],total_statistics['mean_coverage_sum'],total_statistics['mean_coverage_multiple']
	except TypeError: allCon,meanCov,meanMCov = -1,-1,-1
	try: contact_count = kwargs['contact_count']
	except KeyError: contact_count = allCon
	if contact_count: contact_count = int(contact_count) 
	else: contact_count = allCon
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\tmodel %s' % model, logname)
	try: gf.printlog('\t\t\trandomize %s' % random_func.__name__, logname)
	except AttributeError: gf.printlog('\t\t\tno randomization', logname)
	try: gf.printlog('\t\t\tpredict null contact %s' % pred_null_model.__name__, logname)
	except AttributeError: gf.printlog('\t\t\tno prediction null contacts', logname)
	try: gf.printlog('\t\t\tpick contacts %s' % pick_contacts.__name__, logname)
	except AttributeError: gf.printlog('\t\t\tno picking contcts', logname)
	if contact_low and low_res: gf.printlog('\t\t\tadditional low reading with rescale %i %i' % (resolution,low_res), logname)
	else: gf.printlog('\t\t\tno additional low reading', logname)
	if contact_count != allCon: gf.printlog('\t\t\t%i contacts are regressed to %i, %.4f, %.2f' % (allCon,contact_count, 1.*contact_count/allCon,elp), logname)
	else: gf.printlog('\t\t\tno contact regression, %.2f' % (elp), logname)
	out = open(out_name,'w')
	params_model = model, noised, pred_null_model, pab_cont, pab_cov, pab_res, pick_contacts
	params_random = contact_count, allCon, random_func
	pKeys = set([])
	if pointviews: 
		s = ''
		for p in pointviews: s += ('%s %i %i ' % (p[0],p[1],(p[1]+1)))
		gf.printlog('\t\t\tPointviews: %s'% s, logname) 
	else: gf.printlog('\t\t\tNo pointviews', logname) 

	# print(list(LowObjCoorMP.keys())[:10])
	# print(list(LowObjCoorMP2.keys())[:10])
	if LowObjCoorMP:
		low_cl = []
		gf.printlog('\t\t\tstart processing %i real and %i potential interaction by model rescaling' % (len(contact_low),total), logname)
		for i in range(low_lnk):
			for j in range(i*ii,low_lnk2):
				num += 1
				key1,key2 = LowKeys[i],LowKeys2[j]
				#print(LowObjCoorMP[key1],LowObjCoorMP2[key2])
				if set(LowObjCoorMP[key1].keys()) & set(pointviews): pKeys.add(key1)
				elif set(LowObjCoorMP2[key2].keys()) & set(pointviews): pKeys.add(key2)
				else:
					
					low_cl.append( _enumerateContacts(key1+key2, contact_low, cov_low, total_statistics_low, distance_dependence_low, LowObjCoorMP, LowObjCoorMP2, low_res, 
						params_model[:6]) )#, out=True ) )#, model, pred_null_model, pab_cont, pab_cov, pab_res) )
					if num % 2000000 == 0:
						if model == 'easy': print('.',end='')
						try: low_cl = _randomizing(low_cl, allCon, contact_count, random_func)
						except IndexError: print('IndexError')#pass
						cl = _rescaleContacts(low_cl, low_res, ContactsHash, covHash, total_statistics, distance_dependence, ObjCoorMP, ObjCoorMP2, resolution, 
							params_random, params_model) #, out=False)#, model, pred_null_model, pab_cont, pab_cov, pab_res, contact_count, random)
						for c in cl:
							if (ChrIdxs[c[0]] < ChrIdxs[c[2]]) or ((ChrIdxs[c[0]] == ChrIdxs[c[2]]) and (c[1] < c[3])):
								out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
							else:
								out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[2]], c[3], ChrIdxs[c[0]], c[1], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
						low_cl,cl,hash = [],[],{}
					if num % 200000 == 0:
						elp = timeit.default_timer() - start_time
						gf.printlog('\t\t\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
		#print('end1',len(low_cl))
		if len(low_cl) != 0:
			try: low_cl = _randomizing(low_cl, allCon, contact_count, random_func)
			except IndexError: print('IndexError')#pass
			cl = _rescaleContacts(low_cl, low_res, ContactsHash, covHash, total_statistics, distance_dependence, ObjCoorMP, ObjCoorMP2, resolution,
				params_random, params_model) #, out=True)#, model, pred_null_model, pab_cont, pab_cov, pab_res, contact_count, random)
			#print('end2',len(low_cl),cl)
			for c in cl:
				if (ChrIdxs[c[0]] < ChrIdxs[c[2]]) or ((ChrIdxs[c[0]] == ChrIdxs[c[2]]) and (c[1] < c[3])):
					out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
				else:
					out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[2]], c[3], ChrIdxs[c[0]], c[1], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
			low_cl,cl,hash = [],[],{}
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)

	if LowObjCoorMP == False: bKeys,blnk = aKeys,alnk

	if pKeys:
		bKeys = []
		for key in pKeys:
			for i in range(key[1]*low_res//resolution,(key[1]+1)*low_res//resolution): bKeys.append((key[0],i))
		blnk = len(bKeys)
		bKeys.sort()
		aKeys = list(set(aKeys) - set(bKeys))
		aKeys.sort()
		aKeys = bKeys + aKeys
		alnk = len(aKeys)
	
	if untouched:
		aKeys = sorted(untouched)
		alnk = len(aKeys)
	
	if LowObjCoorMP == False or pKeys:
		if pKeys: gf.printlog('\t\t\tPointview edges modelling', logname)
		else: gf.printlog('\t\t\tModelling without rescaling', logname)
		cl = []
		gf.printlog('\t\t\tstart processing %i real and %i potential interaction' % (len(ContactsHash),total), logname)
		#params_model = model, pred_null_model, contact_low, cov_low, low_res
		for i in range(blnk):
			for j in range(i,alnk):
				num += 1
				key1,key2 = bKeys[i],aKeys[j]
				cl.append( _enumerateContacts(key1+key2, ContactsHash, covHash, total_statistics, distance_dependence, ObjCoorMP, ObjCoorMP2, resolution, 
					params_model[:6], out=True) )#, model, pred_null_model, pab_cont, pab_cov, pab_res))
				if num % 2000000 == 0:
					cl = _randomizing(cl, allCon, contact_count, random_func)
					for c in cl:
						if (ChrIdxs[c[0]] < ChrIdxs[c[2]]) or ((ChrIdxs[c[0]] == ChrIdxs[c[2]]) and (c[1] < c[3])):
							out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
						else:
							out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[2]], c[3], ChrIdxs[c[0]], c[1], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
					cl = []
				if num % 200000 == 0:
					elp = timeit.default_timer() - start_time
					gf.printlog('\t\t\t\t %.2f from %i interactions are processed, %.2fs.' % (100.*num/total,total,elp), logname)
		if len(cl) != 0:
			cl = _randomizing(cl, allCon, contact_count, random_func)
			for c in cl: 
				if (ChrIdxs[c[0]] < ChrIdxs[c[2]]) or ((ChrIdxs[c[0]] == ChrIdxs[c[2]]) and (c[1] < c[3])):
					out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[0]], c[1], ChrIdxs[c[2]], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
				else:
					out.write('%s %i %s %i %.8f %.8f %2f %.8f %.8f %.2f %.2f %i %.8f %.2f %.2f\n' % (ChrIdxs[c[2]], c[3], ChrIdxs[c[0]], c[1], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],c[14]) )
			cl = []
			elp = timeit.default_timer() - start_time
			gf.printlog('\t\t\t\t all %i interactions are processed, %.2fs.' % (total,elp), logname)
	out.close()
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t%i interactions are lifovered, %.2fs.' % (total,elp), logname)

def _listing(coor, contactHash, covHash, totalH, distanceH, res, params_list, **kwargs):
	
	c1,b1,c2,b2 = coor
	noised, pred_null_model, pab_cont, pab_cov, pab_res = params_list
	pc,poe = 0,0
	if (c1 != c2): real_dist = -1000
	else: real_dist = abs(b1-b2)
	try: mean_pc = distanceH[real_dist]['mean_contact']
	except KeyError: mean_pc = distanceH[max(distanceH.keys())]['mean_contact']
	if noised: noise = 1
	else: noise = 0
	mcov = covHash[c1,b1] * covHash[c2,b2]
	if ((c1,b1,c2,b2) in contactHash) and (contactHash[c1,b1,c2,b2][0] > noise): pc,poe = contactHash[c1,b1,c2,b2]
	elif ((c2,b2,c1,b1) in contactHash) and (contactHash[c2,b2,c1,b1][0] > noise): pc,poe = contactHash[c2,b2,c1,b1]
	elif pred_null_model and mcov >0:
		try: pc,poe,cov1,cov2 = _predictContacts((c1,b1,c2,b2), real_dist, covHash, totalH, distanceH, res, params_list[1:])
		except KeyError: pass
	else: pass
	_h = (c1,b1,c2,b2,pc,poe,mcov,mean_pc)
	return _h

def _easyRescale(high_cl, high_res, covHash, distanceH, res, params_rescale,**kwargs):
	contact_count, allCon, random_func, pick_probality_func = params_rescale
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
					covHash[c1,ri],covHash[c2,rj]
					if (c1 != c2): distance = -1000
					else: distance = abs(ri-rj)
					try: mean_con = distanceH[distance]['mean_contact']
					except KeyError: mean_con = distanceH[max(distanceH.keys())]['mean_contact']
					_clh.append((c1,ri,c2,rj,1.,1., covHash[c1,ri], covHash[c2,rj], mean_con,distance))
				except KeyError: pass
				except ZeroDivisionError: pass
		try:
			if random_func:
				_clh = np.array(_clh)
				k = _clh[:,6]*_clh[:,7]
				k[k>0] = 1.
				p = pick_probality_func(k*_clh[:,6],k*_clh[:,7],_clh[:,8],contact_distance=_clh[:,9])
				# p = np.apply_along_axis(pick_probality_func,1,[k*_clh[:,6],k*_clh[:,7],_clh[:,8],_clh[:,9]])
				p = p/np.sum(p)
				_clh = _randomizing(_clh, p, count, gf._choice)
			else: _clh = _randomizing(_clh, allCon, contact_count, False)
		except IndexError: pass#
		except ValueError: pass
		_cl.extend(_clh)
		_clh = []
	del _clh
	return _cl

def iContactRegression(covHash,resolution,c1_c2,ChrIdxs,low_ChrSizes,out_name,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: total_statistics = kwargs['total_statistics']
	except KeyError: total_statistics = False
	try: distance_dependence = kwargs['distance_dependence']
	except KeyError: distance_dependence = False
	try: LowObjCoorMP = kwargs['markpoints_low']
	except KeyError: LowObjCoorMP = False
	try: cov_low = kwargs['coverage_low']
	except KeyError:
		cov_low = False
		LowObjCoorMP = False
	try: distance_dependence_low = kwargs['distance_dependence_low']
	except KeyError:
		distance_dependence_low = False
		LowObjCoorMP = False
	try: total_statistics_low = kwargs['total_statistics_low']
	except KeyError: total_statistics_low = False
	try: contact_low = kwargs['contact_low']
	except KeyError:
		contact_low = False
		cov_low = False
		total_statistics_low = False
		distance_dependence_low = False
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
	try: random_func = kwargs['random_func']
	except KeyError: random_func = False
	try: noised = kwargs['noised']
	except KeyError: noised = False
	try: pred_null_model = kwargs['predict_null_contacts']
	except KeyError: pred_null_model = False
	try: pick_contacts = kwargs['pick_contacts']
	except KeyError: pick_contacts = False
	allCon,meanCov,meanMCov = total_statistics['all_contacts'],total_statistics['mean_coverage_sum'],total_statistics['mean_coverage_multiple']
	try: contact_count = kwargs['contact_count']
	except KeyError: contact_count = allCon
	if contact_count == 0: contact_count = allCon
	total = 0
	params_model = noised, pred_null_model, pab_cont, pab_cov, pab_res
	params_rescale = contact_count, allCon, random_func, pick_contacts
	
	for ci in range(len(c1_c2)):
		for cj in range((ci+1),len(c1_c2)):
			chrm1,chrm2 = c1_c2[ci],c1_c2[cj]
			if chrm1 == chrm2: total += (low_ChrSizes[chrm1]+1)*low_ChrSizes[chrm1]//2
			else: total += low_ChrSizes[chrm1]*low_ChrSizes[chrm2]
	num,H,cl = 0, [], []
	start_time = timeit.default_timer()
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t%i contacts are regressed to %i, %.2fs.' % (allCon,contact_count,elp), logname)
	
	chrm1,chrm2 = c1_c2
	full_name = '%s.%s.%s.allCon'% (out_name,chrm1,chrm2)
	gf.printlog('\t\t\t\twrite to %s' % (full_name), logname)
	out = open( full_name,'w')
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t\tcontact_count %s %s chromosome pair, %.2fs.' % (chrm1,chrm2,elp), logname)
	
	if chrm1 == chrm2: k = 1
	else: k = 0
	for i in range(low_ChrSizes[chrm1]):
		j1 = i*k
		for j in range(j1,low_ChrSizes[chrm2]):
			num += 1
			c1,b1,c2,b2 = ChrIdxs[chrm1],i,ChrIdxs[chrm2],j
			H.append(_listing((c1,b1,c2,b2), contact_low, cov_low, total_statistics_low, distance_dependence_low, low_res, params_model))
			if num % 200000 == 0:
				elp = timeit.default_timer() - start_time
				gf.printlog('\t\t\t\t start randomization, %.2fs.' % elp, logname)
				H = _randomizing(H, allCon, contact_count, random_func)
				H = _easyRescale(H, low_res, covHash, distance_dependence, resolution, params_rescale)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t\t\t\twriting, %.2fs.' % (elp), logname)
				for c in H:
					try: out.write('%s\t%i\t%s\t%i\t%.4f\n' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4]))
					except IndexError: pass
				H = []
			if num % 200000 == 0:
				elp = timeit.default_timer() - start_time
				gf.printlog('\t\t\t\t%i interactions from %i are processed, %.2fs.' % (num,total,elp), logname)
				
	if len(H) > 0:
		H = _randomizing(H, allCon, contact_count, random_func)
		H = _easyRescale(H, low_res, covHash, distance_dependence, resolution, params_rescale)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t\t\t\twriting, %.2fs.' % (elp), logname)
		for c in H:
			try: out.write('%s\t%i\t%s\t%i\t%.4f\n' % (ChrIdxs[c[0]], c[1]*resolution, ChrIdxs[c[2]], c[3]*resolution, c[4]))
			except IndexError: pass
	del H
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t\tend contact_count %s %s chromosome pair, %.2fs.' % (chrm1,chrm2,elp), logname)
	out.close()
	H = []
	elp = timeit.default_timer() - start_time
	gf.printlog('\t\t\t%i interactions are randomized, %.2fs.' % (total,elp), logname)

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
