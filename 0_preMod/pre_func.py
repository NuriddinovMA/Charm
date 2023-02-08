import timeit
import sys
import os
try: import numpy as np
except ImportError: print "numpy not found!"
try: import scipy.stats as sc
except ImportError: print "scipy.stats not found!"

def printlog(str, logname):
	print str
	if logname:
		with open(logname, 'a') as log: print >> log, str


def HashTry(Hash, key):
	t = 1
	try: Hash[key]
	except KeyError: t = 0
	return t

def FindKey(Hash, key):
	if key > 0:
		for k in range(key,-1,-1):
			try: 
				Hash[k]
				return k
			except KeyError: pass
 	else: return key

def boolean(x):
	if x == 'False' or x == 'false' or x == 'f' or x == 'F': return False
	elif x == 'True' or x == 'true' or x == 't' or x == 'T': return True
	else: return x

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
		try: ChrSzs[parse[0]] = int(parse[1])/resolution+1
		except IndexError: break
	del lines
	return ChrSzs

def diag_counts(chrSzs,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: pointview = kwargs['pointview']
	except KeyError: pointview = False
	counts = {}
	start_time = timeit.default_timer()
	printlog('\t\tstart diag count',logname)
	if pointview:
		lnc = pointview[2]-pointview[1]
		for l in range(lnc):
			try: counts[l] += (lnc-l)
			except KeyError: counts[l] = (lnc-l)
	else:
		Chrms = chrSzs.keys()
		Chrms.sort()
		lnc = len(Chrms)
		full = (lnc+1)*lnc/2
		for i in range(lnc):
			for l in range(chrSzs[Chrms[i]]):
				try: counts[l] += (chrSzs[Chrms[i]]-l)
				except KeyError: counts[l] = (chrSzs[Chrms[i]]-l)
			for j in range(i+1,lnc):
				l = -1000
				try: counts[l] += chrSzs[Chrms[i]]*chrSzs[Chrms[j]]
				except KeyError: counts[l] = chrSzs[Chrms[i]]*chrSzs[Chrms[j]]
	elp = timeit.default_timer() - start_time
	printlog('\t\tend diag count %.2fs' % elp,logname)
	return counts

def iBinCoverage(path, ChrSzs, resolution, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: l2i = kwargs['chrm_index']
	except KeyError: l2i = False
	try: pointview = kwargs['pointview']
	except KeyError: pointview = False
	try: out = kwargs['out']
	except KeyError: out = False
	try: k1,k2 = kwargs['name_index']
	except KeyError: k1,k2 = 1,2
	try: format = kwargs['format']
	except KeyError: format = 'vp'
	fullBins,CI,Chrms = [],{},ChrSzs.keys()
	Chrms.sort()
	files = os.listdir(path)
	S,SO,N = 0,0,0
	nf,lnf,lnc = 0,len(files),len(Chrms)
	printlog('\t\tstart % i file reading' % lnf,logname)
	start_time = timeit.default_timer()
	if pointview: 
		chrmp,cp1,cp2 = pointview[0],pointview[1],pointview[2]
		CI[chrmp] = [0 for i in range(ChrSzs[chrmp])]
	else: 
		for key in ChrSzs: CI[key] = [[0 for j in range(len(ChrSzs))] for i in range(ChrSzs[key])]
	if format == 'pre': idx = 2,6,8
	else: idx = 0,1,2
	for file in files:
		nf += 1
		printlog('\t\t\tstart reading ' + file,logname)
		name = file.split('.')
		stf = timeit.default_timer()
		if pointview and name[k1] == name[k2] == pointview[0]:
			with open(path + '/' + file,'r') as f: lines = f.readlines()
			for line in lines:
				parse = line.split()
				c1,c2,p = int(parse[idx[0]])/resolution,int(parse[idx[1]])/resolution,float(parse[idx[2]])
				if pointview[1]<=c1<pointview[2] and pointview[1]<=c2<pointview[2]:
					if np.isnan(p): p = 0
					S += p
					CI[name[k1]][c1] += p
					CI[name[k2]][c2] += p
		if pointview == False:
			with open(path + '/' + file,'r') as f: lines = f.readlines()
			for line in lines:
				parse = line.split()
				c1,c2,p = int(parse[idx[0]])/resolution,int(parse[idx[1]])/resolution,float(parse[idx[2]])
				if np.isnan(p): p = 0
				S += p
				CI[name[k1]][c1][l2i[name[k2]]-1] += p
				CI[name[k2]][c2][l2i[name[k1]]-1] += p
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		printlog('\t\t\t...processed %i/%i files %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)

	elp = timeit.default_timer() - start_time
	printlog('\t\tcontacts count %.2f, %.2fs' % (S,elp),logname)
	M,MN = 0,0
	if pointview:
		SO = np.sum(CI[chrmp][cp1:cp2])*np.sum(CI[chrmp][cp1:cp2])
		M = np.nanmean(CI[chrmp][cp1:cp2])
		N = (cp2-cp1)**2
	else:
		for i in range(lnc):
			li = len(CI[Chrms[i]])
			for j in range(i,lnc):
				lj = len(CI[Chrms[j]])
				N += li*lj
				SO += np.sum(CI[Chrms[i]])*np.sum(CI[Chrms[j]])
			M += np.sum(CI[Chrms[i]])
			MN += li
	SO = 1.*SO/N
	M = 1.*M/MN
	elp = timeit.default_timer() - start_time
	printlog('\t\tmean coverage %.2f, calculating %.2fs' % (SO,elp),logname)
	if out:
		fout = open('%s.binCov' % out,'w')
		printlog('\t\twrite to ' + out +'.binCov',logname)
		print >> fout, 'all_contacts: %.2f' % S
		print >> fout, 'mean_coverage: %.2f' % M
		print >> fout, 'mean_coverage_multiple: %.2f' % SO
		if pointview:
			print >> fout, 'chr1 bin1 coverage'
		else:
			chrstr = ''
			for i in range(len(ChrSzs)): chrstr = chrstr + ' ' + l2i[i+1]
			print >> fout, 'chr1 bin1 coverage ' + chrstr
			del chrstr
		if pointview:
			for ki in range(cp1,cp2): print >> fout, '%s %i %f' % (chrmp,ki,CI[chrmp][ki])
		else:
			for i in range(lnc):
				li = len(CI[Chrms[i]])
				for ki in range(li):
					cn = CI[Chrms[i]][ki]
					covstr = str(np.sum(cn))
					for cov in cn: covstr = covstr + ' ' + str(cov)
					print >> fout, '%s %i %s' % (Chrms[i],ki,covstr)
		fout.close()
	elp = timeit.default_timer() - start_time
	printlog('\t\tgenerate count of bin pairs %.2fs' % elp,logname)
	return CI

def iNormContactListing( path, CI, resolution, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: maxcovp = kwargs['max_coverage_percentile']
	except KeyError: maxcovp = 0
	try: mincovp = kwargs['min_coverage_percentile']
	except KeyError: mincovp = 0
	try: msum = kwargs['matrix_sum']
	except KeyError: msum = False
	try: out = kwargs['out']
	except KeyError: out = './..'
	try: k1,k2 = kwargs['name_index']
	except KeyError: k1,k2 = 1,2
	try: format = kwargs['format']
	except KeyError: format = 'vp'
	S = 0
	Vals = []
	printlog('\t\tstart percentile calculating',logname)
	start_time = timeit.default_timer()
	for key in CI:
		for k in CI[key]:
			s = np.sum(k)
			Vals.append(s)
			S += s
	S = S/2
	if mincovp: mincov = np.percentile(Vals,mincovp,interpolation='higher')
	else: mincov = np.min(Vals)
	if maxcovp: maxcov = np.percentile(Vals,maxcovp,interpolation='lower')
	else: maxcov = np.max(Vals)
	elp = timeit.default_timer() - start_time
	printlog('\t\tpercentile calculated %.2fs' % elp,logname)
	Drop = {}
	printlog('\t\tstart drop list creating',logname)
	for key in CI:
		for k in range(len(CI[key])):
			s = np.sum(CI[key][k])
			if mincovp or maxcovp:
				if s < mincov or s > maxcov: 
					S -= s
					Drop[key,k] = True
				else: Drop[key,k] = False
			CI[key][k] = 0
	elp = timeit.default_timer() - start_time
	files = os.listdir(path)
	nf,lnf = 0,len(files)
	printlog('\t\tstart % i file reading %.2fs' % (lnf,elp),logname)
	if format == 'pre': idx = 2,6,8
	else: idx = 0,1,2
	for file in files:
		nf += 1
		printlog('\t\t\tstart reading ' + file,logname)
		name = file.split('.')
		stf = timeit.default_timer()
		with open(path + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			c1,c2,p = int(parse[idx[0]])/resolution,int(parse[idx[1]])/resolution,float(parse[idx[2]])
			try: D1=Drop[name[k1],c1]
			except KeyError: D1 = False
			try: D2=Drop[name[k2],c2]
			except KeyError: D2 = False
			if D1 or D2: pass
			else:
				S += p
				CI[name[k1]][c1] += p
				CI[name[k2]][c2] += p
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		printlog('\t\t\t...processed %i/%i files %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	elp = timeit.default_timer() - start_time
	printlog('\t\tcontacts count %.2f, %.2fs' % (S,elp),logname)
	nf = 0
	if msum: mult=1.*msum/S
	else: mult=1.
	printlog('\t\tcontacts count %.2f, msum=%.2f, %.2fs' % (S,msum,elp),logname)
	for file in files:
		nf += 1
		printlog('\t\t\tstart reading ' + file,logname)
		name = file.split('.')
		stf = timeit.default_timer()
		if format == 'pre': fout = open(out + '/%s.%ikb.norm.%s.%s.pre' % (name[0],resolution/1000,name[k1],name[k2]),'w')
		else: fout = open(out + '/%s.%s.%s.%i.vc' % (name[0],name[k1],name[k2],resolution),'w')
		with open(path + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			c1,c2,p = int(parse[idx[0]])/resolution,int(parse[idx[1]])/resolution,float(parse[idx[2]])
			if np.isnan(p): p = 0
			try: D1=Drop[name[k1],c1]
			except KeyError: D1 = False
			try: D2=Drop[name[k2],c2]
			except KeyError: D2 = False
			if D1 or D2: pass
			else:
				p = mult*p/(CI[name[k1]][c1]*CI[name[k2]][c2])
				if format == 'pre': print >> fout, '%s\t%s\t%i\t%s\t%s\t%s\t%i\t%s\t%f' % (parse[0],parse[1],c1*resolution,parse[3],parse[4],parse[5],c2*resolution,parse[7],p)
				else: print >> fout, c1*resolution,c2*resolution,p
		fout.close()
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		printlog('\t\t\t...processed %i/%i files %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	elp = timeit.default_timer() - start_time
	printlog('\t\tcontacts normed, %.2fs' % elp,logname)
	printlog('\t\tstart normalization writig to ' + out,logname)

def iSemiLow(path,resolution,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	printlog('\tstart low contact analyzing',logname)
	lowCon,lowNorm,CI = {},{},{}
	files = os.listdir(path)
	lnf,nf = len(files),0
	start_time = timeit.default_timer()
	for file in files:
		nf += 1
		printlog('\t\tstart reading '+file,logname)
		name = file.split('.')
		stf = timeit.default_timer()
		with open(path + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			c1,c2,p = int(parse[0])/resolution,int(parse[1])/resolution,float(parse[2])
			if np.isnan(p): p = 0
			try: lowCon[name[1],c1,name[2],c2] += p
			except KeyError: lowCon[name[1],c1,name[2],c2] = p
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		printlog('\t\t\t...end %i/%i files reading %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	Keys = lowCon.keys()
	Keys.sort(key=lambda k: (k[0],k[2],k[1],k[3]))
	ln = len(Keys)
	step = int(ln/10)
	elp = timeit.default_timer() - start_time
	printlog('\tstart semitransformation for %i contats %.2fs' % (ln,elp),logname)
	for i in range(ln):
		k1,c1,k2,c2 = Keys[i]
		try: lowCon[Keys[i]] += lowCon[(k1,(c1+1),k2,c2)]
		except KeyError: pass
		try: lowCon[Keys[i]] += lowCon[(k1,c1,k2,(c2+1))]
		except KeyError: pass
		try: lowCon[Keys[i]] += lowCon[(k1,(c1+1),k2,(c2+1))]
		except KeyError: pass
		if i % step == 0:
			elp = timeit.default_timer() - start_time
			printlog('\t... %i/%i contact transformed end time: %.2fs' % (i,ln,elp),logname)
	elp = timeit.default_timer() - start_time
	printlog('\tbin coverage %.2fs' % elp,logname)
	S = 0
	for key in Keys:
		k1,c1,k2,c2 = key
		try: CI[k1,c1] += lowCon[key]
		except KeyError: CI[k1,c1] = lowCon[key]
		try: CI[k2,c2] += lowCon[key]
		except KeyError: CI[k2,c2] = lowCon[key]
		S += lowCon[key]
	elp = timeit.default_timer() - start_time
	printlog('\tstart coverage normalization %.2fs' % elp,logname)
	for key in lowCon.keys():
		k1,c1,k2,c2 = key
		try: normc = lowCon[key]/(CI[k1,c1]*CI[k2,c2])*S
		except ZeroDivisionError: normc = 0
		try: lowNorm[k1,k2][c1,c2] = lowCon[key],normc
		except KeyError: lowNorm[k1,k2] = { (c1,c2):(lowCon[key],normc) }
		del lowCon[key]
	del CI
	printlog('\t\t...end low contact analyzing %.2fs' % elp,logname)
	return lowNorm

def iDistanceRead(lmax,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: 
		path=kwargs['path']
		resolution=kwargs['resolution']
		hash=False
	except KeyError: path=False
	try: 
		hash=kwargs['hash']
		path=False
	except KeyError: hash=False
	try: bincov = kwargs['coverage']
	except KeyError: bincov=False
	try: pointview = kwargs['pointview']
	except KeyError: pointview = False
	
	iDH = {i:[] for i in range(lmax+1)}
	iDH[-1000] = []
	start_time = timeit.default_timer()
	
	if hash:
		lh = len(hash)
		printlog('\tstart distance read from hash %i chromosome pairs' % lh,logname)
		n = 0
		for chrms in hash:
			k1,k2 = chrms
			for cons in hash[chrms]:
				p=hash[chrms][cons][1]
				c1,c2=cons
				mcov = np.sum(bincov[k1][c1]) * np.sum(bincov[k2][c2])
				if k1 == k2: iDH[abs(c2-c1)].append((p,mcov))
				else: iDH[-1000].append((p,mcov))
				n += 1
				if n % 1000000 == 0:
					elp = timeit.default_timer() - start_time
					printlog('\t... read %i contact time: %.2fs;' % (n, elp),logname)
		elp = timeit.default_timer() - start_time
		printlog('\t\tend distance read from hash, %.2fs' % elp,logname)

	if path:
		files = os.listdir(path)
		lnf,nf = len(files),0
		printlog('\tstart distance read from %i files' % lnf,logname)
		for file in files:
			H = {}
			nf += 1
			stf = timeit.default_timer()
			printlog('\t\tstart reading '+file,logname)
			c1,c2 = file.split('.')[1:3]
			if pointview:
				with open(path+'/'+file,'r') as f: lines = f.readlines()
				if name[1] == name[2] == pointview[0]:
					for line in lines:
						parse = line.split()
						b1,b2,p = int(parse[0])/resolution,int(parse[1])/resolution,float(parse[2])
						if pointview[1]<=c1<pointview[2] and pointview[1]<=c2<pointview[2]:
							if np.isnan(p): p = 0
							try: H[c1,b1,c2,b2] += p
							except KeyError: H[c1,b1,c2,b2] = p
						else: pass
			else:
				with open(path+'/'+file,'r') as f: lines = f.readlines()
				for line in lines:
					parse = line.split()
					b1,b2,p = int(parse[0])/resolution,int(parse[1])/resolution,float(parse[2])
					if np.isnan(p): p = 0
					try: H[c1,b1,c2,b2] += p
					except KeyError: H[c1,b1,c2,b2] = p
				elpf = timeit.default_timer() - stf
				elp = timeit.default_timer() - start_time
			for key in H:
				c1,b1,c2,b2 = key
				p = H[key]
				if c1 == c2: l = abs(b2-b1)
				else: l = -1000
				mcov = np.sum(bincov[c1][b1]) * np.sum(bincov[c2][b2])
				iDH[l].append((p,mcov))
			del H
			elpf = timeit.default_timer() - stf
			elp = timeit.default_timer() - start_time
			printlog('\t\t\t...end %i/%i files reading %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	elp = timeit.default_timer() - start_time
	printlog('\t\t...end distance analyzing %.2fs' % elp,logname)
	return iDH

def iMeanStatistics(contactDistanceHash, count):
	start_time = timeit.default_timer()
	mean = 1.*np.sum(contactDistanceHash,axis=0)[0]/count
	prop = np.mean(contactDistanceHash,axis=0)
	try: prcList = ( np.min(contactDistanceHash,axis=0)[0],mean,np.max(contactDistanceHash,axis=0)[0],prop[0],prop[1])
	except ValueError: prcList = (0,0,0,0,0)
	return prcList

def iMeaner(contactDistanceHash,counts,path,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	start_time = timeit.default_timer()
	Keys = counts.keys()
	Keys.sort()
	ld = len(Keys)
	if ld < 1000: step = int((ld-2)/10)
	else: step = int((ld-2)/100)
	# mxd = max(100,step)
	meanHash = {i:(0,0,0) for i in range(Keys[-1]+1)}
	f = open(path + '.stat','w')
	base_count = counts[0]
	printlog('\tstart statistic for %i distances' % ld,logname)
	print >> f,'distances\tcontact_counts\tdistance_combined\tmin\tmean\tmax'
	for l in Keys:
		if (l % step == 0) or (l < 0): printlog('\t\tstart %i (%i) distance analyze:' % (l,len(contactDistanceHash[l])),logname)
		count,d = counts[l],0,
		try: cdh = [i for i in contactDistanceHash[l]]
		except KeyError: cdh = []
		d = 0
		while (5*count < base_count) or (len(cdh) == 0):
			d += 1
			try:
				cdh.extend(contactDistanceHash[l+d])
				count += counts[l+d]
			except KeyError: pass
			try:
				k = max(0,(l-d))
				cdh.extend(contactDistanceHash[k])
				count += counts[k]
			except KeyError: pass
		# try: del contactDistanceHash[l-d]
		# except KeyError: pass
		cdh.sort()
		meanHash[l] = iMeanStatistics(cdh, count)
		if l >= 0: print >> f, '%i\t%i\t%i' % (l, len(cdh),d),
		else: print >> f, 'interchromosome\t%i\t0' % len(cdh),
		for i in meanHash[l]: print >> f, '%f' % i,
		print >> f, ''
		if (l % step == 0) or (l < 0):
			elp = timeit.default_timer() - start_time
			printlog('\t\t\t...%i distance generate statistic: %.2fs, combine = %i' % (l,elp,d),logname)
	f.close()
	del contactDistanceHash 
	elp = timeit.default_timer() - start_time
	printlog('\t\t...generate distance statistic: %.2fs' % elp,logname)
	return meanHash

def iTotalContactListing( meanHash, binCov, resolution, out, **kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	try: 
		path=kwargs['path']
		hash=False
	except KeyError: path=False
	try: 
		hash=kwargs['hash']
		path=False
	except KeyError: hash=False
	try: pointview = kwargs['pointview']
	except KeyError: pointview = False
	if path:
		files = os.listdir(path)
		n,lnf = 0,len(files)
		printlog('\tStart contact transformation from %i files' % lnf,logname)
		start_time = timeit.default_timer()
		for file in files:
			printlog('\t\topen '+file,logname)
			name = file.split('.')
			H = {}
			if pointview and name[1] == name[2] == pointview[0]:
				with open(path + '/' + file,'r') as f: lines = f.readlines()
				for line in lines:
					parse = line.split()
					c1,c2,p = int(parse[0])/resolution,int(parse[1])/resolution,float(parse[2])
					if pointview[1]<=c1<pointview[2] and pointview[1]<=c2<pointview[2]: 
						try: H[name[1],c1,name[2],c2] += p
						except KeyError: H[name[1],c1,name[2],c2] = p
						n += 1
						if n % 1000000 == 0:
							elp = timeit.default_timer() - start_time
							printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			if pointview == False:
				with open(path + '/' + file,'r') as f: lines = f.readlines()
				for line in lines:
					parse = line.split()
					c1,c2,p = int(parse[0])/resolution,int(parse[1])/resolution,float(parse[2])
					try: H[name[1],c1,name[2],c2] += p
					except KeyError: H[name[1],c1,name[2],c2] = p
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			fout = open('%s.%s.%s.allCon' % (out,name[1],name[2]),'w')
			print >> fout, 'chr1\tbin1\tchr2\tbin2\tcontacts\toe\tcov_mult'
			Keys = H.keys()
			Keys.sort()
			for key in Keys:
				c1,b1,c2,b2 = key
				p = H[key]
				if c1 == c2: l = abs(b2-b1)
				else: l = -1000
				mean = meanHash[l][1]
				mult_cov = np.sum(binCov[c1][b1])*np.sum(binCov[c2][b2])
				try: 
					pm = 1.*p/mean
					print >> fout, '%s %i %s %i %f %f %i' % (c1,b1,c2,b2,p,pm,mult_cov)
				except ZeroDivisionError: printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),lognam)
			fout.close()
			del H
	if hash:
		n,lnc = 0,len(hash)
		printlog('\tstart %i contact transformation from hash' % lnc,logname)
		start_time = timeit.default_timer()
		chrms = hash.keys()
		chrms.sort()
		lnc = len(chrms)
		for i in range(lnc):
			printlog('\t\t processed %i/%i chromosome pairs' % ((i+1),lnc),logname)
			c1,c2 = chrms[i]
			if c1 == c2:
				fout = open('%s.%s.%s.allCon' % (out,c1,c2),'w')
				print >> fout, 'chr1\tbin1\tchr2\tbin2\tcontacts\toe\tcov_mult\tnorm\tnormed_contacts'
				pairs = hash[chrms[i]].keys()
				pairs.sort()
				for pair in pairs:
					b1,b2 = pair
					p1,p2 = hash[chrms[i]][pair]
					l = abs(b2-b1)
					mean = meanHash[l][1]
					mult_cov = np.sum(binCov[c1][b1])*np.sum(binCov[c2][b2])
					try: 
						pm = 1.*p2/mean
						print >> fout, '%s %i %s %i %f %f %i %f' % (c1,b1,c2,b2,p1,pm,mult_cov,p2)
					except ZeroDivisionError: printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			else:
				fout = open('%s.%s.%s.allCon' % (out,c1,c2),'w')
				print >> fout, 'chr1\tbin1\tchr2\tbin2\tcontacts\toe\tcov_mult\tnormed_contacts'
				l = -1000
				mean = meanHash[l][1]
				pairs = hash[chrms[i]].keys()
				pairs.sort()
				for pair in pairs:
					b1,b2 = pair
					p1,p2 = hash[chrms[i]][pair]
					mult_cov = np.sum(binCov[c1][b1])*np.sum(binCov[c2][b2])
					try: 
						pm = 1.*p2/mean
						print >> fout, '%s %i %s %i %f %f %i %f' % (c1,b1,c2,b2,p1,pm,mult_cov,p2)
					except ZeroDivisionError: printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			fout.close()
