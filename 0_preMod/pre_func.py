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
	counts = {}
	Chrms = chrSzs.keys()
	Chrms.sort()
	lnc = len(Chrms)
	full = (lnc+1)*lnc/2
	start_time = timeit.default_timer()
	printlog('\t\tstart diag count',logname)
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

def iBinCoverage(path, ChrSzs, resolution, out,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	fullBins,CI,Chrms = [],{},ChrSzs.keys()
	Chrms.sort()
	files = os.listdir(path)
	S,SO,N = 0,0,0
	nf,lnf,lnc = 0,len(files),len(Chrms)
	printlog('\t\tstart % i file reading' % lnf,logname)
	start_time = timeit.default_timer()
	for key in ChrSzs: CI[key] = [0 for i in range(ChrSzs[key])]
	for file in files:
		nf += 1
		printlog('\t\t\tstart reading ' + file,logname)
		name = file.split('.')
		stf = timeit.default_timer()
		with open(path + '/' + file,'r') as f: lines = f.readlines()
		for line in lines:
			parse = line.split()
			c1,c2,p = int(parse[0]),int(parse[1]),float(parse[2])
			CI[name[1]][c1/resolution] += p
			CI[name[2]][c2/resolution] += p
			S += p
		elpf = timeit.default_timer() - stf
		elp = timeit.default_timer() - start_time
		printlog('\t\t\t...processed %i/%i files %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)
	elp = timeit.default_timer() - start_time
	printlog('\t\tcontacts count %i, %.2fs' % (S,elp),logname)
	for i in range(lnc):
		li = len(CI[Chrms[i]])
		for j in range(i,lnc):
			lj = len(CI[Chrms[j]])
			N += li*lj
			SO += np.sum(CI[Chrms[i]])*np.sum(CI[Chrms[j]])
	SO = 1.*SO/N
	elp = timeit.default_timer() - start_time
	printlog('\t\tmean coverage %.2f, calculating %.2fs' % (SO,elp),logname)
	fout = open(out +'.binCov','w')
	printlog('\t\twrite to ' + out +'.binCov',logname)
	print >> fout, 'all_contacts: %i' % S
	print >> fout, 'mean_coverage: %.2f' % SO
	print >> fout, 'chr1 pos1 coverage'
	for i in range(lnc):
		li = len(CI[Chrms[i]])
		for ki in range(li):
			cn = CI[Chrms[i]][ki]
			print >> fout, '%s %i %i' % (Chrms[i],ki,cn)
	fout.close()
	elp = timeit.default_timer() - start_time
	printlog('\t\tgenerate count of bin pairs %.2fs' % elp,logname)

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
		try: lowNorm[k1,k2][c1,c2] = lowCon[key],lowCon[key]/(CI[k1,c1]*CI[k2,c2])*S
		except KeyError: lowNorm[k1,k2] = { (c1,c2):(lowCon[key],lowCon[key]/(CI[k1,c1]*CI[k2,c2])) }
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
				if k1 == k2: iDH[abs(c2-c1)].append(p)
				else: iDH[-1000].append(p)
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
			nf += 1
			stf = timeit.default_timer()
			printlog('\t\tstart reading '+file,logname)
			name = file.split('.')
			with open(path+'/'+file,'r') as f: lines = f.readlines()
			if name[1] == name[2]:
				for line in lines:
					parse = line.split()
					c1,c2,p = int(parse[0]),int(parse[1]),float(parse[2])
					iDH[abs(c2-c1)/resolution].append(p)
			else:
				for line in lines: 
					parse = line.split()
					p = int(float(parse[2]))
					iDH[-1000].append(p)
			elpf = timeit.default_timer() - stf
			elp = timeit.default_timer() - start_time
			printlog('\t\t\t...end %i/%i files reading %.2fs (%.2fs)' % (nf,lnf,elpf,elp),logname)

	elp = timeit.default_timer() - start_time
	printlog('\t\t...end distance analyzing %.2fs' % elp,logname)
	return iDH

def iMeanStatistics(contactDistanceHash, count):
	start_time = timeit.default_timer()
	mean = 1.*np.sum(contactDistanceHash)/count
	try: prcList = ( np.min(contactDistanceHash),mean,np.max(contactDistanceHash))
	except ValueError: prcList = (0,0,0)
	return prcList

def iMeaner(contactDistanceHash,counts,path,**kwargs):
	try: logname = kwargs['log']
	except KeyError: logname = False
	start_time = timeit.default_timer()
	Keys = counts.keys()
	Keys.sort()
	ld = len(Keys)
	step = int(ld-2/10)
	meanHash = {i:(0,0,0) for i in range(Keys[-1]+1)}
	f = open(path + '.stat','w')
	base_count = counts[0]
	printlog('\tstart statistic for %i distances' % ld,logname)
	print >> f,'distances\tcontact_counts\tdistance_combined\tmin\tmean\tmax'
	for l in Keys:
		if (l % step == 0) or (l < 0): printlog('\t\tstart %i distance analyze:' % l,logname)
		count,d = counts[l],0,
		try: cdh = [i for i in contactDistanceHash[l]]
		except KeyError: cdh = []
		d = 0
		while (5*count < base_count) and (d < 100):
			d += 1
			try:
				cdh.extend(contactDistanceHash[l+d])
				count += counts[l+d]
			except KeyError: pass
			try:
				cdh.extend(contactDistanceHash[l-d])
				count += counts[l-d]
			except KeyError: pass
		try: del contactDistanceHash[l-100]
		except KeyError: pass
		cdh.sort()
		meanHash[l] = iMeanStatistics(cdh, count)
		if l >= 0: print >> f, '%i\t%i\t%i' % (l, len(cdh),d),
		else: print >> f, 'interchromosome\t%i\t0' % len(cdh),
		for i in meanHash[l]: print >> f, '%.2f' % i,
		print >> f, ''
		if (l % step == 0) or (l < 0):
			elp = timeit.default_timer() - start_time
			printlog('\t\t\t...%i distance generate statistic: %.2fs, combine = %i' % (l,elp,d),logname)
	f.close()
	elp = timeit.default_timer() - start_time
	printlog('\t\t...generate distance statistic: %.2fs' % elp,logname)
	return meanHash

def iTotalContactListing( meanHash, resolution, out, **kwargs):
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
	
	if path:
		files = os.listdir(path)
		n,lnf = 0,len(files)
		printlog('\tStart contact transformation from %i files' % lnf,logname)
		start_time = timeit.default_timer()
		for file in files:
			printlog('\t\topen '+file,logname)
			f = open(path + '/' + file,'r')
			lines = f.readlines()
			f.close()
			name = file.split('.')
			if name[1] == name[2]:
				fout = open('%s.%s.%s.allCon' % (out,name[1],name[2]),'w')
				print >> fout, 'chr1\tpos1\tchr2\tpos2\tcontacts\toe\tdistance'
				for line in lines:
					parse = line.split()
					c1,c2,p = int(parse[0])/resolution,int(parse[1])/resolution,float(parse[2])
					l = abs(c2-c1)
					mean = meanHash[l][1]
					try: 
						pm = 1.*p/mean
						print >> fout, '%s %i %s %i %i %.4f %i' % (name[1],c1,name[2],c2,p,pm,l)
					except ZeroDivisionError: printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)

			else:
				fout = open('%s.%s.%s.allCon' % (out,name[1],name[2]),'w')
				print >> fout, 'chr1\tpos1\tchr2\tpos2\tcontacts\toe\tdistance'
				l = -1000
				mean = meanHash[l][1]
				for line in lines:
					parse = line.split()
					c1,c2,p = int(parse[0])/resolution,int(parse[1])/resolution,float(parse[2])
					try: 
						pm = 1.*p/mean
						print >> fout, '%s %i %s %i %i %.4f %i' % (name[1],c1,name[2],c2,p,pm,l)
					except ZeroDivisionError: printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			fout.close()
	if hash:
		n,lnc = 0,len(hash)
		printlog('\tstart %i contact transformation from hash' % lnc,logname)
		start_time = timeit.default_timer()
		chrms = hash.keys()
		chrms.sort()
		lnc = len(chrms)
		for i in range(lnc):
			printlog('\t\t processed %i/%i chromosome pairs' % ((i+1),lnc),logname)
			k1,k2 = chrms[i]
			if k1 == k2:
				fout = open('%s.%s.%s.allCon' % (out,k1,k2),'w')
				print >> fout, 'chr1\tpos1\tchr2\tpos2\tcontacts\tnorm\toe\tdistance'
				pairs = hash[chrms[i]].keys()
				pairs.sort()
				for pair in pairs:
					c1,c2 = pair
					p1,p2 = hash[chrms[i]][pair]
					l = abs(c2-c1)
					mean = meanHash[l][1]
					try: 
						pm = 1.*p2/mean
						print >> fout, '%s %i %s %i %i %.4f %.4f %i' % (k1,c1,k2,c2,p1,p2,pm,l)
					except ZeroDivisionError: printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			else:
				fout = open('%s.%s.%s.allCon' % (out,k1,k2),'w')
				print >> fout, 'chr1\tpos1\tchr2\tpos2\tcontacts\tnorm\toe\tdistance'
				l = -1000
				mean = meanHash[l][1]
				pairs = hash[chrms[i]].keys()
				pairs.sort()
				for pair in pairs:
					c1,c2 = pair
					p1,p2 = hash[chrms[i]][pair]
					try: 
						pm = 1.*p2/mean
						print >> fout, '%s %i %s %i %.i %.4f %.4f %i' % (k1,c1,k2,c2,p1,p2,pm,l)
					except ZeroDivisionError: printlog('\t\tZeroDivisionError %s %s %i %.2f' % (c1,c2,l,meanHash[l]),logname)
					n += 1
					if n % 1000000 == 0:
						elp = timeit.default_timer() - start_time
						printlog('\t... %i contact transformed end time: %.2fs;' % (n, elp),logname)
			fout.close()
