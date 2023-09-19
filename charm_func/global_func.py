import os

def boolean(x):
	if str(x).lower() in ['false','no','n'] : return False
	elif str(x).lower() in ['true','yes','y']: return True
	else: return x

def printlog(str, logname):
	print(str)
	if logname:
		with open(logname, 'a') as log: log.write(str+'\n')

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
		try: ChrSzs[parse[0]] = (int(parse[1])-1)//int(resolution)+1
		except IndexError: break
	del lines
	return ChrSzs

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
	
def colorList():
	color = []
	for i in range(0,256,85):
		for j in range(0,256,85):
			for k in range(0,256,85): color.append('%i,%i,%i' % (i,j,k))
	return color

def alog2(x):
	try: return abs(math.log(x,2))
	except ValueError: return -10