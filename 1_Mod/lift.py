import os
import sys
import timeit
import lift_func as lf

Args = {
	'contact_path':'','contact_files':'','coverage_files':'','distance_files':'','normalization':False,
	'chrom_orders_from':'','chrom_orders_to':'','chosen_chroms_from':False,'chosen_chroms_to':False,
	'remap_files':'','out_path':'','out_name':'',
	'resolution':10000, 'model':'easy','regression':0,'multiples':1.0,'coverage_alignment':False,
	'random':False, 'pointview':False,'predict_null_contacts':False,
	'contact_ab':False,'coverage_ab':False,'resolution_ab':False,
	'contact_coef':False,'coverage_coef':False,'distance_coef':False, 'coef_resolution':False,
	'log_file':False
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		parse = parse[0].strip().split()
		key = parse[0]
		if key in ['chosen_chroms_from','chosen_chroms_to','pointview']: args = parse[2:]
		else: args = parse[2]
		Args[key] = args
	except IndexError: pass

Args['resolution'] = int(Args['resolution'])
Args['regression'] = int(Args['regression'])
Args['multiples'] = float(Args['multiples'])
Args['normalization'] = lf.boolean(Args['normalization'])
random = lf.boolean(Args['random'])
Args['chosen_chroms_to'] = lf.boolean(Args['chosen_chroms_to'])
Args['chosen_chroms_from'] = lf.boolean(Args['chosen_chroms_from'])
Args['pointview'] = lf.boolean(Args['pointview'])
Args['predict_null_contacts'] = lf.boolean(Args['predict_null_contacts'])
Args['contact_coef'] = lf.boolean(Args['contact_coef'])
Args['coverage_coef'] = lf.boolean(Args['coverage_coef'])
Args['distance_coef'] = lf.boolean(Args['distance_coef'])
Args['contact_ab'] = lf.boolean(Args['contact_ab'])
Args['coverage_ab'] = lf.boolean(Args['coverage_ab'])
Args['log_file'] = lf.boolean(Args['log_file'])

try: Args['coef_resolution'] = int(Args['coef_resolution'])
except ValueError: Args['coef_resolution'] = False
try: Args['resolution_ab'] = int(Args['resolution_ab'])/2
except ValueError: Args['resolution_ab'] = False
start_time = timeit.default_timer()

for key in Args: print key, ' = ', Args[key]

try: os.makedirs(Args['out_path'])
except OSError: pass

lf.printlog('log of lift.py script',Args['log_file'])
lf.printlog('Step 0: chromosome indexing...',Args['log_file'])
l2i_from = lf.ChromIndexing(Args['chrom_orders_from'])
print l2i_from
elp = timeit.default_timer() - start_time
lf.printlog('...chromosome indexied, %.2fs' % elp, Args['log_file'])
if Args['chrom_orders_to']:
	l2i_to = lf.ChromIndexing(Args['chrom_orders_to'])
	elp = timeit.default_timer() - start_time
	lf.printlog('...chromosome indexied to, %.2fs' % elp, Args['log_file'])
else: l2i_to = l2i_from
lf.printlog('Step 1: data reading...', Args['log_file'])

fname = Args['contact_path']+'/'+Args['coverage_files']
covHash = lf.readCovHash(fname,l2i_from,log=Args['log_file'])
elp = timeit.default_timer() - start_time
lf.printlog('\t...coverage reading end time %.2fs' % elp, Args['log_file'])
if Args['normalization'] == False: Norm = {}
else: Norm = covHash
if Args['coverage_alignment'] == False: Norm = {}
else:
	fname = Args['contact_path']+'/'+Args['coverage_alignment']
	covHash_2 = lf.readCovHash(fname,l2i_from,log=Args['log_file'])
	Norm = lf.covAlignment(covHash,covHash_2)

fname = Args['contact_path']+'/'+Args['distance_files']
psList = lf.readMeanHash(fname,log=Args['log_file'])
elp = timeit.default_timer() - start_time
lf.printlog('\t...read distance end time %.2fs' % elp, Args['log_file'])

fname = Args['contact_path']+'/'+Args['contact_files']
contactHash = lf.iReadInitialContact(fname,l2i_from,chrms=Args['chosen_chroms_from'],norm=Norm, log=Args['log_file'])
ln = len(contactHash)
elp = timeit.default_timer() - start_time
lf.printlog('\t...contact %i reading end time %.2fs' % (ln,elp), Args['log_file'])


#if Args['pointview']: Args['pointview'] = Args['pointview'][0],l2i_from[Args['pointview'][1]],int(Args['pointview'][2])/Args['resolution'],int(Args['pointview'][3])/Args['resolution']

if Args['contact_coef']:
	if Args['contact_coef'] != Args['contact_files']: scale = 1
	else: scale, Args['coef_resolution'] = 10, 10*Args['resolution']
	fname = Args['contact_path']+'/'+Args['contact_coef']
	contactCoef = lf.iReadInitialContact(fname,l2i_from,chrms=Args['chosen_chroms_from'],scale=scale,log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... multiple coef reading end time %.2fs' % elp, Args['log_file'])
else: contactCoef = False

if Args['coverage_coef']:
	fname = Args['contact_path']+'/'+Args['coverage_coef']
	covCoef = lf.readCovHash(fname,l2i_from,log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... coef coverage reading end time %.2fs' % elp, Args['log_file'])
else: covCoef = False

if Args['distance_coef']:
	fname = Args['contact_path']+'/'+Args['distance_coef']
	psListCoef = lf.readMeanHash(fname,log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t...read coef distance end time %.2fs' % elp, Args['log_file'])
else: psListCoef = False

if Args['contact_ab']:
	fname = Args['contact_path']+'/'+Args['contact_ab']
	contactAB = lf.iReadInitialContact(fname,l2i_from,chrms=Args['chosen_chroms_from'],log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... multiple compartments reading end time %.2fs' % elp, Args['log_file'])
else: contactAB = False

if Args['coverage_ab']:
	fname = Args['contact_path']+'/'+Args['coverage_ab']
	covAB = lf.readCovHash(fname,l2i_from,log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... compartments coverage reading end time %.2fs' % elp, Args['log_file'])
else: covAB = False

elp = timeit.default_timer() - start_time
lf.printlog('...end data reading %.2fs' % elp, Args['log_file'])

lf.printlog('\nStep 2: Reading mark points...', Args['log_file'])
MarkPoints = lf.iReadingMarkPoints(Args['remap_files'],Args['resolution'],chrm_index_from=l2i_from,chrms_from=Args['chosen_chroms_from'],
	chrm_index_to=l2i_to,chrms_to=Args['chosen_chroms_to'], log=Args['log_file'])

if Args['contact_coef']:
	MarkPointsCoef = lf.iReadingMarkPoints(Args['remap_files'],Args['coef_resolution'],chrm_index_from=l2i_from,chrms_from=Args['chosen_chroms_from'],
		chrm_index_to=l2i_to,chrms_to=Args['chosen_chroms_to'], log=Args['log_file'])
else: MarkPointsCoef = False

if MarkPoints: pass
else: 
	print 'ERROR!!! No Markpoints!'
	lf.printlog('ERROR!!! No Markpoints!', Args['log_file'])
	exit()
ln = len(MarkPoints)
elp = timeit.default_timer() - start_time
lf.printlog('...%i mark point readed for %.2f sec' % (ln, elp), Args['log_file'])

lf.printlog('\nStep 3: start contact modeling...', Args['log_file'])
suffix = '%s.%ikb.%s.%s' % (Args['out_name'],Args['resolution']/1000,Args['model'],Args['random'])
try: os.mkdir( '%s/%s'% (Args['out_path'],suffix))
except OSError: pass
out_name = '%s/%s/%s' % (Args['out_path'],suffix,suffix)

if Args['chosen_chroms_from'] == False: Args['chosen_chroms_from'] = l2i_from.keys()
if Args['chosen_chroms_to'] == False: Args['chosen_chroms_to'] = l2i_to.keys()

# if Args['pointview']:
	# lf.iPointviewContact(contactHash, covHash, MarkPoints, l2i_to, Args['pointview'], out_name+'.temp',
		# model=Args['model'], scoring=psList, random=random, regression=Args['regression'], predict=Args['predict'],
		# contact_coef=contactCoef,rescale=rescale,log=Args['log_file'])
# else:
lf.iLiftOverContact(contactHash, covHash, MarkPoints, Args['resolution'], l2i_to, out_name+'.temp',
	model=Args['model'], scoring=psList, random=random, regression=Args['regression'],null_contacts=Args['predict_null_contacts'],
	contact_coef=contactCoef, coverage_coef=covCoef, scoring_coef=psListCoef, markpoints_coef=MarkPointsCoef, coef_resolution=Args['coef_resolution'],
	contact_ab=contactAB, coverage_ab=covAB, resolution_ab=Args['resolution_ab'],
	log=Args['log_file'])

os.system('rm %s/%s/*.liftCon' % (Args['out_path'],suffix))
chrms = list(set(Args['chosen_chroms_to']))
header = "chr1 bin1 chr2 bin2 contact oe mult_cov prev_contact prev_oe prev_mult_cov reality expected normolize_coef balance_coef"
for ci in range(len(chrms)):
	i = chrms[ci]
	for cj in range(len(chrms)):
		j = chrms[cj]
		if l2i_to[i] <= l2i_to[j]:
			try: 
				f = open("%s.%s.%s.liftCon" % (out_name,i,j),'r')
				f.close()
			except IOError:
				with open("%s.%s.%s.liftCon" % (out_name,i,j),'w') as f: print >> f, header
			os.system('grep -E "^%s .* %s " %s.temp >> %s.%s.%s.liftCon' % (i,j,out_name,out_name,i,j))
		else:
			try: 
				f = open("%s.%s.%s.liftCon" % (out_name,j,i),'r')
				f.close()
			except IOError:
				with open("%s.%s.%s.liftCon" % (out_name,j,i),'w') as f: print >> f, header
			os.system('grep -E "^%s .* %s " %s.temp >> %s.%s.%s.liftCon' % (i,j,out_name,out_name,j,i))
os.system('rm %s.temp' % out_name)
elp = timeit.default_timer() - start_time
lf.printlog('end contact writing %.2fs'% elp, Args['log_file'])
exit()
