import os
import sys
import timeit
import lift_func as lf
reload(lf)

Args = {
	'contact_path':'','contact_files':'','coverage_files':'','distance_files':'','normalization':False,
	'chrom_order':'','chosen_chroms':False,'totality':'all', 'out_path':'','out_name':'',
	'resolution':10000, 'model':'easy','regression':0,'multiples':1.0,'coverage_alignment':False,
	'random':False, 'predict':False, 'pointview':False,'null_contacts':False,
	'contact_ab':False,'coverage_ab':False,'ab_resolution':False,
	'contact_coef':False,'coverage_coef':False,'distance_coef':False, 'coef_resolution':False,
	'log_file':False
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		parse = parse[0].strip().split()
		key = parse[0]
		if key == 'chosen_chroms': args = parse[2:]
		else: args = parse[2]
		Args[key] = args
	except IndexError: pass

Args['resolution'] = int(Args['resolution'])
Args['regression'] = int(Args['regression'])
Args['multiples'] = float(Args['multiples'])
Args['normalization'] = lf.boolean(Args['normalization'])
random = lf.boolean(Args['random'])
Args['chosen_chroms'] = lf.boolean(Args['chosen_chroms'])
Args['predict'] = lf.boolean(Args['predict'])
Args['null_contacts'] = lf.boolean(Args['null_contacts'])
Args['contact_coef'] = lf.boolean(Args['contact_coef'])
Args['coverage_coef'] = lf.boolean(Args['coverage_coef'])
Args['distance_coef'] = lf.boolean(Args['distance_coef'])
Args['contact_ab'] = lf.boolean(Args['contact_ab'])
Args['coverage_ab'] = lf.boolean(Args['coverage_ab'])
Args['log_file'] = lf.boolean(Args['log_file'])

try: Args['coef_resolution'] = int(Args['coef_resolution'])
except ValueError: Args['coef_resolution'] = False
try: Args['ab_resolution'] = int(Args['ab_resolution'])/2
except ValueError: Args['ab_resolution'] = False
start_time = timeit.default_timer()

for key in Args: print key, ' = ', Args[key]

try: os.makedirs(Args['out_path'])
except OSError: pass

suffix = '%s.%s' % (Args['contact_files'],Args['regression'])
try: os.makedirs(Args['out_path']+'/'+suffix)
except OSError: pass
try: os.makedirs(Args['out_path']+'/'+suffix +'/' + Args['out_name'])
except OSError: pass

if len(Args['chosen_chroms']) == 1: Args['totality'] = 'intra'
elif len(Args['chosen_chroms']) == 2 and Args['chosen_chroms'][0] == Args['chosen_chroms'][1]: Args['totality'] = 'intra'
elif len(Args['chosen_chroms']) == 2 and Args['chosen_chroms'][0] != Args['chosen_chroms'][1]: Args['totality'] = 'inter'
else: Args['totality'] = 'all'

lf.printlog('Step 0: chromosome indexing...',Args['log_file'])

fname = Args['chrom_orders']
l2i = lf.ChromIndexing(fname)
c2s = lf.ChromSizes(fname,Args['resolution'])
c2s_coef = lf.ChromSizes(fname,Args['coef_resolution'])
elp = timeit.default_timer() - start_time
lf.printlog('...chromosome indexied, %.2fs' % elp, Args['log_file'])


lf.printlog('Step 1: data reading...', Args['log_file'])

fname = Args['contact_path']+'/'+Args['coverage_files']
covHash = lf.readCovHash(fname,l2i,log=Args['log_file'])
elp = timeit.default_timer() - start_time
lf.printlog('\t...coverage reading end time %.2fs' % elp, Args['log_file'])
if Args['normalization'] == False: Norm = {}
else: Norm = covHash
if Args['coverage_alignment'] == False: Norm = {}
else:
	fname = Args['contact_path']+'/'+Args['coverage_alignment']
	covHash_2 = lf.readCovHash(fname,l2i,log=Args['log_file'])
	Norm = lf.covAlignment(covHash,covHash_2)

fname = Args['contact_path']+'/'+Args['distance_files']
psList = lf.readMeanHash(fname,log=Args['log_file'])
elp = timeit.default_timer() - start_time
lf.printlog('\t...read distance end time %.2fs' % elp, Args['log_file'])

fname = Args['contact_path']+'/'+Args['contact_files']
contactHash = lf.iReadInitialContact(fname,l2i,chrms=Args['chosen_chroms'],norm=Norm, log=Args['log_file'],totality=Args['totality'])
ln = len(contactHash)
elp = timeit.default_timer() - start_time
lf.printlog('\t...contact %i reading end time %.2fs' % (ln,elp), Args['log_file'])


#if Args['pointview']: Args['pointview'] = Args['pointview'][0],l2i[Args['pointview'][1]],int(Args['pointview'][2])/Args['resolution'],int(Args['pointview'][3])/Args['resolution']

if Args['contact_coef']:
	if Args['contact_coef'] != Args['contact_files']: scale = 1
	else: scale, Args['coef_resolution'] = 10, 10*Args['resolution']
	fname = Args['contact_path']+'/'+Args['contact_coef']
	contactCoef = lf.iReadInitialContact(fname,l2i,chrms=Args['chosen_chroms'],scale=scale,log=Args['log_file'],totality=Args['totality'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... multiple coef reading end time %.2fs' % elp, Args['log_file'])

else: contactCoef = False

if Args['coverage_coef']:
	fname = Args['contact_path']+'/'+Args['coverage_coef']
	covCoef = lf.readCovHash(fname,l2i,log=Args['log_file'])
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
	contactAB = lf.iReadInitialContact(fname,l2i,chrms=Args['chosen_chroms'],log=Args['log_file'],totality=Args['totality'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... multiple compartments reading end time %.2fs' % elp, Args['log_file'])
else: contactAB = False

if Args['coverage_ab']:
	fname = Args['contact_path']+'/'+Args['coverage_ab']
	covAB = lf.readCovHash(fname,l2i,log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... compartments coverage reading end time %.2fs' % elp, Args['log_file'])
else: covAB = False

elp = timeit.default_timer() - start_time
lf.printlog('...end data reading %.2fs' % elp, Args['log_file'])

lf.printlog('\nStep 2: start contact randomizing...', Args['log_file'])
out_name = '%s/%s/%s/%s.%s' % (Args['out_path'],suffix,Args['out_name'],suffix,Args['out_name'])
lf.iContactRegression( contactHash, covHash, Args['resolution'], Args['chosen_chroms'], l2i, c2s_coef, out_name,
	model=Args['model'], scoring=psList, random=random, regression=Args['regression'], predict=Args['predict'],
	contact_coef=contactCoef, coverage_coef=covCoef, scoring_coef=psListCoef, coef_resolution=Args['coef_resolution'],
	contact_ab=contactAB, coverage_ab=covAB, ab_resolution=Args['ab_resolution'],
	log=Args['log_file'], null_contacts=Args['null_contacts']
	)
elp = timeit.default_timer() - start_time
lf.printlog('\t...end randomize %.2f' % elp, Args['log_file'] )
exit()
