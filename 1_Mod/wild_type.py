import os
import sys
import timeit
import lift_func as lf

Args = {
	'contact_path':'','contact_files':'','cov_files':'','dist_files':'',
	'genome_path':'','chrom_orders':'','chosen_chroms':[],'random':False,
	'out_path':'','out_name':'','resolution':10000, 'regression':0,
	'pointview':False,'contact_coef':{}, 'coef_resolution':False,
	'log_file':False
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		parse = parse[0].strip().split()
		key = parse[0]
		if key == 'chosen_chroms': args = parse[2],parse[3]
		elif key == 'pointview': args = parse[2:]
		else: args = parse[2]
		Args[key] = args
	except IndexError: pass
Args['resolution'] = int(Args['resolution'])
Args['regression'] = int(Args['regression'])
Args['contact_coef'] = lf.boolean(Args['contact_coef'])
random = lf.boolean(Args['random'])
Args['pointview'] = lf.boolean(Args['pointview'])
if Args['pointview']:
	Args['pointview'] = Args['pointview'][0],int(Args['pointview'][1])/Args['resolution'],int(Args['pointview'][2])/Args['resolution']
	Args['chosen_chroms'] = Args['pointview'][0],Args['pointview'][0]
	

Args['log_file'] = lf.boolean(Args['log_file'])
try: 
	Args['coef_resolution'] = int(Args['coef_resolution'])/2
	rescale=(Args['resolution'],Args['coef_resolution'])
except ValueError: 
	Args['coef_resolution'] = False
	rescale=(1,1)
start_time = timeit.default_timer()

for key in Args: print key, ' = ', Args[key]

try: os.makedirs(Args['out_path'])
except OSError: pass
suffix = '%s.%s' % (Args['contact_files'],Args['regression'])
try: os.makedirs(Args['out_path']+'/'+suffix)
except OSError: pass
try: os.makedirs(Args['out_path']+'/'+suffix +'/' + Args['out_name'])
except OSError: pass

lf.printlog('Step 0: chromosome indexing...',Args['log_file'])

fname = Args['genome_path']+'/'+Args['chrom_orders']
l2i = lf.ChromIndexing(fname)
c2s = lf.ChromSizes(fname,Args['resolution'])
elp = timeit.default_timer() - start_time
lf.printlog('...chromosome indexied, %.2fs' % elp, Args['log_file'])

lf.printlog('Step 1: data reading...', Args['log_file'])

fname = Args['contact_path']+'/'+Args['cov_files']
covHash = lf.readCovHash(fname,l2i,log=Args['log_file'])
elp = timeit.default_timer() - start_time
lf.printlog('\t...coverage reading end time %.2fs' % elp, Args['log_file'])

fname = Args['contact_path']+'/'+Args['dist_files']
psList = lf.readMeanHash(fname,log=Args['log_file'])
elp = timeit.default_timer() - start_time
lf.printlog('\t...read distance end time %.2fs' % elp, Args['log_file'])

if Args['contact_coef']:
	fname = Args['contact_path']+'/'+Args['contact_coef']
	contactCoef = lf.iReadInitialContact(fname,l2i,log=Args['log_file'],chrms=Args['chrms'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... multiple  coef reading end time %.2fs' % elp, Args['log_file'])
else: contactCoef = False
 
elp = timeit.default_timer() - start_time
lf.printlog('...end data reading %.2fs' % elp, Args['log_file'])

lf.printlog('\nStep 2: start contact randomizing...', Args['log_file'])
files = os.listdir(Args['contact_path']+'/'+Args['contact_files'])
files.sort()
lnf = len(files)
for i in range(lnf):
	fname = '%s/%s/%s' % (Args['contact_path'],Args['contact_files'],files[i])
	parse = files[i].split('.')
	chrm1, chrm2 = parse[-3], parse[-2]
	if (chrm1, chrm2) == Args['chosen_chroms']:
		lf.printlog('\tstart randomize %s, %i/%s' %(files[i],(i+1),lnf), Args['log_file'] )
		out_name = '%s/%s/%s/%s.%s.%s.%s.pre' % (Args['out_path'],suffix,Args['out_name'],suffix,Args['out_name'],chrm1,chrm2)
		lf.iContactRegression( fname, covHash, l2i, c2s, Args['resolution'], out_name,
			scoring=psList, regression=Args['regression'],random=Args['random'],
			pointview=Args['pointview'],contact_coef=contactCoef, rescale=rescale, log=Args['log_file']
			)
		elp = timeit.default_timer() - start_time
		lf.printlog('\t...end randomize %s, %.2f' %(files[i],elp), Args['log_file'] )
elp = timeit.default_timer() - start_time
lf.printlog('\tend contact writing %.2fs'% elp, Args['log_file'])
exit()
