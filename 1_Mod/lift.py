import os
import sys
import timeit
import lift_func as lf

Args = {
	'contact_path':'','contact_files':'','cov_files':'','dist_files':'',
	'genome_path':'','chrom_orders':'','chosen_chrom':[],
	'remap_path':'','remap_files':'','out_path':'','out_name':'',
	'resolution':10000, 'model':'easy','regression':0,
	'random':False, 'predict':False, 'contact_coef':False, 'coef_resolution':False,
	'log_file':False
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		parse = parse[0].strip().split()
		key = parse[0]
		if key != 'chosen_chrom': args = parse[2]
		else: args = parse[2:]
		Args[key] = args
	except IndexError: pass
Args['resolution'] = int(Args['resolution'])
Args['regression'] = int(Args['regression'])
random = lf.boolean(Args['random'])
Args['predict'] = lf.boolean(Args['predict'])
Args['contact_coef'] = lf.boolean(Args['contact_coef'])
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

lf.printlog('log of lift.py script',Args['log_file'])
lf.printlog('Step 0: chromosome indexing...',Args['log_file'])

fname = Args['genome_path']+'/'+Args['chrom_orders']
l2i = lf.ChromIndexing(fname)
elp = timeit.default_timer() - start_time
lf.printlog('...chromosome indexied, %.2fs' % elp, Args['log_file'])

lf.printlog('Step 1: data reading...', Args['log_file'])

fname = Args['contact_path']+'/'+Args['contact_files']
contactHash = lf.iReadInitialContact(fname,l2i,chrms=Args['chosen_chrom'],log=Args['log_file'])
ln = len(contactHash)
elp = timeit.default_timer() - start_time
lf.printlog('\t...contact %i reading end time %.2fs' % (ln,elp), Args['log_file'])

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
	contactCoef = lf.iReadInitialContact(fname,l2i,chrms=Args['chosen_chrom'],log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	lf.printlog('\t... multiple  coef reading end time %.2fs' % elp, Args['log_file'])
else: contactCoef = False
 
elp = timeit.default_timer() - start_time
lf.printlog('...end data reading %.2fs' % elp, Args['log_file'])

lf.printlog('\nStep 2: Reading mark points...', Args['log_file'])
rname = Args['remap_path']+'/'+Args['remap_files']
MarkPoints = lf.iReadingMarkPoints(rname,Args['resolution'],l2i,chrms=Args['chosen_chrom'], log=Args['log_file'])
ln = len(MarkPoints)
elp = timeit.default_timer() - start_time
lf.printlog('...%i mark point readed for %.2f sec' % (ln, elp), Args['log_file'])

lf.printlog('\nStep 3: start contact modeling...', Args['log_file'])
suffix = '%s.%ikb.%s.%s' % (Args['out_name'],Args['resolution']/1000,Args['model'],Args['random'])
try: os.mkdir( '%s/%s'% (Args['out_path'],suffix))
except OSError: pass
out_name = '%s/%s/%s' % (Args['out_path'],suffix,suffix)
lf.iDifferContact(contactHash, covHash, MarkPoints, l2i, out_name+'.temp',
	model=Args['model'], scoring=psList, random=random, regression=Args['regression'], predict=Args['predict'],
	contact_coef=contactCoef,rescale=rescale,log=Args['log_file'])

for i in Args['chosen_chrom']:
	for j in Args['chosen_chrom']:
		if l2i[i] <= l2i[j]: os.system('grep -E "%s.*%s" %s.temp > %s.%s.%s.liftCon' % (i,j,out_name,out_name,i,j))
os.system('rm %s.temp' % out_name)
elp = timeit.default_timer() - start_time
lf.printlog('\tend contact writing %.2fs'% elp, Args['log_file'])
exit()
