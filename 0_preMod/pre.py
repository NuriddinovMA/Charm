import os
import sys
import timeit
import pre_func as prf
reload(prf)
Args = {
	'chrom_sizes':'','sample_name':'','out_path':'', 'inter': False,
	'contacts':'', 'resolution':10000, 'low_resolution':1000000,
	'log_file':False
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2].strip()
	except IndexError: pass
	else:
		print key, '=', args
		try: Args[key] = int(args)
		except ValueError: Args[key] = args
		except KeyError: pass
Args['inter'] = prf.boolean(Args['inter'])
Args['low_resolution'] = int(Args['low_resolution']/2)
Args['log_file'] = prf.boolean(Args['log_file'])

for key in Args: print key, ' = ', Args[key]

try: os.makedirs(Args['out_path'])
except OSError: pass

suffixH = '/%s.%ikb' % (Args['sample_name'],Args['resolution']/1000)
try: os.makedirs(Args['out_path']+'/'+suffixH)
except OSError: pass
suffixL = '/low.%s.%ikb' % (Args['sample_name'],Args['low_resolution']/500)
try: os.makedirs(Args['out_path']+suffixL)
except OSError: pass

start_time = timeit.default_timer()
with open(Args['log_file'],'w') as log: print >> log, 'log of pre.py script'

prf.printlog('\nStep 0: data preparing...', Args['log_file'])
prf.printlog('\tGenome analysis...', Args['log_file'])

l2i = prf.ChromIndexing(Args['chrom_sizes'])
c2s = prf.ChromSizes(Args['chrom_sizes'],Args['resolution'])
counts = prf.diag_counts(c2s,log=Args['log_file'])
maxd = max(counts.keys())
elp = timeit.default_timer() - start_time
prf.printlog('\t...end genome analysing %.2fs' % elp, Args['log_file'])

out_name = Args['out_path']+'/'+suffixH

prf.printlog('\nStep 1.1: Distance depended statistics...', Args['log_file'])
contactDistanceHash = prf.iDistanceRead(maxd,path=Args['contacts'],resolution=Args['resolution'],log=Args['log_file'])
meanHash = prf.iMeaner( contactDistanceHash, counts, out_name,log=Args['log_file'])
del contactDistanceHash
elp = timeit.default_timer() - start_time
prf.printlog('...distance analyzed for %.2fs' % elp, Args['log_file'])

prf.printlog('\nStep 2: Calculating bin coverage...', Args['log_file'])
prf.iBinCoverage(Args['contacts'],c2s,Args['resolution'],out_name,log=Args['log_file'])
elp = timeit.default_timer() - start_time
prf.printlog('...bin coverage  calculated for %.2fs' % elp, Args['log_file'])
del l2i
out_name = Args['out_path']+'/'+suffixH+'/'+suffixH

prf.printlog('\nStep 3.1: Contact transforming by mean statistic...', Args['log_file'])
prf.iTotalContactListing(meanHash,Args['resolution'],out_name,path=Args['contacts'],log=Args['log_file'])
del meanHash
elp = timeit.default_timer() - start_time
prf.printlog('...contact transformed for %.2fs' % elp, Args['log_file'])


prf.printlog('\nStep 0: data preparing...', Args['log_file'])
prf.printlog('\tLow resolution genome analysis...', Args['log_file'])
c2s_low = prf.ChromSizes(Args['chrom_sizes'],Args['low_resolution'])
counts_low = prf.diag_counts(c2s_low,log=Args['log_file'])
maxd_low = max(counts_low.keys())
elp = timeit.default_timer() - start_time
prf.printlog('\t...end genome analysing %.2fs' % elp, Args['log_file'])

out_name_low =  Args['out_path']+'/'+suffixL

prf.printlog('\nStep 1.2: Distance depended statistics for lower resolution...', Args['log_file'])
lowContacts = prf.iSemiLow(Args['contacts'],Args['low_resolution'],log=Args['log_file'])
contactDistanceHashLow = prf.iDistanceRead(maxd_low,hash=lowContacts,log=Args['log_file'])
meanHashLow = prf.iMeaner( contactDistanceHashLow, counts_low, out_name_low,log=Args['log_file'])
del contactDistanceHashLow
elp = timeit.default_timer() - start_time
prf.printlog('...low resolution analyzed for %.2fs' % elp, Args['log_file'])

out_name_low =  Args['out_path']+'/'+suffixL+'/'+suffixL
prf.printlog('\nStep 3.2: Low contact transforming by mean statistic...', Args['log_file'])
prf.iTotalContactListing(meanHashLow,Args['low_resolution'],out_name_low,hash=lowContacts,log=Args['log_file'])
del meanHashLow
del lowContacts
elp = timeit.default_timer() - start_time
prf.printlog('...low contact transformed time %.2fs' % elp, Args['log_file'])

elp = timeit.default_timer() - start_time
prf.printlog('Full processing for %.2fs' % elp, Args['log_file'])
