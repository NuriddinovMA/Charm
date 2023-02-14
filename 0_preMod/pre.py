import os
import sys
import timeit
import pre_func as prf
reload(prf)
Args = {
	'chrom_sizes':'','sample_name':'','out_path':'', 'inter': False,
	'contacts':'', 'resolution':False, 'low_resolution':False,
	'pointview':False,'drop_bin':False,'log_file':False
}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		if key == 'pointview': args = parse[0].strip().split()[2:]
		else: args = parse[0].strip().split()[2].strip()
	except IndexError: pass
	else:
		print key, '=', args
		try: Args[key] = int(args)
		except ValueError: Args[key] = args
		except TypeError: Args[key] = args
		except KeyError: pass
Args['inter'] = prf.boolean(Args['inter'])
if Args['low_resolution']: Args['low_resolution'] = int(Args['low_resolution']/2)
Args['log_file'] = prf.boolean(Args['log_file'])
Args['pointview'] = prf.boolean(Args['pointview'])
if Args['pointview']: Args['pointview'] = Args['pointview'][0],int(Args['pointview'][1])/Args['resolution'],int(Args['pointview'][2])/Args['resolution']

for key in Args: print key, ' = ', Args[key]

try: os.makedirs(Args['out_path'])
except OSError: pass

try: os.remove(Args['log_file'])
except OSError: pass

start_time = timeit.default_timer()
l2i = prf.ChromIndexing(Args['chrom_sizes'])
if Args['resolution']:
	if Args['pointview']: pointview = Args['pointview'][0],int(Args['pointview'][1])/Args['resolution'],int(Args['pointview'][2])/Args['resolution']
	suffixH = '%s.%ikb' % (Args['sample_name'],Args['resolution']/1000)
	try: os.makedirs(Args['out_path']+'/'+suffixH)
	except OSError: pass
	
	prf.printlog('\nStep 0: data preparing...', Args['log_file'])
	prf.printlog('\tGenome analysis...', Args['log_file'])
	
	c2s = prf.ChromSizes(Args['chrom_sizes'],Args['resolution'])
	counts = prf.diag_counts(c2s,pointview=Args['pointview'],log=Args['log_file'])
	maxd = max(counts.keys())
	elp = timeit.default_timer() - start_time
	prf.printlog('\t...end genome analysing %.2fs' % elp, Args['log_file'])

	out_name = Args['out_path']+'/'+suffixH
	prf.printlog('\nStep 1: Calculating bin coverage...', Args['log_file'])
	binCov=prf.iBinCoverage(Args['contacts'],c2s,Args['resolution'],out=out_name,chrm_index=l2i,pointview=pointview,log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	prf.printlog('...bin coverage calculated for %.2fs' % elp, Args['log_file'])
	print binCov.keys()[:10]
	out_name = Args['out_path']+'/'+suffixH
	prf.printlog('\nStep 2: Distance depended statistics...', Args['log_file'])
	contactDistanceHash = prf.iDistanceRead(maxd,path=Args['contacts'],pointview=pointview,resolution=Args['resolution'],coverage=binCov,log=Args['log_file'])
	meanHash = prf.iMeaner( contactDistanceHash, counts, out_name,log=Args['log_file'])
	del contactDistanceHash
	elp = timeit.default_timer() - start_time
	prf.printlog('...distance analyzed for %.2fs' % elp, Args['log_file'])
	out_name = Args['out_path']+'/'+suffixH+'/'+suffixH
	prf.printlog('\nStep 3.1: Contact transforming by mean statistic...', Args['log_file'])
	prf.iTotalContactListing(meanHash,binCov,Args['resolution'],out_name,pointview=pointview,path=Args['contacts'],log=Args['log_file'])
	del meanHash
	del binCov
	elp = timeit.default_timer() - start_time
	prf.printlog('...contact transformed for %.2fs' % elp, Args['log_file'])

if Args['low_resolution']:
		if Args['pointview']: pointview = Args['pointview'][0],int(Args['pointview'][1])/Args['low_resolution'],int(Args['pointview'][2])/Args['low_resolution']
	suffixL = '/low.%s.%ikb' % (Args['sample_name'],Args['low_resolution']/500)
	try: os.makedirs(Args['out_path']+suffixL)
	except OSError: pass
	prf.printlog('\nStep 0.2: data preparing...', Args['log_file'])
	prf.printlog('\tLow resolution genome analysis...', Args['log_file'])
	c2s_low = prf.ChromSizes(Args['chrom_sizes'],Args['low_resolution'])
	counts_low = prf.diag_counts(c2s_low,log=Args['log_file'])
	maxd_low = max(counts_low.keys())
	elp = timeit.default_timer() - start_time
	prf.printlog('\t...end genome analysing %.2fs' % elp, Args['log_file'])

	out_name_low = Args['out_path']+'/'+suffixL
	prf.printlog('\nStep 1.2: Calculating bin coverage for lower resolution...', Args['log_file'])
	binCovLow=prf.iBinCoverage(Args['contacts'],c2s_low,Args['low_resolution'],out=out_name_low,chrm_index=l2i,pointview=pointview,log=Args['log_file'])
	elp = timeit.default_timer() - start_time
	prf.printlog('...bin coverage calculated for %.2fs' % elp, Args['log_file'])
	prf.printlog('\nStep 2.2: Distance depended statistics for lower resolution...', Args['log_file'])
	lowContacts = prf.iSemiLow(Args['contacts'],Args['low_resolution'],log=Args['log_file'])
	contactDistanceHashLow = prf.iDistanceRead(maxd_low,hash=lowContacts,coverage=binCovLow,log=Args['log_file'])
	meanHashLow = prf.iMeaner( contactDistanceHashLow, counts_low, out_name_low,log=Args['log_file'])
	del contactDistanceHashLow
	elp = timeit.default_timer() - start_time
	prf.printlog('...low resolution analyzed for %.2fs' % elp, Args['log_file'])

	out_name_low = Args['out_path']+'/'+suffixL+'/'+suffixL
	prf.printlog('\nStep 3.2: Low contact transforming by mean statistic...', Args['log_file'])
	prf.iTotalContactListing(meanHashLow,binCovLow,Args['low_resolution'],out_name_low,hash=lowContacts,log=Args['log_file'])
	del meanHashLow
	del binCovLow
	del lowContacts
	elp = timeit.default_timer() - start_time
	prf.printlog('...low contact transformed time %.2fs' % elp, Args['log_file'])

elp = timeit.default_timer() - start_time
prf.printlog('Full processing for %.2fs' % elp, Args['log_file'])
