import os
import sys
import timeit
import lift_func as lf
import importlib as imp
imp.reload(lf)

print( 'Step 1: Sites Reading...' )
start_time = timeit.default_timer()

Args = {
	'mutant_contacts':'','chosen_chroms':False,'wt1_contacts':False,'wt2_contacts':False,
	'chrom_sizes':'', 'resolution':'','out_resolution':1,'pointview': False,'method':'summ',
	'out_path':'','out_name':'','path_to_juicer':'', 'format':'hic'
	}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2:]
		if (key in ['pointview','chosen_chroms']) == False: Args[key] = args[0]
		else: Args[key] = args
	except KeyError: pass
	except IndexError: pass

Args['pointview'] = lf.boolean(Args['pointview'])
Args['chosen_chroms'] = lf.boolean(Args['chosen_chroms'])
for key in Args: print( '\t%s = %s' % (key, Args[key]) )
command = "java -jar %s pre" % Args['path_to_juicer']
command_norm = "java -jar %s addNorm" % Args['path_to_juicer']
print( '\tcommand', command )

RH = {
	2000: "2048000,1024000,512000,256000,128000,64000,32000,16000,8000,4000,2000",
	4000: "2048000,1024000,512000,256000,128000,64000,32000,16000,8000,4000",
	5000: "1000000,500000,250000,100000,50000,25000,10000,5000",
	10000: "1000000,500000,250000,100000,50000,20000,10000",
	1: "1000000,500000,250000,100000,50000,25000,10000",
	20000: "1000000,500000,250000,100000,40000,20000",
	25000: "1000000,500000,250000,100000,50000,25000",
	40000: "1000000,500000,250000,80000,40000",
	50000: "1000000,500000,250000,100000,50000"
}
elp = timeit.default_timer() - start_time
print( 'Step 2: Analyzing', elp )

mut_hic_pre = '%s/%s' % (Args['out_path'],Args['out_name'])
os.system('mkdir ' + Args['out_path'])
os.system('chmod 775 ' + Args['out_path'])
G = Args['chrom_sizes']
Order = lf.ChromIndexing(G)
resolution = int(Args['resolution'])
out_resolution = int(Args['out_resolution'])
if Args['pointview']: Args['pointview'] = Args['pointview'][0],int(Args['pointview'][1]),int(Args['pointview'][2])
elp = timeit.default_timer() - start_time

print( '\tstart generate pre %s, %.2f' % (mut_hic_pre, elp) )

lf.SummingPre(
	Args['mutant_contacts'],Args['wt1_contacts'],Args['wt2_contacts'],Args['out_path'],Args['out_name'],
	out_res=out_resolution, order=Order,format=Args['format']
	)

os.system('rm -r %s/%s' % (Args['out_path'],Args['out_name'] ))
print( '\tpre writing', elp )
if out_resolution > resolution: R = RH[out_resolution ]
else: R = RH[out_resolution ]
if Args['format'] == 'hic': 
	F = '%s/%s.pre' % (Args['out_path'],Args['out_name'] )
	O = '%s/%s.hic' % (Args['out_path'],Args['out_name'] )
	os.system( command + " " + F + " " + O + " " + G + " " + "-r" + " " + R)
else: F = '%s/%s.pre.short' % (Args['out_path'],Args['out_name'] )
os.system('rm ' + F + '.gz')
os.system('gzip ' + F)
os.system('chmod 755 ' + F + '.gz')
os.system('chmod 755 ' + F[:-4] + '.hic')
print( '\tend hic generation %.2f' % elp )
