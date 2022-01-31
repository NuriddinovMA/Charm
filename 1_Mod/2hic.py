import os
import sys
import timeit
import lift_func as lf
reload(lf)

print 'Step 1: Sites Reading...'
start_time = timeit.default_timer()

Args = {
	'contact_path':'','mutant_contacts':'','wt_hic_pre1':'','wt_hic_pre2':'',
	'chrom_path':'','chrom_sizes':'', 'resolution':'',
	'out_path':'','out_name':'','path_to_juicer':'', 'format':'short'
	}

lines = sys.stdin.readlines()
for line in lines:
	parse = line.split('#')
	try:
		key = parse[0].strip().split()[0]
		args = parse[0].strip().split()[2:]
	except IndexError: pass
	else:
		try: Args[key] = args[0]
		except KeyError: pass

for key in Args.keys(): print '\t%s =' % key, Args[key]
command = "java -jar %s pre" % Args['path_to_juicer']
print '\tcommand', command

RH = {
	5000: "1000000,500000,250000,100000,50000,25000,10000,5000",
	10000: "1000000,500000,250000,100000,50000,20000,10000",
	1: "1000000,500000,250000,100000,50000,20000,10000",
	20000: "1000000,500000,250000,100000,40000,20000",
	25000: "1000000,500000,250000,100000,50000,25000",
	40000: "1000000,500000,250000,80000,40000",
	50000: "1000000,500000,250000,100000,50000"
}
elp = timeit.default_timer() - start_time
print 'Step 2: Analyzing', elp

mut_dirs = '%s/%s' % (Args['contact_path'],Args['mutant_contacts'])
mut_hic_pre = '%s/%s' % (Args['out_path'],Args['out_name'])
G = '%s/%s' % (Args['chrom_path'],Args['chrom_sizes'])
Order = lf.ChromIndexing(G)
resolution = int(Args['resolution'])
elp = timeit.default_timer() - start_time
print '\tstart generate mutant pre %s, %.2f' % (mut_hic_pre, elp)
lf.JuiceboxPre(mut_dirs,resolution,Args['out_path'],Args['out_name'])
lf.SummingPre(mut_hic_pre,Args['wt_hic_pre1'],Args['wt_hic_pre2'],Args['out_path'],Args['out_name'],
	order=Order,format=Args['format'])
os.system('rm -r %s/%s' % (Args['out_path'],Args['out_name'] ))
os.system('rm -r %s/%s.summ' % (Args['out_path'],Args['out_name'] ))
print '\tpre writing', elp
R = RH[resolution]
F = '%s/%s.summ.pre.short' % (Args['out_path'],Args['out_name'] )
# O = '%s/%s.summ.hic' % (Args['out_path'],Args['out_name'] )
# os.system( command + " " + F + " " + O + " " + G + " " + "-r" + " " + R + " "+ "-n")
os.system('gzip ' + F)
os.system('chmod 755 ' + F + '.gz')
print '\t%s end hic generation %.2f' % (file, elp)
