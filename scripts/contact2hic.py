import os
import sys
import timeit
import c2h_func as c2h
import global_func as gf

def hic_generate(svs_contacts,wt1_contacts,wt2_contacts,
	chosen_chroms, chrom_sizes, resolution,
	format, path_to_java_dir, path_to_juicertools, hic_resolutions,sim_name,work_dir,log_file
	):
	
	start_time = timeit.default_timer()
	
	out_dir = '%s/out' % work_dir
	mut_hic_pre = '%s/out/%s' % (work_dir,sim_name)
	try: os.makedirs(out_dir)
	except OSError: pass
	resolution = int(resolution)
	path_to_java_dir = gf.boolean(path_to_java_dir)
	hic_resolutions = gf.boolean(hic_resolutions)
	l2i = gf.ChromIndexing(chrom_sizes)
	#if capture: capture = capture[0],int(capture[1]),int(capture[2])
	chosen_chroms = chosen_chroms.split(',')
	svs_contacts = gf.boolean(svs_contacts)
	wt1_contacts = gf.boolean(wt1_contacts)
	wt2_contacts = gf.boolean(wt2_contacts)
	
	elp = timeit.default_timer() - start_time
	gf.printlog( '\t Step 1: summing muntant and wild type contacts, %.2f' % (elp), log_file)

	c2h.SummingPre( svs_contacts, wt1_contacts, wt2_contacts, resolution, sim_name, out_dir, 
		order=l2i, format=format
		)

	elp = timeit.default_timer() - start_time
	gf.printlog( '\tpre writing %.2f' % elp, log_file)
	
	print(format)
	if format == 'hic':
		F = '%s/%s.pre' % ( out_dir, sim_name )
		O = '%s/%s.hic' % ( out_dir, sim_name )
		if path_to_java_dir: command = "%s/java -jar %s pre " % (path_to_java_dir,path_to_juicertools)
		else: command = "java -jar %s pre " % (path_to_juicertools)
		params = "%s %s %s" % (F,O,chrom_sizes)
		if hic_resolutions: resolution_list = ' -r %s' % hic_resolutions
		else: resolution_list= ''
		print(command + params + resolution_list)
		os.system( command + params + resolution_list)
	elif format == 'pre': 
		F = '%s/%s.pre' % ( out_dir, sim_name )
	elif format == 'pre.gz': 
		F = '%s/%s.pre' % ( out_dir, sim_name )
		try: os.remove(F + '.gz')
		except FileNotFoundError: pass
		os.system('gzip ' + F + '.gz')
	else:
		print(format)
		print('Error! Unsupported format, use "hic" or "pre", or "pre.gz" ')
		exit()

	elp = timeit.default_timer() - start_time
	gf.printlog('\tpath to resulted files: %s' % out_dir, log_file)
	gf.printlog('\tend hic generation %.2f' % elp, log_file)
