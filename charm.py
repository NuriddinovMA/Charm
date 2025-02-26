if __name__ == "__main__":
	
	import os
	import shutil
	import timeit
	import argparse
	from configparser import ConfigParser, ExtendedInterpolation
	from charm_func import global_func as gf
	
	try: import numpy as np
	except ModuleNotFoundError:
		print('Lethal Error! NumPy not found!')
		exit()

	parser = argparse.ArgumentParser(description='Parameters for generation of rearrangement maps')
	parser.add_argument( '-i', dest='ini',metavar='ini-file', help='path to ini-file')
	parser.add_argument('-S', dest='stage', metavar='Stage',default='no',help='must be one of "no", "pre+", "SVs+", "sim+", "lift+",  "wt", "hic"')
	parser.add_argument('-g','--global', dest='glo', metavar='global', nargs='?',default=False ,help='must be parameter:value pairs from global section "parameter2:value2"')
	parser.add_argument('-p','--pre', dest='pre', metavar='preprocessing', nargs='?',default=False,help='must be parameter:value pairs from preprocessing section "parameter2:value2"')
	parser.add_argument('-v','--svs', dest='svs', metavar='SVs', nargs='?',default=False, help='must be parameter:value pairs from SVs section "parameter2:value2"')
	parser.add_argument('-s','--sim', dest='sim', metavar='simulation', nargs='?',default=False, help='must be parameter:value pairs from simulation section "parameter2:value2"')
	parser.add_argument('-l','--lift', dest='lift', metavar='liftover', nargs='?',default=False, help='must be parameter:value pairs from liftover section "parameter2:value2"')
	parser.add_argument('-w','--wt', dest='wt', metavar='wild_type', nargs='?',default=False, help='must be parameter:value pairs from wild_type section "parameter2:value2"')
	parser.add_argument('-o','--hic', dest='hic', metavar='hic', nargs='?',default=False, help='must be parameter:value pairs from hic section "parameter2:value2"')
	args = parser.parse_args()
	
	config = ConfigParser(interpolation=ExtendedInterpolation())
	config.read(args.ini)
	commandLineConfig = {}
	
	print(sorted(config.keys()))
	if args.glo:
		print( args.glo )
		commandLineConfig['global'] = {}
		key_value = args.glo.split()
		for kv in key_value: 
			key,value = kv.split('=')
			commandLineConfig['global'][key] = value
			config['global'][key] = value

	if args.stage == 'no': pass
	elif args.stage == 'pre': config['global']['skip_stages'] = "svs,sim,lift,wt,hic"
	elif args.stage == 'pre+': config['global']['skip_stages'] = ""
	elif args.stage == 'SVs': config['global']['skip_stages'] = "pre,sim,lift,wt,hic"
	elif args.stage == 'SVs+': config['global']['skip_stages'] = "pre"
	elif args.stage == 'fast': config['global']['skip_stages'] = "pre,wt"
	elif args.stage == 'sim': config['global']['skip_stages'] = "pre,svs,lift,wt,hic"
	elif args.stage == 'sim+': config['global']['skip_stages'] = "pre,svs"
	elif args.stage == 'evo': config['global']['skip_stages'] = "lift,wt"
	elif args.stage == 'fast_evo': config['global']['skip_stages'] = "pre,lift,wt"
	elif args.stage == 'lift': config['global']['skip_stages'] = "pre,svs,sim,wt,hic"
	elif args.stage == 'lift+': config['global']['skip_stages'] = "pre,svs,sim"
	elif args.stage == 'wt': config['global']['skip_stages'] = "pre,svs,sim,lift,hic"
	elif args.stage == 'wt+': config['global']['skip_stages'] = "pre,svs,sim,lift"
	elif args.stage == 'hic': config['global']['skip_stages'] = "pre,svs,sim,lift,wt"
	else: 
		raise NameError(args.stage,'''is the incorrect stage id!
	Use one from: pre pre+ fast SVs SVs+ sim sim+ lift lift+ wt wt+ hic
	Use no for ignore''')

	try: os.remove(config['global']['log_file'])
	except OSError: pass
	try: os.makedirs(config['global']['work_dir'])
	except FileExistsError: pass
	except KeyError: raise KeyError('The work directory in global section is not defined')

	try: resolution = config['global']['resolution']
	except KeyError: config['global']['resolution'] = 'NO'
	try: resolution_low = config['global']['resolution_low']
	except KeyError: config['global']['resolution_low'] = 'NO'
	try: resolution_pab = config['global']['resolution_pab']
	except KeyError: config['global']['resolution_pab'] = 'NO'
	
	try: global_noised = gf.boolean(config['global']['one_as_null'])
	except KeyError: global_noised = False
	try: coverage_treshold = int(config['global']['coverage_treshold'])
	except KeyError: coverage_treshold = 0
	try: coverage_low_treshold = int(config['global']['coverage_low_treshold'])
	except KeyError: coverage_low_treshold = 0
	try: heterozygous = gf.boolean(config['global']['heterozygous'])
	except KeyError: heterozygous = True
	try:
		path_to_user_func = gf.boolean(config['global']['path_to_user_functions'])
		with open(path_to_user_func,'r') as f: user_file = f.read()
		exec(user_file)
	except KeyError: path_to_user_func = False
	
	try:
		global_count = gf.boolean(config['global']['contact_count'])
		if heterozygous: global_count = str(int(global_count)//2)
	except KeyError: global_count = "NO"
	
	try: skip_stages = set(config['global']['skip_stages'].split(','))
	except KeyError: skip_stages = set([])

	try: cleaning = gf.boolean(config['global']['cleaning'])
	except KeyError: 
		cleaning = True
		config['global']['cleaning'] = "YES"
	try: 
		log_file = config['global']['log_file']
		f = open(log_file,'w')
		f.close()
	except KeyError: log_file = False
		
	start_time = timeit.default_timer()
	
	##############################
	#PREPROCESSING STAGE PRE/PRE+#
	##############################

	if args.pre: 
		print( args.pre )
		commandLineConfig['preprocessing'] = {}
		key_value = args.pre.split()
		for kv in key_value: 
			key,value = kv.split('=')
			commandLineConfig['preprocessing'][key] = value
			config['preprocessing'][key] = value

	try: sim_id = config['preprocessing']['reference_id']
	except KeyError: 
		sim_id = config['global']['reference_id']
		config['preprocessing']['reference_id'] = sim_id
	try: work_dir = config['preprocessing']['work_dir']
	except KeyError:
		work_dir = config['global']['work_dir']
		config['preprocessing']['work_dir'] = work_dir
	try: chrom_sizes = config['preprocessing']['chrom_sizes']
	except KeyError:
		chrom_sizes = config['global']['chrom_sizes']
		config['preprocessing']['chrom_sizes'] = chrom_sizes
	try: resolution = config['preprocessing']['resolution']
	except KeyError: 
		resolution = config['global']['resolution']
		config['preprocessing']['resolution'] = resolution
	try: resolution_low = config['preprocessing']['resolution_low']
	except KeyError:
		resolution_low = config['global']['resolution_low']
		config['preprocessing']['resolution_low'] = resolution_low
	try: resolution_pab = config['preprocessing']['resolution_pab']
	except KeyError:
		resolution_pab = config['global']['resolution_pab']
		config['preprocessing']['resolution_pab'] = resolution_pab
	try: capture = config['preprocessing']['capture']
	except KeyError: capture = False
	try: path_to_hic_map = config['preprocessing']['path_to_hic_map']
	except KeyError: path_to_hic_map = False
	try: normalization = config['preprocessing']['normalization']
	except KeyError: normalization = 'NONE'
	try: path_to_contact_dump = config['preprocessing']['path_to_contact_dump']
	except KeyError: path_to_contact_dump = False
	try: path_to_juicertools = config['preprocessing']['path_to_juicertools']
	except KeyError: 
		try: path_to_juicertools = config['global']['path_to_juicertools']
		except KeyError: path_to_juicertools = False
	try: path_to_java_dir = config['preprocessing']['path_to_java_dir']
	except KeyError:
		try: path_to_java_dir = config['global']['path_to_java_dir']
		except KeyError: path_to_java_dir = ''
	
	try:
		user_coverage_statistic_func_name = gf.boolean(config['preprocessing']['user_coverage_statistic_func_name'])
		try:
			print('''Using user defined function for coverage statistic calculation''')
			exec("user_coverage_statistic_func = %s" % user_coverage_statistic_func_name)
		except NameError:
			raise NameError('''!!!The unknown name of user defined function: %s.''' % user_coverage_statistic_func_name)
	except KeyError: user_coverage_statistic_func = False
	try:
		user_distance_dependent_statistic_func_name = gf.boolean(config['preprocessing']['user_distance_dependent_statistic_func_name'])
		try:
			print('''Using user defined function for distance dependent statistic calculation''')
			exec("user_distance_dependent_statistic_func = %s" % user_distance_dependent_statistic_func_name)
		except NameError:
			raise NameError('''!!!The unknown name of user defined function: %s.''' % user_distance_dependent_statistic_func_name)
	except KeyError: user_distance_dependent_statistic_func = False

	try: log_file = config['preprocessing']['log_file']
	except KeyError:
		log_file = config['global']['log_file']
		config['preprocessing']['log_file'] = log_file
	
	name_res = '%s/pre/%s/%s.%s' % (work_dir,sim_id,sim_id,resolution)
	config['simulation']['resolution'] = resolution
	config['simulation']['contact_dir'] = name_res
	config['simulation']['coverage_file'] = name_res + '.binCov'
	config['simulation']['distance_file'] = name_res + '.stat'
	name_low = '%s/pre/%s/%s.%s' % (work_dir,sim_id,sim_id,resolution_low)
	config['simulation']['resolution_low'] = resolution_low
	config['simulation']['contact_low'] = name_low
	config['simulation']['coverage_low'] = name_low + '.binCov'
	config['simulation']['distance_low'] = name_low + '.stat'
	name_pab,name_cov_pab,name_dist_pab ='','',''
	config['simulation']['resolution_pab'] = resolution_pab
	for pab in resolution_pab.split(','):
		name_pab += '%s/pre/%s/pab.%s.%s,' % (work_dir,sim_id,sim_id,pab)
		name_cov_pab += '%s/pre/%s/pab.%s.%s.binCov,' % (work_dir,sim_id,sim_id,pab)
		name_dist_pab += '%s/pre/%s/pab.%s.%s.stat,' % (work_dir,sim_id,sim_id,pab)
	config['simulation']['contact_pab'] = name_pab[:-1]
	config['simulation']['coverage_pab'] = name_cov_pab[:-1]
	config['simulation']['distance_pab'] = name_dist_pab[:-1]
	config['preprocessing']['reference_id'] = sim_id
	try: config['simulation']['simulation_id']
	except KeyError: config['simulation']['simulation_id'] = sim_id
	
	if 'pre' in skip_stages: gf.printlog('Stage "pre" skipped', log_file)
	else:
		gf.printlog('Stage "pre" - data preprocessing...', log_file)
		from charm_func import pre
		name_res, name_low, name_pab = pre.preprocessing(sim_id, chrom_sizes, resolution, resolution_low, resolution_pab,
			capture, work_dir, path_to_hic_map, normalization, path_to_contact_dump,
			path_to_java_dir, path_to_juicertools, log_file, cleaning,
			user_coverage_statistic_func, user_distance_dependent_statistic_func)
	elp = timeit.default_timer() - start_time
	gf.printlog('... end of stage "pre" %.2f' % elp, log_file)


	############################
	#BILD SV MAP STAGE SVs/SVs+#
	############################
	
	if args.svs: 
		print( args.svs )
		commandLineConfig['SVs'] = {}
		key_value = args.svs.split()
		for kv in key_value: 
			key,value = kv.split('=')
			commandLineConfig['SVs'][key] = value
			config['SVs'][key] = value
	
	if args.stage == 'SVs': stand_alone = True
	else: stand_alone = False
	
	try: chrom_sizes = config['SVs']['chrom_sizes']
	except KeyError: chrom_sizes = config['global']['chrom_sizes']
	try: resolution = config['SVs']['resolution']
	except KeyError: resolution = config['global']['resolution']
	try: work_dir = config['SVs']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
	try: log_file = config['SVs']['log_file']
	except KeyError: log_file = config['global']['log_file']
	
	try: cross_species = gf.boolean(config['SVs']['cross_species'])
	except KeyError: cross_species = False
	
	if cross_species:
		path_to_maf = config['SVs']['path_to_maf_files']
		rname = config['SVs']['reference_id']
		try: chrom_sizes_ref = config['SVs']['chrom_sizes_reference']
		except KeyError: chrom_sizes_ref = chrom_sizes
		qname = config['SVs']['query_id']
		chrom_sizes_qu = config['SVs']['chrom_sizes_query']
	else: 
		path_to_svs_list = config['SVs']['path_to_svs_list']
		sname = config['SVs']['simulation_id']


	if 'svs' in skip_stages: gf.printlog('Stage "SVs" skipped', log_file)
	else:
		from charm_func import sv_maps as sm
		gf.printlog('Stage "SVs" - SV descriptions preparing...', log_file)
		
		if cross_species == False: Map_data = sm.generate_SV_map(path_to_svs_list, sname, resolution, chrom_sizes, work_dir, stand_alone,log_file)
		else: Map_data = sm.generate_SV_map_from_maf(path_to_maf, rname, chrom_sizes_ref, qname, chrom_sizes_qu,  resolution, work_dir, stand_alone, log_file)
		
		if args.stage in ['pre+','SVs+','fast','no','evo+','fast_evo']:
			chosen_chroms,add_pairs,map_SV_from_ref,pointviews,map_SV_to_ref,chrom_sizes_SV = Map_data
			try: config['simulation']['pointviews']
			except KeyError: config['simulation']['pointviews'] = pointviews
			try: config['simulation']['chrom_sizes_to']
			except KeyError: config['simulation']['chrom_sizes_to'] = chrom_sizes_SV
			try: config['simulation']['chosen_chroms_to']
			except KeyError: config['simulation']['chosen_chroms_to'] = chosen_chroms[1].strip()
			try: config['simulation']['map_file']
			except KeyError: config['simulation']['map_file'] = map_SV_from_ref
			try: config['liftover']['chosen_chroms_to']
			except KeyError: config['liftover']['chosen_chroms_to'] = chosen_chroms[0].strip()
			try: config['liftover']['map_file']
			except KeyError: config['liftover']['map_file'] = map_SV_to_ref
			elp = timeit.default_timer() - start_time
			gf.printlog('\tthe rearranged chromosome pairs %s' % chosen_chroms[0], log_file)
			gf.printlog('\tthe resulted chromosome pairs %s' % chosen_chroms[1], log_file)
			if add_pairs[0]:
				try: config['liftover']['add_pairs_to']
				except KeyError:
					config['liftover']['add_pairs_to'] = add_pairs[0].strip()
					gf.printlog('\tthe additional chromosome pairs from %s' % add_pairs[0], log_file)
				try: config['simulation']['add_pairs_to']
				except KeyError:
					config['simulation']['add_pairs_to'] = add_pairs[1].strip()
					gf.printlog('\tthe additional chromosome pairs to %s' % add_pairs[1], log_file)
	elp = timeit.default_timer() - start_time
	gf.printlog('... end of stage "SVs" %.2f' % elp, log_file)

	##################################
	#MUTANT SIMULATION STAGE sim/sim+#
	##################################

	if args.sim:
		print( args.sim )
		commandLineConfig['simulation'] = {}
		key_value = args.sim.split()
		for kv in key_value: 
			key,value = kv.split('=')
			commandLineConfig['simulation'][key] = value
			config['simulation'][key] = value

	try: sim_id = 'in_mut.' + config['simulation']['simulation_id']
	except KeyError:
		config['simulation']['simulation_id'] = config['SVs']['simulation_id']
		sim_id = 'in_mut.' + config['SVs']['simulation_id']
	try: work_dir = config['simulation']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']

	try: resolution = config['simulation']['resolution']
	except KeyError: resolution = config['global']['resolution']
	
	try: contact_dir = config['simulation']['contact_dir']
	except KeyError: 
		config['simulation']['contact_dir'] = '%s/pre/%s/%s.%s' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],resolution)
		contact_dir = config['simulation']['contact_dir']
	try: coverage_file = config['simulation']['coverage_file']
	except KeyError:
		config['simulation']['coverage_file'] = '%s/pre/%s/%s.%s.binCov' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],resolution)
		coverage_file  = config['simulation']['coverage_file']
	try: distance_file = config['simulation']['distance_file']
	except KeyError: 
		config['simulation']['distance_file'] = '%s/pre/%s/%s.%s.stat' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],resolution)
		distance_file = config['simulation']['distance_file']

	try: contact_low = config['simulation']['contact_low']
	except KeyError: contact_low = False
	try: coverage_low = config['simulation']['coverage_low']
	except KeyError: coverage_low = False
	try: distance_low = config['simulation']['distance_low']
	except KeyError: distance_low = False
	
	try:
		resolution_low = config['simulation']['resolution_low']
		if contact_low == False:
			config['simulation']['contact_low'] ='%s/pre/%s/%s.%s' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],resolution_low)
			contact_low = config['simulation']['contact_low']
		if coverage_low == False:
			config['simulation']['coverage_low'] = '%s/pre/%s/%s.%s.binCov' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],resolution_low)
			coverage_low = config['simulation']['coverage_low']
		if distance_low == False:
			config['simulation']['distance_low'] = '%s/pre/%s/%s.%s.stat' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],resolution_low)
			distance_low = config['simulation']['distance_low']
	except KeyError: resolution_low = config['global']['resolution_low']
	
	try: contact_pab = config['simulation']['contact_pab']
	except KeyError: contact_pab = ''
	try: coverage_pab = config['simulation']['coverage_pab']
	except KeyError: coverage_pab = ''
	
	try:
		resolution_pab = config['simulation']['resolution_pab']
		if contact_pab: pass
		else:
			for pab in resolution_pab.split(','):
				contact_pab += '%s/pre/%s/pab.%s.%s,' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],pab)
			config['simulation']['contact_pab'] = contact_pab[:-1]
			contact_pab = config['simulation']['contact_pab']
		if coverage_pab: pass
		else:
			for pab in resolution_pab.split(','):
				coverage_pab +='%s/pre/%s/pab.%s.%s.binCov,' % (config['global']['work_dir'],config['global']['reference_id'],config['global']['reference_id'],pab)
			config['simulation']['coverage_pab'] = coverage_pab[:-1]
			coverage_pab = config['simulation']['coverage_pab']
	except KeyError: resolution_pab = config['global']['resolution_pab']
	
	try: path_to_user_func = gf.boolean(config['simulation']['path_to_user_functions'])
	except KeyError:
		try: path_to_user_func = gf.boolean(config['global']['path_to_user_functions'])
		except KeyError: path_to_user_func == False
	
	if 'sim' in skip_stages and 'lift' in skip_stages : pass
	else:
		try: chrom_sizes_from = config['simulation']['chrom_sizes_from']
		except KeyError: 
			chrom_sizes_from = config['global']['chrom_sizes']
			config['simulation']['chrom_sizes_from'] = chrom_sizes_from

		chrom_sizes_to = config['simulation']['chrom_sizes_to']
		chosen_chroms_to = config['simulation']['chosen_chroms_to'].strip()
		try: add_pairs = gf.boolean(config['simulation']['add_pairs_to'])
		except KeyError: add_pairs = False
	
	model = config['simulation']['model']
	random_func_name = config['simulation']['random']
	try: random_func = gf.define_random_func(random_func_name)
	except NameError:
		try: exec("random_func = %s" % random_func_name)
		except NameError: raise NameError('The function %s not found' % random_func_name)
	
	try: contact_count = config['simulation']['contact_count']
	except KeyError:
		try:
			config['simulation']['contact_count'] = global_count
			contact_count = config['simulation']['contact_count']
		except NameError: raise NameError('The value "contact_count" in section [global] or [simulation] is not defined!')
	
	try: predict_null_contacts_func_name = config['simulation']['predict_null_contacts']
	except KeyError: predict_null_contacts_func_name = 'no'
	try: predict_null_contacts_func = gf.default_predict_function(predict_null_contacts_func_name)
	except NameError:
		try: exec("predict_null_contacts_func = %s" % predict_null_contacts_func_name)
		except NameError: raise NameError('The function %s not found' % predict_null_contacts_func_name)
	
	try: pick_contacts_func_name = config['simulation']['pick_contacts']
	except KeyError: 
		try: pick_contacts_func_name = gf.define_pick_from_predict(predict_null_contacts_func_name)
		except NameError: raise NameError('The parameter "pick_contacts" absent and the default value for "predict_null_contacts" = %s is not found' % predict_null_contacts_func_name)
	try: pick_contacts_func = gf.default_pick_function(pick_contacts_func_name)
	except NameError:
		try: exec("pick_contacts_func = %s" % pick_contacts_func_name)
		except NameError: raise NameError('The function %s not found' % pick_contacts_func_name)

	noised = global_noised
	pair = False
	try: log_file = config['simulation']['log_file']
	except KeyError: log_file = config['global']['log_file']

	sim_dir = '%s/mdl/%s' % (work_dir,sim_id)

	elp = timeit.default_timer() - start_time
	if 'sim' in skip_stages: gf.printlog('Stage "sim" skipped', log_file)
	else:
		from charm_func import sim
		cc_to = [c.split(',') for c in chosen_chroms_to.split(';')]
		cc_ts = set()
		for i in cc_to: cc_ts |= set(i)

		if add_pairs: add_pairs = [c.split(',') for c in add_pairs.split(';')]
		map_file = config['simulation']['map_file']
		pointviews = gf.boolean(config['simulation']['pointviews'])
		
		shutil.rmtree( sim_dir, ignore_errors=True )
		os.makedirs( sim_dir )

		gf.printlog('Stage "sim" - the simulation of contacts in mutant genome...', log_file)
		gf.printlog('\tStep 0: chromosome indexing...',log_file)
		l2i_from = gf.ChromIndexing(chrom_sizes_from)
		if chrom_sizes_to: l2i_to = gf.ChromIndexing(chrom_sizes_to)
		else: l2i_to = l2i_from
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...chromosome indexed, %.2fs' % elp, log_file)

		gf.printlog('\tStep 1: common statistic reading...', log_file)
		contactStatistic = sim.read_Contact_Statistics(
			coverage_file, distance_file,
			coverage_low, distance_low,
			coverage_pab, coverage_treshold, coverage_low_treshold,
			l2i_from, work_dir, log_file
			)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... end data reading %.2fs'% elp, log_file)

		if add_pairs: cct = cc_to + add_pairs
		else: cct = cc_to
		
		lc = len(cct)
		for i in range(lc):
			c1,c2 = cct[i]
			gf.printlog('\tStep %i/%i Simulation contacts on the mutant chromosome pair %s %s %i/%i' % (i,lc+2,c1,c2,i+1,lc), log_file)
			gf.printlog('\tStep %i.1: Reading mark points...' % i, log_file)
			
			if c1 == c2:
				MarkPoints1, MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c1],log_file)
				MarkPoints2, MarkPointsLow2,ccf2 = False, False, ccf1
			elif (c1 in cc_ts) and (c2 in cc_ts):
				MarkPoints1, MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c1],log_file)
				MarkPoints2, MarkPointsLow2,ccf2 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c2],log_file)
			elif (c1 in cc_ts) and (c2 not in cc_ts):
				MarkPoints1, MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c1],log_file)
				if add_pairs:
					from charm_func import sim_func as sf
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution)
					MarkPoints2 = sf.createUntouchedMarkPoints([c2],c2s_to,l2i_to,l2i_from)
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution_low)
					MarkPointsLow2 = sf.createUntouchedMarkPoints([c2],c2s_to,l2i_to,l2i_from)
					ccf2 = [c2]
				else: MarkPoints2,MarkPointsLow2 = False,False
			elif (c1 not in cc_ts) and (c2 in cc_ts):
				MarkPoints1,MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c2],log_file)
				if add_pairs:
					from charm_func import sim_func as sf
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution)
					MarkPoints2 = sf.createUntouchedMarkPoints([c1],c2s_to,l2i_to,l2i_from)
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution_low)
					MarkPointsLow2 = sf.createUntouchedMarkPoints([c1],c2s_to,l2i_to,l2i_from)
					ccf2 = [c1]
				else:  MarkPoints2,MarkPointsLow2 = False,False
			else: pass
			
			chosen_chroms_from = []
			for k1 in range(len(ccf1)):
					for k2 in range(len(ccf2)): chosen_chroms_from += [[ccf1[k1],ccf2[k2]]]
			
			if MarkPoints1:
				gf.printlog('\tchosen_chroms_from pairs %s' % (chosen_chroms_from), log_file)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t...mark point red for %.2f sec' % elp, log_file)
	
				gf.printlog('\tStep %i.2: data reading...' % i, log_file)
				contactData = sim.read_Contact_Data(
					contact_dir, resolution,
					contact_low, resolution_low,
					contact_pab, chosen_chroms_from,
					l2i_from, work_dir, log_file
					)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t... end data reading %.2fs'% elp, log_file)
				
				gf.printlog('\tStep %i.3: The starting of contact simulation...' % i, log_file)
				sim_dir = sim.sv_Simulation(
					contactData+contactStatistic, resolution, resolution_low, resolution_pab, MarkPoints1, MarkPointsLow1,
					l2i_from, l2i_to, (c1,c2), pointviews,
					model, contact_count, random_func, predict_null_contacts_func, pick_contacts_func, noised, 
					add_pairs, MarkPoints2, MarkPointsLow2,
					sim_id, work_dir, log_file, path_to_user_func
					)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t...end of contact simulation %.2fs'% elp, log_file)
			else:
				elp = timeit.default_timer() - start_time
				gf.printlog('\t...no markpoints, no contact simulation %.2fs'% elp, log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('... end of stage "sim" %.2f' % elp, log_file)

	if 'lift' in skip_stages: pass
	else:
		config['liftover']['contact_dir'] = sim_dir
		config['liftover']['distance_file'] = distance_file
		config['liftover']['resolution'] = resolution
		config['liftover']['contact_low'] = str(sim_dir)
		config['liftover']['distance_low'] = str(distance_low)
		config['liftover']['chrom_sizes_from'] = str(chrom_sizes_to)
		config['liftover']['chrom_sizes_to'] = str(chrom_sizes_from)
		config['liftover']['simulation_id'] = config['simulation']['simulation_id']
	
	if 'wt' in skip_stages: pass
	else:
		try: config['wild_type']['simulation_id']
		except KeyError: config['wild_type']['simulation_id'] = '%s.%s' % (config['global']['reference_id'],config['simulation']['predict_null_contacts'])
		try: config['wild_type']['contact_count']
		except KeyError: config['wild_type']['contact_count'] = config['simulation']['contact_count']
		try: config['wild_type']['model']
		except KeyError: config['wild_type']['model'] = model
		try: config['wild_type']['random']
		except KeyError: config['wild_type']['random'] = random_func_name
		try: config['wild_type']['predict_null_contacts']
		except KeyError: config['simulation']['predict_null_contacts'] = predict_null_contacts_func_name
		try: config['wild_type']['pick_contacts']
		except KeyError: config['simulation']['pick_contacts'] = pick_contacts_func_name

	########################################
	#LIFTOVER TO REFERENCE STAGE lift/lift+#
	########################################

	if args.lift: 
		print( args.lift )
		commandLineConfig['liftover'] = {}
		key_value = args.lift.split()
		for kv in key_value: 
			key,value = kv.split('=')
			commandLineConfig['liftover'][key] = value
			config['liftover'][key] = value

	try: sim_id = 'to_ref.' + config['liftover']['simulation_id']
	except KeyError: 
		try: sim_id = 'to_ref.' + config['simulation']['simulation_id']
		except KeyError: sim_id = 'to_ref.' + config['SVs']['simulation_id']
	try: work_dir = config['liftover']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']

	sim_dir = '%s/mdl/%s' % (work_dir,sim_id)
	chosen_chroms_to = 'all'

	if 'lift' in skip_stages: gf.printlog('Stage "lift" skipped', log_file)
	else:
		from charm_func import sim
		
		contact_dir = config['liftover']['contact_dir']
		coverage_file = False
		distance_file = config['liftover']['distance_file']
		try: resolution = config['liftover']['resolution']
		except KeyError: resolution = config['global']['resolution']
		contact_low = config['liftover']['contact_dir']
		coverage_low = False
		distance_low = config['liftover']['distance_low']
		try: resolution_low = config['liftover']['resolution_low']
		except KeyError:
			try: resolution_low = config['global']['resolution_low']
			except KeyError: resolution_low = False
		contact_pab = False
		coverage_pab = False
		distance_pab = False
		resolution_pab = False
		
		chrom_sizes_from = config['liftover']['chrom_sizes_from']
		try: chrom_sizes_to = config['liftover']['chrom_sizes_to']
		except KeyError: chrom_sizes_to = config['global']['chrom_sizes']
		
		path_to_user_func == False
		
		chosen_chroms_to = config['liftover']['chosen_chroms_to']
		cc_to = [c.split(',') for c in chosen_chroms_to.split(';')]
		cc_ts = set()
		for i in cc_to: cc_ts |= set(i)
		
		map_file = config['liftover']['map_file']
		pointviews = False
		model = 'easy'
		random_func = False
		contact_count = False
		predict_null_contacts_func = False
		pick_contacts_func = False
		noised = False
		try: add_pairs = gf.boolean(config['liftover']['add_pairs_to'])
		except KeyError: add_pairs = False
		if add_pairs: add_pairs = [c.split(',') for c in add_pairs.split(';')]
		try: log_file = config['liftover']['log_file']
		except KeyError: log_file = config['global']['log_file']

		shutil.rmtree( sim_dir, ignore_errors=True )
		os.makedirs( sim_dir )
		gf.printlog('Stage "lift" - the contact liftovering to the reference genome...', log_file)
		gf.printlog('\tStep 0: chromosome indexing...',log_file)
		l2i_from = gf.ChromIndexing(chrom_sizes_from)
		if chrom_sizes_to: l2i_to = gf.ChromIndexing(chrom_sizes_to)
		else: l2i_to = l2i_from

		elp = timeit.default_timer() - start_time
		gf.printlog('\t...chromosome indexed, %.2fs' % elp, log_file)
		
		gf.printlog('\tStep 1: common statistic reading...', log_file)
		contactStatistic = sim.read_Contact_Statistics(
			coverage_file, distance_file,
			coverage_low, distance_low,
			coverage_pab, 0, 0,
			l2i_from, work_dir, log_file
			)
		
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... end data reading %.2fs'% elp, log_file)
		if add_pairs: cct = cc_to + add_pairs
		else: cct = cc_to
		
		lc = len(cct)
		for i in range(lc):
			c1,c2 = cct[i]
			gf.printlog('\tLiftover contacts to the chromosome pair %s %s %i/%i' % (c1,c2,i+1,lc), log_file)
			gf.printlog('\tStep 1: Reading mark points...', log_file)
			
			if c1 == c2:
				MarkPoints1, MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c1],log_file)
				MarkPoints2, MarkPointsLow2,ccf2 = False, False, ccf1
			elif (c1 in cc_ts) and (c2 in cc_ts):
				MarkPoints1, MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c1],log_file)
				MarkPoints2, MarkPointsLow2,ccf2 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c2],log_file)
			elif (c1 in cc_ts) and (c2 not in cc_ts):
				MarkPoints1, MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c1],log_file)
				if add_pairs:
					from charm_func import sim_func as sf
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution)
					MarkPoints2 = sf.createUntouchedMarkPoints([c2],c2s_to,l2i_to,l2i_from)
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution_low)
					MarkPointsLow2 = sf.createUntouchedMarkPoints([c2],c2s_to,l2i_to,l2i_from)
					ccf2 = [c2]
				else: MarkPoints2,MarkPointsLow2 = False,False
			elif (c1 not in cc_ts) and (c2 in cc_ts):
				MarkPoints1,MarkPointsLow1,ccf1 = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,l2i_to,[c2],log_file)
				if add_pairs:
					from charm_func import sim_func as sf
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution)
					MarkPoints2 = sf.createUntouchedMarkPoints([c1],c2s_to,l2i_to,l2i_from)
					c2s_to = gf.ChromSizes(chrom_sizes_to,resolution_low)
					MarkPointsLow2 = sf.createUntouchedMarkPoints([c1],c2s_to,l2i_to,l2i_from)
					ccf2 = [c1]
				else:  MarkPoints2,MarkPointsLow2 = False,False
			else: pass
			
			chosen_chroms_from = []
			for k1 in range(len(ccf1)):
					for k2 in range(len(ccf2)): chosen_chroms_from += [[ccf1[k1],ccf2[k2]]]
			
			if MarkPoints1:
				gf.printlog('\tchosen_chroms_from pairs %s' % (chosen_chroms_from), log_file)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t...mark point red for %.2f sec' % elp, log_file)
	
				gf.printlog('\tStep 2: data reading...', log_file)
				contactData = sim.read_Contact_Data(
					contact_dir, resolution,
					contact_low, resolution_low,
					contact_pab, chosen_chroms_from,
					l2i_from, work_dir, log_file
					)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t... end data reading %.2fs'% elp, log_file)
				
				gf.printlog('\tStep 3: The starting of contact liftovering...', log_file)
				sim_dir = sim.sv_Simulation(
					contactData+contactStatistic, resolution, resolution_low, resolution_pab, MarkPoints1, MarkPointsLow1,
					l2i_from, l2i_to, (c1,c2), pointviews,
					model, contact_count, random_func, predict_null_contacts_func, pick_contacts_func, noised,
					add_pairs, MarkPoints2, MarkPointsLow2,
					sim_id, work_dir, log_file, path_to_user_func
					)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t...end of contact liftovering %.2fs'% elp, log_file)
			else:
				elp = timeit.default_timer() - start_time
				gf.printlog('\t...no markpoints, no contact liftovering %.2fs'% elp, log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('... end of stage "lift" %.2f' % elp, log_file)
		if cleaning: shutil.rmtree( contact_dir, ignore_errors=True )
	if 'hic' in skip_stages: pass
	else:
		try: config['hic']['svs_contacts']
		except KeyError: config['hic']['svs_contacts'] = sim_dir
		try: config['hic']['chosen_chroms']
		except KeyError: config['hic']['chosen_chroms'] = chosen_chroms_to
		try: config['hic']['resolution'] 
		except KeyError: config['hic']['resolution'] = resolution
	
	

	################################
	#WILD-TYPE REPLICA STAGE wt/wt+#
	################################
	
	if args.wt: 
		print( args.wt )
		commandLineConfig['wild_type'] = {}
		key_value = args.wt.split()
		for kv in key_value: 
			key,value = kv.split('=')
			commandLineConfig['wild_type'][key] = value
			config['wild_type'][key] = value

	try: sim_id = config['wild_type']['simulation_id']
	except KeyError: 
		try: sim_id = '%s.%s' % (config['global']['reference_id'],config['simulation']['predict_null_contacts'])
		except KeyError: sim_id = config['global']['reference_id']
	try: replica_ids = config['wild_type']['replica_ids']
	except KeyError: replica_ids = '0,1'
	try: work_dir = config['wild_type']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
	try: contact_count = config['wild_type']['contact_count']
	except KeyError: 
		try: contact_count = config['simulation']['contact_count']
		except KeyError: 
			try: contact_count = global_count
			except NameError: raise NameError('The value "contact_count" in sections [global], or [simulation], or [wt] is not defined!')

	if 'wt' in skip_stages: gf.printlog('Stage "wt" skipped', log_file)
	else:
		from charm_func import sim
		
		try: contact_dir = config['wild_type']['contact_dir']
		except KeyError: contact_dir = config['simulation']['contact_dir']
		try: coverage_file = config['wild_type']['coverage_file']
		except KeyError: coverage_file = config['simulation']['coverage_file']
		try: distance_file = config['wild_type']['distance_file']
		except KeyError: distance_file = config['simulation']['distance_file']
		try: resolution = config['wild_type']['resolution']
		except KeyError: resolution = config['global']['resolution']
		try: contact_low = config['wild_type']['contact_low']
		except KeyError: contact_low = config['simulation']['contact_low']
		try: coverage_low = config['wild_type']['coverage_low']
		except KeyError: coverage_low = config['simulation']['coverage_low']
		try: distance_low = config['wild_type']['distance_low']
		except KeyError: distance_low = config['simulation']['distance_low']
		try: resolution_low = config['wild_type']['resolution_low']
		except KeyError: resolution_low = config['global']['resolution_low']
		try: contact_pab = config['wild_type']['contact_pab']
		except KeyError: contact_pab = config['simulation']['contact_pab']
		try: coverage_pab = config['wild_type']['coverage_pab']
		except KeyError: coverage_pab = config['simulation']['coverage_pab']
		try: resolution_pab = config['wild_type']['resolution_pab']
		except KeyError: resolution_pab = config['global']['resolution_pab']
		try: chrom_sizes = config['wild_type']['chrom_sizes_from']
		except KeyError: chrom_sizes = config['global']['chrom_sizes']
		try: chosen_chroms = gf.boolean(config['wild_type']['chosen_chroms'].strip())
		except KeyError: chosen_chroms = config['liftover']['chosen_chroms_to'].strip()
		try: model = config['wild_type']['model']
		except KeyError: model = config['simulation']['model']
		
		try: random_func_name = config['wild_type']['random']
		except KeyError: random_func_name = config['simulation']['random']
		try: random_func = gf.define_random_func(random_func_name)
		except NameError:
			try: exec("random_func = %s" % random_func_name)
			except NameError: raise NameError('The function %s not found' % random_func_name)
		try: predict_null_contacts_func_name = config['simulation']['predict_null_contacts']
		except KeyError: predict_null_contacts_func_name = 'no'
		try: predict_null_contacts_func = gf.default_predict_function(predict_null_contacts_func_name)
		except NameError:
			try: exec("predict_null_contacts_func = %s" % predict_null_contacts_func_name)
			except NameError: raise NameError('The function %s not found' % predict_null_contacts_func_name)
		
		try: pick_contacts_func_name = config['simulation']['pick_contacts']
		except KeyError: 
			try: pick_contacts_func_name = gf.define_pick_from_predict(predict_null_contacts_func_name)
			except NameError: raise NameError('The parameter "pick_contacts" absent and the default value for "predict_null_contacts" = %s is not found' % predict_null_contacts_func_name)
		try: pick_contacts_func = gf.default_pick_function(pick_contacts_func_name)
		except NameError:
			try: exec("pick_contacts_func = %s" % pick_contacts_func_name)
			except NameError: raise NameError('The function %s not found' % pick_contacts_func_name)
		noised = global_noised
		
		try: log_file = config['wild_type']['log_file']
		except KeyError: log_file = config['global']['log_file']
		
		gf.printlog('Stage "wt" - wild-type replica generation - start...',log_file)
		gf.printlog('\tStep 0: chromosome indexing...',log_file)
		l2i = gf.ChromIndexing(chrom_sizes)
		c2s_low = gf.ChromSizes(chrom_sizes,resolution_low)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...chromosome indexed, %.2fs' % elp, log_file)
		
		gf.printlog('\tStart replicas simulation', log_file)
		if replica_ids:
			currentStep,stepCount,currentReplica,replicaCount = 0,0,0,len(replica_ids.split(','))
			for replica_id in replica_ids.split(','):
				chroms,pair_list = [],[]
				if chosen_chroms == 'all' or chosen_chroms == False:
					lnc = len(c2s_low.keys())
					stepCount += lnc*(lnc+1)/2
				else: stepCount += len(chosen_chroms.split(';'))
			for replica_id in replica_ids.split(','):
				currentReplica += 1
				wt_name = '%s/wt/%s/%s/%s' % (work_dir,sim_id,contact_count,replica_id)
				try: os.makedirs( wt_name )
				except FileExistsError: pass
				chroms,pair_list = [],[]
				if chosen_chroms == 'all' or chosen_chroms == False:
					chroms = sorted(c2s_low.keys())
					for ci in range(len(chroms)):
						i = chroms[ci]
						for cj in range(ci,len(chroms)):
							j = chroms[cj]
							if l2i[i] <= l2i[j]: c1_c2 = [chroms[ci],chroms[cj]]
							else: c1_c2 = [chroms[cj],chroms[ci]]
							pair_list.append(c1_c2)
				else: pair_list = [c.split(',') for c in chosen_chroms.split(';')]

				for c1_c2 in pair_list:
					currentStep += 1
					gf.printlog('\t\tStep %i.1/%i: %s data reading...' % (currentStep,stepCount,c1_c2), log_file)
					contactStatistic = sim.read_Contact_Statistics(
						coverage_file, distance_file,
						coverage_low, distance_low,
						coverage_pab, coverage_treshold, coverage_low_treshold,
						l2i, work_dir, log_file
					)
					contactData = sim.read_Contact_Data(
						contact_dir, resolution,
						contact_low, resolution_low,
						contact_pab, [c1_c2],
						l2i, work_dir, log_file
					)
					elp = timeit.default_timer() - start_time
					gf.printlog('\t\t...end data reading, %.2fs' % elp, log_file)
					gf.printlog('\t\tStep %i.2/%i: replicas simulation %s...' % (currentStep,stepCount,c1_c2), log_file)
					sim.wt_Simulation(
							contactData+contactStatistic, resolution, resolution_low, resolution_pab,
							c1_c2, c2s_low, l2i,
							model, contact_count, random_func, predict_null_contacts_func, pick_contacts_func, noised,
							sim_id, replica_id, work_dir, log_file, path_to_user_func
							)
					elp = timeit.default_timer() - start_time
					gf.printlog('\t\tend simulation %i/%i: %.2fs' % (currentStep,stepCount,elp), log_file)
				elp = timeit.default_timer() - start_time
				gf.printlog('\t\tend %s (%i/%i) replica simulation %.2fs' % (replica_id,currentReplica,replicaCount, elp), log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...end of stage "wt" %.2f' % elp, log_file)
		
	wt_path = []
	if replica_ids and ('hic' not in skip_stages):
		for replica_id in replica_ids.split(','):
			wt_path.append('%s/wt/%s/%s/%s' % (work_dir,sim_id,contact_count,replica_id))
		try: config['hic']['wt1_contacts'] 
		except KeyError: config['hic']['wt1_contacts'] = wt_path[0]
		try: config['hic']['wt2_contacts']
		except KeyError: config['hic']['wt2_contacts'] = wt_path[1]

	##################################
	#HIC-MAP GENERATION STAGE hic/hic+#
	##################################

	if args.hic: 
		print( args.hic )
		commandLineConfig['hic'] = {}
		key_value = args.hic.split()
		for kv in key_value: 
			key,value = kv.split('=')
			commandLineConfig['hic'][key] = value
			config['hic'][key] = value

	try: sim_id = config['hic']['simulation_id']
	except KeyError:
		try: sim_id = config['simulation']['simulation_id']
		except KeyError: sim_id = config['SVs']['simulation_id']
	try: work_dir = config['hic']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
	try: resolution = config['hic']['resolution']
	except KeyError: resolution = config['global']['resolution']

	format = config['hic']['format']
	
	try: path_to_java_dir = config['hic']['path_to_java_dir']
	except KeyError:
		try: path_to_java_dir = config['global']['path_to_java_dir']
		except KeyError: path_to_java_dir = ''
	try: path_to_juicertools = config['hic']['path_to_juicertools']
	except KeyError: 
		try: path_to_juicertools = config['global']['path_to_juicertools']
		except KeyError: path_to_juicertools = False
 
	try: hic_resolutions = config['hic']['hic_resolutions']
	except KeyError: hic_resolutions = False
	try: cleaning = gf.boolean(config['hic']['cleaning'])
	except KeyError: cleaning = gf.boolean(config['global']['cleaning'])
	try: log_file = config['hic']['log_file']
	except KeyError: log_file = config['global']['log_file']
	
	if 'hic' in skip_stages: gf.printlog('Stage "hic" skipped', log_file)
	else:
		from charm_func import contact2hic as c2h
		try: svs_contacts = config['hic']['svs_contacts']
		except KeyError: svs_contacts = False
		try: wt1_contacts = config['hic']['wt1_contacts']
		except KeyError: wt1_contacts = False
		if heterozygous:
			try: wt2_contacts = config['hic']['wt2_contacts']
			except KeyError: wt2_contacts = False
		else: wt2_contacts = False
		try: chosen_chroms = gf.boolean(config['hic']['chosen_chroms'])
		except KeyError: chosen_chroms = config['simulation']['chosen_chroms_from']
		try: chrom_sizes = config['hic']['chrom_sizes']
		except KeyError: chrom_sizes = config['global']['chrom_sizes']
	
		gf.printlog('Stage "hic" - hic map generation - start', log_file)
		gf.printlog('\tStep 1: the generation of pre/hic files',log_file)
		c2h.hic_generate(svs_contacts, wt1_contacts, wt2_contacts,
			chosen_chroms, chrom_sizes, resolution, #capture,
			format, path_to_java_dir, path_to_juicertools, hic_resolutions,
			sim_id, work_dir, log_file, cleaning
			)
		elp = timeit.default_timer() - start_time
		gf.printlog('... end "hic" stage, %.2f' % elp, log_file)

	elp = timeit.default_timer() - start_time
	gf.printlog('The End, %.2f' % elp, log_file)
	print(commandLineConfig)