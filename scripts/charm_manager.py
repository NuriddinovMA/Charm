import os
import timeit
import argparse
from configparser import ConfigParser, ExtendedInterpolation
import global_func as gf

parser = argparse.ArgumentParser(description='Parameters for generation of rearrangement maps')
parser.add_argument( '-i', dest='ini',metavar='ini-file', help='path to ini-file')
parser.add_argument('-S', dest='stage', metavar='Stage',default='pre+',help='must be one of "pre+", "SVs+", "simulation+", "liftover+",  "wt", "final"')
args = parser.parse_args()

config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(args.ini)

if args.stage in ['pre','pre+','SVs','SVs+','sim','sim+','lift','lift+','hic']: pass
else:
	print(args.stage,
'''is the incorrect stage id!
Use one from: pre pre+ SVs SVs+ sim sim+ lift lift+ hic''')
	exit()

try: os.remove(config['global']['log_file'])
except OSError: pass

try: skip_stages = config['global']['skip_stages'].split(',')
except KeyError: kip_stages = []

start_time = timeit.default_timer()
if args.stage in ['pre','pre+']:
	
	import pre
	
	try: sim_id = config['preprocessing']['simulation_id']
	except KeyError: 
		sim_id = config['global']['simulation_id']
		config['preprocessing']['simulation_id'] = sim_id
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
	try: path_to_hic = config['preprocessing']['path_to_hic']
	except KeyError: path_to_hic = False
	try: normalization = config['preprocessing']['normalization']
	except KeyError: normalization = 'NONE'
	try: path_to_hic_dump = config['preprocessing']['path_to_hic_dump']
	except KeyError: path_to_hic_dump = False
	try: path_to_juicer = config['preprocessing']['path_to_juicer']
	except KeyError: path_to_juicer = False
	try: path_to_java_dir = config['global']['path_to_java_dir']
	except KeyError: path_to_java_dir = False
	try: log_file = config['preprocessing']['log_file']
	except KeyError:
		log_file = config['global']['log_file']
		config['preprocessing']['log_file'] = log_file
	
	if 'pre' in skip.stages: gf.printlog('Stage "pre" skipped', log_file)
	else:
		gf.printlog('Stage "pre" - data preprocessing...', log_file)
		name_res, name_low, name_pab = pre.preprocessing(sim_id, chrom_sizes, resolution, resolution_low, resolution_pab,
			capture, work_dir, path_to_hic, normalization, path_to_hic_dump, path_to_java_dir, path_to_juicer, log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('... end of stage "pre" %.2f' % elp, log_file)
		if args.stage == 'pre+':
			config['simulation']['contact_dir'] = name_res
			config['simulation']['coverage_file'] = name_res + '.binCov'
			config['simulation']['distance_file'] = name_res + '.stat'
			config['simulation']['contact_low'] = name_low
			config['simulation']['coverage_low'] = name_low + '.binCov'
			config['simulation']['distance_low'] = name_low + '.stat'
			config['simulation']['contact_pab'] = name_pab
			config['simulation']['coverage_pab'] = name_pab + '.binCov'
			config['simulation']['distance_pab'] = name_pab + '.stat'
			config['preprocessing']['simulation_id'] = sim_id
			#config['hic']['capture'] = capture

if args.stage in ['pre+','SVs','SVs+']:

	import sv_maps as sm
	
	if args.stage == 'SVs': stand_alone = True
	else: stand_alone = False
	
	try: chrom_sizes = config['SVs']['chrom_sizes']
	except KeyError: chrom_sizes = config['global']['chrom_sizes']
	try: resolution = config['SVs']['resolution']
	except KeyError: resolution = config['global']['resolution']
	try: work_dir = config['SVs']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
	path_to_svs_list = config['SVs']['path_to_svs_list']
	
	if 'SVs' in skip.stages: gf.printlog('Stage "SVs" skipped', log_file)
	else:
		gf.printlog('Stage "SVs" - SV descriptions preparing...', log_file)
		Map_data = sm.generate_SV_map(chrom_sizes, resolution, path_to_svs_list, work_dir, stand_alone)
		elp = timeit.default_timer() - start_time
		gf.printlog('... end of stage "SVs" %.2f' % elp, log_file)
		if args.stage in ['pre+','SVs+']:
			map_SV_from_ref,chosen_chroms_from,pointviews,map_SV_to_ref,chosen_chroms_to,chrom_sizes_SV = Map_data
			config['simulation']['chosen_chroms_from'] = chosen_chroms_from.strip()
			config['simulation']['pointviews'] = pointviews
			config['simulation']['chrom_sizes_to'] = chrom_sizes_SV
			config['simulation']['chosen_chroms_to'] = chosen_chroms_to.strip()
			config['simulation']['map_file'] = map_SV_from_ref
			config['liftover']['map_file'] = map_SV_to_ref

import sim

if args.stage in ['pre+','SVs+','sim','sim+']:
	
	try: sim_id = 'in_mut.' + config['simulation']['simulation_id']
	except KeyError: sim_id = 'in_mut.' + config['global']['simulation_id']
	try: work_dir = config['simulation']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
	contact_dir = config['simulation']['contact_dir']
	coverage_file = config['simulation']['coverage_file']
	distance_file = config['simulation']['distance_file']
	try: resolution = config['simulation']['resolution']
	except KeyError: resolution = config['global']['resolution']
	contact_low = config['simulation']['contact_low']
	coverage_low = config['simulation']['coverage_low']
	distance_low = config['simulation']['distance_low']
	try: resolution_low = config['simulation']['resolution_low']
	except KeyError: resolution_low = config['global']['resolution_low']
	contact_pab = config['simulation']['contact_pab']
	coverage_pab = config['simulation']['coverage_pab']
	try: resolution_pab = config['simulation']['resolution_pab']
	except KeyError: resolution_pab = config['global']['resolution_pab']
	try: chrom_sizes_from = config['simulation']['chrom_sizes_from']
	except KeyError: 
		chrom_sizes_from = config['global']['chrom_sizes']
		config['simulation']['chrom_sizes_from'] = chrom_sizes_from
	chosen_chroms_from = config['simulation']['chosen_chroms_from'].strip()
	chrom_sizes_to = config['simulation']['chrom_sizes_to']
	chosen_chroms_to = config['simulation']['chosen_chroms_to'].strip()
	map_file = config['simulation']['map_file']
	pointviews = config['simulation']['pointviews']
	model = config['simulation']['model']
	random = config['simulation']['random']
	contact_count = config['simulation']['contact_count']
	predict_null_contacts = config['simulation']['predict_null_contacts']
	try: log_file = config['simulation']['log_file']
	except KeyError: log_file = config['global']['log_file']

	elp = timeit.default_timer() - start_time
	if 'sim' in skip.stages: gf.printlog('Stage "sim" skipped', log_file)
	else:
		gf.printlog('Stage "sim" - the simulation of contacts in mutant genome...', log_file)
		gf.printlog('\tStep 0: chromosome indexing...',log_file)
		l2i_from = gf.ChromIndexing(chrom_sizes_from)
		if chrom_sizes_to: l2i_to = gf.ChromIndexing(chrom_sizes_to)
		else: l2i_to = l2i_from
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...chromosome indexed, %.2fs' % elp, log_file)

		gf.printlog('\tStep 1: data reading...', log_file)
		contactData = sim.read_Contact_Data(
			contact_dir, coverage_file, distance_file, resolution,
			contact_low, coverage_low, distance_low, resolution_low,
			contact_pab, coverage_pab, chosen_chroms_from,
			l2i_from, work_dir, log_file
			)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... end data reading %.2fs'% elp, log_file)
		gf.printlog('\tStep 2: Reading mark points...', log_file)
		MarkPoints,MarkPointsLow = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,chosen_chroms_from,l2i_to,chosen_chroms_to,log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...%i mark point red for %.2f sec' % (len(MarkPoints), elp), log_file)
		
		gf.printlog('\tStep 3: The starting of contact simulation...', log_file)
		sim_dir = sim.sv_Simulation(
			contactData, resolution, resolution_low, resolution_pab, MarkPoints, MarkPointsLow,
			l2i_from, chosen_chroms_from, l2i_to, chosen_chroms_to, pointviews,
			model, contact_count, random, predict_null_contacts,
			sim_id, work_dir, log_file
			)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... end of contact simulation %.2fs'% elp, log_file)
		
		if args.stage in ['pre+','sim+']:
			config['liftover']['contact_dir'] = sim_dir
			config['liftover']['coverage_file'] = name_res + '.binCov'
			config['liftover']['distance_file'] = name_res + '.stat'
			config['liftover']['resolution'] = resolution
			config['liftover']['contact_low'] = sim_dir
			config['liftover']['distance_low'] = name_low + '.stat'
			config['liftover']['resolution_low'] = resolution_low
			config['liftover']['chrom_sizes_from'] = chrom_sizes_to
			config['liftover']['chosen_chroms_from'] = chosen_chroms_to
			config['liftover']['chrom_sizes_to'] = chrom_sizes_from
			config['liftover']['chosen_chroms_to'] = chosen_chroms_from

if args.stage in ['pre+','SVs+','sim+','lift','lift+']:
	
	try: sim_id = 'to_ref.' + config['liftover']['simulation_id']
	except KeyError: sim_id = 'to_ref.' + config['global']['simulation_id']
	try: work_dir = config['liftover']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
	contact_dir = config['liftover']['contact_dir']
	coverage_file = config['liftover']['coverage_file']
	distance_file = config['liftover']['distance_file']
	try: resolution = config['liftover']['resolution']
	except KeyError: resolution = config['global']['resolution']
	contact_low = config['liftover']['contact_dir']
	coverage_low = False
	distance_low = config['liftover']['distance_low']
	try: resolution_low = config['liftover']['resolution_low']
	except KeyError: resolution_low = config['global']['resolution_low']
	contact_pab = False
	coverage_pab = False
	distance_pab = False
	resolution_pab = False
	chrom_sizes_from = config['liftover']['chrom_sizes_from']
	chosen_chroms_from = config['liftover']['chosen_chroms_from']
	try: chrom_sizes_to = config['liftover']['chrom_sizes_to']
	except KeyError: chrom_sizes_to = config['global']['chrom_sizes']
	chosen_chroms_to = config['liftover']['chosen_chroms_to']
	map_file = config['liftover']['map_file']
	pointviews = False
	model = 'easy'
	random = 'norandom'
	contact_count = False
	predict_null_contacts = False
	try: log_file = config['liftover']['log_file']
	except KeyError: log_file = config['global']['log_file']
	
	if 'lift' in skip.stages: gf.printlog('Stage "lift" skipped', log_file)
	else:
		gf.printlog('Stage "lift" - the contact liftovering to the reference genome...', log_file)
		gf.printlog('\tStep 0: chromosome indexing...',log_file)
		l2i_from = gf.ChromIndexing(chrom_sizes_from)
		if chrom_sizes_to: l2i_to = gf.ChromIndexing(chrom_sizes_to)
		else: l2i_to = l2i_from

		elp = timeit.default_timer() - start_time
		gf.printlog('\t...chromosome indexed, %.2fs' % elp, log_file)

		gf.printlog('\tStep 1: data reading...', log_file)
		contactData = sim.read_Contact_Data(
			contact_dir, coverage_file, distance_file, resolution,
			contact_low, coverage_low, distance_low, resolution_low,
			contact_pab, coverage_pab, chosen_chroms_from,
			l2i_from, work_dir, log_file
		)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t... end data reading %.2fs'% elp, log_file)
		gf.printlog('\tStep 2: Reading mark points...', log_file)
		MarkPoints,MarkPointsLow = sim.read_RearMap(map_file,resolution,resolution_low,l2i_from,chosen_chroms_from,l2i_to,chosen_chroms_to,log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...%i mark point readed for %.2f sec' % (len(MarkPoints), elp), log_file)
		
		gf.printlog('\tStep 3: The starting of contact liftovering...', log_file)
		sim_dir = sim.sv_Simulation(
			contactData, resolution, resolution_low, resolution_pab, MarkPoints, MarkPointsLow,
			l2i_from, chosen_chroms_from, l2i_to, chosen_chroms_to, pointviews,
			model, contact_count, random, predict_null_contacts,
			sim_id, work_dir, log_file
			)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...end of contact liftover and simulation %.2fs'% elp, log_file)
		if args.stage in ['pre+','sim+','lift+']: 
			config['hic']['svs_contacts'] = sim_dir
			config['hic']['chosen_chroms'] = chosen_chroms_to
			config['hic']['resolution'] = resolution
			
if args.stage in ['pre+','sim+','lift','lift+','wt']:

	try: sim_id = config['wild_type']['simulation_id']
	except KeyError: sim_id = config['global']['simulation_id']
	try: replica_ids = config['wild_type']['replica_ids']
	except KeyError: replica_ids = '0,1'
	try: work_dir = config['wild_type']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
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
	except KeyError: config['simulation']['coverage_low']
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
	try: chosen_chroms = config['wild_type']['chosen_chroms'].strip()
	except KeyError: chosen_chroms = config['simulation']['chosen_chroms_from'].strip()
	try: model = config['wild_type']['model']
	except KeyError: model = config['simulation']['model']
	try: random = config['wild_type']['random']
	except KeyError: random = config['simulation']['random']
	try: contact_count = config['wild_type']['contact_count']
	except KeyError: contact_count = config['simulation']['contact_count']
	try: predict_null_contacts = config['wild_type']['predict_null_contacts']
	except KeyError: predict_null_contacts = config['simulation']['predict_null_contacts']
	try: log_file = config['wild_type']['log_file']
	except KeyError: log_file = config['global']['log_file']
	
	resolution,resolution_low,resolution_pab = int(resolution),int(resolution_low),int(resolution_pab)
	if 'wt' in skip.stages: gf.printlog('Stage "wt" skipped', log_file)
	else:
		gf.printlog('Stage "wt" - wild-type replica generation - start...',log_file)
		gf.printlog('\tStep 0: chromosome indexing...',log_file)
		l2i = gf.ChromIndexing(chrom_sizes)
		c2s_low = gf.ChromSizes(chrom_sizes,resolution_low)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...chromosome indexed, %.2fs' % elp, log_file)
		
		gf.printlog('\tStep 1: data reading...', log_file)
		contactData = sim.read_Contact_Data(
			contact_dir, coverage_file, distance_file, resolution,
			contact_low, coverage_low, distance_low, resolution_low,
			contact_pab, coverage_pab, chosen_chroms,
			l2i, work_dir, log_file
		)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...end data reading, %.2fs' % elp, log_file)
		gf.printlog('\tStep 2: replicas simulation...', log_file)
		wt_path = []
		if replica_ids:
			for replica_id in replica_ids.split(','):
				wt_path.append( 
					sim.wt_Simulation(
					contactData, resolution, resolution_low, resolution_pab,
					chosen_chroms, c2s_low, l2i,
					model, contact_count, random, predict_null_contacts,
					sim_id, replica_id, work_dir, log_file
					)
				)
			elp = timeit.default_timer() - start_time
			gf.printlog('end %s replica simulation %.2fs' % (replica_id, elp), log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...end replicas simulation %.2fs'% elp, log_file)
		if args.stage in ['pre+','SVs+','sim+','lifover+']:
			config['hic']['wt1_contacts'] = 'NO'
			config['hic']['wt2_contacts'] = 'NO'
			if replica_ids:
				config['hic']['wt1_contacts'] = wt_path[0]
				try: config['hic']['wt2_contacts'] = wt_path[1]
				except IndexError: pass

if args.stage in ['pre+','SVs+','sim+','lift','lifover+','hic']:
	import contact2hic as c2h
	try: sim_id = config['hic']['simulation_id']
	except KeyError: sim_id = config['global']['simulation_id']
	try: work_dir = config['hic']['work_dir']
	except KeyError: work_dir = config['global']['work_dir']
	try: resolution = config['hic']['resolution']
	except KeyError: resolution = config['global']['resolution']
	try: svs_contacts = config['hic']['svs_contacts']
	except KeyError: svs_contacts = False
	try: wt1_contacts = config['hic']['wt1_contacts']
	except KeyError: wt1_contacts = False
	try: wt2_contacts = config['hic']['wt2_contacts']
	except KeyError: wt2_contacts = False
	chosen_chroms = config['hic']['chosen_chroms']
	try: chrom_sizes = config['hic']['chrom_sizes']
	except KeyError: chrom_sizes = config['global']['chrom_sizes']
	#capture =
	format = config['hic']['format']
	try: path_to_java_dir = config['hic']['path_to_java_dir']
	except KeyError: path_to_java_dir = config['global']['path_to_java_dir']
	try: path_to_juicertools = config['hic']['path_to_juicertools']
	except KeyError: path_to_juicertools = config['global']['path_to_juicertools']
	hic_resolutions = config['hic']['hic_resolutions'] 
	try: log_file = config['hic']['log_file']
	except KeyError: log_file = config['global']['log_file']
	if 'hic' in skip.stages: gf.printlog('Stage "hic" skipped', log_file)
	else:
		gf.printlog('Stage "hic" - hic map generation - start', log_file)
		gf.printlog('\tStep 1: the generation of pre/hic files',log_file)
		c2h.hic_generate(svs_contacts,wt1_contacts,wt2_contacts,
			chosen_chroms, chrom_sizes, resolution, #capture,
			format, path_to_java_dir, path_to_juicertools, hic_resolutions,sim_id,work_dir,log_file
			)
		elp = timeit.default_timer() - start_time
		gf.printlog('... end "hic" stage, %.2f' % elp, log_file)

elp = timeit.default_timer() - start_time
gf.printlog('The End, %.2f' % elp, log_file)
