import os
import sys
import timeit
import pre_func as prf
import global_func as gf

def preprocessing(sim_name, chrom_sizes, resolution, resolution_low, resolution_pab,
	capture, work_dir, path_to_hic, norm, path_to_hic_dump, path_to_java_dir, path_to_juicer, log_file
	):
	
	start_time = timeit.default_timer()
	l2i = gf.ChromIndexing(chrom_sizes)
	capture = gf.boolean(capture)
	if resolution:
		resolution = int(resolution)
		gf.printlog('\nStatistic for %sbp resolution' % resolution, log_file)
		if capture: locus = capture[0],int(capture[1])/resolution,int(capture[2])/resolution
		suffixH = '%s.%ikb' % (sim_name,resolution/1000)
		try: os.makedirs(work_dir+'/pre/'+suffixH)
		except OSError: pass
		
		gf.printlog('\nStep 0: data preparing...', log_file)
		c2s = gf.ChromSizes(chrom_sizes,resolution)
		if path_to_hic and path_to_juicer and path_to_hic_dump == False:
			gf.printlog('\tDump contactcs from hic-map', log_file)
			chr_num = len(c2s)
			path_to_hic_dump = '%s/bcm/%s' % (work_dir,suffixH)
			try: os.makedirs(path_to_hic_dump)
			except OSError: pass
			command = "%s/java -jar %s dump observed %s %s %s %s BP %i %s/%s.%i.%s.%s.%s"
			for i in range(1,chr_num+1):
				for j in range(i,chr_num+1):
					gf.printlog(command % (path_to_java_dir, path_to_juicer, norm, path_to_hic, l2i[i],l2i[j],resolution,path_to_hic_dump,sim_name,resolution,l2i[i],l2i[j],norm) , log_file)
					os.system(command % (path_to_java_dir, path_to_juicer, norm, path_to_hic, l2i[i],l2i[j],resolution,path_to_hic_dump,sim_name,resolution,l2i[i],l2i[j],norm) )
			elp = timeit.default_timer() - start_time
			gf.printlog('\t...end dumping %.2fs' % elp, log_file)
		
		gf.printlog('\tGenome analysis...', log_file)
		counts = prf.diag_counts(c2s,capture=capture,log=log_file)
		maxd = max(counts.keys())
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...end genome analysing %.2fs' % elp, log_file)

		out_name = work_dir+'/pre/'+suffixH
		out_name_res = work_dir+'/pre/'+suffixH
		gf.printlog('\nStep 1: Calculating bin coverage...', log_file)
		binCov=prf.iBinCoverage(path_to_hic_dump,c2s,resolution,out=out_name,chrm_index=l2i,capture=capture,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('...bin coverage calculated for %.2fs' % elp, log_file)
		
		gf.printlog('\nStep 2: Distance depended statistics...', log_file)
		contactDistanceHash = prf.iDistanceRead(maxd,path=path_to_hic_dump,capture=capture,resolution=resolution,coverage=binCov,log=log_file)
		meanHash = prf.iMeaner( contactDistanceHash, counts, out_name,log=log_file)
		del contactDistanceHash
		elp = timeit.default_timer() - start_time
		gf.printlog('...distance analyzed for %.2fs' % elp, log_file)
		out_name = work_dir+'/pre/'+suffixH+'/'+suffixH
		
		gf.printlog('\nStep 3: Contact transforming by mean statistic...', log_file)
		prf.iTotalContactListing(meanHash,binCov,resolution,out_name,capture=capture,path=path_to_hic_dump,log=log_file)
		del meanHash
		del binCov
		elp = timeit.default_timer() - start_time
		gf.printlog('...contact transformed for %.2fs' % elp, log_file)

	if resolution_low and resolution:
		resolution_low = int(resolution_low)
		gf.printlog('\nStatistic for %sbp resoluion' % resolution_low, log_file)
		if capture: locus = capture[0],int(capture[1])/resolution_low,int(capture[2])/resolution_low
		suffixH = '%s.%ikb' % (sim_name,resolution_low/1000)
		try: os.makedirs(work_dir+'/pre/'+suffixH)
		except OSError: pass
		
		gf.printlog('\nStep 0.1: data preparing...', log_file)
		gf.printlog('\tGenome analysis...', log_file)
		
		c2s = gf.ChromSizes(chrom_sizes,resolution_low)
		counts = prf.diag_counts(c2s,capture=capture,log=log_file)
		maxd = max(counts.keys())
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...end genome analysing %.2fs' % elp, log_file)

		out_name = work_dir+'/pre/'+suffixH
		out_name_low = work_dir+'/pre/'+suffixH
		gf.printlog('\nStep 1.1: Calculating bin coverage...', log_file)
		binCov=prf.iBinCoverage(path_to_hic_dump,c2s,resolution_low,out=out_name,chrm_index=l2i,capture=capture,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('...bin coverage calculated for %.2fs' % elp, log_file)

		gf.printlog('\nStep 2.1: Distance depended statistics...', log_file)
		contactDistanceHash = prf.iDistanceRead(maxd,path=path_to_hic_dump,capture=capture,resolution=resolution_low,coverage=binCov,log=log_file)
		meanHash = prf.iMeaner( contactDistanceHash, counts, out_name,log=log_file)
		del contactDistanceHash
		elp = timeit.default_timer() - start_time
		
		gf.printlog('...distance analyzed for %.2fs' % elp, log_file)
		out_name = work_dir+'/pre/'+suffixH+'/'+suffixH
		gf.printlog('\nStep 3.1: Contact transforming by mean statistic...', log_file)
		prf.iTotalContactListing(meanHash,binCov,resolution_low,out_name,capture=capture,path=path_to_hic_dump,log=log_file)
		del meanHash
		del binCov
		elp = timeit.default_timer() - start_time
		gf.printlog('...contact transformed for %.2fs' % elp, log_file)

	if resolution_pab and resolution:
		resolution_pab = int(resolution_pab)
		gf.printlog('\nStatistic for pseudo AB-compartment %sbp resoluion' % resolution_pab, log_file)
		if capture: locus = capture[0],int(capture[1])/resolution_pab,int(capture[2])/resolution_pab
		suffixL = '/pab.%s.%ikb' % (sim_name,resolution_pab/1000)
		resolution_pab = resolution_pab//2
		try: os.makedirs(work_dir+'/pre/'+suffixL)
		except OSError: pass
		
		gf.printlog('\nStep 0.2: data preparing...', log_file)
		gf.printlog('\tPseudocompartment genome analysis...', log_file)
		c2s_ab = gf.ChromSizes(chrom_sizes,resolution_pab)
		counts_ab = prf.diag_counts(c2s_ab,log=log_file)
		maxd_ab = max(counts_ab.keys())
		elp = timeit.default_timer() - start_time
		gf.printlog('\t...end genome analysing %.2fs' % elp, log_file)

		out_name = work_dir+'/pre/'+suffixL
		out_name_pab = work_dir+'/pre/'+suffixL
		gf.printlog('\nStep 1.2: Calculating bin coverage for pseudocompartment resolution...', log_file)
		binCovAB=prf.iBinCoverage(path_to_hic_dump,c2s_ab,resolution_pab,out=out_name,chrm_index=l2i,capture=capture,log=log_file)
		elp = timeit.default_timer() - start_time
		gf.printlog('...bin coverage calculated for %.2fs' % elp, log_file)
		
		gf.printlog('\nStep 2.2: Distance depended statistics for pseudocompartment resolution...', log_file)
		abContacts = prf.iPsuedoAB(path_to_hic_dump,resolution_pab,log=log_file)
		contactDistanceHashAB = prf.iDistanceRead(maxd_ab,hash=abContacts,coverage=binCovAB,log=log_file)
		meanHashAB = prf.iMeaner( contactDistanceHashAB, counts_ab, out_name,log=log_file)
		del contactDistanceHashAB
		elp = timeit.default_timer() - start_time
		gf.printlog('...pseudocompartment resolution analyzed for %.2fs' % elp, log_file)

		out_name = work_dir+'/pre/'+suffixL+'/'+suffixL
		gf.printlog('\nStep 3.2: Pseudocompartments contact transforming by mean statistic...', log_file)
		prf.iTotalContactListing(meanHashAB,binCovAB,resolution_pab,out_name,hash=abContacts,log=log_file)
		del meanHashAB
		del binCovAB
		del abContacts
		elp = timeit.default_timer() - start_time
		gf.printlog('...Pseudocompartments transformed time %.2fs' % elp, log_file)

	elp = timeit.default_timer() - start_time
	gf.printlog('Full processing for %.2fs' % elp, log_file)
	return out_name_res,out_name_low,out_name_pab
