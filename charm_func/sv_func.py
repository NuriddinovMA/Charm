def markgenerate(MT,resolution,name,output_dir):

	f1 = open('%s/%s.%s.%i.mark' % (output_dir, name[0],name[1],resolution), 'w')
	f2 = open('%s/%s.%s.%i.mark' % (output_dir, name[1],name[0],resolution), 'w')
	f1.write('chrm_ref\tpos_ref_1\tpos_ref_2\tchrm_mut\tpos_mut_1\tpos_mut_2\tpoint number\n')
	f2.write('chrm_mut\tpos_mut_1\tpos_mut_2\tchrm_ref\tpos_ref_1\tpos_ref_2\tpoint number\n')
	r,k = resolution//2,0
	for key in MT:
		for i in range(len(MT[key])):
			k += 1
			chrm = key
			if key[-1] == '|': chrm = key[:-1]
			f1.write('%s\t%i\t%i\t%s\t%i\t%i\t%i\n' % (MT[key][i][0], MT[key][i][1]*resolution+r-50,MT[key][i][1]*resolution+r+50, chrm, i*resolution+r-50,i*resolution+r+50, k) )
			f2.write('%s\t%i\t%i\t%s\t%i\t%i\t%i\n' % (chrm, i*resolution+r-50,i*resolution+r+50, MT[key][i][0], MT[key][i][1]*resolution+r-50,MT[key][i][1]*resolution+r+50, k) )
	f1.close()
	f2.close()

