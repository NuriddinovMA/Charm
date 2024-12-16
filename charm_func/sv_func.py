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
	
def maf2mark(maf_path, resolution, rname, qname, chromSizes_ref, chromSizes_qu, output_dir):
	import os
	try:
		maf_list = os.listdir(maf_path)
		mafs = [f'{maf_path}/{maf}' for maf in maf_list]
	except OSError: mafs = [maf_path]
	
	Maf_seq,Syn_blocks = [],[]
	for maf in mafs:
		with open(maf,'r') as m: lines = m.readlines()
		for i in range(len(lines)):
			if lines[i][0] == 'a':
				c1,p11,ln = lines[i+1].split()[1:4]
				p11,ln = int(p11),int(ln)
				p12 = p11 + ln
				c2,p21,ln,dir,end = lines[i+2].split()[1:6]
				p21,ln,end = int(p21),int(ln),int(end)
				dir = int(dir+'1')
				if dir > 0: p22 = p21 + ln
				else: 
					p22 = end - p21
					p21 = p22 - ln
				Maf_seq.append((c1,p11,p12,c2,p21,p22,dir))
		del lines
	Maf_seq.sort()

	# f = open('%s/%s.%s.%i.synm' % (output_dir, name[0],name[1],resolution), 'w')
	# for c1,p11,p12,c2,p21,p22,dir in Maf_seq: f.write(f'{c1} {p11} {p12} {c2} {p21} {p22} {dir}\n')
	# f.close()
	
	Syn_blocks.append(Maf_seq[0])
	for c1,p11,p12,c2,p21,p22,dir in Maf_seq[1:]:
		cs1,ps11,ps12,cs2,ps21,ps22,dirs = Syn_blocks[-1]
		gap1 = p11 - ps11
		if dir > 0: gap2 = p21-ps22
		else: gap2 = ps21 - p22
		if gap1 < 100: gap1 = 100
		if gap2 < 100: gap2 = 100
		if c1 == cs1 and c2 == cs2 and dir == dirs and gap1 < 1000000 and gap2 < 1000000 and gap1/gap2 > 0.5:
			if dir > 0: Syn_blocks[-1] = cs1,ps11,p12,cs2,ps21,p22,dirs
			else: Syn_blocks[-1] = cs1,ps11,p12,cs2,p21,ps22,dirs
		else: Syn_blocks.append((c1,p11,p12,c2,p21,p22,dir))
	
	f = open('%s/%s.%s.%i.synb' % (output_dir, rname,qname,resolution), 'w')
	for c1,p11,p12,c2,p21,p22,dir in Syn_blocks: f.write(f'{c1} {p11} {p12} {c2} {p21} {p22}\n')
	f.close()
	
	Mark = []
	res = resolution//2
	n = 0
	chrm_from = set([])
	chrm_to = set([])
	for i in range(len(Maf_seq)):
		c1,p11,p12,c2,p21,p22,dir = Maf_seq[i]
		cs1,ps11,ps12,cs2,ps21,ps22,dirs = Maf_seq[i-1]
		chrm_from.add(c1)
		chrm_to.add(c2)
		gap1 = p11 - ps12
		if dir > 0: gap2 = p21-ps22
		else: gap2 = ps21 - p22
		if gap1 < 100: gap1 = 100
		if gap2 < 100: gap2 = 100
		if i > 0 and c1 == cs1 and c2 == cs2 and dir == dirs and resolution//2 < gap1 < 1000000 and resolution//2 < gap2 < 1000000 and ( ( (gap1/gap2 > 0.5) and (gap2/gap1 > 0.5) ) or ( abs(gap1-gap2) < 100000) ):
			step = gap1//resolution + 2
			ln1,ln2 = gap1//step,gap2//step
			for j in range(1,step):
				n += 1
				if dir > 0: Mark.append((cs1,ps12+j*ln1-100,ps12+j*ln1+100,cs2,ps22+j*ln2-100,ps22+j*ln2+100,dirs, n, 'fill'))
				else: Mark.append((cs1,ps12+j*ln1-100,ps12+j*ln1+100,cs2,ps21-j*ln2-100,ps21-j*ln2+100,dirs, n, 'fill'))
		elif i > 0 and p11 > resolution and ps21 > resolution and p21 > resolution:
			if c1 == cs1 and gap1 > resolution*2 : step = 3
			else: step = 2
			for j in range(1,step):
				n += 1
				if dirs > 0: Mark.append((cs1,ps12+j*res-100,ps12+j*res+100,cs2,ps22+j*res-100,ps22+j*res+100,dirs, n, 'extend'))
				else: Mark.append((cs1,ps12+j*res-100,ps12+j*res+100,cs2,ps21-j*res-100,ps21-j*res+100,dirs, n, 'extend'))
			for j in range(step-1,0,-1):
				n += 1
				if dir > 0: Mark.append((c1,p11-j*res-100,p11-j*res+100,c2,p21-j*res-100,p21-j*res+100,dir, n, 'extend'))
				else: Mark.append((c1,p11-j*res-100,p11-j*res+100,c2,p22+j*res-100,p22+j*res+100,dir, n, 'extend'))
		elif i == 0 and p11 > resolution and p21 > resolution:
			step = 3
			for j in range(step-1,0,-1):
				n += 1
				if dir > 0: Mark.append((c1,p11-j*res-100,p11-j*res+100,c2,p21-j*res-100,p21-j*res+100,dir, n, 'extend'))
				else: Mark.append((c1,p11-j*res-100,p11-j*res+100,c2,p22+j*res-100,p22+j*res+100,dir, n, 'extend'))
		else: pass
		
		step = (p12 - p11)//resolution + 2
		ln1,ln2 = (p12-p11)//step,(p22-p21)//step
		for j in range(1,step): 
			n += 1
			if dir > 0: Mark.append((c1,p11+j*ln1-100,p11+j*ln1+100,c2,p21+j*ln2-100,p21+j*ln2+100,dir, n, 'direct'))
			else: Mark.append((c1,p11+j*ln1-100,p11+j*ln1+100,c2,p22-j*ln2-100,p22-j*ln2+100,dir, n, 'direct'))
	
	for j in range(1,3):
		n += 1
		if dir > 0: Mark.append((c1,p11+j*res-100,p11+j*res+100,c2,p21+j*res-100,p21+j*res+100,dir, n, 'extend'))
		else: Mark.append((c1,p11+j*res-100,p11+j*res+100,c2,p22-j*res-100,p22-j*res+100,dir, n, 'extend'))
	
	map_name = '%s/%s.%s.%i.mark' % (output_dir, rname, qname, resolution)
	f = open( map_name, 'w')
	for c1,p11,p12,c2,p21,p22,dir,number,type in Mark: 
		if 0 < p11 and 0 < p21 and p12 < chromSizes_ref[c1] and p12 < chromSizes_qu[c2]: f.write(f'{c1} {p11} {p12} {c2} {p21} {p22} {dir} {number} {type}\n')
	f.close()
	
	
	cf, ct = '', ''
	for chr in chromSizes_ref: cf += f'{chr} '
	for chr in chromSizes_qu: ct += f'{chr} '
	
	return map_name,cf,ct

def maf2synblocks(maf_path, resolution, rname, qname, output_dir):
	import os
	try:
		maf_list = os.listdir(maf_path)
		mafs = [f'{maf_path}/{maf}' for maf in maf_list]
	except OSError: mafs = [maf_path]
	
	Maf_seq,Syn_blocks = [],[]
	for maf in mafs:
		with open(maf,'r') as m: lines = m.readlines()
		for i in range(len(lines)):
			if lines[i][0] == 'a':
				c1,p11,ln = lines[i+1].split()[1:4]
				p11,ln = int(p11),int(ln)
				p12 = p11 + ln
				c2,p21,ln,dir,end = lines[i+2].split()[1:6]
				p21,ln,end = int(p21),int(ln),int(end)
				dir = int(dir+'1')
				if dir > 0: p22 = p21 + ln
				else: 
					p22 = end - p21
					p21 = p22 - ln
				Maf_seq.append((c1,p11,p12,c2,p21,p22,dir))
		del lines
	Maf_seq.sort()

	f = open('%s/%s.%s.%i.synm' % (output_dir, rname,qname,resolution), 'w')
	for c1,p11,p12,c2,p21,p22,dir in Maf_seq: f.write(f'{c1} {p11} {p12} {c2} {p21} {p22} {dir}\n')
	f.close()
	
	for i in range(3):
		Syn_blocks = [Maf_seq[0]]
		for c1,p11,p12,c2,p21,p22,dir in Maf_seq[1:]:
			cs1,ps11,ps12,cs2,ps21,ps22,dirs = Syn_blocks[-1]
			gap1 = p11 - ps12
			if dir > 0: gap2 = p21-ps22
			else: gap2 = ps21 - p22
			if gap1 < 100: gap1 = 100
			if gap2 < 100: gap2 = 100
			if c1 == cs1 and c2 == cs2 and dir == dirs and gap1 < 1000000 and gap2 < 1000000 and ( ( (gap1/gap2 > 0.5) and (gap2/gap1 > 0.5) ) or ( abs(gap1-gap2) < 100000 ) ):
				if dir > 0: Syn_blocks[-1] = cs1,ps11,p12,cs2,ps21,p22,dirs
				else: Syn_blocks[-1] = cs1,ps11,p12,cs2,p21,ps22,dirs
			else: Syn_blocks.append((c1,p11,p12,c2,p21,p22,dir))
		Maf_seq = []
		for c1,p11,p12,c2,p21,p22,dir in Syn_blocks: 
			if p12-p11 > 50000: Maf_seq.append((c1,p11,p12,c2,p21,p22,dir))
		Maf_seq.sort()
	f = open('%s/%s.%s.%i.synb' % (output_dir, rname,qname,resolution), 'w')
	for c1,p11,p12,c2,p21,p22,dir in Syn_blocks: 
		if p12-p11 > 50000: f.write(f'{c1} {p11} {p12} {c2} {p21} {p22} {dir}\n')
	f.close()



