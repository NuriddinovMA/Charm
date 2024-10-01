import os
import numpy as np

def boolean(x):
	if str(x).lower() in ['false','no','n'] : return False
	elif str(x).lower() in ['true','yes','y']: return True
	else: return x

def printlog(str, logname):
	print(str)
	if logname:
		with open(logname, 'a') as log: log.write(str+'\n')

def ChromIndexing(path):
	ChrInd = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)):
		parse = lines[i].split()
		try: 
			ChrInd[parse[0]] = i+1
			ChrInd[i+1] = parse[0]
		except IndexError: break
	del lines
	return ChrInd

def ChromSizes(path,resolution):
	ChrSzs = {}
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)):
		parse = lines[i].split()
		try: ChrSzs[parse[0]] = int(parse[1])//int(resolution)+1
		except IndexError: break
	del lines
	return ChrSzs

def HashTry(Hash, key):
	t = True
	try: Hash[key]
	except KeyError: t = False
	return t

def FindKey(Hash, key):
	if key != 'interchromosome':
		for k in range(key,-1,-1):
			try: 
				Hash[k]
				return k
			except KeyError: pass
	else: return key

def boolean(x):
	
	if x in ['False', 'false', 'no','No','NO']: var = False
	elif x in [ 'True', 'true','yes','Yes','YES' ]: var = True
	else: var = x
	try:
		if var[0] in ['False', 'false', 'no','No','NO']: return False
		elif var[0] in [ 'True', 'true','yes','Yes','YES' ]: return True
		else:var = x
	except: pass
	return var
	
def colorList():
	color = []
	for i in range(0,256,85):
		for j in range(0,256,85):
			for k in range(0,256,85): color.append('%i,%i,%i' % (i,j,k))
	return color

def _pick_sq(cov1, cov2, distance_coef, **kwargs):
	return np.sqrt(cov1*cov2)*distance_coef
def _pick_mult(cov1, cov2, distance_coef, **kwargs):
	return 1.*cov1*cov2*distance_coef
def _pick_mixsq(cov1, cov2, distance_coef, **kwargs):
	contact_distance = kwargs['contact_distance']
	h = np.sign(contact_distance)
	h[h == 0] = 1
	return ((cov1+cov2)*distance_coef*(h+1) - 1.*np.sqrt(cov1*cov2)*distance_coef*(h-1))/2
def _pick_mixed(cov1, cov2, distance_coef, **kwargs):
	contact_distance = kwargs['contact_distance']
	h = np.sign(contact_distance)
	h[h == 0] = 1
	return ((cov1+cov2)*distance_coef*(h+1) - 1.*cov1*cov2*distance_coef*(h-1))/2
	#if contact_distance == -1000: return cov1*cov2*distance_coef
	#else: return 1.*cov1*cov2*distance_coef
	
def _pick_sum(cov1,cov2,distance_coef, **kwargs):
	return 1.*(cov1+cov2)*distance_coef
def _pick_ps(cov1, cov2, distance_coef, **kwargs):
	return 1.*distance_coef
def _pick_default(cov1, cov2, distance_coef, **kwargs):
	return 1
	
def _predict_model_pts(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = 1.
	pc = poe*distance_dependent_coef[contact_distance]['mean_contact']
	return pc,poe
def _predict_model_pts_ab(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = oe_AB
	pc = poe*distance_dependent_coef[contact_distance]['mean_contact']
	return pc,poe
def _predict_model_cov_mult_f(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = oe_AB*(cov1*cov2)/distance_independent_coef['mean_coverage_multiple']
	pc = poe*distance_dependent_coef[contact_distance]['mean_contact']
	return pc,poe
def _predict_model_cov_sq_f(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = oe_AB*np.sqrt((cov1*cov2)/distance_independent_coef['mean_coverage_multiple'])
	pc = poe*distance_dependent_coef[contact_distance]['mean_contact']
	return pc,poe
def _predict_model_cov_sum_f(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = oe_AB*(cov1+cov2)/(2*distance_independent_coef['mean_coverage_sum'])
	pc = poe*distance_dependent_coef[contact_distance]['mean_contact']
	return pc,poe
def _predict_model_cov_mult_f1(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = oe_AB*(cov1*cov2)/distance_dependent_coef[contact_distance]['mean_coverage_multiple']
	pc = poe*distance_dependent_coef[contact_distance]['no_null_mean']
	return pc,poe
def _predict_model_cov_sq_f1(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = oe_AB*np.sqrt((cov1*cov2)/distance_dependent_coef[contact_distance]['mean_coverage_multiple'])
	pc = poe*distance_dependent_coef[contact_distance]['no_null_mean']
	return pc,poe
def _predict_model_cov_sum_f1(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	poe = oe_AB*(cov1+cov2)/(2*distance_independent_coef['mean_coverage_sum'])
	pc = poe*distance_dependent_coef[contact_distance]['no_null_mean']
	return pc,poe
def _predict_model_cov_mixed_f(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	if contact_distance == -1000: poe = oe_AB*cov1*cov2/distance_independent_coef['mean_coverage_multiple']
	else: poe = oe_AB*(cov1+cov2)/(2*distance_independent_coef['mean_coverage_sum'])
	pc = poe*distance_dependent_coef[contact_distance]['mean_contact']
	return pc,poe
def _predict_model_cov_mixsq_f(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	if contact_distance == -1000: poe = oe_AB*np.sqrt((cov1*cov2)/distance_independent_coef['mean_coverage_multiple'])
	else: poe = oe_AB*(cov1+cov2)/(2*distance_independent_coef['mean_coverage_sum'])
	pc = poe*distance_dependent_coef[contact_distance]['mean_contact']
	return pc,poe
def _predict_model_cov_mixed_f1(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	if contact_distance == -1000: poe = oe_AB*cov1*cov2/distance_dependent_coef[contact_distance]['mean_coverage_multiple']
	else: poe = oe_AB*(cov1+cov2)/(2*distance_independent_coef['mean_coverage_sum'])
	pc = poe*distance_dependent_coef[contact_distance]['no_null_mean']
	return pc,poe
def _predict_model_cov_mixsq_f1(cov1, cov2, cont_AB, oe_AB, contact_distance, distance_dependent_coef, distance_independent_coef):
	if contact_distance == -1000: poe = oe_AB*np.sqrt((cov1*cov2)/distance_dependent_coef[contact_distance]['mean_coverage_multiple'])
	else: poe = oe_AB*(cov1+cov2)/(2*distance_independent_coef['mean_coverage_sum'])
	pc = poe*distance_dependent_coef[contact_distance]['no_null_mean']
	return pc,poe

def default_pick_function(pick_func_name):
	if pick_func_name == 'multiply': pick_func = _pick_mult
	elif pick_func_name == 'square_root': pick_func = _pick_sq
	elif pick_func_name == 'mixed': pick_func = _pick_mixed
	elif pick_func_name == 'mixsq': pick_func = _pick_mixsq
	elif pick_func_name == 'sum': pick_func = _pick_sum
	elif pick_func_name == 'distance_mean': pick_func = _pick_ps
	elif pick_func_name == 'default': pick_func = _pick_default
	else: raise NameError('The function %s is not default' % pick_func_name)
	return pick_func

def default_predict_function(pred_null_model_name):
	pred_null_model_name = boolean(pred_null_model_name)
	if pred_null_model_name == 'cov_mult_f': pred_null_model_func = _predict_model_cov_mult_f
	elif pred_null_model_name == 'cov_sq_f': pred_null_model_func = _predict_model_cov_sq_f
	elif pred_null_model_name == 'cov_sum_f': pred_null_model_func = _predict_model_cov_sum_f
	elif pred_null_model_name == 'cov_mult_f1': pred_null_model_func = _predict_model_cov_mult_f1
	elif pred_null_model_name == 'cov_sq_f1': pred_null_model_func = _predict_model_cov_sq_f1
	elif pred_null_model_name == 'cov_sum_f1': pred_null_model_func = _predict_model_cov_sum_f1
	elif pred_null_model_name == 'cov_mixed_f': pred_null_model_func = _predict_model_cov_mixed_f
	elif pred_null_model_name == 'cov_mixsq_f': pred_null_model_func = _predict_model_cov_mixsq_f
	elif pred_null_model_name == 'cov_mixed_f1': pred_null_model_func = _predict_model_cov_mixed_f1
	elif pred_null_model_name == 'cov_mixsq_f1': pred_null_model_func = _predict_model_cov_mixsq_f1
	elif pred_null_model_name == 'pts_ab': pred_null_model_func = _predict_model_pts_ab
	elif pred_null_model_name == 'pts': pred_null_model_func = _predict_model_pts
	elif boolean(pred_null_model_name) == False: pred_null_model_func = False
	else: raise NameError('The function %s is not default' % predict_null_model_name)
	return pred_null_model_func

def define_pick_from_predict(pred_null_model_name):
	if pred_null_model_name == 'cov_mult_f': 'multiply'
	elif pred_null_model_name == 'cov_sq_f': pick_func_name = 'square_root'
	elif pred_null_model_name == 'cov_sum_f': pick_func_name = 'sum'
	elif pred_null_model_name == 'cov_mult_f1': pick_func_name = 'multiply'
	elif pred_null_model_name == 'cov_sq_f1': pick_func_name = 'square_root'
	elif pred_null_model_name == 'cov_sum_f1': pick_func_name = 'sum'
	elif pred_null_model_name == 'cov_mixed_f': pick_func_name = 'mixed'
	elif pred_null_model_name == 'cov_mixsq_f': pick_func_name = 'mixsq'
	elif pred_null_model_name == 'cov_mixed_f1': pick_func_name = 'mixed'
	elif pred_null_model_name == 'cov_mixsq_f1': pick_func_name = 'mixsq'
	elif pred_null_model_name == 'pts_ab': pick_func_name = 'distance_mean'
	elif pred_null_model_name == 'pts': pick_func_name = 'distance_mean'
	elif boolean(pred_null_model_name) == False: pick_func_name = 'distance_mean'
	else: raise NameError('The function %s is not default' % predict_null_model_name)
	return pick_func_name

def _binomial(contact_array,old_counts,new_counts):
	new_array = np.random.binomial(new_counts,contact_array/old_counts)
	# except ValueError: #pass
		# print('VE1',h[contact_array<0])
		# print('VE2',h[contact_array>1])
		# print('VE3',h[(contact_array<1) & (contact_array>0)])
	return new_array
def _normal(contact_array,old_counts,new_counts):
	loc = np.sqrt(contact_array,old_counts)
	loc[loc == 0] = 1
	new_array = np.random.normal(contact_array,loc)
	return new_array
def _hypergeometric(contact_array,old_counts,new_counts):
	ngood = contact_array
	ngood[ngood < 1] = 1
	nsample = new_counts
	nbad = old_counts - ngood
	nsample,ngood,nbad = np.int32(np.round(nsample)),np.int32(np.round(ngood)),np.int32(np.round(nbad))
	new_array = np.random.hypergeometric(ngood,nbad,nsample)
	return new_array
def _round(contact_array,old_counts,new_counts):
	ngood = contact_array
	ngood[ngood < 1] = 1
	nsample = new_counts
	nbad = old_counts - ngood
	nsample,ngood,nbad = np.int32(np.round(nsample)),np.int32(np.round(ngood)),np.int32(np.round(nbad))
	new_array = np.round(1.*new_counts*contact_array/old_counts,decimals=8)
	return new_array
def _choice(contact_array,old_counts,new_counts):
	lnc = len(contact_array)
	new_array = np.zeros(lnc)
	idx = np.random.choice(lnc,size=int(new_counts),p=old_counts)
	# print('new1',idx)
	# print(new_counts)
	# print(old_counts)
	# print(new_array)
	# print(new_array[idx])
	idx,repeat = np.unique(idx,return_counts=True,axis=0)
	# print('new',new_array,repeat)
	new_array[idx] = repeat
	# print('new',new_array,repeat)
	return new_array


def define_random_func(random_func_name):
	if random_func_name == 'binomial': random_func = _binomial
	elif random_func_name == 'normal': random_func = _normal
	elif random_func_name == 'hypergeometric': random_func = _hypergeometric
	elif random_func_name == 'round': random_func = _round
	elif random_func_name == 'choice': random_func = _choice
	elif boolean(random_func_name) == False: random_func = False
	else: raise NameError('''!!!The unknown name of randomizing function!!!''')
	return random_func