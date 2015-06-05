from __future__ import division, with_statement, print_function, unicode_literals, absolute_import

lastmodified = "11 May 2015"

def memorypercent():
	try: import psutil
	except: return 0.0
	mem = 0.0
	for p in psutil.get_process_list():
		try:mem += p.get_memory_percent()
		except: pass
	return mem
	

def rank(sortby, handleties=0):
	""" return rank for each value in sortby, in the same order """
	zippedin = list(zip(sortby, list(range(len(sortby)))))
	if handleties:
		zippedin.sort()
		ranks = list(range(len(zippedin)))
		lastsameindex = -1
		lastsamevalue = None
		for ii in range(len(zippedin)):
			if zippedin[ii][0] != lastsamevalue:
				if lastsameindex != -1:
					targetrank = sum(ranks[lastsameindex:ii])/float(ii-lastsameindex)
					for jj in range(lastsameindex, ii):
						ranks[jj] = targetrank
				lastsameindex = ii
				lastsamevalue = zippedin[ii][0]
		if lastsameindex != -1:
			ii = len(zippedin)
			targetrank = sum(ranks[lastsameindex:ii])/float(ii-lastsameindex)
			for jj in range(lastsameindex, ii):
				ranks[jj] = targetrank
	else:
		import random
		random.shuffle(zippedin)
		zippedin.sort(key=lambda o:o[0])
		ranks = list(range(len(zippedin)))
	zippedout = list(zip([z[1] for z in zippedin], ranks))
	zippedout.sort()
	return [z[1] for z in zippedout]

def globalFDR(pvalues):
	""" returns list of global FDR values (corrected p values) in same order as input list
	    Uses Benjamini-Hochberg method """
	pnum = len(pvalues)
	pi = 1
	zippedin = list(zip(pvalues, list(range(pnum))))
	zippedin.sort()
	fdr = [z[0] for z in zippedin]
	for rank in range(pnum-1,0,-1):
		index = rank-1
		fdr[index] = min(fdr[index+1], pi*zippedin[index][0]*pnum/rank)
	zippedout = list(zip([z[1] for z in zippedin], fdr))
	zippedout.sort()
	return [z[1] for z in zippedout]
	
def esttrue(pvalues, minp=0.5, maxp=0.9):
	""" crudely estimate true number of non-null-hypothesis-followers from p-value distribution
	    will underestimate, unsure if it works """
	undermin = sum([1.0 for v in pvalues if v <= minp])
	undermax = sum([1.0 for v in pvalues if v <= maxp])
	return 1-(undermax-undermin)/undermax/(maxp-minp)*maxp
	
def PtoZ(p):
	from scipy.special import erfinv
	from math import sqrt
	return sqrt(2.0)*erfinv(2.0*p-1.0)
	
def ZtoP(Z):
	from scipy.special import erf
	from math import sqrt
	return 0.5*(1.0+erf(Z/sqrt(2.0)))

def combinedP(pvalues, weights=None):
	""" takes arrays of p-values (preferably 1-sided) and weights (sample sizes if equal variance), returns p-value (1-sided if 1-sided input) """
	from math import sqrt
	Zs = [PtoZ(p) for p in pvalues]
	if weights is None:
		Zcombined = sum([Z for Z in Zs])/sqrt(len(Zs))
	else:
		Zcombined = sum([w*Z for w,Z in zip(weights, Zs)])/sqrt(sum([w**2 for w in weights]))	
	return ZtoP(Zcombined)

def Ztest(values, popstd, popmean=0):
	""" two-tailed test if values are from a normal distributed with standard deviation popstd and mean popmean """
	import numpy, math
	Z = math.sqrt(len(values))*(numpy.mean(values)-popmean)/popstd
	return 2.0*ZtoP(-abs(Z))

def ftest(*args):
	""" input: lists of counts, output: p-value """
	import ctools
	if len(args) == 4 and not hasattr(args[0], '__iter__'):
		return ctools.ftest(args[0], args[1], args[2], args[3])
	else:
		array = list(args[0])
		dim = len(array)
		for a in args[1:]:
			if len(a) != dim: raise ValueError
			array += list(a)
		return ctools.ftest3(array, dim)

def bootstrap(func, datatuple, controls=1000, confidence=0.95, nullval=0, processes=1, give_median=False):
	if processes == 1:
		arr = _bootstrap_loop(func, datatuple, controls)
	else:
		import multiprocessing
		arr=[]
		jobs=[]
		cLeft = controls
		q = multiprocessing.Queue()
		for pi in range(processes, 0, -1):
			jobs.append(multiprocessing.Process(target=_bootstrap_loop, args=(func, datatuple, cLeft//pi, q)))
			cLeft -= cLeft//pi
		for job in jobs: job.start()	
		for job in jobs: arr.extend(q.get())
	if hasattr(arr[0], '__iter__'):
		parr = []
		minvarr = []
		maxvarr = []
		medvarr = []
		for a in map(list,list(zip(*arr))):
			p, minv, maxv, medv = _bootstrap_points(a, confidence, controls, nullval)
			minvarr.append(minv)
			maxvarr.append(maxv)
			medvarr.append(medv)
			parr.append(p)
		if give_median:
			return tuple(parr), tuple(minvarr), tuple(maxvarr), tuple(medvarr), tuple(func(*datatuple))
		else:
			return tuple(parr), tuple(minvarr), tuple(maxvarr)
	else:
		p, minv, maxv, medv = _bootstrap_points(arr, confidence, controls, nullval)
		if give_median:
			return p, minv, maxv, medv, func(*datatuple)
		else:
			return p, minv, maxv

def _bootstrap_points(arr, confidence, controls, nullval):
	if hasattr(confidence, '__iter__'):
		minv, maxv = list(zip(*[pointinarray(arr, [(1-c)/2, (1+c)/2]) for c in confidence]))
		medv = pointinarray(arr, 0.5)
	else:
		minv, maxv, medv = pointinarray(arr, [(1-confidence)/2, (1+confidence)/2, 0.5])
	p = sum(v <= nullval for v in arr)/float(controls)
	return p, minv, maxv, medv

def _bootstrap_loop(func, datatuple, controls, q=None):
	import numpy
	arr = []
	for ci in range(controls):
		resamplings = tuple([data[i] for i in numpy.random.randint(0,len(data),len(data))] for data in datatuple)
		arr.append(func(*resamplings))
	if q is None:
		return arr
	else:
		q.put(arr)

def p_from_distr(vals, distr):
	# same as return [sum(d < v for d in distr)/float(len(distr)) for v in vals], but faster
	envals = sorted((v,i) for i,v in enumerate(vals))
	p_arr = []
	distr = sorted(distr)
	ldistr = float(len(distr))
	r = 0
	p_arr = [-1 for v in vals]
	for v,i in envals:
		try:
			while distr[r] < v:
				r += 1
		except IndexError: pass
		p_arr[i] = 1 - r/ldistr
		if p_arr[i] == 0:
			p_arr[i] = 1/ldistr
	return p_arr

def variance_shrinkage_t_test(values1, values2, params=(0.5, 0.9), permutN=100):
	"""
	values1, values2 = NxM matrices, where N=number of replicates, M=number of tests
	params = (fraction from shrinkage, quantile to use for shrinkage) (default 0.5, 0.9)
	permutN = as high as possible (but makes it take more time)
	"""
	import numpy, math
	means1 = numpy.mean(values1, axis=1)
	means2 = numpy.mean(values2, axis=1)
	def new_var(values):
		variance = numpy.var(values, axis=1)
		c = pointinarray(variance, params[1])
		return variance * (1-params[0]) + c * params[0], c
	def new_var_perm(values, c):
		variance = numpy.var(values, axis=1)
		return variance * (1-params[0]) + c * params[0]
	def t_like_values(val1, val2, combined_variance):
		# ignore the factor 2/sqrt(2/n), since it's constant
		return (numpy.mean(val1, axis=1)-numpy.mean(val2, axis=1))/combined_variance
	var1, c1 = new_var(values1)
	var2, c2 = new_var(values2)
	t_like_real = t_like_values(values1, values2, var1+var2)
	stackval = numpy.hstack((values1, values2))
	n1 = len(values1[0])
	nboth = len(stackval[0])
	def permut_t():
		i_p = numpy.random.permutation(nboth)
		val1 = numpy.hstack(tuple(stackval[:,i:i+1] for i in i_p[:n1]))
		val2 = numpy.hstack(tuple(stackval[:,i:i+1] for i in i_p[n1:]))
		return t_like_values(val1, val2, new_var_perm(val1, c1) + new_var_perm(val2, c2))
	t_like_distr = [t for i in range(permutN) for t in permut_t()]
	return [1-2*abs(0.5-p) for p in p_from_distr(t_like_real, t_like_distr)]

def strip_end_zeros(number):
	string = str(number)
	while string.endswith('0'): string = string[:-1]
	if string.endswith('.'): string = string[:-1]
	return string

class Parsed_rpkms(dict):
	def __init__(self, infiles, counts):
		# contains the table of values and names as dictionary with 'symbols', 'IDs' or sample name a key
		self.samples = []
		self.filenames = infiles
		self.allmappedreads = []
		self.normalizationreads = []
		self.is_counts = counts
		self.symbol_to_index = dict()
		self.ID_to_index = dict()

def writeexpr(filename, rpkm_expr, counts_expr=None, samples=None, row_indices=None, extra_comment_lines=[]):
	import sys, time
	if samples is None: samples = rpkm_expr.samples
	if row_indices is None: row_indices = range(len(rpkm_expr['symbols']))
	with open(filename, 'w') as outfh:
		print(join('#samples', samples), file=outfh)
		totalreadsD = dict(zip(rpkm_expr.samples, rpkm_expr.allmappedreads))
		normreadsD = dict(zip(rpkm_expr.samples, rpkm_expr.normalizationreads))
		print(join('#allmappedreads', [totalreadsD.get(s, 0) for s in samples]), file=outfh)
		print(join('#normalizationreads', [normreadsD.get(s, 0) for s in samples]), file=outfh)
		print(join('#arguments', ' '.join(sys.argv), 'time: '+time.asctime()), file=outfh)
		for line in extra_comment_lines:
			if not line[0] == '#': line = '#' + line
			line = line.rstrip('\r\n')
			print(line, file=outfh)
		for i in row_indices:
			symbol = rpkm_expr['symbols'][i]
			ID = rpkm_expr['IDs'][i]
			values_rpkm = (rpkm_expr[s][i] for s in samples)
			if counts_expr is None:
				print(join(symbol, ID, values_rpkm), file=outfh)
			else:
				values_reads = (counts_expr[s][i] for s in samples)
				print(join(symbol, ID, values_rpkm, map(strip_end_zeros, values_reads)), file=outfh)

def loadexpr(infiles, counts=False):
	"""
	loads from output of rpkmforgenes.py
	"""
	if isinstance(infiles, str): infiles = [infiles]
	values = Parsed_rpkms(infiles, counts)
	numsymbols = None
	for infile in infiles:
		with open(infile, 'r') as infh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				if p[0] == '#samples':
					samples = p[1:]
					values.update(dict((n,[]) for n in samples))
					values['symbols'] = []
					values['IDs'] = []
					indexstart = 2+len(samples) if counts else 2
					values.samples.extend(samples)
				elif p[0] == '#allmappedreads':
					values.allmappedreads.extend([float(v) for v in p[1:]])
				elif p[0] == '#normalizationreads':
					values.normalizationreads.extend([float(v) for v in p[1:]])
				elif line.startswith('#'):
					continue
				else:
					for s,v in zip(samples, p[indexstart:]):
						try:
							values[s].append(None if v == '-1' else float(v))
						except:
							import sys
							#print('Problem line:', repr(line), file=sys.stderr)
							raise
					values['symbols'].append(p[0])
					values['IDs'].append(p[1])
		if not (numsymbols is None or numsymbols == len(values['symbols'])):
			raise Exception('Mismatch in number of gene symbols between files')
		numsymbols = len(values['symbols'])
	
	# prepare some dictionaries
	values.symbol_to_index = dict((s, i) for i, S in enumerate(values['symbols']) for s in S.split('+'))
	values.symbol_to_index.update(dict((S, i) for i, S in enumerate(values['symbols'])))
	values.ID_to_index = dict((s, i) for i, S in enumerate(values['IDs']) for s in S.split('+'))
	values.ID_to_index.update(dict((S, i) for i, S in enumerate(values['IDs'])))
	
	return values

def getsequence(chromosome, start, end, genomedir, filesuffix=".fa"):
	""" returns nucleotide sequence string, for 0-based inclusive to exclusive interval """
	import os
	chromosomefile = os.path.join(genomedir, chromosome + filesuffix)
	cfileh = open(chromosomefile, "r")
	if start < 0: start = 0
	
	global chromosomefile_infodict
	try: chromosomefile_infodict
	except: chromosomefile_infodict = {}
	try: offset, linelength, seqlength = chromosomefile_infodict[chromosomefile]
	except:
		line1 = cfileh.readline(1000)
		if len(line1) < 1000 and line1[0] == '>': offset = len(line1)
		else:
			cfileh.seek(0)
			offset = 0
		line2 = cfileh.readline(1000)
		if len(line2) < 1000:
			linelength = len(line2)
			seqlength = len(line2.rstrip())
		else:
			linelength = 0
			seqlength = 0
		chromosomefile_infodict[chromosomefile] = offset, linelength, seqlength
	if linelength == 0:
		startfilepos = start + offset
		endfilepos = end + offset
	else:
		startfilepos = offset + (start // seqlength)*linelength + (start % seqlength)
		endfilepos = offset + (end // seqlength)*linelength + (end % seqlength)	
	cfileh.seek(startfilepos, 0)
	seq = ''.join(cfileh.read(endfilepos-startfilepos).split())
	cfileh.close()
	return seq

def reverseDNA(seq_in):
	""" returns nucleotide sequence string """
	sequencetools_reverseDNAdict = {"A":"T", "C":"G", "G":"C", "T":"A", "R":"Y","Y":"R","K":"M","M":"K","S":"S","W":"W","B":"V","D":"H","H":"D","V":"B","N":"N", "a":"t","c":"g","g":"c","t":"a","n":"n","\n":"\n"}
	seq_out = ""
	for bi in range(len(seq_in)):
		seq_out = sequencetools_reverseDNAdict[seq_in[bi]] + seq_out
	return seq_out
	
def expandsequence(sequence):
	return expandsequences([sequence])
	
def expandsequences(sequences):
	""" returns array of sequences """
	expanddict = {"A":"A", "C":"C", "G":"G", "T":"T", "R":"GA", "Y":"TC", "K":"GT", "M":"AC", "S":"GC", "W":"AT", "B":"GTC", "D":"GAT", "H":"ACT", "V":"GCA", "N":"ACGT"}
	seq_out = []
	for seqin in sequences:
		seq_l = [""]
		for bi in range(len(seqin)):
			expanded = expanddict[seqin[bi]]
			oldsequences = seq_l
			seq_l = []
			for ei in range(len(expanded)):
				for seq in oldsequences:
					seq_l.append(seq + expanded[ei])
		seq_out += seq_l
	return seq_out
	
def tocolour(seq):
	""" returns colourspace sequence string """
	plainseq = seq.upper()
	colourdict = {"AA":"0", "CC":"0", "GG":"0", "TT":"0", "CA":"1", "AC":"1", "GT":"1", "TG":"1", "GA":"2", "AG":"2", "TC":"2", "CT":"2", "TA":"3", "AT":"3", "CG":"3", "GC":"3"}
	colourseq = ""
	for pos in range(len(plainseq)-1):
		try:
			colour = colourdict[plainseq[pos:pos+2]]
		except:
			colour = "."
		colourseq += colour
	return plainseq[0] + colourseq
	
def fromcolour(seq):
	colour0 = {"A":"A", "C":"C", "G":"G", "T":"T"}
	colour1 = {"A":"C", "C":"A", "G":"T", "T":"G"}
	colour2 = {"G":"A", "A":"G", "T":"C", "C":"T"}
	colour3 = {"T":"A", "A":"T", "C":"G", "G":"C"}
	
	cs_sequence = seq
	pos = 1
	prevbase = cs_sequence[0]
	seq = prevbase
	while pos < len(cs_sequence):
		colour = cs_sequence[pos]
		if colour == "0":
			nextbase = colour0[prevbase]
		elif colour == "1":
			nextbase = colour1[prevbase]
		elif colour == "2":
			nextbase = colour2[prevbase]
		elif colour == "3":
			nextbase = colour3[prevbase]
		else:
			seq += "N"
			return seq
		seq += nextbase
		prevbase = nextbase
		pos += 1
	return seq
	
def loadlist(infile, index=None, func=None, ignore='#', ignorelines=0):
	""" returns array of strings """
	infileh = open(infile, "r")
	for i in range(ignorelines):
		infileh.readline()
	outarray = [l.rstrip() for l in infileh.readlines()]
	infileh.close()
	if ignore is not None:
		outarray = [l for l in outarray if not l.startswith(ignore)]
	if index is not None:
		outarray = [l.split("\t")[index] for l in outarray]	
	if func is not None:
		outarray = list(map(func, outarray))
	return outarray
	
def printlist(outfile, inlist, method="w"):
	outfileh = open(outfile, method)
	for string in inlist:
		print(string, file=outfileh)
	outfileh.close()
	
def histogramheights(array, start, end, step, cumulative=0, fractions=False):
	""" returns 2 arrays: x and y """
	binpositions = []
	binheights = []
	pos = start
	arraylen = 0
	while pos <= end:
		binpositions.append(pos)
		binheights.append(0)
		pos += step
	for element in array:
		bin = int((element-start)/step)
		arraylen += 1
		if cumulative > 0:
			for bi in range(bin, len(binheights)):
				try:
					assert bi >= 0
					binheights[bi] += 1
				except: pass
		elif cumulative < 0:
			for bi in range(0, bin+1):
				try:
					assert bi >= 0
					binheights[bi] += 1
				except: pass
		else:
			try:
				assert bin >= 0
				binheights[bin] += 1
			except: pass
	if fractions:
		return (binpositions, [h/arraylen for h in binheights])
	else:
		return (binpositions, binheights)
	
def bin(array, start, end, step, cumulative=0, fractions=False):
	return histogramheights(array, start, end, step, cumulative, fractions)
	
def mixcolours(colours, weights):
	cout = "#"
	for si in [1,3,5]:	# red, green, blue
		cvals = [int(c[si:si+2], 16) for c in colours]
		mixed = max(0,min(255,int(0.5+sum([cvals[i]*weights[i] for i in range(len(cvals))]))))
		outstr = "%X" % mixed
		if len(outstr) == 1: outstr = "0" + outstr
		cout += outstr
	return cout

def rainbowmix(fraction, stops=['#fe0000', '#00fe00', '#0000fe', '#000000', 'f0f0f0', '#fe0000', '#0000fe', '#f0f0f0', '#00fe00', '#000000']):
	n = len(stops)
	mf = 1-(fraction*n)%1
	lstop = int(fraction*n)
	return mixcolours([stops[lstop%n], stops[(lstop+1)%n]], [mf, 1-mf])

def randomcolour():
	import random
	def hexconv(n):
		r = hex(n)[2:]
		return r if len(r) == 2 else '0'+r
	return '#' + ''.join([hexconv(random.randint(0,255)) for c in 'rgb'])

def chisquare(observed, total):
	""" returns fold enrichment, p-value """
	from scipy import stats
	import numpy
	expected = [float(v)*sum(observed)/sum(total) for v in total]
	obs = numpy.asarray(observed)
	obs = obs.astype(float)
	exp = numpy.asarray(expected)
	exp = exp.astype(float)
	fc = [observed[ii]/expected[ii] if expected[ii] != 0 else 0.0 for ii in range(len(observed))]
	return fc, stats.chisquare(obs, exp)[1]

def chisquaretable(a, b, c, d):
	fc, p = chisquare([a,b], [c+a,d+b])
	return fc[0], p

def permutationtest(func, a, b, controls=1000):
	""" returns 2-sided p-value for func(a) = func(b) """
	import random
	arr = [0 for i in range(controls+1)]
	arr[0] = abs(func(a) - func(b))
	c = a + b
	lena = len(a)
	for i in range(1, controls+1):
		random.shuffle(c)
		arr[i] = abs(func(c[:lena]) - func(c[lena:]))
	r = rank(arr)
	return 1-float(r[0])/(controls+1)

def pointinarray(array, quantile):
	#sortedarray = array[:]
	#sortedarray.sort()
	import numpy
	sortedarray = numpy.array(sorted(array))
	def _point(q):
		indexlow = int(q*(len(sortedarray)-1))
		fl = q*(len(sortedarray)-1)-indexlow
		if indexlow == len(sortedarray) - 1:
			return float(array[indexlow])
		else:
			return float(sortedarray[indexlow]*(1-fl)+sortedarray[indexlow+1]*fl)
	if isinstance(quantile, list):
		return list(map(_point, quantile))
	else:
		return _point(quantile)

def permutationtest_confint(func, a, b, confidence=0.95, controls=1000):
	""" returns func(a)-func(b), min, max, assumes equal distribution except func (eg mean) """
	import random
	pvalue = confidence
	if pvalue > 1: pvalue /= 100.0
	if pvalue > 0.5: pvalue = 1-pvalue
	arr = [0 for i in range(controls)]
	arr_0 = func(a) - func(b)
	c = a + b
	lena = len(a)
	for i in range(controls):
		random.shuffle(c)
		arr[i] = func(c[:lena]) - func(c[lena:])
	return arr_0+pointinarray(arr,pvalue/2.0), arr_0, arr_0+pointinarray(arr,1-pvalue/2.0)

def adjWald(a, b, confidence=0.95):
	""" calculate confidence interval for a/b, where a and b are integers """
	# http://measuringux.com/AdjustedWald.htm
	n = float(b)
	p = a/n
	z2 = PtoZ((1+confidence)/2)**2
	padj = (n*p + z2/2)/(n+z2)
	nadj = n + z2
	d = (z2*padj*(1-padj)/nadj)**0.5
	return padj-d, padj+d

def bootstrap_confint(func, a, b=None, confidence=0.95, resamplings=1000):
	""" returns min, max; dose not assume similar distribution, requires >20 values in a and b """
	import random
	pvalue = confidence
	if pvalue > 1: pvalue /= 100.0
	if pvalue > 0.5: pvalue = 1-pvalue
	arr = []
	if b is None:
		for i in range(resamplings):
			rs_a = [random.choice(a) for j in range(len(a))]
			arr.append(func(rs_a))
	else:
		for i in range(resamplings):
			rs_a = [random.choice(a) for j in range(len(a))]
			rs_b = [random.choice(b) for j in range(len(b))]
			arr.append(func(rs_a)-func(rs_b))
	return pointinarray(arr,pvalue/2.0), pointinarray(arr,1-pvalue/2.0)

def violin_plot(ax,data,pos, bp=False):
    '''
    create violin plots on an axis
    run with e.g violin_plot(pylab.axes(), [[3,4,5],[7,8]], [0, 1])
    '''
    # from http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
    from matplotlib.patches import Rectangle
    from scipy.stats import gaussian_kde
    from numpy.random import normal
    from numpy import arange
    
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/500.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
    if bp:
        #ax.boxplot(data,notch=1,positions=pos,vert=1)
        boxplotborders = [pointinarray(D, [0.25,0.5,0.75]) for D in data]
        for x, borders in zip(pos, boxplotborders):
            ax.add_patch(Rectangle((x-0.05, borders[0]), 0.1, borders[2]-borders[0], linewidth=0, facecolor='k'))
        ax.plot(pos, [d[1] for d in boxplotborders], 'wo')

def GOgenelist(GO_file, term, shortened=0):
	files = {"BP":"/home/danielr/ChIP-seq-Sox3/perGOcat/data/BP_goterm.txt","MF":"/home/danielr/ChIP-seq-Sox3/perGOcat/data/MF_goterm.txt", "CC":"/home/danielr/ChIP-seq-Sox3/perGOcat/data/CC_goterm.txt"}
	try: GO_file = files[GO_file]
	except: pass
	foundcat  = 0
	GO_fileh = open(GO_file, "r")
	for line in GO_fileh:
		p = line[:-1].split("\t")
		if p[0] == term or (shortened and term.upper().replace(" ","_") in p[0].upper().replace(" ","_")):
			GO_fileh.close()
			return p[1].split(";")
	GO_fileh.close()
	raise UserWarning("Did not find term " + term + " in " + GO_file)

def loadmotif(infile, trimstart=0, trimend=0):
	from TAMO import MotifTools
	lines = loadlist(infile)
	if lines[0] == "A\tC\tG\tT":
		ma = []
		for l in lines[1:]:
			p = l.split("\t")
			ma.append({'A':float(p[0]), 'C':float(p[1]), 'G':float(p[2]), 'T':float(p[3])})
		if trimend == 0: ma = ma[trimstart:]
		else: ma = ma[trimstart:-trimend]
		return MotifTools.Motif_from_counts(ma)
	elif lines[0][0] in 'ACGT':
		if trimend == 0: lines = lines[trimstart:]
		else: lines = lines[trimstart:-trimend]
		return MotifTools.Motif(lines)
	else:
		na = []
		for line in lines:
			na.append(list(map(int, line.split())))
		ma = []
		for i in range(len(na[0])):
			ma.append({'A':na[0][i], 'C':na[1][i], 'G':na[2][i], 'T':na[3][i]})
		return MotifTools.Motif_from_counts(ma)

def join(*args, **kwargs):
	""" returns tab-separated string """
	try: sep = kwargs["sep"]
	except: sep = "\t"
	array = []
	for a in args:
		if hasattr(a, '__iter__') and not (isinstance(a, str) or isinstance(a, unicode)): array.extend(a)
		else: array.append(a)
	return sep.join(map(str, array))

def split(line, sep='\t'):
	return line.rstrip('\r\n').split(sep)

def splitlines(infile, ignore='', sep='\t'):
	if infile.endswith('.gz'):
		import gzip
		infh = gzip.open(infile, 'r')
	else:
		infh = open(infile, 'rU')
	try:
		for line in infh:
			if ignore and line.startswith(ignore): continue
			yield line.rstrip('\r\n').split(sep)
	finally:
		infh.close()

class Cregion:
	allchromosomes = {}
	indexdict = {} # inverse of allchromosomes
	allwindows = []
	WL = 3000
	
	def __init__(self, chromosome, start, end=None, strand='?'):
		self.start = start
		if end == None: self.end = start
		else: self.end = end
		try: self.chrindex = Cregion.allchromosomes[chromosome+strand]
		except KeyError:
			self.chrindex = len(Cregion.allchromosomes)
			Cregion.allchromosomes[chromosome+strand] = self.chrindex
			Cregion.indexdict[self.chrindex] = chromosome+strand
			Cregion.allwindows.append([])
	
	def addtowindows(self):
		# adds instance to Cregion.allwindows
		wchr = Cregion.allwindows[self.chrindex]
		if len(wchr) <= self.end//Cregion.WL: wchr.extend([[] for i in range(1+self.end//Cregion.WL-len(wchr))])
		for wi in range(self.start//Cregion.WL, self.end//Cregion.WL+1):
			wchr[wi].append(self)
		
	def getwindow(self):
		# returns list of Cregion instances which could overlap
		wchr = Cregion.allwindows[self.chrindex]
		s = min(len(wchr), self.start//Cregion.WL)
		e = min(len(wchr), self.end//Cregion.WL+1)
		return list(set([v for l in wchr[s:e] for v in l])) # flattens wchr[s:e], removes duplicates
		
	def overlaps(self, other):
		return self.start <= other.start < self.end or other.start <= self.start < other.end
	
	def overlapping(self):
		return [r for r in self.getwindow() if r.overlaps(self)]
	
	def getchromosome(self):
		return Cregion.indexdict[self.chrindex][:-1]
	
	def startingwithin(self):
		# returns list of Cregion instances whose start coordinate is within the region
		return [r for r in self.getwindow() if self.start <= r.start < self.end]
	
	def getstrand(self):
		strand = Cregion.indexdict[self.chrindex][-1]
		if strand == "?": raise Exception("No strand given")
		return strand
	
	def __repr__(self):
		return self.name(1, 0)
	
	def name(self, start_add=0, end_add=0):
		try:
			strand = self.getstrand()
		except:
			return self.getchromosome()+":"+str(self.start+start_add)+"-"+str(self.end+end_add)
		else:
			return self.getchromosome()+":"+str(self.start+start_add)+"-"+str(self.end+end_add)+":"+strand

	@staticmethod
	def clearwindows(new_windowsize=None):
		if new_windowsize is not None: Cregion.WL = new_windowsize
		Cregion.allwindows = [[] for c in Cregion.allchromosomes]
		
	@staticmethod
	def overlappingpoint(chromosome, pos, strand='?'):
		try:
			wchr = Cregion.allwindows[Cregion.allchromosomes[chromosome+strand]]
		except KeyError:
			return []
		s = pos//Cregion.WL
		try:
			return [r for r in wchr[s] if r.start <= pos < r.end]
		except IndexError:
			return []
	
	@staticmethod
	def closesttopoint(chromosome, pos, strand='?', mindist=0, maxdist=1e30, check_forward=True, check_backward=True):
		try:
			wchr = Cregion.allwindows[Cregion.allchromosomes[chromosome+strand]]
		except KeyError:
			return []
		s = pos//Cregion.WL
		windist = 0
		candidates = set()
		
		def closesttopoint_distance(r):
			if r.start <= pos <= r.end: return 0
			return min(abs(r.start-pos), abs(r.end-pos))
		
		while windist < 5+maxdist/Cregion.WL:
			new_candidates = set()
			if check_backward:
				try: new_candidates |= set(r for r in wchr[s-windist])
				except IndexError:pass
			if check_forward:
				try: new_candidates |= set(r for r in wchr[s+windist])
				except IndexError:pass
			new_candidates = set(r for r in new_candidates if mindist <= closesttopoint_distance(r) < maxdist)
			if candidates:
				candidates |= new_candidates
				break
			candidates |= new_candidates
			windist += 1
		
		if not candidates: return []
		closest_dist = min(closesttopoint_distance(r) for r in candidates)
		return [r for r in candidates if closesttopoint_distance(r) == closest_dist]

def flag(flags, default=None, array=None):
	if array == None:
		import sys
		argv = sys.argv
	else: argv = array
	if "--help" in argv: raise UserWarning("--help in arguments")
	for f in str(flags).split("/"):
		if f[0] != "-":
			try: return flagpos(int(f), None, argv)
			except: pass
		try:
			index = argv.index(f)
		except ValueError: continue
		try:
			return argv[index+1]
		except IndexError: raise UserWarning("No value after " + f + " in arguments")
	if default == None: raise UserWarning(str(flags) + " not in arguments")
	return default
	
def flagarray(flags, default=None, array=None):
	if array == None:
		import sys
		argv = sys.argv
	else: argv = array
	if "--help" in argv: raise UserWarning("--help in arguments")
	index = None
	for f in str(flags).split("/"):
		if f[0] != "-":
			v = []
			i = 0
			try:
				while 1:
					v.append(flagpos(int(f)+i, None, argv))
					i += 1
			except: pass
			if len(v) > 0: return v
		try:
			index = argv.index(f)
		except ValueError: continue
	if index == None:
		if default == None: raise UserWarning(str(flags) + " not in arguments")
		return default
	nextflag = index+1
	while nextflag < len(argv):
		if argv[nextflag][0] == "-" and len(argv[nextflag]) > 1 and argv[nextflag][1] not in "0123456789.": break
		nextflag += 1
	return argv[index+1:nextflag]
	
def ifflag(flags, present=1, absent=0, array=None):
	if array == None:
		import sys
		argv = sys.argv
	else: argv = array
	if "--help" in argv: raise UserWarning("--help in arguments")
	for f in str(flags).split("/"):
		if f in argv: return present
	return absent
	
def flagpos(index, default=None, array=None):
	if array == None:
		import sys
		argv = sys.argv
	else: argv = array
	if "--help" in argv: raise UserWarning("--help in arguments")
	try:
		value = argv[index]
		for v in argv[1:index+1]:
			if len(v) > 1 and v[0] == "-" and v[1] not in "0123456789.": raise IndexError
	except:
		if default == None:
			raise UserWarning("Argument " + str(index) + " is missing")
		else: return default
	return value
