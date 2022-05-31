#!/usr/bin/python
"""
Calculates gene expression from a read mapping file
"""

lastmodified = "26 May 2022"
author = "Daniel Ramskold"
#4 june 2010: now assumes sam and gff files use 1-based coordinates, not 0-based, -exonnorm is default
#6 june 2010: partially rewrote code for overlapping exons from different isoforms, swapped ID and symbol fields for some annotation types
#10 june 2010: corrected gff end coordinates
#17 june 2010: removed exons of other isoforms from -intronsinstead, normalise by only mRNA by default
#28 june 2010: bed12 format as annotation, -flat option for collapsing before RPKM calculation, put strand in unused gene name column if -strand
#3 july 2010: changed output header to give number of reads normalised by regardless of normalisation option
#8 july 2010: -bedann for 3 columns gives 1-based chrX:start:end but no change for -bed3ann, -bothends map end-1 position instead of end, added -n option, more stringent check (strand, length) to filter with -diffreads
#28 july 2010: added -unique option (removal of multimapping reads)
#6 aug 2010: -u option can call bigWigSummary on bigWig files
#11 aug 2010: rewrote parsing of sam files to make -bothends and -diffreads work better, removed the need for -unique option
#19 aug 2010: changed gene/ID fields with -bedann
#23 aug 2010: changed ouput of 0 rpkm to 0.0 rpkm
#10 sep 2010: indexing exons instead of reads
#17 nov 2010: reenabled old sam format input function (now with -samse), added bam file support
#9 dec 2010: allowed -u to read from files with FASTA header and fixed line wrap
#5 jan 2011: looks for .bigWig file suffix, fixed reporting for -readcount -nocollapse
#20 jan 2011: changed GTF format ID to transcript ID
#24 jan 2011: added -randomreads
#9 feb 2011: speed optimizations for -randomreads, fixed error reporting at 0 total reads
#7 mars 2011: debugging of isoform deconvolution
#12 apr 2011: added -namecollapse option
#3 may 2011: added -dynuniq
#8 aug 2011: let -dynuniq/-ulen be followed by a read length value
#31 aug 2011: ignores lines starting with # in annotation file
#26 sep 2011: changed output sorting, -sortpos for old sorting
#27 sep 2011: added -rmnameoverlap
#30 sep 2011: improved read length auto-detection for -ulen, added -table
#4 oct 2011: altered -rmnameoverlap
#19 oct 2011: added -quite, -p, made auto-detect for -ulen handle multiple read lengths
#20 oct 2011: -ulen not needed to run with Helena's files
#3 nov 2011: added detection of file suffix for -ulen
#4 nov 2011: fixed bug with -p combined wth -forcedtotal
#29 nov 2011: -exportann
#24 apr 2012: debugging of isoform deconvolution
#30 may 2012: added -bamse flag, also works for sam files
#8 jun 2012: renamed -bamse to -bamu, fixed a crash with unknown query sequence in bam files
#25 jun 2012: paired-end support with -bamu, -samse, made -bamu default for bam/sam, made -fulltranscript default instead of -no3utr, added last modified date to #arguments line
#27 jul 2012: added -minqual
#17 aug 2012: changed -readcount to not counting reads for completely non-unique exons
#20 aug 2012: -forcetotal changes
#22 aug 2012: added -maxNM, -exportann now allowed for when using multiple input files
#11-13 sep 2012: .unique20-255.btxt suffix added for -u
#19 sep 2012: catch runtime warnings from matrix determinant
#21 sep 2012: more fixes for rare much-too-large rpkm values
#19 oct 2012: made -strand work with -u
#26 oct 2012: added -bothendsceil
#5 Nov 2012: some more error reporting, made * a valid CIGAR string for SAM files
#6 Nov 2012: considers reads with empty CIGAR in sam and bam files unmapped
#15 Jan 2013: debugged -bothends for -samu/-bamu, added -midread
#6 Feb 2013: bug fix for paired-end bam/sam that caused it to ignore all reads
#8 Feb 2013: added -ensgtfann and -norandom
#12 Feb 2013: swapped the strand of the 2nd read for paired-end sequencing
#11 Apr 2013: added -readpresent (from Ramu Chenna)
#3 May 2013: added -rmregions
#7 May 2013: bug fix for -rmregions used with -intronsinstead/-rmnameoverlap
#5 Aug 2013: fixed a bug for -bamu and -mapends where it double-counted paired-end reads
#7 Feb 2014: Added early detection for -p and -table clash
#13 Nov 2014: added -addchr
#10 Dec 2014: backup code for old version of numpy lacking numpy.linalg.slogdet
#8 jan 2015: added *_alt to removed chromosomes unless -keephap is used
#23 Jan 2015: made -maxNM also check nM flag (used by rna-star)
#13 Mar 2015: added hidden flag -blocksize
#7 Jun 2018: added flag -skipiffewreads
#9 Aug 2018: added chrUn_ to chromosome name partial matches that are removed from annotaion by -norandom (fix for hg38), modified how -rnnameoverlap interact with gene model collapse (can be returned to old behaviour with -limitcollapse2), added -namesum
#16 Aug 2018: warned for -midread
#17 aug 2018: removed a print statement
#14 Sep 2018: added -tpm
#19 sep 2018: readded 16-17 aug changes that had disappeared by accidental forking
#20 nov 2018: fixed a bug where -namesum would push RPKM values into the reads section
#4 apr 2022: added -tmm
#26 May 2022: fixed crash when using -tmm
#31 May 2022: added -csv, fixed zero-stripping TMM count export bug, added .gz/.bz2 support here and there

from numpy import matrix, linalg
import numpy, sys, time, os, subprocess, math
from collections import defaultdict

uniqueposdir = "none"

if False: # check if print statement works, SyntaxError if it doesn't
	print 'This program does not work under python 3, run it in python 2.5/2.6/2.7'

def join(*args, **kwargs):
	""" returns tab-separated string """
	try: sep = kwargs["sep"]
	except: sep = "\t"
	array = []
	for a in args:
		iterable = hasattr(a, '__iter__')
		if iterable and isinstance(a, str): iterable = False
		try:
			if iterable and isinstance(a, unicode): iterable = False
		except NameError: pass
		if iterable: array.extend(a)
		else: array.append(a)
	return sep.join(map(str, array))

def strip_end_zeros(number):
	string = str(number)
	if '.' in string or ',' in string:
		while string.endswith('0'): string = string[:-1]
		if string.endswith('.') or string.endswith(','): string = string[:-1]
	return string

def quote_if_comma_within(text):
	if ',' in text or '\n' in text or '"' in text: text = '"'+text+'"'
	return text

def writeexpr_csv(filename, rpkm_expr, samples=None, row_indices=None):
	import sys, time
	if samples is None: samples = rpkm_expr.samples
	if row_indices is None: row_indices = range(len(rpkm_expr['symbols']))
	with openfile(filename, 'w') as outfh:
		_print(join('symbol',[quote_if_comma_within(s) for s in samples], sep=','), file=outfh)
		for i in row_indices:
			symbol = rpkm_expr['symbols'][i]
			values_rpkm = (quote_if_comma_within(strip_end_zeros(rpkm_expr[s][i])) for s in samples)
			_print(join(quote_if_comma_within(symbol), values_rpkm, sep=','), file=outfh)

def getpotentialargument(flag):
	# the argument after the flag
	flagindex = sys.argv.index(flag)
	if flagindex + 1 == len(sys.argv) or sys.argv[flagindex+1][0] == '-': return None
	return sys.argv[flagindex+1]

def convert_to_csv(outfile, readcountcsvfile):
	expr_rpkms = loadexpr(outfile, False)
	if readcountcsvfile is not None: expr_counts = loadexpr(outfile, True)
	writeexpr_csv(outfile, expr_rpkms)
	if readcountcsvfile is not None: writeexpr_csv(readcountcsvfile, expr_counts)

def run_tmm(infile, outfile, copy_counts=True, run_on_counts=False, ref_samples=None):
	def log2(v):
		return math.log(v, 2)
	def M(gi, Y_k, Y_r, N_k, N_r):
		return log2((float(Y_k[gi])/N_k)/(float(Y_r[gi])/N_r))
	def M_logdivlog(gi, Y_k, Y_r, N_k, N_r):
		return log2(float(Y_k[gi])/N_k)/log2(float(Y_r[gi])/N_r)
	def w(gi, Y_k, Y_r, N_k, N_r):
		return float(N_k - Y_k[gi])/N_k/Y_k[gi] + float(N_r - Y_r[gi])/N_r/Y_r[gi]
	def A(gi, Y_k, Y_r, N_k, N_r):
		return 0.5*log2(float(Y_k[gi])/N_k * Y_r[gi]/N_r) if Y_k[gi] > 0 else -10000
	expr_in = loadexpr(infile, counts=run_on_counts)
	
	ref_samples = expr_in.samples if ref_samples is None else ref_samples
	Y_r = [numpy.mean([expr_in[s][gi] for s in ref_samples]) for gi in range(len(expr_in['symbols']))]
	N_r = sum(Y_r)
	
	expr_out = Parsed_rpkms([], False)
	normalization_factors = []
	failed = False
	for s in expr_in.samples:
		Y_k = expr_in[s]
		N_k = sum(Y_k)
		nonzero = [gi for gi in range(len(expr_in['symbols'])) if Y_k[gi] > 0 and Y_r[gi] > 0]
		A_distr = sorted((A(gi, Y_k, Y_r, N_k, N_r), gi) for gi in nonzero)
		M_distr = sorted((M(gi, Y_k, Y_r, N_k, N_r), gi) for gi in nonzero)
		
		Gstar = set(gi for A_val,gi in A_distr[int(0.05*len(A_distr)):-int(0.05*len(A_distr))]) & set(gi for M_val,gi in M_distr[int(0.3*len(M_distr)):-int(0.3*len(M_distr))])
		
		if len(nonzero) == 0: f_k = 1
		else:
			try: log2TMM = float(sum(w(gi, Y_k, Y_r, N_k, N_r) * M(gi, Y_k, Y_r, N_k, N_r) for gi in Gstar))/sum(w(gi, Y_k, Y_r, N_k, N_r) for gi in Gstar)
			except ZeroDivisionError:
				if vocal: print "TMM normalisation failed!"
				failed = True
				f_k = 1
			else:
				f_k = 2**log2TMM # multipy non-reference by this value
		
		
		expr_out[s] = [Y_k[gi]/f_k for gi in range(len(expr_in['symbols']))]
		normalization_factors.append(f_k)
	expr_out.allmappedreads = expr_in.allmappedreads
	expr_out.normalizationreads = [nr*nf for nr, nf in zip(expr_in.normalizationreads, normalization_factors)]
	expr_out.samples = expr_in.samples
	expr_out['symbols'] = expr_in['symbols']
	expr_out['IDs'] = expr_in['IDs']
	
	writeexpr(outfile, expr_out, counts_expr=(loadexpr(infile, counts=True) if copy_counts else None))
	return failed

def loadexpr(infiles, counts=False):
	"""
	loads from output of rpkmforgenes.py
	"""
	if isinstance(infiles, str): infiles = [infiles]
	values = Parsed_rpkms(infiles, counts)
	numsymbols = None
	for infile in infiles:
		samples = None
		with openfile(infile, 'r') as infh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				if p[0] == '#samples':
					samples = p[1:]
					values.update(dict((n,[]) for n in samples))
					values['symbols'] = []
					values['IDs'] = []
					indexstart = 2+len(samples) if counts else 2
					values.samples.extend(samples)
				elif samples is None and p[0] == 'gene_name' and p[1] == 'geneID':
					samples = p[2:]
					if samples[:len(samples)//2] == samples[len(samples)//2:]:
						samples = samples[:len(samples)//2]
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
						if v == '-1': break
						else: values[s].append(float(v))
					else:
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
	
	def to_dataframe(self, indexname='symbols'):
		import pandas
		return pandas.DataFrame(dict((s,self[s]) for s in self.samples), index=(None if indexname is None else self[indexname]))

def _print(txt, **kwargs):
	print >>kwargs['file'], txt

def writeexpr(filename, rpkm_expr, counts_expr=None, samples=None, row_indices=None, extra_comment_lines=[]):
	import sys, time
	if samples is None: samples = rpkm_expr.samples
	if row_indices is None: row_indices = range(len(rpkm_expr['symbols']))
	with openfile(filename, 'w') as outfh:
		_print(join('#samples', samples), file=outfh)
		totalreadsD = dict(zip(rpkm_expr.samples, rpkm_expr.allmappedreads))
		normreadsD = dict(zip(rpkm_expr.samples, rpkm_expr.normalizationreads))
		if rpkm_expr.allmappedreads != []:
			_print(join('#allmappedreads', [totalreadsD.get(s, 0) for s in samples]), file=outfh)
		if rpkm_expr.normalizationreads != []:
			_print(join('#normalizationreads', [normreadsD.get(s, 0) for s in samples]), file=outfh)
		_print(join('#arguments', ' '.join(sys.argv), 'time: '+time.asctime()), file=outfh)
		for line in extra_comment_lines:
			if not line[0] == '#': line = '#' + line
			line = line.rstrip('\r\n')
			_print(line, file=outfh)
		for i in row_indices:
			symbol = rpkm_expr['symbols'][i]
			ID = rpkm_expr['IDs'][i]
			values_rpkm = (rpkm_expr[s][i] for s in samples)
			if counts_expr is None:
				_print(join(symbol, ID, values_rpkm), file=outfh)
			else:
				values_reads = (counts_expr[s][i] for s in samples)
				_print(join(symbol, ID, values_rpkm, map(strip_end_zeros, values_reads)), file=outfh)

class Cexon:
	def __init__(self, start, end, normalise_with):
		self.start = start
		self.end = end
		self.reads = 0
		self.transcripts = []
		self.readspersample = []
		self.forbidden = False # only for normalisation, hide from expression calculation and output
		self.normalise_with = normalise_with # include in exon normalisation
		self.length = self.end - self.start
	
	def calclength(self, chromosome, readlength_dict, filesuffix):
		if ONLYUNIQUEPOS:
			if USESTRANDINFO:
				chromosome = chromosome[:-1]
			if uniquenessfiletype == 'fasta':
				self.setuniquelength(os.path.join(uniqueposdir, chromosome + ".fa"))
			elif bigWigSummary_path:
				self.uniquelengthfromBigWig(chromosome)
			else:
				self.setuniquelength_H(os.path.join(uniqueposdir, chromosome + filesuffix), readlength_dict)
		else:
			self.length = self.end - self.start

	def exonstring(self):
		outstr = "(" + str(self.start) + "," + str(self.end)
		for tx in self.transcripts:
			outstr += "," + tx.ID
		outstr += ") " + str(self.reads) + " " + str(self.length)
		return outstr
	
	def setuniquelength_H(self, chromosomefile, readlength_dict):
		try:
			cfileh = openfile(chromosomefile, 'rb')
		except:
			global warnedchromosomes
			if not chromosomefile in warnedchromosomes:
				warnedchromosomes.append(chromosomefile)
				if vocal: print "Warning: Did not find", chromosomefile
			self.length = self.end - self.start
			return
		self.length = 0
		for readlength, weight in readlength_dict.items():
			cfileh.seek(self.start, 0)
			sequence = cfileh.read(self.end-self.start)
			self.length += sum(l != '\0' and ord(l) <= readlength for l in sequence) * weight
	
	def setuniquelength(self, chromosomefile):
		try:
			cfileh = openfile(chromosomefile)
		except:
			global warnedchromosomes
			if not chromosomefile in warnedchromosomes:
				warnedchromosomes.append(chromosomefile)
				print "Warning: Did not find", chromosomefile # if a SyntaxError is here, check the python version, use version 2.5/2.6/2.7
			self.length = self.end - self.start
			return
		
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
			startfilepos = self.start + offset
			endfilepos = self.end + offset
		else:
			startfilepos = offset + (self.start // seqlength)*linelength + (self.start % seqlength)
			endfilepos = offset + (self.end // seqlength)*linelength + (self.end % seqlength)	
		cfileh.seek(startfilepos, 0)
		sequence = cfileh.read(endfilepos-startfilepos)
		self.length = sequence.count('A')+sequence.count('C')+sequence.count('G')+sequence.count('T')	# upper-case means unique
	
	def uniquelengthfromBigWig(self, chromosome):
		try:
			proc = subprocess.Popen([bigWigSummary_path, uniqueposdir, chromosome, str(self.start), str(self.end), '1'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		except KeyboardInterrupt:
			raise
		except:
			global bigwig_norun
			try: bigwig_norun
			except:
				if vocal: print 'Error: Failed running bigWigSummary, required by the -u option for bigWig files. If it was\'t aborted by the user, perhaps its permission to execute has not been set (\'chmod +x '+ bigWigSummary_path + '\' on Unix-like systems)'
				sys.exit(1)
				bigwig_norun = 1
			self.length = self.end - self.start
			return
		ret, err = proc.communicate()
		if proc.returncode:
			if err.startswith('no data'):
				self.length = 0
			else:
				if vocal:
					print "Warning: bigWigSummary returned an error:"
					print err
				self.length = self.end - self.start
		else:
			self.length = float(ret)*(self.end - self.start)
			

class Ctranscript:
	def __init__(self, ID, sortnum):
		self.ID = ID
		self.expression = []
		self.reads = []
		self.overlaptx = None
		self.genename = "?"
		self.sortnum = sortnum

class Cgene:
	def __init__(self, name):
		self.name = name
		self.transcripts = []
		self.exons = []
		self.chromosome = ""
		self.overlapsets = None

	def exonstring(self):
		outstr = self.name
		for exon in self.exons:
			outstr += "\t" + exon.exonstring()
		return outstr

class Cunit:
	# after collapse of genes
	def __init__(self, sortnum):
		self.name1 = "."
		self.name2 = "."
		self.rpkms = []
		self.reads = []
		self.sortnum = sortnum

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

def split_by_rmregions(exon, gene):
	if exon.start >= exon.end: return []
	r = Cregion(gene.chromosome, exon.start, exon.end)
	rmregions = r.overlapping()
	if not rmregions: return [exon]
	rmregion = rmregions[0]
	rightexon = Cexon(rmregion.end, exon.end, exon.normalise_with)
	rightexon.transcripts = exon.transcripts
	rightexon.forbidden = exon.forbidden
	exon.end = rmregion.start
	return split_by_rmregions(exon, gene) + split_by_rmregions(rightexon, gene)

def contractarrays(Ll, Rl):
	row = 0
	while row < len(Rl):
		zeropattern = Ll[row]
		(Ll, Rl) = contractarrays_inner(zeropattern, Ll, Rl, row)
		row += 1
	return (Ll, Rl)

def contractarrays_inner(zeropattern, Ll, Rl, firstmatchingrow):
	currentrow = firstmatchingrow+1
	while currentrow < len(Ll):
		matches = 1
		for place in range(len(zeropattern)):
			if (zeropattern[place] and not Ll[currentrow][place]) or \
			 (Ll[currentrow][place] and not zeropattern[place]):
				matches = 0
		if matches:
			for place in range(len(zeropattern)):
				Ll[firstmatchingrow][place] += Ll[currentrow][place]
			Rl[firstmatchingrow][0] += Rl[currentrow][0]
			del Ll[currentrow]
			del Rl[currentrow]
		else:
			currentrow += 1
	return (Ll, Rl)
		
def rowincludes(inclusiverow, includedrow):
	# checks if all elements of includedrows are found within inclusiverow
	for element in range(len(inclusiverow)):
		if includedrow[element] and not inclusiverow[element]:
			return 0
	return 1
	
def elementsinrow(row):
	count = 0
	for element in range(len(row)):
		if row[element]:
			count += 1
	return count
			
def removecommonrow(Ll, Rl):
	currentrow = 0
	while currentrow < len(Ll):
		matches = 1
		for place in range(len(Ll[currentrow])):
			if not Ll[currentrow][place]:
				matches = 0
		if matches:
			del Ll[currentrow]
			del Rl[currentrow]
		else:
			currentrow += 1
	return (Ll, Rl)
				
def removecolumn(Ll, column):
	currentrow = 0
	while currentrow < len(Ll):
		del Ll[currentrow][column]
		currentrow += 1		
	return Ll

def removezeros(Ll, Rl):
	currentrow = 0
	while currentrow < len(Ll):
		matches = 1
		for place in range(len(Ll[currentrow])):
			if Ll[currentrow][place]:
				matches = 0
		if matches:
			del Ll[currentrow]
			del Rl[currentrow]
		else:
			currentrow += 1
	return (Ll, Rl)

def gtf_field(p, tags):
	l1 = [v.strip() for v in p[-1].split(';')]
	for tag in tags:
		for f in l1:
			if f.startswith(tag):
				return f.split(' "')[-1].split('"')[0]
	return '.'

def fromannotationline(line, l_annotationtype=None):
	if l_annotationtype is None: l_annotationtype = annotationtype
	inferred_strand = False
	if l_annotationtype == 0:
		# from refGene.txt
		p = line.rstrip('\r\n').split("\t")
		exonstarts = [int(f) for f in p[9].split(",")[:-1]]	# start positions for exons for the gene
		exonends = [int(f) for f in p[10].split(",")[:-1]]
		ID = p[1]
		chromosome = p[2]
		genename = p[12]
		strand = p[3]
		cdsstart = min(int(p[7]), int(p[6]))
		cdsend = max(int(p[7]), int(p[6]))
	elif l_annotationtype == 1:
		# from knownGene.txt
		p = line.rstrip('\r\n').split("\t")
		exonstarts = [int(f) for f in p[8].split(",")[:-1]]
		exonends = [int(f) for f in p[9].split(",")[:-1]]
		ID = p[11]
		chromosome = p[1]
		genename = p[0]
		strand = p[2]
		cdsstart = min(int(p[5]), int(p[6]))
		cdsend = max(int(p[5]), int(p[6]))
	elif l_annotationtype == 2:
		# from ensGene.txt or sibGene.txt
		p = line.rstrip('\r\n').split("\t")
		ID = p[1]
		strand = p[3]
		chromosome = p[2]
		genename = p[12]
		cdsstart = min(int(p[7]), int(p[6]))
		cdsend = max(int(p[7]), int(p[6]))
		exonstarts = [int(f) for f in p[9].split(",")[:-1]]
		exonends = [int(f) for f in p[10].split(",")[:-1]]
	elif l_annotationtype == 3:
		# 3 column bed file
		p = line.rstrip('\r\n').split("\t")
		exonstarts = [int(p[1])]
		exonends = [int(p[2])]
		cdsstart = 0
		cdsend = 0
		chromosome = p[0]
		genename = None
		ID = p[0]+":"+p[1]+"-"+p[2]
		strand = "+"
		inferred_strand = True
	elif l_annotationtype == 4:
		# 6 column bed file
		p = line.rstrip('\r\n').split("\t")
		exonstarts = [int(p[1])]
		exonends = [int(p[2])]
		cdsstart = 0
		cdsend = 0
		chromosome = p[0]
		genename = None
		ID = p[3]
		if len(p) > 5: strand = p[5]
		else:
			strand = "+"
			inferred_strand = True
	elif l_annotationtype == 5:
		# gtf file format
		p = line.rstrip('\r\n').split("\t")
		exonstarts = [int(p[3])-1]
		exonends = [int(p[4])]
		cdsstart = 0
		cdsend = 0
		chromosome = p[0]
		transInfo = p[8].split(';')
		transIDs = transInfo[1].lstrip('transcript_id ').strip('"').replace(' ', '').split(',')
		ID = '+'.join(transIDs)
		genename = p[0]+":"+p[3]+"-"+p[4]
		strand = p[6]
		if strand == '.': strand = '+'
	elif l_annotationtype == 6:
		# up to 12 column bed file
		p = line.rstrip('\r\n').split("\t")
		chromosome = p[0]
		absstart = int(p[1])
		try: absend = int(p[2])
		except: absend = absstart
		try: strand = p[5]
		except:
			strand = "+"
			inferred_strand = True
		try:	cdsstart = int(p[6]); cdsend = int(p[7])
		except: cdsstart = 0; cdsend = 0
		try: blocksizes = map(int, p[10].split(',')) ; blockstarts = map(int, p[11].split(','))
		except: blocksizes = [absend - absstart]; blockstarts = [0]
		exonstarts = [absstart + rs for rs in blockstarts]
		exonends = [es + size for es,size in zip(exonstarts, blocksizes)]
		try: genename = p[3]
		except: 
			genename = None
		ID = str(chromosome)+":"+str(min(exonstarts)+1)+"-"+str(max(exonends)+1)
		if len(blockstarts) != len(blocksizes) or max(exonends) != absend:
			if vocal: print "Warning: Block structure for " + ID + " is malformed"
	elif l_annotationtype == 7:
		# ensembl gtf
		exonstarts = []
		exonends = []
		cdspoints = []
		for exonline, chrom, counter in line:
			p = exonline.rstrip('\r\n').split('#')[0].split("\t")
			pos1 = int(p[3])-1
			pos2 = int(p[4])-1
			chromosome = chrom
			ID = gtf_field(p, ['transcript_name', 'transcript_id'])
			genename = gtf_field(p, ['gene_name', 'gene_id'])
			strand = p[6]
			exontype = p[2]
			if exontype in ('start_codon', 'stop_codon'):
				if strand == '+':
					cdspoints.append(pos1)
				else:
					cdspoints.append(pos2)
			else:
				exonstarts.append(pos1)
				exonends.append(pos2)
			if exontype == 'CDS':
				cdspoints.extend([pos1,pos2])
		if cdspoints:
			cdsstart = min(cdspoints)
			cdsend = max(cdspoints)
		else:
			cdsstart = 0
			cdsend = 0
		
	if USESTRANDINFO:
		if genename is None: genename = strand
		if swapstrands:
			if strand == "+": chromosome += "-"
			else: chromosome += "+"
		else: 
			chromosome += strand	# so e.g. "chr1+" is different from "chr1-"
	else:
		if genename is None: genename = '.'
		inferred_strand = False
	return (chromosome, strand, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand)
		
class Cline:
	def __init__(self, line, counter):
		global allchromosomes, allchromosomes_dict
		
		self.line = line
		(chromosome, direction, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand) = fromannotationline(line)
		try: self.chromosome = allchromosomes_dict[chromosome]
		except:
			allchromosomes.append(chromosome)
			self.chromosome = allchromosomes.index(chromosome)
			allchromosomes_dict[chromosome] = self.chromosome
		self.start = min(exonstarts)
		self.end = max(exonends)
		self.strand = direction
		self.sortnum = counter
		
	def __cmp__(self, other):
		if self.chromosome != other.chromosome:
			return cmp(self.chromosome, other.chromosome)
		if self.start > other.end:
			return 1
		if other.start > self.end:
			return -1
		return cmp(self.start, other.start)

def openfile(filename, mode='r'):
	if filename.endswith(".gz"):
		import gzip
		fileh = gzip.open(filename, mode)
	elif filename.endswith(".bz2"):
		import bz2
		fileh = bz2.BZ2File(filename, mode)
	else: fileh = open(filename,mode)
	return fileh	

def addread(chromosome, start, end, strand, halfweight=False):
	if addchr: chromosome = 'chr'+chromosome
	try: chrID = allchromosomes_dict[chromosome]
	except:
		return -1, None
	endpos = 1-end if strand == '-' else end-1
	if halfweight:
		return chrID, (start, endpos, None, None)
	else:
		return chrID, (start, endpos)

def addread_tuple(chromosome, postuple):
	if addchr: chromosome = 'chr'+chromosome
	try: chrID = allchromosomes_dict[chromosome]
	except:
		return -1, None
	return chrID, postuple

def readsfrombedtabfile(filename, justreadnum, readlength_dict):
	try: 
		if removemultimappers:
			tmpfile = sorttotmpfile(filename, 3, ignoredprefix='track', sep='\t')
			fileh = open(tmpfile,'r')
		else:
			fileh = openfile(filename)
	except IOError:
		if vocal: print "Warning: No such file", filename
		return
	try:
		if not fileh.readline().startswith("track"): fileh.seek(0)
		while 1:
			line = fileh.readline()
			if not line: break
			if justreadnum: yield 1; continue
			p = line.rstrip('\r\n').split("\t")
			if minqual and  int(p[4]) < minqual: continue # too low score
			if readlength_dict is not None:
				if len(p) <= 10:
					readlength_dict[int(p[2])-int(p[1])] += 1
				else:
					readlength_dict[sum(map(int, p[10].split(',')))] += 1
			chromosome = p[0]
			try: strand = p[5]
			except: strand = '.'
			if USESTRANDINFO: chromosome += strand
			yield addread(chromosome, int(p[1]), int(p[2]), strand)
	finally:
		fileh.close()
		if removemultimappers:
			os.remove(tmpfile)

def readsfrombedtabfile_rs(filename, justreadnum, readlength_dict):
	if removemultimappers and vocal: print "Warning: The -unique option was ignored, unsupported for this format"
	try: fileh = openfile(filename)
	except:
		if vocal: print "Warning: No such file", filename
		return
	try:
		if not fileh.readline().startswith("track"): fileh.seek(0)
		while 1:
			line = fileh.readline()
			if not line: break
			p = line.rstrip('\r\n').split("\t")
			if minqual and  int(p[4]) < minqual: continue # too low score
			if readlength_dict is not None and (len(p) <= 3 or not 'jxn' in p[3]):
				readlength_dict[int(p[2])-int(p[1])] += 1
			chromosome = p[0]
			try: strand = p[5]
			except: strand = '.'
			if USESTRANDINFO: chromosome += strand
			for nb_reads in xrange(int(p[4])):
				yield addread(chromosome, int(p[1]), int(p[2]), strand)
	finally:
		fileh.close()
		if removemultimappers:
			os.remove(tmpfile)

def readsfrombedspacefile(filename, justreadnum, readlength_dict):
	try: 
		if removemultimappers:
			tmpfile = sorttotmpfile(filename, 3, ignoredprefix='track', sep=' ')
			fileh = open(tmpfile,'r')
		else:
			fileh = openfile(filename)
	except IOError:
		if vocal: print "Warning: No such file", filename
		return

	try:
		if not fileh.readline().startswith("track"): fileh.seek(0)
		while 1:
			line = fileh.readline()
			if not line: break
			if justreadnum: yield 1; continue
			p = line.rstrip('\r\n').split()
			if minqual and  int(p[4]) < minqual: continue # too low score
			if readlength_dict is not None:
				if len(p) <= 10:
					readlength_dict[int(p[2])-int(p[1])] += 1
				else:
					readlength_dict[sum(map(int, p[10].split(',')))] += 1
			chromosome = p[0]
			try: strand = p[5]
			except: strand = '.'
			if USESTRANDINFO: chromosome += strand
			yield addread(chromosome, int(p[1]), int(p[2]), strand)
	finally:
		fileh.close()

def readsfromgtffile(filename, justreadnum, readlength_dict):
	if removemultimappers and vocal: print "Warning: The -unique option was ignored, unsupported for this format"
	try: fileh = openfile(filename)
	except:
		if vocal: print "Warning: No such file", filename
		return
	try:
		if not fileh.readline().startswith("track"): fileh.seek(0)
		while 1:
			line = fileh.readline()
			if not line: break
			if justreadnum: yield 1; continue
			p = line.rstrip('\r\n').split("\t")
			if minqual and  int(p[5]) < minqual: continue # too low score
			if readlength_dict is not None:
				readlength_dict[int(p[4])-(int(p[3])-1)] += 1
			chromosome = p[0]
			try: strand = p[6]
			except: strand = '.'
			if USESTRANDINFO: chromosome += strand
			yield addread(chromosome, int(p[3])-1, int(p[4]), strand)
	finally:
		fileh.close()

def splitcigar(string):
	numarr = ['']
	symbolarr = []
	wasnum = 1
	for letter in string:
		if letter.isdigit():
			if wasnum:
				numarr[-1] += letter
			else:
				numarr.append(letter)
			wasnum = 1
		else:
			if not wasnum:
				symbolarr[-1] += letter
			else:
				symbolarr.append(letter)
			wasnum = 0
	if '' in numarr:
		if vocal: print 'Warning: strange CIGAR field in SAM file:', string
		return [0 if v == '' else int(v) for v in numarr], symbolarr
	return [int(v) for v in numarr], symbolarr

def rowparse(row):
	chromosome = row[2]
	start = row[1]-1
	endincl = start + row[4]-1
	if row[3] & 16:
		if USESTRANDINFO:
			chromosome += '+' if row[3] & 0x80 else "-"
		endincl = -endincl
	else:
		if USESTRANDINFO:	
			chromosome += '-' if row[3] & 0x80 else "+"
	return chromosome, start, endincl

def process_rowinfo(rowinfo, justreadnum):
	lastname = ''
	while 1:
		try: row = rowinfo.pop()
		except IndexError: break
		try: nrow = rowinfo[-1]
		except IndexError: pass
		else:
			if row[0] == nrow[0]:	# same name
				try:
					if rowinfo[-2][0] == row[0]:
						# multimapper
						try:
							while rowinfo[-1][0] == row[0]:
								rowinfo.pop()
						except IndexError: pass
						continue
				except IndexError: pass
				if row[3]&0x40 or row[3]&0x80:
					# paired end reads
					if row[4]&0x80:
						chromosome_1, start_1, endincl_1 = rowparse(nrow)
						chromosome_2, start_2, endincl_2 = rowparse(row)
					else:
						chromosome_1, start_1, endincl_1 = rowparse(row)
						chromosome_2, start_2, endincl_2 = rowparse(nrow)
					yield addread_tuple(chromosome_1, (start_1, endincl_1, start_2, endincl_1))
					continue
				# multimapper
				try:
					while rowinfo[-1][0] == row[0]:
						rowinfo.pop()
				except IndexError: pass
				continue
		chromosome_1, start_1, endincl_1 = rowparse(row)
		yield addread_tuple(chromosome_1, (start_1, endincl_1))

def readsfromsam(filename, justreadnum, readlength_dict):
	global saved_rowinfo
	try:
		if justreadnum: raise NameError
		rowinfo, fn = saved_rowinfo
		del saved_rowinfo
		assert fn == filename
	except NameError:
		try: 
			fileh = openfile(filename)
		except IOError:
			if vocal: print "Warning: No such file", filename
			return
		
		has_warned_cigar = False
		rowinfo = []
		for line in fileh:
			try:
				if line.startswith("@"): continue	# header
				p = line.rstrip('\r\n').split("\t")
				bits = int(p[1])
				if bits & 4: continue	# unmapped
				if p[5] == '*': continue # unmapped
				if minqual and int(p[4]) < minqual: continue # too low mapping quality
				if maxNM >= 0 and any((field.startswith('NM:i:') or field.startswith('nM:i:')) and int(field[5:]) > maxNM for field in p[11:]): continue # too many mismatches/indels
				if readlength_dict is not None:
					readlength_dict[len(p[9])] += 1
				chromosome = p[2]
				name = p[0]
				if p[5] == '*':
					if vocal and not has_warned_cigar: print 'Warning: line in SAM file has empty (*) CIGAR field but is not flagged as unmapped:', line.rstrip('\r\n')
					length = 1
					has_warned_cigar = True
				else:
					length = sum([l for l,s in zip(*splitcigar(p[5])) if s in ["M","D","X","=","N"]])
				rowinfo.append((name, int(p[3]), chromosome, bits, length))
			except:
				if vocal: print 'Error: could not recognize/parse this line as SAM format:', line.rstrip('\r\n')
				raise
		rowinfo.sort()
		import copy
		if justreadnum: saved_rowinfo = copy.copy(rowinfo), filename
		fileh.close()
		
	return process_rowinfo(rowinfo, justreadnum)

def pysamopen(filename):
	try:
		import pysam
	except:
		if vocal: print "Error: tried to read BAM file, but missing the pysam module"
		exit(1)
	try:
		return pysam.Samfile(filename, 'rb') # open as bam, gives ValueError for sam file
	except ValueError:
		return pysam.Samfile(filename, 'r') # open as sam

def readsfrombam(filename, justreadnum, readlength_dict):
	global saved_rowinfo
	try:
		if justreadnum: raise NameError
		rowinfo, fn = saved_rowinfo
		del saved_rowinfo
		assert fn == filename
	except NameError:
		try: 
			fileh = pysamopen(filename)
		except IOError:
			if vocal: print "Warning: No such file", filename
			return
		
		rowinfo = []
		for read in fileh:
			if read.flag & 4: continue	# unmapped
			if read.cigar is None: continue #unmapped, http://seqanswers.com/forums/showthread.php?t=3551
			if minqual and read.mapq < minqual: continue # too low mapping quality
			if maxNM >= 0 and any(((key=='NM' or key == 'nM') and val > maxNM) for key, val in read.tags): continue # too many mismatches/indels
			try:
				chromosome = fileh.getrname(read.rname)
			except ValueError:
				# probably the result of mapping to chrUn_random, chr9_random etc but then not mentioning theses reference sequences in the header
				continue
			if readlength_dict is not None:
				readlength_dict[read.rlen] += 1
			name = read.qname
			length = sum(l for o,l in read.cigar if o in [0,2,3]) #0=M,1=I,2=D,3=N
			rowinfo.append((name, read.pos+1, chromosome, read.flag, length))
		rowinfo.sort()
		import copy
		if justreadnum: saved_rowinfo = copy.copy(rowinfo), filename
		fileh.close()
		
	return process_rowinfo(rowinfo, justreadnum)

def readsfrombam_se(filename, justreadnum, readlength_dict):
	if removemultimappers and vocal: print "Warning: The -unique option was ignored, unsupported for this format"
	
	try: 
		fileh = pysamopen(filename)
	except IOError:
		if vocal: print "Warning: No such file", filename
		return
	
	try:
		for read in fileh:
			if read.flag & 4: continue	# unmapped
			if read.cigar is None: continue #unmapped
			if minqual and read.mapq < minqual: continue # too low mapping quality
			if maxNM >= 0 and any(((key=='NM' or key == 'nM') and val > maxNM) for key, val in read.tags): continue # too many mismatches/indels
			if (not mapends) and (read.flag & 0x80) and (read.flag & 0x2):
				# is last segment and the other segment is unmapped
				continue
			try:
				chromosome = fileh.getrname(read.rname)
			except ValueError:
				# probably the result of mapping to chrUn_random, chr9_random etc but then not mentioning theses reference sequences in the header
				continue
			if readlength_dict is not None:
				readlength_dict[read.rlen] += 1
			length = sum(l for o,l in read.cigar if o in [0,2,3]) #0=M,1=I,2=D,3=N
			strand = '-' if read.flag & 16 else '+'
			if USESTRANDINFO:
				if mapends and read.flag & 0x80:
					# swap strand for the 2nd read for paired-end sequencing
					chromosome += '-' if strand == '+' else '+'
				else:
					chromosome += strand
			start = read.pos
			yield addread(chromosome, start, start+length, strand, mapends and (read.flag & 2))
	finally:
		fileh.close()
	
	

def readsfromsam_se(filename, justreadnum, readlength_dict):
	try: 
		if removemultimappers:
			tmpfile = sorttotmpfile(filename, 0, ignoredprefix='@', sep='\t')
			fileh = open(tmpfile,'r')
		else:
			fileh = openfile(filename)
	except IOError:
		if vocal: print "Warning: No such file", filename
		return
	
	has_warned_cigar = False
	line = None
	try:
		while 1:
			line = fileh.readline()
			if not line: break
			if line.startswith("@"): continue	# header
			p = line.rstrip('\r\n').split("\t")
			read_flag = int(p[1])
			if read_flag & 4: continue	# unmapped
			if p[5] == '*':
				if vocal and not has_warned_cigar: print 'Warning: line in SAM file has empty (*) CIGAR field, will be considered as unmapped:', line.rstrip('\r\n')
				length = 1
				has_warned_cigar = True
				continue # unmapped
			if minqual and int(p[4]) < minqual: continue # too low mapping quality
			if maxNM >= 0 and any((field.startswith('NM:i:') or field.startswith('nM:i:')) and int(field[5:]) > maxNM for field in p[11:]): continue # too many mismatches/indels
			if (not mapends) and read_flag & 0x80 and read_flag & 2:
				# is last segment and the other segment is mapped
				continue
			if readlength_dict is not None:
				readlength_dict[len(p[9])] += 1
			if justreadnum: yield 1; continue
			chromosome = p[2]
			if int(p[1]) & 16: strand = "-"
			else: strand = "+"
			if USESTRANDINFO:
				if mapends and read_flag & 0x80:
					# swap strand for the 2nd read for paired-end sequencing
					chromosome += '-' if strand == '+' else '+'
				else:
					chromosome += strand
			
			cigarints, cigarsymbols = splitcigar(p[5])
			length = sum([l for l,s in zip(cigarints, cigarsymbols) if s in ["M","D","X","=","N"]])
			
			start = int(p[3])-1
			yield addread(chromosome, start, start+length, strand, mapends and (read_flag & 0x2))
	except:
		if line is not None:
			if vocal: print 'Error: could not recognize/parse this line as SAM format:', line.rstrip('\r\n')
		raise
	finally:
		fileh.close()
		if removemultimappers:
			os.remove(tmpfile)
	return

def readsfrombowtieoutput(filename, justreadnum, readlength_dict):
	try: 
		if removemultimappers:
			tmpfile = sorttotmpfile(filename, 0, sep='\t')
			fileh = open(tmpfile,'r')
		else:
			fileh = openfile(filename)
	except IOError:
		if vocal: print "Warning: No such file", filename
		return
	try:
		while 1:
			line = fileh.readline()
			if not line: break
			p = line.rstrip('\r\n').split("\t")
			if readlength_dict is not None:
				readlength_dict[len(p[4])] += 1
			chromosome = p[2]
			try: strand = p[1]
			except: strand = '.'
			if USESTRANDINFO: chromosome += strand
			yield addread(chromosome, int(p[3]), int(p[3])+len(p[4]), strand)
	finally:
		fileh.close()
		if removemultimappers:
			os.remove(tmpfile)
	return 

def sorttotmpfile(filename, idindex, ignoredprefix='', sep='\t'):
	# sort file, remove duplicates, export to temporary file
	# returned file needs later removal with os.remove()
	import tempfile
	outf = tempfile.mkstemp(suffix='_rpkmforgenes')[1]
	try:
		if vocal: print "Sorting to ", outf
		infh = openfile(filename, 'r')
		lines = []
		for line in infh:
			if ignoredprefix and line.startswith(ignoredprefix): continue
			p = line.rstrip('\r\n').split('\t')
			if ignoredprefix=='@' and int(p[2])&0x80: continue	# remove 2nd in paired read
			if len(p) <= idindex:
				lines.append(('', line))
			else:
				lines.append((p[idindex],line))
		infh.close()
		lines.sort()
		outfh = open(outf, 'w')
		for i in range(len(lines)):
			this = lines[i][0]
			if i > 0:
				neighbour = lines[i-1][0]
				if neighbour and neighbour == this: continue
			try:
				neighbour = lines[i+1][0]
			except IndexError: pass
			else:
				if neighbour and neighbour == this: continue
			outfh.write(lines[i][1])
		outfh.close()
	except:
		try: os.remove(outf)
		except: pass
		raise
	return outf

def mappos(chrindex, readpos):
	try: window = allwindows[chrindex][readpos//WINDOWSIZE]
	except IndexError:
		return []
	return [exon for exon in window if exon.start <= readpos < exon.end]
	
def make_overlapsets(currentgene, genecollapse, limitcollapse, limitcollapse2, add_overlapsets = None):
	currentgene.overlapsets = []
	currentgene_exons = [exon for exon in currentgene.exons if not exon.forbidden]
	# generate the sets of overlapping transcripts
	# firstly see which ones are connected by overlap
	for fromtx in currentgene.transcripts:	# fromtx = transcript others are connected from, the one to add to
		fromtx.overlaptx = [fromtx]
		for exon in currentgene_exons:
			if fromtx in exon.transcripts:
				for totx in exon.transcripts:	# totx = transcript others are connected to, the one to be added
					if not totx in fromtx.overlaptx:
						if genecollapse and fromtx.genename != totx.genename: continue
						fromtx.overlaptx.append(totx)
	
	# in the middle, add in genes who lost. or didn't lose, internal conection due to -rmnameoverlap
	# should have helped, doesn't seem to though
	if (add_overlapsets is not None) and not limitcollapse2:
		for overlapset in add_overlapsets:
			for txi, fromtx in enumerate(overlapset):
				for totx in overlapset[txi+1:]:
					if totx not in fromtx.overlaptx:
						fromtx.overlaptx.append(totx)
					if fromtx not in totx.overlaptx:
						totx.overlaptx.append(fromtx)
	
	# secondly follow the overlap to see indirect overlap
	donetx = []
	for fromtx in currentgene.transcripts:
		if fromtx in donetx: continue
		currentoverlaptx = []
		for overtx in fromtx.overlaptx:
			currentoverlaptx.append(overtx)
			donetx.append(overtx)
		hasadded = 1
		while hasadded:
			hasadded = 0
			for overtx in currentoverlaptx:
				if limitcollapse: break
				for totx in overtx.overlaptx:
					if totx in donetx: continue
					currentoverlaptx.append(totx)
					donetx.append(totx)
					hasadded = 1
		currentgene.overlapsets.append(currentoverlaptx)
	
	if sum([len(s) for s in currentgene.overlapsets]) != len(currentgene.transcripts):
		if vocal: print "Warning at gene model collapse. Info:", currentgene.name, [tx.ID for tx in currentgene.transcripts], [[tx.ID for tx in s] for s in currentgene.overlapsets]

def merge(o_infiles, o_outfile):
	# parse files
	numgenes = None
	header = ['#samples', '#allmappedreads', '#normalizationreads', '#arguments']
	header[3] += '\t' + ' '.join(sys.argv) + '\ttime: ' + time.asctime()
	genelines = []
	readlines = []
	for inf in o_infiles:
		lnumgenes = 0
		with openfile(inf, 'r') as infh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				if p[0] == '#samples':
					header[0] += ''.join(['\t' + v for v in p[1:]])
					lnumsamples = len(p) - 1
				elif p[0] == '#allmappedreads':
					header[1] += ''.join(['\t' + v for v in p[1:]])
				elif p[0] in ('#genemappedreads', '#normalizationreads'):
					header[2] += ''.join(['\t' + v for v in p[1:]])
				elif line[0] != '#':
					if numgenes is None:
						genelines.append('\t'.join(p[:2]))
						readlines.append('')
					genelines[lnumgenes] += ''.join('\t'+s for s in p[2:2+lnumsamples])
					readlines[lnumgenes] += ''.join('\t'+s for s in p[2+lnumsamples:2+2*lnumsamples])
					lnumgenes += 1
			if numgenes is None:
				numgenes = lnumgenes
			elif numgenes != lnumgenes:
				raise Exception('Unequal number of genes between files (try running without the -p option)')
			
	# write output
	with openfile(o_outfile, 'w') as outfh:
		for line in header:
			print >>outfh, line
		for linenum in xrange(numgenes):
			print >>outfh, genelines[linenum] + readlines[linenum]
			
def rmflag(argv, flag):
	# remove flag and list of arguments that follow the flag
	if flag not in argv: return
	flagindex = argv.index(flag)
	index = flagindex + 1
	del_index = [flagindex]
	while index < len(argv) and argv[index][0] != "-":
		del_index.append(index)
		index += 1
	for index in reversed(del_index):
		del argv[index]

def getarguments(flag):
	# list of arguments that follow the flag
	flagindex = sys.argv.index(flag)
	ret = []
	index = flagindex + 1
	while index < len(sys.argv) and sys.argv[index][0] != "-":
		ret.append(sys.argv[index])
		index += 1
	return ret
		
def getargument(flag):
	# the argument after the flag
	flagindex = sys.argv.index(flag)
	return sys.argv[flagindex+1]
	
def testargumentflag(flag):
	# if the flag is among the arguments
	if flag in sys.argv: return 1
	else: return 0

def main():
	global MAXGENES
	global MAXREADS
	global ONLYUNIQUEPOS
	global ONLYTXUNIQUEEXONS
	global INTRONSINSTEAD 
	global NODECONVOLUTION 
	global COLLAPSEGENES 
	global FULLTRANSCRIPTS
	global USESTRANDINFO, WINDOWSIZE
	global allchromosomes
	global allchromosomes_dict
	global annotationtype
	global uniqueposdir
	global swapstrands
	global mapends, midread
	global removemultimappers, minqual, maxNM
	global bigWigSummary_path, uniquenessfiletype
	global allwindows
	global vocal, addchr

	# interpret program arguments
	if len(sys.argv) < 2 or testargumentflag("--help") or testargumentflag("-h"):
		print "Non-optional arguments:"
		print " -o followed by output file"
		print " -i followed by list of input files (by default, guesses format from file extension)"
		print " -a followed by annotation file"
		print "Gene model-related options:"
		print " -u followed by a bigWig file, alternatively a directory for files for non-unique positions (lower case for nonunique k-mers (where k is the read length), upper case for unique; filenames are e.g. chr1.fa, can also be chr1_unique20-255.btxt etc"
		print " -no3utr to remove 3'UTRs"
		print " -fulltranscript to not remove 3'UTRs (default)"
		print " -maxlength followed by a distance to cut each transcript from the 3' end, from 5' if negative (never seems to give better values)"
		print " -maxgenes limit how many genes expression is calculated for (for testing purposes)"
		#print " -limitcollapse to not consider indirect transcript overlap"
		print " -namecollapse to only consider overlap between isoform with the same gene identifier (shaky)"
		print " -nocollapse to get isoform expressions (shaky)"
		print " -nooverlap to ignore that transcripts can overlap (will count some reads several times)"
		print " -rmnameoverlap to ignore regions shared my multiple genes (seems to work well)"
		print " -rmregions followed by a bed file of regions which should be removed from genes"
		print " -flat to flatten all isoforms to one gene model (likely to give too low RPKM values)"
		print " -namesum to sum together all entries with the same gene name even if they are apart on the genome"
		print " -txunique to ignore regions shared by multiple gene isoforms"
		print " -onlycoding to ignore noncoding transcripts"
		print " -swapstrands to make reads on + strand map to genes on - and vice versa (and sets -strand)"
		print " -introns gives gene expression from introns rather than exons (also removes exons of other isoforms)"
		print " -keephap to not remove haplotype chromosome (_hap/_alt) annotation"
		print " -norandom to remove genes on unplaced contigs"
		print "Annotation file formats:"
		print " -genePred if annotation file uses format of refGene.txt etc (default if cannot guess from file name suffix)"
		print " -bedann tab-separated 0-based bed file, chromosome start end and 9 optional fields"
		print " -ensgtfann for gtf format gene annotation (slow, generally better to convert to bed or genepred)"
		print " -addchr to change a read aligning to 1 to chr1, 2 to chr2 etc"
		print "Input formats (flags generally not needed, will detect from file ending):"
		print " -bed tab separated bed file (default if cannot guess from file name suffix)"
		print " -bedcompacted bed file with number of reads in the score column"
		print " -bedspace space separated bed file"
		print " -bowtie the default output format of bowtie"
#		print " -sam SAM format, and to remove non-unique hits" # broken read counting
#		print " -bam BAM format, and to remove non-unique hits" # broken read counting
		print " -samse SAM format, uniquely mapped reads (default for SAM)"
		print " -bamu BAM or SAM format, uniquely mapped reads (faster than -samse, default for BAM)"
		print " -gff GFF file with reads, no groups"
		print "Normalisation options:"
		print " -mRNAnorm to normalize by the number of reads matching mRNA exons (default)"
		print " -exonnorm to normalize by the number of reads matching exons, including ncRNA"
		print " -allmapnorm to normalize by the total number of mapped reads (default if annotation contains no mRNA)"
		print " -forcedtotal followed by a number of reads for each sample to set a constant to normalise by"
		print " -tpm to calculate tags per million instead of RPKM"
		print " -tmm to run TMM normalisation on top of RPKM (or TPM if -tmp is used) normalisation"
		print "Output format options:"
		print " -readcount to add the number of reads to the output"
		print " -csv [read_count_file_out] to output in csv format, optionally with a second output file with read counts in"
		print " -table another output format"
		print " -sortpos for output sorted by genome position"
		print " -exportann followed by a filename to write which exons have been used, also prints exon read counts for the last input file"
		print " -readpresent - to suppress zero count entries - Ramu"
		print "Read-related arguments:"
		print " -strand to use strand information of reads"
		print " -bothends to also map the end positions to genes, each end counted as 0.5 (or 0.25 for paired-end reads)"
		print " -bothendsceil to set -bothends but round the read count upward"
		print " -midread to use middle of the read as read position  (buggy for RNA)"
		print " -diffreads to count only one read if several have the same position, strand and length (use with -bam or -sam if paired-end; samtools rmdup is generally better)"
		print " -maxreads followed by maximum number of reads to be used"
		print " -randomreads to make -maxreads pick reads at random"
		print " -skipiffewreads to ignore samples with fewer than -maxreads reads in the input"
		print " -minqual followed by an integer, to restrict reads to minimum this mapping quality (for sam, bam) or score (for bed, gff), default use all"
		print " -maxNM followed by an integer, to restrict reads to maximum this edit distance (nM or NM flag in sam, bam), default use all"
		print " -addchr to add the text 'chr' to read chromosome field, to resolve Ensembl genome/UCSC annotation mismatch"
		print "Other optional arguments:"
		print " -n followed by list of sample names (input file names are otherwise used)"
		print " -p followed by number of files to process in parallel"
		print " -quite to skip progress messages and warnings"
		print " -h to print this message and quit"
		print "Special output values:"
		print " 0 gene has no reads, -1 gene has no exons"
		print " otherwise the output is in reads per kilobase and million mappable reads (or rather FPKM for paired-end reads)"
		print "Output:"
		print " gene symbol -tab- ID -tab- RPKM values [-tab- read count]"
		return 0
	try: outfile = getargument("-o")
	except:
		print "No outfile (-o)"
		return 1
	try: infiles = getarguments("-i")
	except:
		print "No infile (-i)"
		return 1
	try: annotationfile = getargument("-a")
	except:
		print "No annotationfile (-a)"
		return 1
	try:
		names = getarguments("-n")
		if len(names) != len(infiles):
			print "Not the same number of sample names as files with reads"
			return 1
	except: names = [filename.rsplit("/")[-1] for filename in infiles]
	compressionsuffixes = ['gz','bz2']
	vocal = not testargumentflag("-quite")
	use_tmm = testargumentflag("-tmm")
	if testargumentflag("-forcedtotal"):
		forcedtotal = [int(v) for v in getarguments("-forcedtotal")]
		if len(forcedtotal) != len(infiles):
			print "Error: -forcedtotal has fewer or more numbers than -i has files"
			return 1
	elif testargumentflag("-forcetotal"):
		forcedtotal = [int(v) for v in getarguments("-forcetotal")]
		if len(forcedtotal) != len(infiles):
			print "Error: -forcetotal has fewer or more numbers than -i has files"
			return 1
	else: forcedtotal = 0
	if testargumentflag("-exportann"):
		exportannotation = getargument("-exportann")
	else: exportannotation = 0
	
	if testargumentflag("-readcount"): addreadcount = 1
	else: addreadcount = 0
	if testargumentflag("-readpresent") and testargumentflag('-p') and vocal:
		print 'Warning: -readpresent and -p do not work together'
	if testargumentflag("-readpresent") and testargumentflag('-table'):
		print 'Error: -table and -p do not work together'
		sys.exit(1)
	elif testargumentflag("-readpresent") and len(infiles)>1 and vocal:
		print 'Warning: -readpresent will filter based on the first file (%s)'%names[0]
	OUTPUT_CSV = testargumentflag('-csv')
	if OUTPUT_CSV:
		READCOUNTCSV = getpotentialargument("-csv")
		if READCOUNTCSV is not None:
			addreadcount = 1
	if testargumentflag('-p'):
		# multi-processing mode
		# this instance of the program will be master and not perform calculations, only merging
		
		import tempfile, subprocess
		out = []
		prlist = []
		samplenames = names[:]
		inputfiles = infiles[:]
		if forcedtotal:
			ft_copy = forcedtotal[:]
		processes = int(getargument('-p'))
		restargs = sys.argv[:]
		rmflag(restargs, '-p')
		rmflag(restargs, '-n')
		rmflag(restargs, '-i')
		rmflag(restargs, '-o')
		rmflag(restargs, '-forcedtotal')
		rmflag(restargs, '-forcetotal')
		rmflag(restargs, '-exportann')
		rmflag(restargs, '-tmm')
		rmflag(restargs, '-csv')
		if addreadcount and not testargumentflag("-readcount"): restargs.append("-readcount")
		try:
			# launch child instances
			for procleft in range(processes, 0, -1):
				numfiles = len(inputfiles) // procleft
				if numfiles == 0: continue
				out.append(tempfile.mkstemp(suffix='_rpkmforgenes_out.txt')[1])
				cmd = [sys.executable] + restargs + ['-o', out[-1], '-i']
				names_child = []
				ft_child = []
				for i in range(numfiles):
					cmd.append(inputfiles.pop(0))
					names_child.append(samplenames.pop(0))
					if forcedtotal: ft_child.append(str(ft_copy.pop(0)))
				cmd += ['-n'] + names_child
				if exportannotation and len(inputfiles) == 0:
					# put -exportann back in for the last batch of files
					cmd += ['-exportann', exportannotation]
				if forcedtotal: cmd += ['-forcedtotal'] + ft_child
				prlist.append(subprocess.Popen(cmd))
		
			for pr in prlist:
				pr.wait()
				if pr.returncode: raise Exception, 'A child process did not finish normally, look somewhere further up to find its error message'
			
			# merge the output
			merge(out, outfile)
			if vocal: print 'Merged to', outfile
		
		finally:
			# delete temporary files
			for tmpfile in out:
				try: os.remove(tmpfile)
				except: pass
		
		if use_tmm:
			ret = run_tmm(outfile, outfile, copy_counts=addreadcount)
			if ret: return ret
			
		if OUTPUT_CSV:
			convert_to_csv(outfile, READCOUNTCSV)
		return 0
	
	if testargumentflag("-genePred") or testargumentflag("-refseq"): annotationtype = 0
	elif testargumentflag("-ucsc"): annotationtype = 1
	elif testargumentflag("-ensembl"): annotationtype = 2
	elif testargumentflag("-bedann"): annotationtype = 6
	elif testargumentflag("-bed3ann"): annotationtype = 3
	elif testargumentflag("-bed6ann"): annotationtype = 4
	elif testargumentflag("-gtfann") or testargumentflag("-gffann"): annotationtype = 5
	elif testargumentflag("-ensgtfann"): annotationtype = 7
	else:
		# if no annotation type given, guess from file suffix
		anndict = {'bed':6, 'gff':5, 'gtf':5, 'txt':0}
		try:
			annpre, ending = annotationfile.rsplit('.',1)
			if ending in compressionsuffixes: ending = annpre.rsplit('.',1)[1]
			annotationtype = anndict[ending]
		except:
			annotationtype = 0
	
	if testargumentflag("-blocksize"): WINDOWSIZE = int(getargument("-blocksize"))
	else: WINDOWSIZE = 5000
	if testargumentflag("-diffreads"): usediffreads = 1
	else: usediffreads = 0
	if testargumentflag("-minqual"): minqual = int(getargument("-minqual"))
	else: minqual = 0
	if testargumentflag("-maxNM"): maxNM = int(getargument("-maxNM"))
	else: maxNM = -1
	if testargumentflag("-bed"): infileformat = "bed"
	elif testargumentflag("-bedcompacted"): infileformat = "bedcompacted"
	elif testargumentflag("-bedspace"): infileformat = "bedspace"
	elif testargumentflag("-gff") or testargumentflag("-gtf"): infileformat = "gtf"
	elif testargumentflag("-sam"): infileformat = "sam"
	elif testargumentflag("-bowtie"): infileformat = "bowtie"
	elif testargumentflag("-samse"): infileformat = "samse"
	elif testargumentflag("-samu"): infileformat = "samse"
	elif testargumentflag("-bamu"): infileformat = "bamse"
	elif testargumentflag("-bamse"): infileformat = "bamse"
	elif testargumentflag("-bam"): infileformat = "bam"
	else: infileformat = "guess"
		
	if testargumentflag("-nocollapse"): COLLAPSEGENES = 0
	else: COLLAPSEGENES = 1
	if testargumentflag("-limitcollapse"): limitcollapse = 1
	else: limitcollapse = 0
	if testargumentflag("-limitcollapse2"): limitcollapse2 = 1   # like -limitcollapse, this is for legacy behaviour
	else: limitcollapse2 = 0
	genecollapse = testargumentflag("-namecollapse")
	if testargumentflag("-fulltranscript"): FULLTRANSCRIPTS = 1
	elif testargumentflag("-no3utr"): FULLTRANSCRIPTS = 0
	else: FULLTRANSCRIPTS = 1
	bigWigSummary_path = 0
	if testargumentflag("-dynuniq"):
		uniquenessfiletype = 'Helena'
		arr = getarguments('-dynuniq')
		if len(arr) == 0:
			global_readlength = None
		else:
			global_readlength = int(arr[0])
	elif testargumentflag("-ulen"):
		uniquenessfiletype = 'Helena'
		arr = getarguments('-ulen')
		if len(arr) == 0:
			global_readlength = None
		else:
			global_readlength = int(arr[0])
	else:
		uniquenessfiletype = 'fasta'
		global_readlength = 0
	if testargumentflag("-u"):
		uniqueposdir = getargument("-u")
		if uniqueposdir.endswith('.bw') or uniqueposdir.endswith('.bw.gz') or uniqueposdir.endswith('.bw.bz2') or uniqueposdir.endswith('.bigWig') or uniqueposdir.endswith('.bigWig.gz') or uniqueposdir.endswith('.bigWig.bz2'):
			uniquenessfiletype = 'bigWig'
			bigWigSummary_path = os.path.split(sys.argv[0])[0]
			if bigWigSummary_path == '':
				bigWigSummary_path = './'
			bigWigSummary_path = os.path.join(bigWigSummary_path, 'bigWigSummary')
			if not os.path.isfile(bigWigSummary_path):
				print 'Error: ' + bigWigSummary_path + ' does not exist, needed to read bigWig file specified by -u option, try downloading from http://hgdownload.cse.ucsc.edu/admin/exe/'
				return 1
		elif uniquenessfiletype == 'fasta':
			# detect if it uses variable length uniqueness files
			if any(f.endswith('.btxt') for f in os.listdir(uniqueposdir)):
				global_readlength = None
				uniquenessfiletype = 'Helena'
		ONLYUNIQUEPOS = 1
	else: ONLYUNIQUEPOS = 0
	if testargumentflag("-maxgenes"): MAXGENES = int(getargument("-maxgenes"))
	else: MAXGENES = 0
	if testargumentflag("-maxreads"):
		MAXREADS = int(getargument("-maxreads"))
		randomreads = testargumentflag("-randomreads")
		skipiffewreads = testargumentflag("-skipiffewreads")
	else: MAXREADS = 0; randomreads = 0; skipiffewreads = 0
	if testargumentflag("-table"): outputformat = 'table'
	else: outputformat = 'v2'
	if testargumentflag("-readpresent"): readpresent = True 
	else: readpresent = False
	if testargumentflag("-intronsinstead") or testargumentflag("-introns"): INTRONSINSTEAD = 1
	else: INTRONSINSTEAD = 0
	if testargumentflag("-strand"): USESTRANDINFO = 1
	else: USESTRANDINFO = 0
	if testargumentflag("-swapstrands"):
		swapstrands = 1
		USESTRANDINFO = 1
	else: swapstrands = 0
	removemultimappers = testargumentflag("-unique")
	if testargumentflag("-maxlength"): maxgenelength = int(getargument("-maxlength"))
	else: maxgenelength = 0
	if testargumentflag("-onlycoding"): useonlycoding = 1
	else: useonlycoding = 0
	if testargumentflag("-sortpos"): originalorder = 0
	else: originalorder = 1
	if testargumentflag("-mRNAnorm") or testargumentflag("-exonnorm"): genereadnormalise = 1
	elif testargumentflag("-allmapnorm"): genereadnormalise = 0
	else: genereadnormalise = 1
	exonnorm_include_ncRNA = testargumentflag("-exonnorm")
	if testargumentflag("-nooverlap") or testargumentflag("-nodeconvolution"): NODECONVOLUTION = 1
	else: NODECONVOLUTION = 0
	if testargumentflag("-txunique"): ONLYTXUNIQUEEXONS = 1
	else: ONLYTXUNIQUEEXONS = 0
	if testargumentflag("-rmnameoverlap"): removenameoverlap = 1
	else: removenameoverlap = 0
	if testargumentflag("-bothendsceil"):
		mapends = 1
		mapends_ceil = 1
	elif testargumentflag("-bothends"):
		mapends = 1
		mapends_ceil = 0
	else:
		mapends = 0
		mapends_ceil = 0
	midread = testargumentflag("-midread")
	flattengenes = testargumentflag("-flat") or testargumentflag("-tpm") # auto-set of calculating TPM
	use_tpm = testargumentflag("-tpm")
	fullgenenamesum = testargumentflag("-namesum")
	ignoredchrfragments = []
	if not testargumentflag("-keephap"):
		ignoredchrfragments.append("_hap")
		ignoredchrfragments.append("_alt")
		keephap = False
	else:
		keephap = True
	if testargumentflag("-norandom"):
		ignoredchrfragments.append("_random")
		ignoredchrfragments.append("chrUn_")
		keeprandom = False
	else:
		keeprandom = True
	if testargumentflag("-addchr"): addchr = 1
	else: addchr = 0
	allchromosomes = []
	usedreadspersample = []
	totalreadspersample = []
	normreadfactors = []
	allchromosomes_dict = {}
	
	if annotationtype == 7: # -ensgtfann
		chrom_translation = {'X':'chrX', 'Y':'chrY', 'MT':'chrM'}
		for i in range(1,1000):
			chrom_translation[str(i)] = 'chr%d'%i
		from collections import defaultdict
		lines_per_ID = defaultdict(list)
	
	annotationfileh = openfile(annotationfile)
	counter = 0
	annlines = []
	for line in annotationfileh:
		if line[0] == '#': continue
		if MAXGENES and counter > MAXGENES: break
		if annotationtype in [3,4,6] and line.startswith("track"): continue # bed
		counter += 1
		if annotationtype == 7: # -ensgtfann
			p_ann = line.rstrip('\r\n').split('\t')
			ID = gtf_field(p_ann, ['transcript_name', 'transcript_id'])
			chrom = chrom_translation.get(p_ann[0], p_ann[0])
			if '_MHC_' in chrom:
				if not keephap: continue # throw out haplotype chromosomes
			elif not keeprandom:
				if p_ann[0] not in chrom_translation:
					continue # throw out _random chromsosomes / unplaced contigs
			lines_per_ID[(chrom, ID)].append((line, chrom, counter))
		else:
			newline = Cline(line, counter)
			annlines.append(newline)
	annotationfileh.close()
	if annotationtype == 7: # -ensgtfann
		annlines = [Cline(L, min(e[-1] for e in L)) for L in lines_per_ID.values()]
	annlines.sort()
	if vocal: print "Has sorted entries from gene annotation file"

	global warnedchromosomes
	warnedchromosomes = []

	genes = []
	currentgene = Cgene("none")
	oldend = -1000
	annotation_includes_mRNA = 0
	for annline in annlines:
		line = annline.line
		
		# file format-specific code
		(chromosome, direction, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand) = fromannotationline(line)
		
		# remove _hap and maybe _random chromosomes
		disallowedchr = 0
		for a in ignoredchrfragments:
			if a in chromosome: disallowedchr = 1
		if disallowedchr: continue

		# if -onlycoding argument then skip any noncoding transcripts
		if useonlycoding and cdsstart == cdsend: continue

		# initiate the transcript and gene
		transcript = Ctranscript(ID, annline.sortnum)
		if oldend < min(exonstarts) or chromosome != currentgene.chromosome or NODECONVOLUTION:
			currentgene = Cgene(genename)
			genes.append(currentgene)
			currentgene.chromosome = chromosome
			oldend = max(exonends)
		else:
			oldend = max([max(exonends), oldend])
		currentgene.transcripts.append(transcript)
		transcript.genename = genename
		
		if maxgenelength > 0:
			# trim from 3' end
			totaltxlength = 0
			if direction == "+":
				for exoni in range(len(exonstarts)-1, -1, -1):
					if totaltxlength >= maxgenelength:
						exonstarts = exonstarts[exoni+1:]
						exonends = exonends[exoni+1:]
						break
					else:
						totaltxlength += exonends[exoni] - exonstarts[exoni]
						if totaltxlength > maxgenelength:
							exonstarts[exoni] += totaltxlength - maxgenelength
			else:
				for exoni in range(len(exonstarts)):
					if totaltxlength >= maxgenelength:
						exonstarts = exonstarts[:exoni]
						exonends = exonends[:exoni]
						break
					else:
						totaltxlength += exonends[exoni] - exonstarts[exoni]
						if totaltxlength > maxgenelength:
							exonends[exoni] -= totaltxlength - maxgenelength
		elif maxgenelength < 0:
			# trim from 5' end
			totaltxlength = 0
			if direction == "-":
				for exoni in range(len(exonstarts)-1, -1, -1):
					if totaltxlength >= -maxgenelength:
						exonstarts = exonstarts[exoni+1:]
						exonends = exonends[exoni+1:]
						break
					else:
						totaltxlength += exonends[exoni] - exonstarts[exoni]
						if totaltxlength > -maxgenelength:
							exonstarts[exoni] += totaltxlength + maxgenelength
			else:
				for exoni in range(len(exonstarts)):
					if totaltxlength >= -maxgenelength:
						exonstarts = exonstarts[:exoni]
						exonends = exonends[:exoni]
						break
					else:
						totaltxlength += exonends[exoni] - exonstarts[exoni]
						if totaltxlength > -maxgenelength:
							exonends[exoni] -= totaltxlength + maxgenelength

		for exoni in range(len(exonstarts)):
			if cdsstart != cdsend and not FULLTRANSCRIPTS:
				# remove 3'UTR
				if direction == "+" and exonends[exoni] <= cdsend:
					newexon = Cexon(exonstarts[exoni], exonends[exoni], True)
				elif direction == "+" and exonstarts[exoni] < cdsend and exonends[exoni] > cdsend:
					newexon = Cexon(exonstarts[exoni], cdsend, True)
				elif direction == "-" and exonstarts[exoni] >= cdsstart:
					newexon = Cexon(exonstarts[exoni], exonends[exoni], True)
				elif direction == "-" and exonstarts[exoni] < cdsstart and exonends[exoni] > cdsstart:
					newexon = Cexon(cdsstart, exonends[exoni], True)
				else:
					continue
			elif cdsstart != cdsend or exonnorm_include_ncRNA:
				newexon = Cexon(exonstarts[exoni], exonends[exoni], True)
			else:
				newexon = Cexon(exonstarts[exoni], exonends[exoni], False)
			if cdsstart != cdsend:
				annotation_includes_mRNA = 1

			if INTRONSINSTEAD:
				# add exon to forbidden set
				newexon.forbidden=True

			currentgene.exons.append(newexon)
			newexon.transcripts.append(transcript)
			
		if INTRONSINSTEAD:
			for exoni in range(len(exonstarts)):		
				# add intron to allowed set
				if exoni + 1 == len(exonstarts):
					continue
				newexon = Cexon(exonends[exoni], exonstarts[exoni+1], False)
				currentgene.exons.append(newexon)
				newexon.transcripts.append(transcript)
				
	annotationfileh.close()
	if vocal: print "Has read positions from gene annotation file"
	
	uniquepos_filesuffix = None
	if ONLYUNIQUEPOS and uniquenessfiletype == 'Helena':
		# find the right value for uniquepos_filesuffix
		for suffix in ('.genomeGeneLevelMerge.unique20-255.btxt','.unique20-255.btxt', '_unique20-255.btxt'):
			if any(f.endswith(suffix) for f in os.listdir(uniqueposdir)):
				uniquepos_filesuffix = suffix
				break
		else:
			for chromosome in sorted(list(allchromosomes), key=lambda s: -len(s)):
				if USESTRANDINFO: chromosome = chromosome[:-1]
				for f in os.listdir(uniqueposdir):
					if f.startswith(chromosome) and f.endswith('.btxt'):
						uniquepos_filesuffix = f[len(chromosome):]
						break
				if uniquepos_filesuffix is not None: break
			if vocal: print 'Guessed -u file suffix:', uniquepos_filesuffix
	
	if not annotation_includes_mRNA and not (testargumentflag("-mRNAnorm") or testargumentflag("-exonnorm")):
		genereadnormalise = 0
	
	# read file with problem regions to remove from annotated genes
	if testargumentflag("-rmregions"):
		removalfile = getargument("-rmregions")
		with openfile(removalfile, 'r') as infh:
			for line in infh:
				chromosome, direction, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand = fromannotationline(line, 6) # load bed format
				for start, end in zip(exonstarts, exonends):
					r = Cregion(chromosome, start, end)
					r.addtowindows()
				
				# if -strand and -rmregion file doesn't specify strand, remove from both strands
				strand1, strand2 = ('-','+') if swapstrands else ('+','-')
				if inferred_strand:
					if chromosome[-1] != strand1:
						if vocal: print "Warning: something wrong in annotation parsing"
						continue
					for start, end in zip(exonstarts, exonends):
						r = Cregion(chromosome[:-1]+strand2, start, end)
						r.addtowindows()
	
	# decompress bigwig file if compressed
	has_decompressed_bigwig = 0
	if bigWigSummary_path and (uniqueposdir.endswith(n) for n in ('.bw.gz','.bw.bz2', '.bigWig.gz', '.bigWig.bz2')):
		import tempfile
		f = tempfile.mkstemp(suffix='_rpkmforgenes')[1]
		h = open(f,'w')
		zh = openfile(uniqueposdir, 'r')
		h.write(zh.read())
		h.close()
		zh.close()
		uniqueposdir = f
		has_decompressed_bigwig = 1
	try:
		# sort out isoform overlaps
		for gene in genes:
			newexons = []
			# remove from a list of problem regions
			for ei, exon in enumerate(gene.exons):
				overlap = Cregion.overlappingpoint(gene.chromosome, exon.end)
				if overlap:
					exon.end = min(r.start for r in overlap)
				overlap = Cregion.overlappingpoint(gene.chromosome, exon.start)
				if overlap:
					exon.end = max(r.end for r in overlap)
				if exon.start >= exon.end:
					continue
				newexons_here = split_by_rmregions(exon, gene)[1:] # [1:] so it doesn't add 'exon' which is already in the list
				newexons.extend(newexons_here)
			gene.exons.extend(newexons) 
			
			#gene.exons = [exon for exon in gene.exons if exon.start < exon.end]
			
			# generate possible exons
			borders = list(set([exon.start for exon in gene.exons] + [exon.end for exon in gene.exons]))
			borders.sort()
			possibleexons = [Cexon(start, end, False) for start, end in zip(borders[:-1], borders[1:])]
		
			# populate exons
			for e_exon in gene.exons:
				for p_exon in possibleexons:
					# possibleexons shorter or as long as gene.exons elements at this stage, so little checks
					if e_exon.start <= p_exon.start < e_exon.end:
						p_exon.transcripts.extend(e_exon.transcripts)
						p_exon.normalise_with = p_exon.normalise_with or e_exon.normalise_with
						p_exon.forbidden = p_exon.forbidden or e_exon.forbidden
		
			# remove constitutive introns
			gene.exons = [exon for exon in possibleexons if len(exon.transcripts) > 0]
			
			if removenameoverlap:
				make_overlapsets(gene, genecollapse, limitcollapse, limitcollapse2, None)
				for overlapset in gene.overlapsets:
					num_exons = 0
					removal_set = []
					for exon in gene.exons:
						if any(tx in overlapset for tx in exon.transcripts):
							num_exons += 1
							if len(set(tx.genename for tx in exon.transcripts)) > 1:
								removal_set.append(exon)
					if len(removal_set) != num_exons:
						for exon in removal_set:
							exon.forbidden = True # side effect: the normalizationreads value will be higher than the sum of genes' read counts, by ~1%
			
			# forbidden exons do not belong to transcripts
			for exon in gene.exons:
				if exon.forbidden:
					exon.transcripts = []
			
			if removenameoverlap:
				singlegeneoverlaps = []
				for overlapset in gene.overlapsets:
					for genename in set(tx.genename for tx in overlapset):
						transcripts = [tx for tx in overlapset if tx.genename == genename]
						if len(transcripts) >= 2:
								singlegeneoverlaps.append(transcripts)
			else:
				singlegeneoverlaps = None
			
			make_overlapsets(gene, genecollapse, limitcollapse, limitcollapse2, singlegeneoverlaps)
			
		# generate overlapsets
		if COLLAPSEGENES or flattengenes:
			for gene in genes:
				make_overlapsets(gene, genecollapse, limitcollapse, limitcollapse2, singlegeneoverlaps)
			
			# flatten to one isoform if requested
			if flattengenes:
				for currentgene in genes:
					newtxset = []
					for overlapset in currentgene.overlapsets:
						containedexons = []
						containedgenenames = []
						for tx in overlapset:
							if tx.genename not in containedgenenames:
								containedgenenames.append(tx.genename)
						containedgenenames.sort()
						name1 = "+".join(containedgenenames)
						name2 = "+".join([tx.ID for tx in overlapset])
						flattx = Ctranscript(name2, min(tx.sortnum for tx in overlapset))
						flattx.genename = name1
						newtxset.append(flattx)
						for exon in currentgene.exons:
							for tx in exon.transcripts:
								if tx in overlapset:
									containedexons.append(exon)
									break
						for exon in containedexons:
							exon.transcripts.append(flattx)
					for exon in currentgene.exons:
						exon.transcripts = [tx for tx in exon.transcripts if tx in newtxset]
					currentgene.transcripts = newtxset
	finally:
		if has_decompressed_bigwig:
			try: os.remove(uniqueposdir)
			except: pass
	
	allwindows = [[] for chromosome in allchromosomes]
	for gene in genes:
		for exon in gene.exons:
			chrID = allchromosomes_dict[gene.chromosome]
			w_start = exon.start//WINDOWSIZE
			w_end = exon.end//WINDOWSIZE
			lw_chr = len(allwindows[chrID])
			if lw_chr <= w_end:
				allwindows[chrID].extend([] for i in xrange(w_end-lw_chr+1))
			for w_i in range(w_start, w_end+1):
				allwindows[chrID][w_i].append(exon)
	if vocal: print "Has indexed exons"
	
	has_not_calc_ulen = True
	
	# get the reads
	if skipiffewreads: infiles_i_skip = []
	for infile_i, infile in enumerate(infiles):
		from collections import defaultdict
		readlength_dict = defaultdict(int) if global_readlength is None else None
		for gene in genes:
			for exon in gene.exons:
				exon.differentreadsarray = []
				exon.reads = 0
		
		if infileformat == "guess":
			# if no input type given, guess from file suffix
			anndict = {'bed':'bed', 'gff':'gtf', 'gtf':'gtf', 'bam':'bamse'}
			try:
				annpre, ending = infile.rsplit('.',1)
				if ending in compressionsuffixes: ending = annpre.rsplit('.',1)[1]
				if ending == 'sam':
					try:
						import pysam
					except ImportError:
						infileformat_l = 'samse'
					else:
						infileformat_l = 'samse' #before: 'bamse'
				else:
					infileformat_l = anndict[ending]
			except:
				infileformat_l = "bed"
		else:
			infileformat_l = infileformat
		
		if infileformat_l == "bed":
			readgen = readsfrombedtabfile
		elif infileformat_l == "bedspace":
			readgen = readsfrombedspacefile
		elif infileformat_l == "bedcompacted":
			readgen = readsfrombedtabfile_rs
		elif infileformat_l == "gtf":
			readgen = readsfromgtffile
		elif infileformat_l == "sam":
			readgen = readsfromsam
		elif infileformat_l == "samse":
			readgen = readsfromsam_se
		elif infileformat_l == "bamse":
			readgen = readsfrombam_se
		elif infileformat_l == "bam":
			readgen = readsfrombam
		elif infileformat_l == "bowtie":
			readgen = readsfrombowtieoutput
		else:
			if vocal: print "Error: unrecognized format:", infileformat
			return 1
		
		usedreads = 0
		totalreads = 0
		
		if randomreads or skipiffewreads:
			l_totalreads = sum(1 for x in readgen(infile, True, None))
			if skipiffewreads and l_totalreads < MAXREADS:
				infiles_i_skip.append(infile_i)
				continue
			if randomreads:
				inclusionlist = tuple(v>=MAXREADS for v in numpy.random.permutation(l_totalreads))
		for chrID, read_tuple in readgen(infile, False, readlength_dict):
			totalreads += 1
			if randomreads:
				if inclusionlist[totalreads-1]: continue
			elif totalreads == MAXREADS: break
			if chrID == -1: continue
			
			if midread:
				if len(read_tuple) == 2:
					read_tuple = ((abs(read_tuple[0])+abs(read_tuple[1]))//2,)
				elif len(read_tuple) == 4:
					if read_tuple[2] is None:
						read_tuple = ((abs(read_tuple[0])+abs(read_tuple[1]))//2, None)
					else:
						read_tuple = ((abs(read_tuple[0])+abs(read_tuple[1]))//2, (abs(read_tuple[2])+abs(read_tuple[3]))//2)
				else:
					if vocal: print "Warning: problem with read_tuple length and -midread"
			
			if mapends:
				exons = []
				for pos in read_tuple:
					if pos is None: continue
					exons.extend(mappos(chrID, abs(pos)))
			else:
				exons = mappos(chrID, read_tuple[0])
				
			if usediffreads and len(exons) > 0:
				if read_tuple in exons[0].differentreadsarray:
					continue
				else:
					exons[0].differentreadsarray.append(read_tuple)
			
			for exon in exons:
				if mapends:
					exon.reads += 1.0/len(read_tuple)
				else:
					exon.reads += 1
		if randomreads: totalreads = min(totalreads, MAXREADS)   # changed 7 Jun 2018
		if vocal: print "Has compared reads to exons"
		
		if has_not_calc_ulen and global_readlength is not None:
			# fixed read length
			calcuniqlength(genes, {global_readlength: 1}, uniquepos_filesuffix)
			has_not_calc_ulen = False
		elif global_readlength is None:
			# multiple read lengths determined by the sample
			# different for different samples
			calcuniqlength(genes, readlength_dict, uniquepos_filesuffix)
		
		# RPKM calculation
		for gene in genes:
			for exon in gene.exons:
				exon.readspersample.append(exon.reads)
				if exon.normalise_with: usedreads += exon.reads
			
			if sum(not exon.forbidden for exon in gene.exons) == 0:
				for tx in gene.transcripts:
					tx.expression.append(-1)
					tx.reads.append(0)
				continue
			for exon in gene.exons:
				if exon.reads > 0 and not exon.forbidden: break
			else: # no reads found
				for tx in gene.transcripts:
					tx.expression.append(0.0)
				continue
				
			# deconvolute gene isoforms
			#  L is the lengths
			#  R is the reads
			#  D is the densities
			#  L*D = R, D=(LT*L)^-1 * LT * R, where ^-1 is inversion and T is transposition

			# first prepare the arrays
			Ll = []	# list that is to become the matrix L
			for exon in gene.exons:
				if exon.forbidden: continue
				Lc = []	# column of L
				for tx in gene.transcripts:
					if tx in exon.transcripts:
						Lc.append(float(exon.length))
					else:
						Lc.append(0.0)
				Ll.append(Lc)
			Rl = []	# list that is to become the matrix R
			for exon in gene.exons:
				if exon.forbidden: continue
				Rl.append([float(exon.reads)])
			
			Dindices = range(len(gene.transcripts))	# which element in gene.transcripts that an element in D corresponds to
			redofornegative = 1
			counter = 0
			
			if ONLYTXUNIQUEEXONS:
				# only use exons that are only found in one isoform
				row = 0
				while row < len(Rl):
					if elementsinrow(Ll[row]) != 1:
						del Ll[row]
						del Rl[row]
					else:
						row += 1
				if len(Rl) == 0:
					redofornegative = 0
					Dindices = []
			
			while redofornegative:			
				# contract the arrays
				(Ll, Rl) = contractarrays(Ll, Rl)
					
				if len(Ll) == 0:
					Dindices = []
					break
					
				# find isoforms that have no reads and remove them
				readthreshold = 0
				zeroreadcolumns = []
				for column in range(len(Ll[0])):
					zeroreadcolumns.append(0)
				row = 0
				while row < len(Rl):
					nonzeros = 0
					crow = 0
					while crow < len(Rl):
						if rowincludes(Ll[row], Ll[crow]):
							if Rl[crow][0] > readthreshold:
								nonzeros += 1
						crow += 1
					if nonzeros == 0:
						for element in range(len(Ll[row])):
							if Ll[row][element]:
								zeroreadcolumns[element] = 1
					row += 1
				count = 0	# to go through the rows with the fewest exons first
				while count < len(Rl):
					row = 0
					while row < len(Rl):
						if Rl[row][0] > readthreshold and elementsinrow(Ll[row]) == count and \
						 rowincludes(zeroreadcolumns, Ll[row]):
							for element in range(len(Ll[row])):
								if Ll[row][element]:
									zeroreadcolumns[element] = 0
						row += 1		
					count += 1			
				columnsremoved = 0	
				
				for column in range(len(Ll[0])):
					if zeroreadcolumns[column] and len(Ll[0]) > 1:
						Ll = removecolumn(Ll, column-columnsremoved)
						del Dindices[column-columnsremoved]
						columnsremoved += 1	
				if columnsremoved:
					(Ll, Rl) = contractarrays(Ll, Rl)
					if len(Ll) > 1:
						(Ll, Rl) = removezeros(Ll, Rl)
			
				# create the matrices and do the matrix calculations
				L = matrix(Ll)
				R = matrix(Rl)
				LTL = L.transpose() * L
				inverseturn = 0
				
				# these three lines fixed a problem 24 Apr 2012 with massively too large RPKM values for certain samples and genes
				# that problem appears to result from a singular matrix which contains rounding errors and therefore didn't make linalg.inv raise an exception, but produce invLTL with lots of big values
				# this code forces the algorithm to treat borderline genes as giving a potentially singular matrix, instead of testing for exception with linalg.inv
				# the number is set to be large enough in relation to the numbers in LTL to not be rounded away
				# side-effect is an 10^-6 error on average for RPKM values
				# 
				try: detsign, logdet = linalg.slogdet(LTL)
				except AttributeError:
					determinant = numpy.linalg.det(LTL)
					if determinant == 0:
						detsign = 0
						logdet = -100
					else:
						detsign = 1 if determinant > 0 else -1
						logdet = numpy.log(abs(determinant))
				if detsign == 0 or logdet - 2*math.log(numpy.mean(abs(LTL))) < -4:
					for ii in range(len(LTL)):
						LTL[ii, ii] += numpy.mean(abs(LTL))*0.0000001
				
				while 0 <= inverseturn <= 1:
					try:
						invLTL = linalg.inv(LTL)
					except:
						# this exception should no longer get thrown
						if inverseturn == 1:
							if vocal: print "Could not take the inverse of LTL for", gene.name
						for ii in range(len(LTL)):
							LTL[ii, ii] += 1	# bit cheating, but this ends up with the same expression estimate for both tx, though calculated with a slightly longer length
						inverseturn += 1
					else:
						D = invLTL * L.transpose() * R
						inverseturn += 1000
				
				# check if any density is lower than zero, end the loop if none
				redofornegative = 0
				negatives = []
				for index in range(len(Dindices)):
					if D[index, 0] <= 0:
						negatives.append(index)
						redofornegative = 1
						
				# remove columns that give negative values from the matrix
				# (those densities will instead be set to zero)
				if redofornegative:
					columnsremoved = 0
					for negativepos in negatives:
						Ll = removecolumn(Ll, negativepos-columnsremoved)
						del Dindices[negativepos-columnsremoved]
						columnsremoved += 1
						
					(Ll, Rl) = removezeros(Ll, Rl)
				counter += 1
				if counter == 1000:
					if vocal: print "Warning: high counter in negloop"
			
			# calculate expression from densities
			Danswers = numpy.zeros(len(gene.transcripts))	# so that found to be negative are set to zero
			for index in range(len(Dindices)):
				Danswers[Dindices[index]] = D[index, 0]
			for index in range(len(gene.transcripts)):
				tx = gene.transcripts[index]
				density = Danswers[index]
				tx.expression.append(1e9 * density / totalreads)
				txreads = 0
				for exon in gene.exons:
					if exon.forbidden: continue
					if ONLYTXUNIQUEEXONS and len(exon.transcripts) > 1:
						continue
					if tx in exon.transcripts:
						txreads +=exon.reads
				tx.reads.append(txreads)
		
		if vocal: print "Has calculated expression values for", infile
		if genereadnormalise and not forcedtotal:
			if vocal: print "Reads for normalisation:", usedreads
			if usedreads == 0:
				if vocal: print "Warning: no reads for normalisation, will use _1 read_ as normalisation factor"
				usedreads = 1
		
		usedreadspersample.append(usedreads)
		totalreadspersample.append(totalreads)
	
	if skipiffewreads:
		infiles = [f for i,f in enumerate(infiles) if i not in infiles_i_skip]
		names = [f for i,f in enumerate(names) if i not in infiles_i_skip]
	
	tunits = []
	if COLLAPSEGENES and not flattengenes:
		for currentgene in genes:
			for overlapset in currentgene.overlapsets:
				containedgenenames = []
				tunit = Cunit(min(tx.sortnum for tx in overlapset))
				tunits.append(tunit)
				for tx in overlapset:
					if tx.genename not in containedgenenames:
						containedgenenames.append(tx.genename)
				containedgenenames.sort()
				tunit.name1 = "+".join(containedgenenames)
				tunit.name2 = "+".join([tx.ID for tx in overlapset])
				tunit.rpkms = [sum([f.expression[t] for f in overlapset]) for t in range(len(tx.expression))]
				tunit.reads = [0.0 for rpkm in tunit.rpkms]
				for exon in currentgene.exons:
					if exon.length == 0: continue # don't count reads that didn't affect the rpkm value, can occur if non-unique exon with reads extending into intron
					if exon.forbidden: continue #doesn't do anything
					if sum([1 for tx in exon.transcripts if tx in overlapset]):
						for ii in range(len(tunit.rpkms)):
							tunit.reads[ii] += exon.readspersample[ii]
	else:	
		for gene in genes:
			for tx in gene.transcripts:
				tunit = Cunit(tx.sortnum)
				tunit.name1 = tx.genename
				tunit.name2 = tx.ID
				tunit.rpkms = tx.expression
				tunit.reads = [0.0 for rpkm in tunit.rpkms]
				for exon in gene.exons:
					if exon.length == 0: continue # don't count reads that didn't affect the rpkm value, can occur if non-unique exon with reads extending into intron
					if exon.forbidden: continue #doesn't do anything
					if tx in exon.transcripts:
						for ii in range(len(tunit.rpkms)):
							tunit.reads[ii] += exon.readspersample[ii]
				if use_tpm:
				    tunit.rpkms = [reads*1e6/totalreads for reads, totalreads in zip(tunit.reads, totalreadspersample)]
				tunits.append(tunit)

	# change normalisation if requested
	for si in range(len(infiles)):
		normreadfactors.append(totalreadspersample[si])
		if forcedtotal or genereadnormalise:
			if forcedtotal:
				samplefactor = totalreadspersample[si]/float(forcedtotal[si])
				normreadfactors[si] = forcedtotal[si]
			elif genereadnormalise:
				samplefactor = totalreadspersample[si]/float(usedreadspersample[si])
				normreadfactors[si] = usedreadspersample[si]
			for tu in tunits:
				if tu.rpkms[si] >= 0:
					tu.rpkms[si] *= samplefactor
					if vocal and tu.rpkms[si] > 1e6:
						print 'Warning: %s has a very high RPKM value of %f'%(tu.name1, tu.rpkms[si])
	
	# sum up expression for each gene name
	if fullgenenamesum:
		tunits.sort(key=lambda unit: unit.name1)
		tui = 1
		while tui < len(tunits):
			if tunits[tui].name1 == tunits[tui-1].name1:
				tunits[tui-1].name2 += '++' + tunits[tui].name2
				tunits[tui-1].rpkms = [r1+r2 for r1,r2 in zip(tunits[tui-1].rpkms, tunits[tui].rpkms)]
				tunits[tui-1].reads = [r1+r2 for r1,r2 in zip(tunits[tui-1].reads, tunits[tui].reads)]
				tunits[tui-1].sortnum = min(tunits[tui-1].sortnum, tunits[tui].sortnum)
				del tunits[tui]
			else:
				tui += 1
	
	# output to file
	outfileh = openfile(outfile, "w")
	if outputformat == 'v2':
		print >>outfileh, "\t".join(["#samples"]+names)
		print >>outfileh, "\t".join(["#allmappedreads"] + [str(v) for v in totalreadspersample])
		if forcedtotal:
			print >>outfileh, "\t".join(["#normalizationreads"] + [str(v) for v in forcedtotal])
		else:
			print >>outfileh, "\t".join(["#normalizationreads"] + [str(v) for v in normreadfactors])
		print >>outfileh, "#arguments\t"+" ".join(sys.argv) + "\ttime:" + time.asctime() + '\tversion:'+lastmodified
	
		if originalorder:
			tunits.sort(key=lambda unit: unit.sortnum)
	
		for tu in tunits:
			if readpresent:
				if tu.reads[0] <= 0: continue # Ramu, skip 0 count reads
			if addreadcount:
				if mapends_ceil:
					print >>outfileh, '\t'.join([tu.name1, tu.name2] + [str(v) for v in tu.rpkms] + [str(int(math.ceil(v))) for v in tu.reads])
				else:
					print >>outfileh, '\t'.join([tu.name1, tu.name2] + [str(v) for v in tu.rpkms] + [str(v)[:-2] if str(v).endswith('.0') else str(v) for v in tu.reads])
			else:
				print >>outfileh, '\t'.join([tu.name1, tu.name2] + [str(v) for v in tu.rpkms])
	elif outputformat == 'table':
		print >>outfileh, "\t".join(["Gene"]+names)
		if originalorder:
			tunits.sort(key=lambda unit: unit.sortnum)
		for tu in tunits:
			if any(lvl < 0 for lvl in tu.rpkms): continue
			if tu.name2 in ('.','+','-',''):
				name = tu.name1
			elif tu.name1 in ('.','+','-',''):
				name = tu.name2
			else:
				name = tu.name1.replace('+',' ')+' '+tu.name2
			if addreadcount:
				if mapends_ceil:
					print >>outfileh, '\t'.join([name] + [str(int(math.ceil(v))) for v in tu.reads])
				else:
					print >>outfileh, '\t'.join([name] + [str(v)[:-2] if str(v).endswith('.0') else str(v) for v in tu.reads])
			else:
				print >>outfileh, '\t'.join([name] + [str(v) for v in tu.rpkms])
	outfileh.close()
	if vocal: print "Has written to", outfile
	
	# print information about the gene annotation used to file
	if exportannotation:
		outfileh = openfile(exportannotation, 'w')
		print >>outfileh, "#numbers for " + infiles[-1]
		print >>outfileh, '#genename\tchromosome\texons'
		print >>outfileh, '#exon:\tstart\tend\tlength\treads\ttranscript_IDs'
		for gene in genes:
			exons = [exon for exon in gene.exons if not exon.forbidden]
			print >>outfileh, gene.name + '\t' + gene.chromosome +'\t' + str(len(exons))
			for exon in exons:
				print >>outfileh, '\t'.join(map(str,['exon:', exon.start, exon.end, exon.length, exon.reads] + [tx.ID+'|'+tx.genename for tx in exon.transcripts]))
	
	if use_tmm:
		ret = run_tmm(outfile, outfile, copy_counts=addreadcount)
		if ret: return ret
		
	if OUTPUT_CSV:
		convert_to_csv(outfile, READCOUNTCSV)
	
	return 0

def calcuniqlength(genes, readlength_dict, filesuffix):
	readlength_dict_weights = dict((r,float(n)/sum(readlength_dict.values())) for r,n in readlength_dict.items())
	for gene in genes:
		# forbidden exons do not belong to transcripts
		# calculate length for the rest
		for exon in gene.exons:
			if not exon.forbidden:
				exon.calclength(gene.chromosome, readlength_dict_weights, filesuffix)

if '__main__' == __name__:
	ret = main()
	if ret:	sys.exit(ret)
