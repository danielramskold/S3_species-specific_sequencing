import sys

def getargument(flag, default=None):
	testhelpargument()
	try:
		flagindex = sys.argv.index(flag)
	except:
		if default is None: raise Exception, flag + " is not in argument"
		else: return default
	try:
		return sys.argv[flagindex+1]
	except:
		raise Exception, flag + " is not followed by a value"

def stringfromflag(flags, default=None):
	testhelpargument()
	for flag in flags:
		try:
			return getargument(flag, None)
		except: continue
	if default is None: raise Exception, '/'.join(flags) + " is not in argument"
	return default

def flag(flag, default=None):
	return stringfromflag(flag.split("/"), default)
	
def flagarray(flag, default=None):
	return arrayfromflag(flag.split("/"), default)

def ifflag(flag, ifpresent=1, ifabsent=0):
	return setifflag(flag.split("/"), ifpresent, ifabsent)

def intfromflag(flags, default=None):
	testhelpargument()
	for flag in flags:
		try:
			return int(getargument(flag, None))
		except: continue
	if default is None: raise Exception, '/'.join(flags) + " is not in argument"
	return default
	
def floatfromflag(flags, default=None):
	testhelpargument()
	for flag in flags:
		try:
			return float(getargument(flag, None))
		except: continue
	if default is None: raise Exception, '/'.join(flags) + " is not in argument"
	return default

def arrayfromflag(flags, default=None):
	testhelpargument()
	for flag in flags:
		try:
			return getargumentarray(flag, None)
		except: continue
	if default is None: raise Exception, '/'.join(flags) + " is not in argument"
	return default
	
def setifflag(flags, ifpresent=1, ifabsent=0):
	testhelpargument()
	for flag in flags:
		if flag in sys.argv: return ifpresent
	return ifabsent

def setifargument(flag, ifpresent=1, ifabsent=0):
	testhelpargument()
	if flag in sys.argv: return ifpresent
	else: return ifabsent
	
def getargumentarray(flag, default=None):
	testhelpargument()
	try:
		flagindex = sys.argv.index(flag)
	except:
		if default is None: raise Exception, flag + " is not in argument"
		else: return default
	nextflag = flagindex+1
	while nextflag < len(sys.argv):
		if sys.argv[nextflag][0] == "-" and len(sys.argv[nextflag]) > 1 and sys.argv[nextflag][1] not in "0123456789.": break
		nextflag += 1
	return sys.argv[flagindex+1:nextflag]
	
def testhelpargument():
	if "-h" in sys.argv:
		raise Exception, "-h is among arguments"
	if "--help" in sys.argv:
		raise Exception, "--help is among arguments"
		
def helpflag():
	testhelpargument()


if '__main__' == __name__:
	try:
		testhelpargument()
		infile = getargument("-i")
		outfiles = getargumentarray("-o", 0)
		giveID = setifargument("-ID")
		threshold = float(getargument("-t"))
		multiplegenes = setifargument("-m")
	except Exception, string:
		print string
		print "Arguments: -i <file in> [-o <genelist files out ('skip' for no output)>] -t <threshold rpkm> [-ID (for ID instead of symbol)] [-m (for print all names of a gene)]"
		exit()
	
	infileh = open(infile, "r")
	samples = infileh.readline()[:-1].split("\t")
	if samples[0][0] == "#":
		samples = samples[1:]
	
	if outfiles and len(samples) != len(outfiles):
		print "There should be as many output files as samples ("+str(len(samples))+")"
		infileh.close()
		exit()
	
	genelists = []
	for s in samples:
		genelists.append([])
	
	for line in infileh:
		if line[0] == '#': continue
		p = line[:-1].split("\t")
		
		for si in range(len(samples)):
			if float(p[2+si]) > threshold:
				if giveID: genestxt = p[1]
				else: genestxt = p[0]
				genes = genestxt.split("+")
				if not multiplegenes: genes = [genes[0]]
				for gene in genes:
					genelists[si].append(gene)
	infileh.close()
	
	if outfiles:
		for si in range(len(outfiles)):
			if outfiles[si] == "skip": continue
			outfileh = open(outfiles[si], "w")
			for gene in genelists[si]:
				print >>outfileh, gene
			outfileh.close()
		print "Done"
	else:
		for si in range(len(samples)):
			print str(len(genelists[si]))+ "\t"+samples[si]
