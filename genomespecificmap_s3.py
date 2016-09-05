import subprocess, os, argparse, tempfile, sys
from collections import defaultdict

# 21 May 2015: first version
# 24 Aug 2016: made it work for paire-end reads
# 5 Sep 2016: added -pe option

version = '5 Sep 2016'

def open_as_sam(filename):
	if filename.endswith('.sam'):
		return open(filename, 'rU')
	else:
		pr = subprocess.Popen(['samtools', 'view', '-h', filename], stdout=subprocess.PIPE)
		return pr.stdout

def write_as_bam(filename, fai=None):
	if fai is None:
		pr = subprocess.Popen(['samtools', 'view', '-bSo', filename, '-'], stdin=subprocess.PIPE)
	else:
		pr = subprocess.Popen(['samtools', 'view', '-bt', fai, '-o', filename, '-'], stdin=subprocess.PIPE)
	return pr.stdin

def sort_and_index(filename):
	subprocess.check_call(['samtools', 'sort', filename, filename[:-4]])
	subprocess.check_call(['samtools', 'index', filename])

def twocomesfirst(num):
	return -1000 if num == 2 else num

def has_least_mismatch(mismatches_per_file, max_mismatch_num_unique, min_mismatch_diff, max_mismatch_num_shared):
	# can return None, 's1', 's2' etc, where None means insufficiently dissimilar
	
	if not mismatches_per_file: return None
	
	# prefer in order 2 reads > 1 read > 3+ reads > 0 reads, then by sum of mismatches
	label_and_mm = sorted((twocomesfirst(len(e[1])), sum(mm for mm,line in e[1]), e[0]) for e in mismatches_per_file.items())
	#ordered = sorted(mismatches_per_file.items(), key = lambda e: twocomesfirst(len(e[1])), sum(mm for mm,line in e[1]))
	label_and_mm = [(numreads_2as0, mmsum, label) for numreads_2as0, mmsum, label in label_and_mm]
	
	if label_and_mm[0][0] >= 3: return None
	
	if not label_and_mm[0][1] <= max_mismatch_num_unique: return None
	
	if len(label_and_mm) == 1 or (label_and_mm[1][1] - label_and_mm[0][1] >= min_mismatch_diff and label_and_mm[1][1] > max_mismatch_num_shared):
		return label_and_mm[0][2]
	return None

def remove_genome_junction_redundant_mappingreads(merged_samfile, tmp_file, max_mismatch_num_unique, max_mismatch_num_shared, min_mismatch_diff):
	
	# sort by readname
	subprocess.check_call(['sort', merged_samfile,'-o', tmp_file])
	reads = []
	reads2 = []
	lastreadname = None
	if o.pairedend:
		mismatches_per_file = defaultdict(list)
	
	with open(tmp_file) as infh:
		for line in infh:
			if line.startswith('@'): continue
			parts = line.strip().split("\t")
			readname = parts[0].split()[0]
			flags = parts[11:]
			mismatches = 0
			for flag in flags:
				if flag.startswith("nM:i"): # STAR
					mismatches = int(flag.split(":")[-1])
					break
			else:
				for flag in flags:
					if flag.startswith("NM:i"): # bowtie
						mismatches = int(flag.split(":")[-1])
						break
			
			if readname != lastreadname:
				# process and write unique hits (if any)
				if o.pairedend:
					if has_least_mismatch(mismatches_per_file, max_mismatch_num_unique, min_mismatch_diff, max_mismatch_num_shared) == 's1':
						for mm, read in sorted(reads):
							if read.endswith('s1\n'):
								yield read
								break
						for mm, read in sorted(reads2):
							if read.endswith('s1\n'):
								yield read
								break
				else:
				
				
					if len(reads)>0:
						reads.sort()
						if reads[0][0] <= max_mismatch_num_unique and (len(reads) == 1 or (reads[1][0] - reads[0][0] >= min_mismatch_diff and reads[1][0] > max_mismatch_num_shared)):
							# matches with <= x mismatches and other read has >= y mismatches or is unmapped
							yield reads[0][-1]
					if len(reads2)>0:
						reads2.sort()
						if reads2[0][0] <= max_mismatch_num_unique and (len(reads2) == 1 or (reads2[1][0] - reads2[0][0] >= min_mismatch_diff and reads2[1][0] > max_mismatch_num_shared)):
							# matches with <= x mismatches and other read has >= y mismatches or is unmapped
							yield reads2[0][-1]
				reads = []
				reads2 = []
				if o.pairedend:
					mismatches_per_file = defaultdict(list)
			if int(parts[1]) & 0x80:
				reads2.append((mismatches, line))
			else:
				reads.append((mismatches, line))
			if o.pairedend:
				mismatches_per_file[parts[-1]].append((mismatches, line))
			lastreadname = readname
		
		if len(reads)>0:
			reads.sort()
			if reads[0][0] <= max_mismatch_num_unique and (len(reads) == 1 or (reads[1][0] - reads[0][0] >= min_mismatch_diff and reads[1][0] > max_mismatch_num_shared)):
				# matches with <= x mismatches and other read has >= y mismatches or is unmapped
				yield reads[0][-1]
		if len(reads2)>0:
			reads2.sort()
			if reads2[0][0] <= max_mismatch_num_unique and (len(reads2) == 1 or (reads2[1][0] - reads2[0][0] >= min_mismatch_diff and reads2[1][0] > max_mismatch_num_shared)):
				# matches with <= x mismatches and other read has >= y mismatches or is unmapped
				yield reads2[0][-1]
	

def single_run(o):
	o.bam_out = o.strain1_sam_out.endswith('.bam')
	header = []
	tmpfile_merged = tempfile.mkstemp(suffix='_genomespec_m.sam')[1]
	tmpfile_sorted = tempfile.mkstemp(suffix='_genomespec_s.sam')[1]
	try:
		with open(tmpfile_merged, 'w') as outfh:
			infh = open_as_sam(o.unique_strain1_bam_in)
			for line in infh:
				if line.startswith('@'): continue
				if o.minqual:
					parts = line.strip().split("\t")
					if int(parts[4]) < o.minqual: continue
				print >>outfh, line[:-1]+ '\ts1'
			infh.close()
			for fi, inf in enumerate(o.merged_strain2_bam_in):
				infh = open_as_sam(inf)
				for line in infh:
					if line.startswith('@'):
						header.append(line)
					else:
						print >>outfh, line[:-1]+ '\ts%d'%(fi+2)
				infh.close()
		
		if o.bam_out:
			outfh = write_as_bam(o.strain1_sam_out, fai=o.fai)
		else:
			outfh = open(o.strain1_sam_out, 'w')
		for line in header:
			outfh.write(line)
		for line in remove_genome_junction_redundant_mappingreads(tmpfile_merged, tmpfile_sorted, o.max_mismatch_num_unique, o.max_mismatch_num_shared, o.min_mismatch_diff):
			samline = line[:-4]
			strain = line[-3:-1]
			if strain == 's1':
				print >>outfh, samline
		if o.bam_out:
			outfh.close()
		else:
			outfh.close()
		
	finally:
		os.remove(tmpfile_merged)
		os.remove(tmpfile_sorted)

def batch_run(o):
	import taskmanager
	tasks = taskmanager.Tasklist(12, True)
	
	for d in os.listdir(o.unique_strain1_bam_in):
		files = ['-u', os.path.join(o.unique_strain1_bam_in, d, d+'_unique.bam')]
		for merged_strain2_bam_in in o.merged_strain2_bam_in:
			if o.bowtie:
				files += ['-m', os.path.join(merged_strain2_bam_in, d, d+'_merged.sam')]
			else:
				files += ['-m', os.path.join(merged_strain2_bam_in, d, 'Aligned.out.bam')]
		if o.batch_one_folder:
			files += ['-o', os.path.join(o.strain1_sam_out, d+o.batch_suffix+'.sam')]
		else:
			try: os.mkdir(os.path.join(o.strain1_sam_out, d))
			except: pass
			files += ['-o', os.path.join(o.strain1_sam_out, d, d+o.batch_suffix+'.sam')]
		if os.path.exists(files[-1]):
			print 'Skipping', d
			continue
		options = ['-Mu', str(o.max_mismatch_num_unique), '-Ms', str(o.max_mismatch_num_shared), '-D', str(o.min_mismatch_diff), '--minqual', str(o.minqual)]
		if o.pretend: print ' '.join(['python', sys.argv[0]]+files+options)
		else: tasks.add(subprocess.call, (['python', sys.argv[0]]+files+options,), sample=d)

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('-i', '--in', dest='unique_strain1_bam_in', help='input sam or bam file for species1', required=True, metavar='in_species1.sam/bam')
	opts.add_argument('-f', '--filter', dest='merged_strain2_bam_in', nargs='+', help='input sam or bam file for species2', required=True, metavar='in_otherspecies.sam/bam')
	opts.add_argument('-o', '--out', dest='strain1_sam_out', help='output sam or bam file, reads for species1', required=True, metavar='out_species1.sam/bam')
	opts.add_argument('-D', '--minmismatchdiff', dest='min_mismatch_diff', type=int, help='minumim number of mismatch different between the species\' alignments (default: 2, or 4 if -pe)', metavar='min_mismatch_difference')
	opts.add_argument('-Mu', '--maxmismatches-unique', dest='max_mismatch_num_unique', type=int, help='maximum number of mismatches to species1 (default: 1, or 2 if -pe)', metavar='max_species1_mismatches')
	opts.add_argument('-Ms', '--maxmismatches-shared', dest='max_mismatch_num_shared', type=int, help='read is removed if it has up to this many mismatches to species2 (default: 1, or 2 if -pe)', metavar='max_otherspecies_mismatches_filter')
	#opts.add_argument('-b', '--batchmode', dest='batch_mode', action='store_true', help='run several files, specify -u, -m, -o as folders (e.g. bowtie_hg19) instead of sam files')
	#opts.add_argument('-s', '--batchsuffix', dest='batch_suffix', default='')
	#opts.add_argument('-S', '--batch1folder', dest='batch_one_folder', action='store_true')
	opts.add_argument('-q', '--minqual', dest='minqual', type=int, default=255, metavar='min_mapping_qual')
	#opts.add_argument('--pretend', action='store_true')
	#opts.add_argument('--bowtie', action='store_true')
	#opts.add_argument('--bam_out', action='store_true')
	opts.add_argument('--fai', help='.fai file from samtools faidx, for bam output if the input sam file is missing the header', metavar='species1.fa.fai')
	opts.add_argument('-pe', '--paired-end', dest='pairedend', action='store_true', help='to consider the number of mismatches together for both reads of a pair')
	o = opts.parse_args()
	
	if o.max_mismatch_num_unique is None:
		o.max_mismatch_num_unique = 2 if o.pairedend else 1
	if o.max_mismatch_num_shared is None:
		o.max_mismatch_num_shared = 2 if o.pairedend else 1
	if o.min_mismatch_diff is None:
		o.min_mismatch_diff = 4 if o.pairedend else 2
	
	
	single_run(o)
	#if o.batch_mode:
	#	batch_run(o)
	#else:
	#	if o.pretend: raise Exception
	#	single_run(o)
	

