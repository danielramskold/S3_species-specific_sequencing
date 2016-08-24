from __future__ import division
import argparse, dr_tools, numpy, math

def log2(v):
	return math.log(v, 2)

def M(gi, Y_k, Y_r, N_k, N_r):
	return log2((Y_k[gi]/N_k)/(Y_r[gi]/N_r))

def M_logdivlog(gi, Y_k, Y_r, N_k, N_r):
	return log2(Y_k[gi]/N_k)/log2(Y_r[gi]/N_r)

def w(gi, Y_k, Y_r, N_k, N_r):
	return (N_k - Y_k[gi])/N_k/Y_k[gi] + (N_r - Y_r[gi])/N_r/Y_r[gi]

def A(gi, Y_k, Y_r, N_k, N_r):
	return 0.5*log2(Y_k[gi]/N_k * Y_r[gi]/N_r) if Y_k[gi] > 0 else -10000

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('outfile')
	parser.add_argument('--ref_samples', nargs='+', metavar='samplename')
	parser.add_argument('--copy_counts', action='store_true', help='does not work with stdin as input')
	parser.add_argument('--run_on_counts', action='store_true')
	o = parser.parse_args()
	
	expr_in = dr_tools.loadexpr(o.infile, counts=o.run_on_counts)
	
	ref_samples = expr_in.samples if o.ref_samples is None else o.ref_samples
	Y_r = [numpy.mean([expr_in[s][gi] for s in ref_samples]) for gi in range(len(expr_in['symbols']))]
	N_r = sum(Y_r)
	
	expr_out = dr_tools.Parsed_rpkms([], False)
	normalization_factors = []
	
	for s in expr_in.samples:
		Y_k = expr_in[s]
		N_k = sum(Y_k)
		nonzero = [gi for gi in range(len(expr_in['symbols'])) if Y_k[gi] > 0 and Y_r[gi] > 0]
		A_distr = sorted((A(gi, Y_k, Y_r, N_k, N_r), gi) for gi in nonzero)
		M_distr = sorted((M(gi, Y_k, Y_r, N_k, N_r), gi) for gi in nonzero)
		
		Gstar = set(gi for A_val,gi in A_distr[int(0.05*len(A_distr)):-int(0.05*len(A_distr))]) & set(gi for M_val,gi in M_distr[int(0.3*len(M_distr)):-int(0.3*len(M_distr))])
		
		if len(nonzero) == 0: f_k = 1
		else:
			log2TMM = sum(w(gi, Y_k, Y_r, N_k, N_r) * M(gi, Y_k, Y_r, N_k, N_r) for gi in Gstar)/sum(w(gi, Y_k, Y_r, N_k, N_r) for gi in Gstar)
			f_k = 2**log2TMM # multipy non-reference by this value
		
		#print s, f_k
		
		expr_out[s] = [Y_k[gi]*f_k for gi in range(len(expr_in['symbols']))]
		normalization_factors.append(f_k)
	expr_out.allmappedreads = expr_in.allmappedreads
	expr_out.normalizationreads = expr_in.normalizationreads
	expr_out.samples = expr_in.samples
	expr_out['symbols'] = expr_in['symbols']
	expr_out['IDs'] = expr_in['IDs']
	
	dr_tools.writeexpr(o.outfile, expr_out, counts_expr=(dr_tools.loadexpr(o.infile, counts=True) if o.copy_counts else None), extra_comment_lines=[dr_tools.join('#TMM_normalization_factors', normalization_factors)])
	

