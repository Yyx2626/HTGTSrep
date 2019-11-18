
#### Usage: python3 this.py <F|B> <should_annotate_prob> <codon_usage.txt> <input_sequences.txt>
#### Options:
####    <F|B>  forward (from nt to aa) or backward (from aa to any possible nt)
####    <should_annotate_prob>  True or False
#### Output: STDOUT
####   if <should_annotate_prob>==False, one column: output sequence
####   if <should_annotate_prob>==True, two columns: output sequence and 
####      posterior ratio of input sequence to the max posteior of backtranslated sequence, for forward translation
####      posterior probability of each backtranslated sequence, for backward translation


import sys

FB, should_annotate_prob, codon_filename, input_filename = sys.argv[1:]

def str2bool(x):
	try:
		x = int(x)
		return x > 0
	except ValueError:
		return x.lower() in ('true', 't', 'yes', 'y')

def regularize_FB(FB):
	FB = FB.upper()
	if FB[0] == 'F':
		FB = 'F'
	elif FB[0] == 'B':
		FB = 'B'
	else:
		print('Error: cannot recognize FB = {}'.format(FB), file=sys.stderr)
		sys.exit(-1)
	return FB
	
FB = regularize_FB(FB)
should_annotate_prob = str2bool(should_annotate_prob)


def read_codon_file(codon_filename):
	import re
	one_codon_pattern = re.compile('([ACGTU]{3})\s+([\w*])\s+([0-9.]+)\s+([0-9.]+)\s+\(\s*(\d+)\)')
	
	triplet2aa = {}
	# triplet2aa['TTT'] = 'F'
	aa2triplet = {}
	# aa2triplet['F'] = ['TTT', 'TTC']
	triplet2freq = {}
	# triplet2freq['TTT'] = 17.2/1000

	with open(codon_filename, 'r') as fin:
		for line in fin:
			if line[0] == '#':   # comment
				continue
			line = line.rstrip()
			for one_match in one_codon_pattern.finditer(line):
				# [triplet] [amino acid] [fraction] [frequency: per thousand] ([number])
				triplet, aa, fraction, frequency_per1000, number = one_match.groups()
#				print('{} {}'.format(triplet, aa))   # for debug
				
				triplet = triplet.replace('U', 'T')
				triplet2aa[triplet] = aa
				if aa not in aa2triplet:
					aa2triplet[aa] = []
				aa2triplet[aa].append(triplet)
				triplet2freq[triplet] = float(frequency_per1000) / 1000
	
	triplet_posterior_given_aa = {}
	# triplet_posterior_given_aa['F']['TTT'] = 17.2 / (17.2 + 21.8) = 0.441
	triplet2posterior = {}
	# triplet_posterior['TTT'] = 17.2 / (17.2 + 21.8) = 0.441
	for aa in aa2triplet:
		triplet_posterior_given_aa[aa] = {}
		now_freq_sum = 0
		for triplet in aa2triplet[aa]:
			now_freq_sum += triplet2freq[triplet]
		for triplet in aa2triplet[aa]:
			triplet_posterior_given_aa[aa][triplet] = triplet2freq[triplet] / now_freq_sum
			triplet2posterior[triplet] = triplet2freq[triplet] / now_freq_sum
	
	return triplet2aa, aa2triplet, triplet2freq, triplet_posterior_given_aa, triplet2posterior
	
triplet2aa, aa2triplet, triplet2freq, triplet_posterior_given_aa, triplet2posterior = read_codon_file(codon_filename)
#print(triplet2aa['TAG'])   # for debug
#print(triplet_posterior_given_aa['*'])   # for debug
#print(triplet2posterior['TAG'])   # for debug

complementNt = {
'A' : 'T',
'C' : 'G',
'G' : 'C',
'T' : 'A',
'N' : 'N'
}
def get_reverseComplement(input_nt_seq):
	ans = []
	for x in reversed(input_nt_seq):
		if x in complementNt:
			ans.append(complementNt[x])
		else:
			print('Warning: cannot recognize nt = {} in get_reverseComplement()'.format(x), file=sys.stderr)
			ans.append(x)
	return ''.join(ans)

def translate(input_seq, should_annotate_prob, frame_shift=+1, fout=None):
	# frame_shift should be +1~+3 or -1~-3, 0 will be convert to +1
	if frame_shift == 0:
		frame_shift = 1
	elif frame_shift < 0:
		input_seq = get_reverseComplement(input_seq)
		frame_shift = -frame_shift
	input_seq = input_seq.replace('U', 'T')
	
	ans = []
	input_seq_len = len(input_seq)
	for i in range(frame_shift-1, input_seq_len, 3):
		now_triplet = input_seq[i:i+3]
		if now_triplet in triplet2aa:
			ans.append(triplet2aa[now_triplet])
		else:
			print('Warning: cannot recognize triplet = {} in translate()'.format(now_triplet), file=sys.stderr)
			ans.append('?')
	ans = ''.join(ans)
	
	if should_annotate_prob:
		tmp_list = backtranslate(ans, True)
		max_posterior = -1
		input_posteior = -1
		for x in tmp_list:
			seq, posterior = x.split('\t')
			posterior = float(posterior)
			if seq == input_seq:
				input_posteior = posterior
			if posterior > max_posterior:
				max_posterior = posterior
		ans += '\t' + str(input_posteior / max_posterior)
	
	if fout is None:
		return ans
	else:
		print(ans, file=fout)
		return
	
def backtranslate(input_seq, should_annotate_prob, fout=None):
	triplet_list = []
	for aa in input_seq:
		if aa in aa2triplet:
			triplet_list.append(aa2triplet[aa])
		else:
			print('Warning: cannot recognize aa = {} in backtranslate()'.format(aa), file=sys.stderr)
			triplet_list.append(['???'])
	
	ans = []
	import itertools
	for one_choice in itertools.product(*triplet_list):
		nt_seq = ''.join(one_choice)
		if should_annotate_prob:
			posterior = 1
			for one_triplet in one_choice:
				if one_triplet in triplet2posterior:
					posterior *= triplet2posterior[one_triplet]
			tmp = nt_seq + '\t' + str(posterior)
			if fout is None:
				ans.append(tmp)
			else:
				print(tmp, file=fout)
		else:
			if fout is None:
				ans.append(nt_seq)
			else:
				print(nt_seq, file=fout)
	if fout is None:
		return ans
	else:
		return



with open(input_filename, 'r') as fin:
	for line in fin:
		line = line.rstrip()
		if FB == 'F':
			translate(line, should_annotate_prob, fout=sys.stdout)
		elif FB == 'B':
			backtranslate(line, should_annotate_prob, fout=sys.stdout)
	