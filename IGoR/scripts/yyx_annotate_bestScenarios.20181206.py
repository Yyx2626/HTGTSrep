
#### Usage: python3 this.py
####    <reference_dir> <indexed_sequences.csv> <parms.txt> <marginals.txt> <best_scenarios_counts.csv> [Pgen_counts.csv] [should_sort_idx_output]
#### Options: [should_sort_idx_output] (default: True)
#### Input: <best_scenarios_counts.csv>  15 columns
####  Note: <reference_dir> is actually not used so far
#### Output: STDOUT  .csv (';'-separated, idx-sorted)
####            add 3+4+8+8+6(+1)=29(30) columns
####              3(16-18th) V_name,D_name,J_name
####              4(19-22th) P(V),P(J|V),P(D|V,J),P(V,D,J)
####              2*4=8(23-30th) V3_delLen, P(V3_delLen|V), D5_delLen, P(.|D), D3_delLen, P(,|D), J5_delLen, P(.|J)
####              4*2=8(31-38th) VD_insLen, P(VD_insLen), VD_insNt, P(VD_insNt), DJ_...
####              1(39th) P_recomb = P(V,D,J)*P(V3_delLen|V)*...*P(VD_insLen)*P(VD_InsNt)*...
####              1(40th) P_mismatch = error_rate ^ num_mismatch
####              1(41th) P_scenario = P_recomb * P_mismatch
####              1(42th) P_posterior = P_scenario / P_read
####              1(43th) P_read = sum P_scenario
####              1(44th) P_gen  = sum P_recomb
####              1(45th) Pgen (if [Pgen_counts.csv] is provided)
####           each best_scenario followed by 6 lines for sequence alignment: read, match, VDJjoin, V, D, J

import sys
import numpy as np

ref_dirPath, idx_seq_filename, parms_filename, marginals_filename, bestScenarios_filename = sys.argv[1:6]
Pgen_filename = None
if len(sys.argv) > 6:
	Pgen_filename = sys.argv[6]

def str2bool(x):
	try:
		x = int(x)
		return x > 0
	except ValueError:
		return x.lower() in ('true', 't', 'yes', 'y')

should_sort_idx_output = True
if len(sys.argv) > 7:
	should_sort_idx_output = str2bool(sys.argv[7])
	

# for debug
print('ref_dirPath = {}'.format(ref_dirPath), file=sys.stderr)
print('idx_seq_filename = {}'.format(idx_seq_filename), file=sys.stderr)
print('parms_filename = {}'.format(parms_filename), file=sys.stderr)
print('marginals_filename = {}'.format(marginals_filename), file=sys.stderr)
print('bestScenarios_filename = {}'.format(bestScenarios_filename), file=sys.stderr)
print('Pgen_filename = {}'.format(Pgen_filename), file=sys.stderr)



def str_left(str, Len, fill=''):
	if Len > len(str):
		return str + fill * (Len - len(str))
	else:
		return str[0:Len]
	
def str_right(str, Len, fill=''):
	if Len > len(str):
		return fill * (Len - len(str)) + str
	else:
		return str[(len(str)-Len):]

def str_trim_left(str, Len):
	return str[Len:]

def str_trim_right(str, Len):
	if Len > len(str):
		return ''
	else:
		return str[0:(len(str)-Len)]

def str_left_match(str, target):
	return str_left(str, len(target)) == target

def str_right_match(str, target):
	return str_right(str, len(target)) == target


def read_parms_Event_list(parms_filename):
	parms_event = {}
	# parms_event[nick] = {key_idx:name, ..., name:int(idx)}
	VDJ_alleles = {'V':{}, 'D':{}, 'J':{}}
	# VDJ_alleles['V'/'D'/'J'] = {key_idx:seq, name:seq}
	error_rate = -1
	
	with open(parms_filename, 'r') as fin:
		part = ''
		subpart = ''
		type, gene, side, prio, nick = [''] * 5
		for line in fin:
			line = line.rstrip()
			if line[0] == '@':
				part = line[1:]
				if part != 'Event_list':
					pass   # do nothing
#					break   # end of reading
			elif line[0] == '#':
				if part == 'Event_list':
					subpart = line[1:]
					type, gene, side, prio, nick = subpart.split(';')
					if nick not in parms_event:
						parms_event[nick] = {}
				elif part == 'ErrorRate':
					subpart = line[1:]
			elif line[0] == '%':
				if part == 'Event_list':
					F = line[1:].split(';')
					name = F[0]
					idx = F[len(F)-1]
					parms_event[nick]['key_' + idx] = name
					parms_event[nick][name] = int(idx)
				if nick == 'v_choice' or nick == 'j_choice' or nick == 'd_gene':
					VDJ = nick[0].upper()
					seq = F[1]
					VDJ_alleles[VDJ]['key_' + idx] = seq
					VDJ_alleles[VDJ][name] = seq
			else:
				if part == 'ErrorRate' and subpart == 'SingleErrorRate':
					error_rate = float(line)
					subpart = 'SingleErrorRate_alreadyRead'
					
	return parms_event, VDJ_alleles, error_rate
				
def read_marginals(marginals_filename):
	marginals = {}
	# marginals[nick] = np.array()

	with open(marginals_filename, 'r') as fin:
		for line in fin:
			line = line.rstrip()
			if line[0] == '@':
				nick = line[1:]
			elif line[0] == '$':
				dim = line[5:(len(line)-1)]
				dim = list(map(int, dim.split(',')))
				marginals[nick] = np.zeros(dim)
			elif line[0] == '#':
				namedDim = line[2:(len(line)-1)]
				if namedDim == '':
					headIdxes = []
				else:
					namedDim = namedDim.split('],[')
					headIdxes = [int(x.split(',')[1]) for x in namedDim]
			elif line[0] == '%':
				vec = list(map(float, line[1:].split(',')))
				if len(headIdxes) == 0:
					marginals[nick] = np.array(vec)
					if str_right_match(nick, '_dinucl'):
						marginals[nick] = np.reshape(vec, (4,4))
				elif len(headIdxes) == 1:
					marginals[nick][headIdxes[0]] = vec
				elif len(headIdxes) == 2:
					marginals[nick][headIdxes[0]][headIdxes[1]] = vec
				else:
					print('Warning: dim > 3 ?!', file=sys.stderr)
	
	return marginals


VDJ2nick = {
	'V' : 'v_choice',
	'J' : 'j_choice',
	'D' : 'd_gene'
}
VDJdel2fullPrefix = {
	'V3' : 'Deletion_V_gene_Three_prime',
	'D5' : 'Deletion_D_gene_Five_prime',
	'D3' : 'Deletion_D_gene_Three_prime',
	'J5' : 'Deletion_J_gene_Five_prime'
}
ntCode2Nt = {
	'0' : 'A',
	'1' : 'C',
	'2' : 'G',
	'3' : 'T'
}
nt2intCode = {
	'A' : 0,
	'C' : 1,
	'G' : 2,
	'T' : 3
}
def getNtFromNtCode(ntCode):
	if ntCode in ntCode2Nt:
		return ntCode2Nt[ntCode]
	else:
		print('Warning: unknown ntCode={}, so I return ?'.format(ntCode), file=sys.stderr)
		return '?'

def get_VDJidx(VDJname, VDJ):
	VDJidx = -1
	if VDJname in parms_event[VDJ2nick[VDJ]]:
		VDJidx = parms_event[VDJ2nick[VDJ]][VDJname]
	else:
		print('Warning: cannot find {} {}'.format(VDJ, VDJname), file=sys.stderr)
	return VDJidx

def get_P_VDJchoice(Vidx, Didx, Jidx):
	ans = [-1] * 4
	# ans = [P(V), P(J|V), P(D|V,J), P(V,D,J)]
	if Vidx < marginals['v_choice'].shape[0]:
		ans[0] = marginals['v_choice'][Vidx]
		if Jidx < marginals['j_choice'][Vidx].shape[0]:
			ans[1] = marginals['j_choice'][Vidx][Jidx]
			if Didx < marginals['d_gene'][Vidx][Jidx].shape[0]:
				ans[2] = marginals['d_gene'][Vidx][Jidx][Didx]
				ans[3] = ans[0] * ans[1] * ans[2]
	return ans


def fields_to_field2idx(fields):
	ans = {}
	for i, x in enumerate(fields):
		ans[x] = i
	return ans

def fields_prefix_match_idx(query, fields):
	query_len = len(query)
	ans = []
	for i in range(len(fields)):
		if fields[i][0:query_len] == query:
			ans.append(i)
	return ans
	
def read_Pgen(Pgen_filename):
	ans = {}
	# ans[idx] = Pgen_estimate
	
	is_headline = True
	fields = []
	field2idx = {}
	with open(Pgen_filename, 'r') as fin:
		for line in fin:
			line = line.rstrip()
			F = line.split(';')
			if is_headline:
				is_headline = False
				fields = F
				field2idx = fields_to_field2idx(F)
			else:
				if 'seq_index' not in field2idx:
					print('Warning: I cannot find seq_idx column in Pgen file {}'.format(Pgen_filename), file=sys.stderr)
				if 'Pgen_estimate' not in field2idx:
					print('Warning: I cannot find Pgen_estimate column in Pgen file {}'.format(Pgen_filename), file=sys.stderr)
#				idx, Pgen = F[0:2]
				idx = F[field2idx['seq_index']]
				Pgen = F[field2idx['Pgen_estimate']]
				ans[idx] = Pgen

	return ans

def read_idx_seq(idx_seq_filename):
	input_idx2seq = {}
	is_headline = True
	fields = []
	field2idx = {}
	with open(idx_seq_filename, 'r') as fin:
		for line in fin:
			line = line.rstrip()
			F = line.split(';')
			if is_headline:
				is_headline = False
				fields = F
				field2idx = fields_to_field2idx(F)
			else:
				if 'seq_index' not in field2idx:
					print('Warning: I cannot find seq_idx column in idx_seq file {}'.format(idx_seq_filename), file=sys.stderr)
				if 'sequence' not in field2idx:
					print('Warning: I cannot find sequence column in idx_seq file {}'.format(idx_seq_filename), file=sys.stderr)
					
				idx = F[field2idx['seq_index']]
				seq = F[field2idx['sequence']]
				input_idx2seq['key_' + idx] = seq
	return input_idx2seq
	
	
def get_P_delLen(VDJidx, VDJdelLenIdx):
	VDJdelLen = {}
	P_VDJdel = {}
	for VDJdel in ('V3', 'D5', 'D3', 'J5'):
		nick = VDJdel[0].lower() + '_' + VDJdel[1] + '_del'
		nowLenIdx = VDJdelLenIdx[VDJdel]
		P_VDJdel[VDJdel] = marginals[nick][VDJidx[VDJdel[0]]][nowLenIdx]
		if VDJdel == 'D3':
			# P(D3|D,D5)
			P_VDJdel[VDJdel] = marginals[nick][VDJidx[VDJdel[0]]][VDJdelLenIdx['D5']][nowLenIdx]
			# I have checked and confirmed this extraction is correct
		VDJdelLen[VDJdel] = int(parms_event[nick]['key_' + str(nowLenIdx)])
#		## for simple convenience
#		if VDJdelLen[VDJdel] < 0:
#			VDJdelLen[VDJdel] = 0
	return P_VDJdel, VDJdelLen

def get_P_insLen(VDJinsLenIdx):
	VDJinsLen = {}
	P_VDJins = {}
	for VDJins in ('VD', 'DJ'):
		nick = VDJins.lower() + '_ins'
		nowLenIdx = VDJinsLenIdx[VDJins]
		P_VDJins[VDJins] = marginals[nick][nowLenIdx]
		VDJinsLen[VDJins] = int(parms_event[nick]['key_' + str(nowLenIdx)])
	return P_VDJins, VDJinsLen

def get_P_insDinucl(VDJinsNt, initNt={}):
	P_VDJinsNt = {}
	for VDJins in ('VD', 'DJ'):
		nick = VDJins.lower() + '_dinucl'
		ntVec = VDJinsNt[VDJins]
		ans = 1
		for i in range(len(ntVec)):
			prevNt = None
			if i > 0:
				prevNt = ntVec[i-1]
			elif VDJins in initNt:
				prevNt = initNt[VDJins]
			nextNt = ntVec[i]
#			print('{} {}'.format(prevNt, nextNt))   # for debug
			if prevNt is None or prevNt == '':
				print('Warning: no initNt, so I skip the first base transition (when calling get_P_insDinucl)', file=sys.stderr)
				continue   # skip the first nt if initNt == None (prevNt==None)
			if prevNt == '?' or nextNt == '?':
				print('Warning: prevNt == "?" or nextNt == "?", so I skip this base transition (when calling get_P_insDinucl)', file=sys.stderr)
				continue   # skip the first nt if prevNt == "?" or nextNt == "?"
			prevNtCode = nt2intCode[prevNt]
			nextNtCode = nt2intCode[nextNt]
			ans *= marginals[nick][prevNtCode, nextNtCode]
		P_VDJinsNt[VDJins] = ans
	return P_VDJinsNt

def get_match_mismatch_str(query, subjt, caseSensitive=False):
	if not caseSensitive:
		query = query.upper()
		subjt = subjt.upper()
		
	from itertools import zip_longest
	ans = []
	for a, b in zip_longest(query, subjt):
		if a is None:
			break
		if b is None:
			ans.append(a)
		elif a==b:
			ans.append('.')
		else:
			ans.append(a)
	return ''.join(ans)
	
def find_max_match_shift(query, subjt):
	max_match = -1
	max_match_shift = 0
	for shift in range(0, -len(query), -1):
		num_match = get_match_mismatch_str(query, ' '*(-shift) + subjt).count('.')
#		print('{} {}'.format(shift, num_match))   # for debug
		if num_match > max_match:
			max_match_shift = shift
			max_match = num_match
	for shift in range(1, len(subjt)):
		num_match = get_match_mismatch_str(' '*(shift) + query, subjt).count('.')
#		print('{} {}'.format(shift, num_match))   # for debug
		if num_match > max_match:
			max_match_shift = shift
			max_match = num_match
	return max_match_shift

def continuous_LCS_DP(query, subjt, caseSensitive=False):
	if not caseSensitive:
		query = query.upper()
		subjt = subjt.upper()
		
	len1 = len(query)
	len2 = len(subjt)
	max_LCS_len = 0
	max_len_i, max_len_j = -1, -1
	DPmat = np.zeros((len1+1, len2+1), dtype='int')
	for i in range(1, len1+1):
		for j in range(1, len2+1):
			if query[i-1] == subjt[j-1]:
				DPmat[i,j] = DPmat[i-1,j-1] + 1
				if DPmat[i,j] > max_LCS_len:
					max_LCS_len = DPmat[i,j]
					max_len_i, max_len_j = i, j
	
	return max_LCS_len, max_len_i, max_len_j
	# max_match_shift = j - i

def generate_scenario_alignment_str(input_read_seq, VDJidx, VDJdelLen, VDJinsNt, VDJinsLen={}, output_prefix=''):
	VDJseq = {}
	for VDJ in ('V', 'D', 'J'):
		VDJseq[VDJ] = VDJ_alleles[VDJ]['key_' + str(VDJidx[VDJ])]
	myVDJinsLen = {}
	for VDJins in ('VD', 'DJ'):
		if VDJins in VDJinsLen:
			if VDJinsLen[VDJins] != len(VDJinsNt[VDJins]):
				print('Warning: len(VDJinsNt["{}"]={} != VDJinsLen={}'.format(VDJins, len(VDJinsNt[VDJins]), VDJinsLen[VDJins]), file=sys.stderr)
		myVDJinsLen[VDJins] = len(VDJinsNt[VDJins])
	
	# del
#	Vseq = str_trim_right(VDJseq['V'], VDJdelLen['V3'])
#	Dseq = str_trim_right(str_trim_left(VDJseq['D'], VDJdelLen['D5']), VDJdelLen['D3'])
#	Jseq = str_trim_left(VDJseq['J'], VDJdelLen['J5'])
	if VDJdelLen['V3'] >= 0:
		Vseq = str_trim_right(VDJseq['V'], VDJdelLen['V3'])
	else:
		Vseq = VDJseq['V'] + '?' * - VDJdelLen['V3']
	if VDJdelLen['D5'] >= 0:
		Dseq = str_trim_left(VDJseq['D'], VDJdelLen['D5'])
	else:
		Dseq = '?' * - VDJdelLen['D5'] + VDJseq['D']
	if VDJdelLen['D3'] >= 0:
		Dseq = str_trim_right(Dseq, VDJdelLen['D3'])
	else:
		Dseq = Dseq + '?' * - VDJdelLen['D3'] 
	if VDJdelLen['J5'] >= 0:
		Jseq = str_trim_left(VDJseq['J'], VDJdelLen['J5'])
	else:
		Jseq = '?' * - VDJdelLen['J5'] + VDJseq['J']
	## 2018-12-04 Note: I found that IGoR DJins dinuc is probably from J to D (3'->5' direction, + strand)
	VDJjoinSeq = Vseq.upper() + VDJinsNt['VD'].lower() + Dseq.upper() + ''.join(reversed(VDJinsNt['DJ'])).lower() + Jseq.upper()
	VDJseqStart = {'V' : 0}
	VDJseqStart['D'] = VDJseqStart['V'] + len(Vseq) + myVDJinsLen['VD'] - VDJdelLen['D5']
	VDJseqStart['J'] = VDJseqStart['V'] + len(Vseq) + myVDJinsLen['VD'] + len(Dseq) + myVDJinsLen['DJ'] - VDJdelLen['J5']
	
#	myShift = find_max_match_shift(input_read_seq, VDJjoinSeq)
	tmp = continuous_LCS_DP(input_read_seq, VDJjoinSeq)
	myShift = tmp[2] - tmp[1]
	
	ans = output_prefix
	if myShift >= 0:
		ans += ' ' * myShift + input_read_seq
		ans += '\n' + output_prefix + get_match_mismatch_str(' ' * myShift + input_read_seq, VDJjoinSeq)
	else:
		ans += input_read_seq
		ans += '\n' + output_prefix + get_match_mismatch_str(input_read_seq, ' ' * (-myShift) + VDJjoinSeq)
		output_prefix += ' ' * (-myShift)
	ans += '\n' + output_prefix + VDJjoinSeq
	ans += '\n' + output_prefix + ' '*VDJseqStart['V'] + Vseq.upper() + str_right(VDJseq['V'], VDJdelLen['V3']).lower()
	ans += '\n' + output_prefix + ' '*VDJseqStart['D'] + str_left(VDJseq['D'], VDJdelLen['D5']).lower() + Dseq.upper() + str_right(VDJseq['D'], VDJdelLen['D3']).lower()
	ans += '\n' + output_prefix + ' '*VDJseqStart['J'] + str_left(VDJseq['J'], VDJdelLen['J5']).lower() + Jseq.upper()
	return ans, VDJjoinSeq, Vseq, Dseq, Jseq
	
	
def read_bestScenarios(bestScenarios_filename):
	contents = {}
	input_idx_order = []
	# contents[idx][rank] = F(list of the line)

	is_headline = True
	fields = []
	field2idx = {}
	VDJchoice_colidx = {'V':-1, 'D':-1, 'J':-1}
	VDJdel_colidx = {'V':-1, 'D':-1, 'J':-1}
	VDJins_colidx = {'V':-1, 'D':-1, 'J':-1}
	VDJinsNt_colidx = {'V':-1, 'D':-1, 'J':-1}
	seq_idx_colidx, scenario_rank_colidx, mismatch_colidx = -1, -1, -1
	with open(bestScenarios_filename, 'r') as fin:
		for line in fin:
			line = line.rstrip()
			F = line.split(';')
			if is_headline:
				is_headline = False
				fields = F
				field2idx = fields_to_field2idx(F)
				
				for VDJ in ('V','D','J'):
					VDJchoice_colidx[VDJ] = fields_prefix_match_idx('GeneChoice_'+VDJ+'_gene_', fields)
					if len(VDJchoice_colidx[VDJ]) <= 0:
						print('Error: cannot find GeneChoice_{}_gene_ column in bestScenarios input file {}'.format(VDJ, bestScenarios_filename), file=sys.stderr)
						return None
					if len(VDJchoice_colidx[VDJ]) > 1:
						print('Warning: multiple GeneChoice_{}_gene_ columns in bestScenarios input file {}. I just use the left-most column'.format(VDJ, bestScenarios_filename), file=sys.stderr)
					VDJchoice_colidx[VDJ] = VDJchoice_colidx[VDJ][0]
				for VDJdel in ('V3', 'D5', 'D3', 'J5'):
					VDJdel_colidx[VDJdel] = fields_prefix_match_idx(VDJdel2fullPrefix[VDJdel], fields)
					if len(VDJdel_colidx[VDJdel]) <= 0:
						print('Error: cannot find ' + VDJdel2fullPrefix[VDJdel] + ' column in bestScenarios input file {}'.format(VDJ, bestScenarios_filename), file=sys.stderr)
						return None
					if len(VDJdel_colidx[VDJdel]) > 1:
						print('Warning: multiple ' + VDJdel2fullPrefix[VDJdel] + ' columns in bestScenarios input file {}. I just use the left-most column'.format(VDJ, bestScenarios_filename), file=sys.stderr)
					VDJdel_colidx[VDJdel] = VDJdel_colidx[VDJdel][0]
				for VDJins in ('VD', 'DJ'):
					VDJins_colidx[VDJins] = fields_prefix_match_idx('Insertion_'+VDJins+'_gene', fields)
					if len(VDJins_colidx[VDJins]) <= 0:
						print('Error: cannot find Insertion_{}_gene column in bestScenarios input file {}'.format(VDJins, bestScenarios_filename), file=sys.stderr)
						return None
					if len(VDJins_colidx[VDJins]) > 1:
						print('Warning: multiple Insertion_{}_gene columns in bestScenarios input file {}. I just use the left-most column'.format(VDJins, bestScenarios_filename), file=sys.stderr)
					VDJins_colidx[VDJins] = VDJins_colidx[VDJins][0]
				
					VDJinsNt_colidx[VDJins] = fields_prefix_match_idx('DinucMarkov_'+VDJins+'_gene', fields)
					if len(VDJinsNt_colidx[VDJins]) <= 0:
						print('Error: cannot find DinucMarkov_{}_gene column in bestScenarios input file {}'.format(VDJins, bestScenarios_filename), file=sys.stderr)
						return None
					if len(VDJinsNt_colidx[VDJins]) > 1:
						print('Warning: multiple DinucMarkov_{}_gene columns in bestScenarios input file {}. I just use the left-most column'.format(VDJins, bestScenarios_filename), file=sys.stderr)
					VDJinsNt_colidx[VDJins] = VDJinsNt_colidx[VDJins][0]
				
				if 'seq_index' in field2idx:
					seq_idx_colidx = field2idx['seq_index']
				else:
					print('Error: cannot find Mismatches column in bestScenarios input file {}'.format(bestScenarios_filename), file=sys.stderr)
					return None
				if 'scenario_rank' in field2idx:
					scenario_rank_colidx = field2idx['scenario_rank']
				else:
					print('Error: cannot find Mismatches column in bestScenarios input file {}'.format(bestScenarios_filename), file=sys.stderr)
					return None
				if 'Mismatches' in field2idx:
					mismatch_colidx = field2idx['Mismatches']
				else:
					print('Warning: cannot find Mismatches column in bestScenarios input file {}'.format(bestScenarios_filename), file=sys.stderr)
			else:
#				idx, rank = F[0:2]
				idx = F[seq_idx_colidx]
				rank = F[scenario_rank_colidx]
				input_read_seq = ''
				if 'key_' + idx in input_idx2seq:
					input_read_seq = input_idx2seq['key_' + idx]
				VDJidx = {'V':-1, 'D':-1, 'J':-1}
				VDJname = {'V':'?','D':'?','J':'?'}
				VDJdelLenIdx = {'V3':-1, 'D5':-1, 'D3':-1, 'J5':-1}
				VDJinsLenIdx = {'VD':-1, 'DJ':-1}
				VDJinsNt = {'VD':'?', 'DJ':'?'}
				for VDJ in ('V','D','J'):
					VDJidx[VDJ] = F[VDJchoice_colidx[VDJ]]
					VDJidx[VDJ] = VDJidx[VDJ][1:(len(VDJidx[VDJ])-1)]   # remove outer ()
					if 'key_' + VDJidx[VDJ] in parms_event[VDJ2nick[VDJ]]:
						VDJname[VDJ] = parms_event[VDJ2nick[VDJ]]['key_' + VDJidx[VDJ]]
					VDJidx[VDJ] = int(VDJidx[VDJ])
				for VDJdel in ('V3', 'D5', 'D3', 'J5'):
					VDJdelLenIdx[VDJdel] = F[VDJdel_colidx[VDJdel]]
					VDJdelLenIdx[VDJdel] = int(VDJdelLenIdx[VDJdel][1:(len(VDJdelLenIdx[VDJdel])-1)])   # remove outer ()
				for VDJins in ('VD', 'DJ'):
					VDJinsLenIdx[VDJins] = F[VDJins_colidx[VDJins]]
					VDJinsLenIdx[VDJins] = int(VDJinsLenIdx[VDJins][1:(len(VDJinsLenIdx[VDJins])-1)])   # remove outer ()
					VDJinsNt[VDJins] = F[VDJinsNt_colidx[VDJins]]
					VDJinsNt[VDJins] = VDJinsNt[VDJins][1:(len(VDJinsNt[VDJins])-1)]   # remove outer ()
#					print(VDJins)   # for debug
#					print('"' + VDJinsNt[VDJins] + '"')
#					print(F)   # for debug
					if len(VDJinsNt[VDJins]) > 0:
						VDJinsNt[VDJins] = ''.join([ getNtFromNtCode(x) for x in VDJinsNt[VDJins].split(',') ])
					else:
						VDJinsNt[VDJins] = ''
				mismatches = []
				if mismatch_colidx >= 0:
					mismatches = F[mismatch_colidx]
					mismatches = mismatches[1:(len(mismatches)-1)]   # remove outer ()
					mismatches = mismatches.split(',')
				if len(mismatches) == 1 and mismatches[0] == '':
					mismatches = []
					
####            add 3+4+1(+1)=8(9) columns
####              3 V_name,D_name,J_name
				F.append(VDJname['V'])
				F.append(VDJname['D'])
				F.append(VDJname['J'])
####              4 P(V),P(J|V),P(D|V,J),P(V,D,J)
#				print(VDJidx)   # for debug
#				print(get_P_VDJchoice(VDJidx['V'], VDJidx['D'], VDJidx['J']))
				P_VDJchoice = get_P_VDJchoice(VDJidx['V'], VDJidx['D'], VDJidx['J'])
				F.extend(list(map(str, P_VDJchoice)))
				P_scenario = P_VDJchoice[3]
####              2*4=8 V3_delLen, P(V3_delLen|V), D5_delLen, P(.|D), D3_delLen, P(,|D, D5_delLen), J5_delLen, P(.|J)
				P_VDJdel, VDJdelLen = get_P_delLen(VDJidx, VDJdelLenIdx)
#				print(VDJidx)   # for debug
#				print(VDJdelLenIdx)   # for debug
#				print(VDJdelLen)   # for debug
#				print(P_VDJdel)   # for debug
				
				for VDJdel in ('V3', 'D5', 'D3', 'J5'):
					F.extend(list(map(str, (VDJdelLen[VDJdel], P_VDJdel[VDJdel]))))
					P_scenario *= P_VDJdel[VDJdel]
####              4*2=8 VD_insLen, P(VD_insLen), VD_insNt, P(VD_insNt), DJ_...
				P_VDJins, VDJinsLen = get_P_insLen(VDJinsLenIdx)
				alignStr, VDJjoinSeq, Vseq, Dseq, Jseq = generate_scenario_alignment_str(input_read_seq, VDJidx, VDJdelLen, VDJinsNt, VDJinsLen)
				initNt = {'VD':str_right(Vseq,1), 'DJ':str_right(Dseq,1)}
				if initNt['DJ'] == '':
					initNt['DJ'] = str_right(VDJinsNt['VD'], 1)
				P_VDJinsNt = get_P_insDinucl(VDJinsNt, initNt)
				for VDJins in ('VD', 'DJ'):
					F.extend(list(map(str, (VDJinsLen[VDJins], P_VDJins[VDJins], VDJinsNt[VDJins], P_VDJinsNt[VDJins]))))
					P_scenario *= P_VDJins[VDJins] * P_VDJins[VDJins] * P_VDJinsNt[VDJins]
####              1 P_mismatch = error_rate ^ num_mismatch
				P_mismatch = (error_rate/3) ** len(mismatches) * (1-error_rate) ** len(input_read_seq)
#				print('{} {} {}'.format(mismatches, len(mismatches), P_mismatch))   # for debug
				P_recomb = P_scenario
				P_scenario *= P_mismatch
				F.append(str(P_recomb))
				F.append(str(P_mismatch))
####              1 P_scenario = P(V,D,J)*P(V3_delLen|V)*...*P(VD_insLen)*P(VD_InsNt)*...
				F.append(str(P_scenario))
#				for i, x in enumerate(F):   # for debug
#					print('{} {}'.format(i, x))
#				break
				
				F.append(alignStr)   # push, wait for pop
				
####              1 P(V,D,J)_weighted_mean (weighted by scenario_proba_cond_seq)
####              1 P_scenario_weighted_mean (weighted by scenario_proba_cond_seq)
####              1 Pgen (if [Pgen_counts.tsv] is provided)

				if idx not in contents:
					contents[idx] = {}
					input_idx_order.append(idx)
				contents[idx][rank] = F

	output_fields = fields
	output_fields.extend(['V_name', 'D_name', 'J_name', 'P_V', 'P_J_given_V', 'P_D_given_V_J', 'P_VDJ'])
	for VDJdel in ('V3', 'D5', 'D3', 'J5'):
		output_fields.extend([VDJdel + '_delLen', 'P_' + VDJdel + '_delLen_given_' + VDJdel[0]])
		if VDJdel == 'D3':
			output_fields[len(output_fields)-1] = 'P_D3_delLen_given_D_and_D5_delLen'
	for VDJins in ('VD', 'DJ'):
		output_fields.extend([VDJins + '_insLen', 'P_' + VDJins + '_insLen', VDJins + 'insNt', 'P_' + VDJins + 'insNt'])
	output_fields.append('P_recomb')
	output_fields.append('P_mismatch')
	output_fields.append('P_scenario')
#	for i, x in enumerate(output_fields):   # for debug
#		print('{} {}'.format(i, x))
	return {'fields':output_fields, 'contents':contents, 'input_idx_order':input_idx_order}




print('Now read two model files ...', file=sys.stderr)
parms_event, VDJ_alleles, error_rate = read_parms_Event_list(parms_filename)
marginals = read_marginals(marginals_filename)
input_idx2seq = read_idx_seq(idx_seq_filename)

if Pgen_filename is not None:
	print('Now read Pgen file ...', file=sys.stderr)
	Pgen_hash = read_Pgen(Pgen_filename)

print('Now read best_scenariors file ...', file=sys.stderr)
bestScenario_hash = read_bestScenarios(bestScenarios_filename)

####              1(39th) P_recomb = P(V,D,J)*P(V3_delLen|V)*...*P(VD_insLen)*P(VD_InsNt)*...
####              1(40th) P_mismatch = error_rate ^ num_mismatch
####              1(41th) P_scenario = P_recomb * P_mismatch
####              1(42th) P_posterior = P_scenario / P_read
####              1(43th) P_read = sum P_scenario
####              1(44th) P_gen  = sum P_recomb
####              1(45th) Pgen (if [Pgen_counts.csv] is provided)
output_fields = bestScenario_hash['fields']
field2idx = fields_to_field2idx(output_fields)
output_fields.append('P_posterior')
output_fields.append('P_read')
output_fields.append('P_gen')
if Pgen_filename is not None:
	output_fields.append('Pgen_estimate')
print(';'.join(output_fields))
## sort, 2-pass for weighted mean, and output
if should_sort_idx_output:
	idxes = sorted(bestScenario_hash['contents'].keys(), key=int)
else:
	idxes = bestScenario_hash['input_idx_order']
for idx in idxes:
	this_seqIdx_Fs = bestScenario_hash['contents'][idx].values()
#	print(this_seqIdx_Fs)   # for debug
	num_sce = len(this_seqIdx_Fs)
#	sce_probs = [ F[field2idx['scenario_proba_cond_seq']]  for F in this_seqIdx_Fs ]
#	sce_probs = list(map(float, sce_probs))
#	sum_sce_prob = sum(sce_probs)
	
	P_scenarios = [ F[field2idx['P_scenario']]  for F in this_seqIdx_Fs ]
	P_scenarios = list(map(float, P_scenarios))
	P_read = sum(P_scenarios)
	P_recombs = [ F[field2idx['P_recomb']]  for F in this_seqIdx_Fs ]
	P_recombs = list(map(float, P_recombs))
#	P_posteriors = [x/P_read for x in P_scenarios]
	P_gen = sum(P_recombs)
	
	for rank in sorted(bestScenario_hash['contents'][idx].keys(), key=int):
		F = bestScenario_hash['contents'][idx][rank]
		alignStr = F[len(F)-1]
		F = F[0:(len(F)-1)]   # pop
		
		F.append(str(float(F[field2idx['P_scenario']]) / P_read))   # P_posterior
		F.append(str(P_read))
		F.append(str(P_gen))
		if Pgen_filename is not None:
			if idx in Pgen_hash:
				F.append(Pgen_hash[idx])
			else:
				F.append('-1')
#		print(F)   # for debug
		print(';'.join(F))
		print(alignStr)
