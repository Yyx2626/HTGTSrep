#!/usr/bin/env python3
usage = 'Usage: cat <input> | python3 this.py <show|onlyshow|skip|retrieve> [column_pattern_1] [column_pattern_2] ...'

import sys, re


def match_any_patterns(patterns, one_query_string):
	for now_pattern in patterns:
		if now_pattern.search(one_query_string):
			return True
	return False

def which(bool_vec):
	ans = []
	for i in range(len(bool_vec)):
		if bool_vec[i]:
			ans.append(i)
	return ans


empty_line_pattern = re.compile('^\s*$')

if len(sys.argv) < 2:
	print(usage, file=sys.stderr)
	sys.exit(-1)

mode = sys.argv[1]
if mode not in ['show', 'skip', 'retrieve', 'onlyshow']:
	print('Error: unrecognized mode="{}", which should be either "show" or "skip" or "retrieve"'.format(mode))
	sys.exit(-2)

column_pattern_strings = sys.argv[2:]
pattern_len = len(column_pattern_strings)
column_patterns = [re.compile(x) for x in column_pattern_strings]


is_headline = True
NC = -1
column_match_idxes = []
column_not_match_idxes = []
for line in sys.stdin:
	if empty_line_pattern.search(line):
		continue
	line = line.rstrip()
	F = line.split('\t')
	if is_headline:
		is_headline = False
		NC = len(F)
		column_match_bool_vec = [match_any_patterns(column_patterns, x) for x in F]
		column_not_match_bool_vec = [not x for x in column_match_bool_vec]
		column_match_idxes = which(column_match_bool_vec)
		column_not_match_idxes = which(column_not_match_bool_vec)
		
		if mode == 'show' or mode == 'onlyshow':
			if pattern_len > 0:
				print('# {} matched columns:'.format(len(column_match_idxes)))
				for j in column_match_idxes:
					print('{}\t{}'.format(j+1, F[j]))
				if mode == 'show':
					print('# {} unmatched columns:'.format(len(column_not_match_idxes)))
					for j in column_not_match_idxes:
						print('{}\t{}'.format(j+1, F[j]))
			else:
				print('# {} columns:'.format(NC))
				for j in range(NC):
					print('{}\t{}'.format(j+1, F[j]))
			break
	if mode == 'skip':
		print('\t'.join([F[j] for j in column_not_match_idxes]))
	if mode == 'retrieve':
		print('\t'.join([F[j] for j in column_match_idxes]))

