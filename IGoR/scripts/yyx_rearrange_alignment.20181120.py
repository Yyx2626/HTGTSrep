
#### Usage: cat <_alignments.csv> | python3 this.py <idx_shift> <idx_start> <idx_end> [should_headline]
#### Options:
####   [should_headline] (default: True = true, T, yes, y, >0)
#### Input : STDIN (_alignments.csv, separated by ';', automatically judging headline
#### Output: STDOUT


import sys


### reference: https://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python
def str2bool(x):
	try:
		x = int(x)
		return x > 0
	except ValueError:
		return x.lower() in ('true', 't', 'yes', 'y')

def is_int(x):
	try:
		x = int(x)
		return True
	except ValueError:
		return False


idx_shift , idx_start, idx_end = sys.argv[1:4]
idx_shift , idx_start, idx_end = map(int, [idx_shift , idx_start, idx_end])
#print((idx_shift , idx_start, idx_end))
should_headline = True
if len(sys.argv) > 4:
	should_headline = str2bool(sys.argv[4])

early_stop_ahead = 100000   # hard-coded early stopping

if False:   ## all store then output code
	## first read in and store in memory
	storage = {}
	is_headline = True
	for line in sys.stdin:
		line = line.rstrip()
		F = line.split(';')
		if is_headline:
			is_headline = False
			## automatically judge input headline
			if not is_int(F[0]):   # input has headline
				if should_headline:
					print(line)
				continue
			else:   # input no headline
				if should_headline:
					print('Warning: seems no headline in input but request headline in output, so I will output default _alignments.csv headline.', file=sys.stderr)
					print('seq_index;gene_name;score;offset;insertions;deletions;mismatches;length;5_p_align_offset;3_p_align_offset')
		if not is_int(F[0]):
			# strange line format, skip it
			print('Warning: idx cannot be converted to integer, so I will skip it. line={}'.format(line), file=sys.stderr)
			continue
		idx = int(F[0])
		idx += idx_shift;
	#	print(idx)
		if idx < idx_start:
			continue
		if idx > idx_end + early_stop_ahead:   # hard-coded early stopping
			break
		if idx > idx_end:
			continue
	
		F[0] = str(idx)
		line = ';'.join(F)
		if str(idx) not in storage:
			storage[str(idx)] = []
		storage[str(idx)].append(line)

	now_output_idx = idx_start
	for idx in sorted(storage.keys(), key=int):
		idx = int(idx)
		if idx < idx_start:
			continue
		if idx > idx_end:
			break
		if idx < now_output_idx:
			print('Warning: idx={} < now_output_idx={}'.format(idx, now_output_idx), file=sys.stderr)
			continue
		while now_output_idx <= idx_end and now_output_idx < idx:
			now_output_idx += 1
			if now_output_idx % 10000 == 0:
				print('Now {}'.format(now_output_idx), file=sys.stderr)
			if idx == now_output_idx:
				print('\n'.join(storage[str(now_output_idx)]))
				already[str(now_output_idx)] = 1
				break

if True:   ## original code
	storage = {}
	already = {}
	now_output_idx = idx_start
	is_headline = True
	for line in sys.stdin:
		if now_output_idx > idx_end:  break
		line = line.rstrip()
		F = line.split(';')
		if is_headline:
			is_headline = False
			## automatically judge input headline
			if not is_int(F[0]):   # input has headline
				if should_headline:
					print(line)
				continue
			else:   # input no headline
				if should_headline:
					print('Warning: seems no headline in input but request headline in output, so I will output default _alignments.csv headline.', file=sys.stderr)
					print('seq_index;gene_name;score;offset;insertions;deletions;mismatches;length;5_p_align_offset;3_p_align_offset')
		if not is_int(F[0]):
			# strange line format, skip it
			print('Warning: idx cannot be converted to integer, so I will skip it. line={}'.format(line), file=sys.stderr)
			continue
		idx = int(F[0])
		idx += idx_shift;
	#	print(idx)
		if idx < idx_start:
			continue
		if idx > idx_end + early_stop_ahead:   # hard-coded early stopping
			break
		if idx > idx_end:
			continue
	
		F[0] = str(idx)
		line = ';'.join(F)
		if idx == now_output_idx:
			print(line)
			already[str(idx)] = 1
		elif idx < now_output_idx:
			print('Error: idx={} < now_output_idx={}'.format(idx, now_output_idx), file=sys.stderr)
			exit(-1)
		else:   # idx >= now_output_idx + 1
			if str(now_output_idx) in already:
				while str(now_output_idx + 1) in storage:
					now_output_idx += 1;
					if now_output_idx > idx_end:  break
					if now_output_idx % 10000 == 0:
						print('Now {}'.format(now_output_idx), file=sys.stderr)
					print('\n'.join(storage[str(now_output_idx)]))
					del storage[str(now_output_idx)]
					already[str(now_output_idx)] = 1
			if str(now_output_idx) in already and idx == now_output_idx + 1:
				now_output_idx += 1
				if now_output_idx > idx_end:  break
				if now_output_idx % 10000 == 0:
					print('Now {}'.format(now_output_idx), file=sys.stderr)
				print(line)
				already[str(idx)] = 1
			else:
				if idx <= now_output_idx:
					print('Error: idx={} <= now_output_idx={}'.format(idx, now_output_idx), file=sys.stderr)
					exit(-2)
				else:   # idx > now_output_idx + 1
					if str(idx) not in storage:
						storage[str(idx)] = []
				storage[str(idx)].append(line)
	if len(storage.keys()) > 0:
		for idx in sorted(storage.keys(), key=int):
			idx = int(idx)
			if idx > idx_end:
				break
			if idx < now_output_idx:
				print('Warning: idx={} < now_output_idx={}'.format(idx, now_output_idx), file=sys.stderr)
				continue
			while now_output_idx <= idx_end and now_output_idx < idx:
				now_output_idx += 1
				if now_output_idx % 10000 == 0:
					print('Now {}'.format(now_output_idx), file=sys.stderr)
				if idx == now_output_idx:
					print('\n'.join(storage[str(now_output_idx)]))
					already[str(now_output_idx)] = 1
					break
