
#### Usage: python3 this.py <name1> <file1> <name2> <file2> [...]
#### Input: <file1> <file2> ...   (.csv, separated by ';')
####        should have headline, with field_name including 'seq_index' (unsorted)
#### Output: STDOUT  (.tsv, separated by '\t')
####         add <name?>_ prefix to each field name, except key 'seq_index' (sorted)

import sys

argv = sys.argv[1:]
if len(argv) < 4:
	print('Error: no enough command-line arguments; please check the usage by head this.py', file=sys.stderr)
	sys.exit(-1)

storage = {}
# storage[seq_idx][name] = str

## read in and store in storage (hash)
names = []
name2fieldnames = {}
name2emptyStr = {}
while len(argv) >= 2:
	name, filename = argv[0:2]
	print('Now process '+name+' = '+filename, file=sys.stderr)
	argv = argv[2:]
	
	with open(filename, 'r') as fin:
		name2fieldnames[name] = []
		field2idx = {}
		is_headline = True
		for line in fin:
			line = line.rstrip()
			F = line.split(';')
			if is_headline:
				is_headline = False
				for i, x in enumerate(F):
					field2idx[x] = i
				if 'seq_index' not in field2idx:
					print('Warning: cannot find seq_index column in file ' + filename + ', so I skip it.', file=sys.stderr)
					break
				name2fieldnames[name] = [x for x in F if x != 'seq_index']
				name2emptyStr[name] = '\t'.join([''] * len(name2fieldnames[name]))
				names.append(name)
				continue
			idx = F[field2idx['seq_index']]
			G = []
			for fieldname in name2fieldnames[name]:
				if field2idx[fieldname] < len(F):
					G.append(F[field2idx[fieldname]])
				else:
					G.append('')
			remain = '\t'.join(G)
			if idx not in storage:
				storage[idx] = {}
			storage[idx][name] = remain

## sort idx and output
print('Now sort and output', file=sys.stderr)
output_fieldnames = ['seq_index']
output_fieldnames.extend([name+'_'+fieldname for name in names for fieldname in name2fieldnames[name]])
print('\t'.join(output_fieldnames))
for idx in sorted(storage.keys(), key=int):
	strs = [idx]
	for name in names:
		if name in storage[idx]:
			strs.append(storage[idx][name])
		else:
			strs.append(name2emptyStr[name])
	print('\t'.join(strs))
