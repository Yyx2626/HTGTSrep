
#### Usage: cat <model_parms.txt> | python3 this.py
#### Input: STDIN
#### Output: STDOUT

import sys

def sort_and_output(items):
	if len(items) <= 0:
		return
	hash = {}
	for line in items:
		F = line.split(';')
		idx = F[len(F)-1]
#		idx = int(idx.rstrip())
		hash[idx] = line
	for idx in sorted(hash.keys(), key=int):
		print(hash[idx])

should_sort = False
items = []
for line in sys.stdin:
	line = line.rstrip()
	if line[0] == '@':
		if line == '@Event_list':
			should_sort = True
		else:
#			sort_and_output(items)
#			items = []
			should_sort = False
	if should_sort:
		if line[0] == '%':
			items.append(line)
		else:
			sort_and_output(items)
			items = []
			print(line)
	else:
		sort_and_output(items)
		items = []
		print(line)
sort_and_output(items)
items = []

