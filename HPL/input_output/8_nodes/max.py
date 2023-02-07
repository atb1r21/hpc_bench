import os,sys

raw_lst = []

with open(sys.argv[1], 'r') as infile:
	for line in infile:
		line = line.strip()
		columns = line.split()
		if '+0' in line:
			raw_lst.append(str(columns[6]))
processed = []

for x in raw_lst:
	values = x.split('e')
	base = float(values[0])
	if values[1] == '+02':
		base = base * 100
	elif values[1] == '+03':
		base = base * 1000
	else:
		base = base * 10000
	processed.append(base)

maximum = max(processed)

print('HIGHEST GFLOPS')

print(" ")

print(maximum)

