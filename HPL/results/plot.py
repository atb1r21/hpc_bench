import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

node_lst, perf_lst = [],[]

with open(sys.argv[1],'r') as infile:
	for line in infile:
		if not line.startswith("Number"):
			line = line.strip()
			if len(line) > 0:
				columns  = line.split()
				node_lst.append(int(columns[0]))
				perf_lst.append(float(columns[1]))

sns.set_style('darkgrid')



plt.plot(node_lst, perf_lst)
plt.xlabel("Number of nodes")
plt.ylabel("Performance (Gflops)")
plt.title("HPL scaling benchmarks (more Gflops is faster)")
plt.tight_layout(pad=5)
plt.savefig(str(sys.argv[2]))

