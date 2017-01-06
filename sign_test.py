# Performs sign test using a .gmt file and an ase file (containing gene names and ase values)

import scipy.stats as stats
import sys
CUTOFF = 150

ase, categories, outfile = sys.argv[1:]

top_genes = []
top_values = []
bottom_genes = []
bottom_values = []

with open(ase, 'r') as ase:
	ase.readline()

	counter = 0
	for line in ase:
		line = line.strip().split('\t')
		gene = line[0]
		ase_value = line[-1]
		if ase_value == "NA":
			continue
		else:
			ase_value = float(ase_value)

		if counter < CUTOFF:
			top_genes.append(gene) 
			top_values.append(ase_value)
			bottom_genes.append(gene) 
			bottom_values.append(ase_value)

		else:
			if ase_value > min(top_values):
				del top_genes[top_values.index(min(top_values))]
				top_genes.append(gene)
				top_values.remove(min(top_values))
				top_values.append(ase_value)
				
			elif ase_value < max(bottom_values):
				del bottom_genes[bottom_values.index(max(bottom_values))]
				bottom_genes.append(gene)
				bottom_values.remove(max(bottom_values))
				bottom_values.append(ase_value)	

		counter += 1

out = open(outfile, 'w')

with open(categories, 'r') as categories:
	for line in categories:
		line = line.strip().split('\t')
		category = line[0]
		genes = line[2:]

		bottom_yes = 0
		bottom_no = 0
		top_yes = 0
		top_no = 0

		for i in top_genes:
			if i in genes:
				top_yes += 1
			else:
				top_no += 1
		for i in bottom_genes:
			if i in genes:
				bottom_yes += 1
			else:
				bottom_no += 1

		oddsratio, pvalue = stats.fisher_exact([[bottom_yes, top_yes], [bottom_no, top_no]])
		outlist = [category, str(oddsratio), str(pvalue)]
		output = '\t'.join(outlist)
		out.write(output + '\n')

out.close()	
