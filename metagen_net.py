# Script to generate networks from metagenomics data
# Criteria: Nodes = Function from a certain sample; Edges = functions with similar abundance where the dozens are the same (i.e 0.92 and 0.97)
# Example: node1 (Function A from Sample X) and node2 (Function B from Sample Y), where both nodes have abundance of 0.57 and 0.52. Therefore, they are connected.
# Goal of network: Show functions with similar abundance across different samples
# Input: transcript abundance spreadsheet

import sys

raw = sys.argv[1] # Excel table tab delimited with functions, samples and abundance
fh = open(raw)

# Output files for GePhi
nodes = open('Nodes_Metagenome.csv','w')
edges = open('Edges_Metagenome.csv','w')

nodes.write('Id\tLabel\tDescription\tSample\tSample_equal\tSize\n') # Label = function; Size = abundance
edges.write('Source\tTarget\tType\n')

links = {} # Dictionary with nodes and their size to make the edges file. key = node id; value = size

fh.readline()
ids = 0
for x in fh:
	if 'E+' in x: continue
	if 'E-' in x: continue
	y = x.strip().split('\t')
	label = y[3].upper() # Function
	description = y[1].upper()
	avg = (float(y[-1])+float(y[-2])+float(y[-3])+float(y[-4])+float(y[-5])+float(y[-6])+float(y[-7])+float(y[-8])+float(y[-9]))/9.0
	if avg < 0.2: continue
	size1b = float(y[-9])#*100
	ids += 1
	out1b = '\t'.join([str(ids),label,description,'KMC_1b1','KMC_1',str(size1b)])
	nodes.write(out1b+'\n')
	links[str(ids)] = size1b
	
	size1c = float(y[-8])#*100
	ids += 1
	out1c = '\t'.join([str(ids),label,description,'KMC_1c1','KMC_1',str(size1c)])
	nodes.write(out1c+'\n')
	links[str(ids)] = size1c
	
	size2b = float(y[-7])#*100
	ids += 1
	out2b = '\t'.join([str(ids),label,description,'KMC_2b1','KMC_2',str(size2b)])
	nodes.write(out2b+'\n')
	links[str(ids)] = size2b
	
	size2c = float(y[-6])#*100
	ids += 1
	out2c = '\t'.join([str(ids),label,description,'KMC_2c1','KMC_2',str(size2c)])
	nodes.write(out2c+'\n')
	links[str(ids)] = size2c
	
	size3b = float(y[-5])#*100
	ids += 1
	out3b = '\t'.join([str(ids),label,description,'KMC_3b1','KMC_3',str(size3b)])
	nodes.write(out3b+'\n')
	links[str(ids)] = size3b
	
	size3c = float(y[-4])#*100
	ids += 1
	out3c = '\t'.join([str(ids),label,description,'KMC_3c1','KMC_3',str(size3c)])
	nodes.write(out3c+'\n')
	links[str(ids)] = size3c
	
	size4b = float(y[-3])#*100
	ids += 1
	out4b = '\t'.join([str(ids),label,description,'KMC_4b1','KMC_4',str(size4b)])
	nodes.write(out4b+'\n')
	links[str(ids)] = size4b
	
	size4c = float(y[-2])#*100
	ids += 1
	out4c = '\t'.join([str(ids),label,description,'KMC_4c1','KMC_4',str(size4c)])
	nodes.write(out4c+'\n')
	links[str(ids)] = size4c
	
	size5 = float(y[-1])#*100
	ids += 1
	out5 = '\t'.join([str(ids),label,description,'KMC_5_1','KMC_5',str(size5)])
	nodes.write(out5+'\n')
	links[str(ids)] = size5
nodes.close()

for i in links:
	source = i
	value = str(links[i])
	for i2 in links:
		target = i2
		if source == target: continue
		value2 = str(links[i2])
		if value[2] == value2[2]:
			edges.write('\t'.join([i,i2,'Undirected'])+'\n')
edges.close()