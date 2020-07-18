#!/usr/bin/python3.7
#		 Computing statistics on species tree branch supports
#
#		 Da Rocha Coimbra, Nilson Antonio
#		 nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#		 Ouangraoua, Aida
#		 aida.ouangraoua@usherbrooke.ca
#
#		 2020
#

#from Utils import *
from ete3 import Tree
import decimal
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np

###############################################################################

# Directory for statistics files

###############################################################################

for pe in ["20_", "40_", "60_", "80_", "100_"]:
	files=[]
	sub_path = "100_0"
	species_tree_LT = "Astral/test_data/"+ pe + "SpeciesTreeEnsembl.tre"
	species_tree_LTBD = "Astral/test_data/"+ pe + "SpeciesTreeReconstructEnsembl.tre"
	species_tree_IsoSel = "Astral/test_data/"+ pe + "SpeciesTreeIsoSel.tre"

	files.append(species_tree_LT)
	files.append(species_tree_LTBD)
	files.append(species_tree_IsoSel)

	listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		
	support_values_per_method = []

	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))

		species_tree_SHT = "Astral/test_data/" + pe + sub_path + "_SpeciesTreeOrthoGroup.tre"
		
		files.append(species_tree_SHT)


	for f in files:
		support_values = []
		p = f.split("/")[-1]
		print(f)
		t_astral = Tree(f, quoted_node_names=False)
		outpufile = open("Astral/test_data/support_" + p, "w")
		for node in t_astral.traverse("postorder"):
			if( not node.is_leaf()):
				leaves = []
				for leaf in node:
					leaves.append(leaf.name)
				outpufile.write(str(node.support)+"\n")
				support_values.append(node.support)
		support_values_per_method.append(support_values)		
		outpufile.close()

	
	c = "black"
	all_stat = support_values_per_method
	plt.boxplot(all_stat,
					vert=True,  # vertical box alignment
					patch_artist=True,  # fill with color
					capprops=dict(color=c),
					whiskerprops=dict(color=c),
					flierprops=dict(color=c, markeredgecolor=c),
					medianprops=dict(color=c),
					labels=["LT", "LT_db", "IsoSel", "SHT "+r'$\alpha=0.0$', "SHT "+r'$\alpha=0.2$', "SHT "+r'$\alpha=0.4$', "SHT "+r'$\alpha=0.5$', "SHT "+r'$\alpha=0.6$', "SHT "+r'$\alpha=0.8$', "SHT "+r'$\alpha=1.0$'])

	plt.ylim(0, 1)
	#plt.xlim(0, len(all_stat)*2)
	plt.tight_layout()
	plt.xticks(fontsize=80, rotation=50)
	plt.yticks(fontsize=50)
	plt.ylabel('Quartet support values of branches', fontsize=50)
	plt.xlabel(pe.split("_")[0] + "% of gene trees used", fontsize=70)
	plt.legend(prop={'size': 60})	
	plt.show()



###############################################################################

# Generate statistics on common clades between trees

###############################################################################
print("\nGenerating statistics on common clades between trees...")
all_files = []
for pe in ["20_", "40_", "60_", "80_", "100_"]:
	files=[]
	sub_path = "100_0"
	species_tree_LT = "Astral/test_data/"+ pe + "SpeciesTreeEnsembl.tre"
	species_tree_LTBD = "Astral/test_data/"+ pe + "SpeciesTreeReconstructEnsembl.tre"
	species_tree_IsoSel = "Astral/test_data/"+ pe + "SpeciesTreeIsoSel.tre"
	
	files.append(species_tree_LT)
	files.append(species_tree_LTBD)
	files.append(species_tree_IsoSel)

	listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		

	support_values_per_method = []

	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))

		species_tree_SHT = "Astral/test_data/" + pe + sub_path + "_SpeciesTreeOrthoGroup.tre"
	
		files.append(species_tree_SHT)
	all_files.append(files)
print(all_files)


common_clades_stat_file = "Astral/test_data/common_clades.stats"

stat = open(common_clades_stat_file,"w")
line = ""

all_percentage = []
for files in all_files:
	local_per = []
	for i in range(len(files)):
		filei = files[i]	

		ti = Tree(filei)

		filej = "Astral/test_data/originalEnsemblSpeciesTree.nw"

		tj = Tree(filej)
		
		leaves_i = ti.get_leaf_names()
		leaves_j = tj.get_leaf_names()

		for node in ti.traverse("postorder"):
			if node.is_leaf():
				if not(node.name in leaves_j):					
					node.detach()

		for node in tj.traverse("postorder"):
			if node.is_leaf():
				if not(node.name in leaves_i):					
					node.detach()

		nb_ti = 0
		nb_tj = 0
		nb_common_clade = 0
		for node in ti.traverse("postorder"):
			if(not node.is_leaf()):
				nb_ti += 1
				leaves = []

				for leaf in node:
					name = leaf.name
					leaves.append(name.lower())

				if len(leaves) != 0:
					if(tj.check_monophyly(values=leaves, target_attr="name", ignore_missing=True)[0]):
						nb_common_clade +=1
		for node in tj.traverse("postorder"):
			if( not node.is_leaf()):
				nb_tj += 1
		percent_common_clade = 100.0*nb_common_clade/min(nb_ti,nb_tj)
		print(percent_common_clade)
		local_per.append(int(round(percent_common_clade,0)))
	all_percentage.append(local_per)


percentages_per_method = [[]]*10	
for i in range(len(percentages_per_method)):
	l = []
	for j in range(len(all_percentage)):
		l.append(all_percentage[j][i])
	percentages_per_method[i] = l

barWidth = 0.3

r1 = 2*np.arange(len(all_percentage[0]))
r2 = [x + barWidth for x in r1]
r3 = [x + 2*barWidth for x in r1]
r4 = [x + 3*barWidth for x in r1]
r5 = [x + 4*barWidth for x in r1]

rects1 = plt.bar(r1, all_percentage[0], width = barWidth, color = '#4285F4', edgecolor = '#4285F4', capsize=30, label='20%')

rects2 = plt.bar(r2, all_percentage[1], width = barWidth, color = '#EA4335', edgecolor = '#EA4335', capsize=30, label='40%')

rects3 = plt.bar(r3, all_percentage[2], width = barWidth, color = '#FBBC03', edgecolor = '#FBBC03', capsize=30, label='60%')

rects4 = plt.bar(r4, all_percentage[3], width = barWidth, color = '#6b392b', edgecolor = '#6b392b', capsize=30, label='80%')

rects5 = plt.bar(r5, all_percentage[4], width = barWidth, color = '#0b0812', edgecolor = '#6b392b', capsize=30, label='100%')


position = [2*r + 2*barWidth for r in range(len(all_percentage[0]))]


plt.xticks(position, ['LT', 'LT_db', 'IsoSel', r'$SHT \alpha = 0.0$', r'$SHT \alpha = 0.2$', r'$SHT \alpha = 0.4$', r'$SHT \alpha = 0.5$', r'$SHT \alpha = 0.6$', r'$SHT \alpha = 0.8$', r'$SHT \alpha = 1.0$'])

plt.xticks(fontsize=40, rotation=40)
plt.yticks(fontsize=60)
plt.legend(prop={'size': 35})
plt.legend(loc="best", borderaxespad=0., prop={'size': 35}, ncol=5)
#plt.xlabel('xlabel', fontsize=22)
plt.ylabel('Percentages of common clades', fontsize=45)
#plt.savefig("RI.png")
# Show graphic
def autolabel(rects):
	for rect in rects:
	    height = rect.get_height()        
	    plt.text(rect.get_x() + rect.get_width()/2., 0.99*height,
	             str(height),
	            ha='center', va='bottom',  fontsize=30)


autolabel(rects1)                
autolabel(rects2)                
autolabel(rects3)                
autolabel(rects4)
autolabel(rects5)

plt.show()

