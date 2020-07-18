 #!/usr/bin/python
import sys
import requests
import os
from os import getcwd
import argparse
from ete3 import Tree
import glob
import re
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
import pandas as pd
import random, string
from Bio.Align.Applications import MafftCommandline
#from StringIO import StringIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from copy import deepcopy, copy
import traceback
import numpy as np
from scipy.stats.stats import pearsonr  
import math
import matplotlib.pyplot as plt
import pandas as pd

def splided_gene_family(id):	
	sourceToTarget = "initialSource/" + id + "_initialsource2target.txt"
	file = open(sourceToTarget, "r")
	dict_ = {}	
	lines = file.readlines()
	
	flag = False
	for l in lines:
		l = l.replace("\n", "")
		p = l.split(" ")
		g_id = p[1]
		c_id = p[0]
		if g_id in dict_.keys():
			dict_[g_id] +=1
			flag = True
		else:
			dict_[g_id] =1

	return flag

def ratio(all=0):	
	SHT_similaritry = [[],[],[],[],[],[],[]]
	LT_similaritry = [[],[],[],[],[],[],[]]
	Ration_SHT_LT_similaritry = [[],[],[],[],[],[],[]]
	for x in glob.glob("clusters_100_0/fuzzyCMeans/*_cluster_1_root.nhx"): 		
	    try:
		    x2 =x
		    file = open(x, "r")
		    x = str(x)
		    x= x.split("/")
		    x = x[-1]		
		    x = x.split("_")[0]
		    if not(splided_gene_family(x)):
		    	continue
		    mapping_species_file = open("tmp/"+x+"_mapping_species.txt", "r")
		    mapping_genes_file = open("tmp/"+x+"_mapping_transcript.txt", "r")

		    mapping_transcript_gene_protein_from_emsenbl_file = open("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt", "r")

		    trees_build_by_orthoGroup_nw = x2
		    mapping_genes = {}
		    lines = mapping_genes_file.readlines()	
		    #print(lines)
		    for line in lines:
			    line = line.replace("\n", "")
			    line = line.split("\t")
			    mapping_genes[line[0]] = [line[1], line[2], line[3]]
			

		    file = open(trees_build_by_orthoGroup_nw, "r")
		    tree_str = file.read()
		    tree_str = tree_str.replace("\n", "")
		    tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
		    tree = Tree(tree_str)
		    leavesNames = []
		    SHT_CDS_dict = {}
		    SHT_CDS = []

		    for node in tree.traverse("postorder"):
			    if node.is_leaf():
				    name = node.name
				    SHT_CDS_dict[mapping_genes[name][1]] = mapping_genes[name][0]
				    name =  mapping_genes[name][0]
				    SHT_CDS.append(name)

		    #print("isoselSeq/" + x + "_root.nhx")
		    file = open("isoselSeq/" + x + "_root.nhx" , "r")
		    tree_str = file.read()
		    tree_str = tree_str.replace("\n", "")
		    tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
		    tree = Tree(tree_str)

		    leavesNames = []
		    LT_CDS = []
		    LT_CDS_dict = {}
		    genes_used_by_LT = []

		    #print(SHT_CDS)
		    #print(mapping_genes)
		    for node in tree.traverse("postorder"):
			    if node.is_leaf():
				    name = node.name
				    if name in mapping_genes.keys():

				        var = mapping_genes[name][1]
				        LT_CDS_dict[mapping_genes[name][1]] = mapping_genes[name][0]

				        name =  mapping_genes[name][0]

				        LT_CDS.append(name)

				        genes_used_by_LT.append(var)



		    #print(LT_CDS)
		    #print(SHT_CDS)
		    for k, v in SHT_CDS_dict.items():
			    if k not in genes_used_by_LT:
				    SHT_CDS.remove(v)

		    values = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1]
		    if all == 1 and len(LT_CDS)==len(SHT_CDS):
			    for alpha in values:
				    sim_LT, sim_SHT = mean_simalarity(x, LT_CDS, SHT_CDS, alpha)

				    SHT_similaritry[values.index(alpha)].append(sim_SHT)
				    LT_similaritry[values.index(alpha)].append(sim_LT)
				    Ration_SHT_LT_similaritry[values.index(alpha)].append(SHT_similaritry/LT_similaritry)
				    #print(x, len(LT_CDS), SHT_similaritry/LT_similaritry)

		    elif all ==0 and len(LT_CDS)==len(SHT_CDS):
			    if (set(SHT_CDS) != set(LT_CDS)):
				    for alpha in values:
					    sim_LT, sim_SHT = mean_simalarity(x, LT_CDS, SHT_CDS, alpha)

					    SHT_similaritry[values.index(alpha)].append(sim_SHT)
					    LT_similaritry[values.index(alpha)].append(sim_LT)
					    Ration_SHT_LT_similaritry[values.index(alpha)].append(sim_SHT/sim_LT)
	    except Exception as e:
	        #print(e)
	        pass
	#print(SHT_similaritry)
	#print(LT_similaritry)
	c = "black"
	plt.boxplot(Ration_SHT_LT_similaritry,
	                     vert=True,  # vertical box alignment
	                     #patch_artist=True,  # fill with color
	                      boxprops=dict(color=c),
                          capprops=dict(color=c),
                          whiskerprops=dict(color=c),
                          flierprops=dict(color=c, markeredgecolor=c),
                          medianprops=dict(color=c),
	                     labels=[r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])  # will be used to label x-ticks
	plt.xticks(fontsize=70, rotation=50)
	plt.yticks(fontsize=70)
	plt.ylim(0, 5)
	plt.ylabel('Ratio of mean similarity', fontsize=60)
	plt.legend(prop={'size': 60})	
	plt.show()

	
	
	mean_std_SHT = [[np.mean(x), np.std(x)] for x in SHT_similaritry]
	mean_std_LT = [[np.mean(x), np.std(x)] for x in LT_similaritry]
	mean_std_ratio = [[np.mean(x), np.std(x)] for x in Ration_SHT_LT_similaritry]
	print(mean_std_SHT)
	print(mean_std_LT)
	print(mean_std_ratio)

	mean_SHT = [item[0] for item in mean_std_SHT]
	std_SHT = [item[1] for item in mean_std_SHT]

	mean_LT = [item[0] for item in mean_std_LT]
	std_LT = [item[1] for item in mean_std_LT]	

	mean_ratio = [item[0] for item in mean_std_ratio]
	std_ratio = [item[1] for item in mean_std_ratio]	

	print(mean_SHT)
	print(std_SHT)
	print(mean_LT)
	print(std_LT)
	print(mean_ratio)	
	print(std_ratio)
def retrieveName(matrix):
	names = []
	for val in range(len(matrix)): 
		names.append(matrix[val][0].replace("\n",""))
	return names

'''
Reconstruct the matrix from the file, return a matrix with only the values and a dataframe with the header
'''
def retrieveValues(matrix, size):
	values = np.eye(size, size+1)
	for i, v in enumerate(matrix): 
		for j, o in enumerate(v): 
			if type(o) == str  : 
				continue 
			else:
				values[i][j] = o
	
	values = np.delete(values,0,1)
	return values


def mean_simalarity(x, LT_CDS, SHT_CDS, alpha):
	struc_file = "similarityScores/microalignment_"+x+"_score.csv"
	seq_file = "similaritySeqScores/"+x+"_microalignment.csv"
	matrixStruct = pd.read_csv(struc_file)		
	sizestruc = matrixStruct.shape[0]
	matrixStructValues = matrixStruct.values
	namesStruct = retrieveName(matrixStructValues)
	#namesStruct = [x.split("_")[0] for x in namesStruct]
	matrixStructValues = retrieveValues(matrixStructValues, sizestruc)	 


	matrixSeq = pd.read_csv(seq_file, sep="\t")
	sizeseq = matrixSeq.shape[0]
	matrixSeqValues = matrixSeq.values
	namesSeq = retrieveName(matrixSeqValues)
	matrixSeqValues = retrieveValues(matrixSeqValues, sizeseq) 
	LT_sim = 0.0
	SHT_sim = 0.0
	cmpt1 = 0
	#print(LT_CDS)
	#print(SHT_CDS)
	for i in range(len(namesStruct)):
		for j in range(i):
			#print(namesStruct[i].split("_")[0], namesStruct[j])
			if (namesStruct[i].split("_")[0] in LT_CDS) and  (namesStruct[j].split("_")[0] in LT_CDS):
				cmpt1 +=1
				stuct_sim = float(matrixStructValues[i,j])
				seq_sim = float(getValueSeqMatrix(i,j,matrixSeqValues, namesStruct, namesSeq))
				LT_sim += stuct_sim*alpha + (1-alpha)*seq_sim
	cmpt2 = 0
	#print("______________________________--")
	for i in range(len(namesStruct)):
		for j in range(i):
			#print(namesStruct[i], namesStruct[j])
			if (namesStruct[i].split("_")[0] in SHT_CDS) and  (namesStruct[j].split("_")[0] in SHT_CDS):
				cmpt2 +=1
				stuct_sim = float(matrixStructValues[i,j])
				seq_sim = float(getValueSeqMatrix(i,j,matrixSeqValues, namesStruct, namesSeq))
				SHT_sim += stuct_sim*alpha + (1-alpha)*seq_sim

	if cmpt1 ==0 or cmpt2 == 0:
		#exit()
		return -1, -1
	else:
		return LT_sim/cmpt1, SHT_sim/cmpt2

def getValueSeqMatrix(position1, position2, matrixSequence, names, header):
	namesonly = []

	for nameandgene in names:
		tmp = nameandgene.split("_")
		tmp_name = ""
		
		if len(tmp)==2:
			namesonly.append(tmp[0])
		else:
			del tmp[-1]
			tmp_name = "_".join(tmp)		
			namesonly.append(tmp_name)


	for i, position in enumerate(header):
		if position == namesonly[position1]:
			pos1 = i

	for j,pos in enumerate(header): 
		if pos == namesonly[position2]:
			pos2 = j


	return matrixSequence[pos1, pos2]


def correlation_matrix():
	correlations = []
	correlations2 = []
	for x in glob.glob("similaritySeqScores/*_microalignment.csv"): 
		seq_file =x
		file = open(x, "r")
		x = str(x)
		x= x.split("/")
		x = x[-1]		
		x = x.split("_")[0]
		print(x)
		struc_file = "similarityScores/microalignment_" + x + "_score.csv"
		if os.path.exists(struc_file):
			matrixStruct = pd.read_csv(struc_file)		
			sizestruc = matrixStruct.shape[0]
			matrixStructValues = matrixStruct.values
			namesStruct = retrieveName(matrixStructValues)
			#namesStruct = [x.split("_")[0] for x in namesStruct]
			matrixStructValues = retrieveValues(matrixStructValues, sizestruc)	 


			matrixSeq = pd.read_csv(seq_file, sep="\t")
			sizeseq = matrixSeq.shape[0]
			matrixSeqValues = matrixSeq.values
			namesSeq = retrieveName(matrixSeqValues)
			matrixSeqValues = retrieveValues(matrixSeqValues, sizeseq) 
			struc_values = []
			seq_values = []
			flag = False
			for e in namesStruct:
				if len(e.split("_")) > 2:
					flag = True
			if flag == False:				
				for i in range(len(namesStruct)):
					for j in range(i):
						stuct_sim = float(matrixStructValues[i,j])
						seq_sim = float(getValueSeqMatrix(i,j,matrixSeqValues, namesStruct, namesSeq))
						struc_values.append(stuct_sim)
						seq_values.append(seq_sim)
				correlation = pearsonr(struc_values, seq_values)
				
				correlation2 = np.corrcoef(struc_values,seq_values)

				if np.isnan(correlation[0]) or correlation[0]<0:
					#print(struc_values, seq_values)
					pass
				else:
					correlations.append(correlation[0])

				if np.isnan(correlation2[1][0]) or correlation2[1][0]<0.2:
					#print(struc_values, seq_values)
					pass
				else:
					correlations2.append(correlation2[1][0])

	
	#print(correlations)
	print(np.mean(correlations))
	print(np.median(correlations))
	
	#print(correlations2)
	print(np.mean(correlations2))
	print(np.median(correlations2))
	print(np.min(correlations2), np.max(correlations2))
	plt.xticks(fontsize=25)
	plt.yticks(fontsize=25)
	plt.legend(prop={'size': 25})

	df = pd.DataFrame({'Correlation':correlations2})
	df.boxplot()
	plt.show()



ratio()		

