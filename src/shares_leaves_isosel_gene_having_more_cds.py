import os
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob
# libraries
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import random, string
from Bio.SeqRecord import SeqRecord

def splided_gene_family(id):	
	sourceToTarget = "initialSource/" + id + "_initialsource2target.txt"
	file = open(sourceToTarget, "r")
	dict_ = {}	
	lines = file.readlines()

	flag = False
	if len(lines)<3:
		return False
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

#share_leaves_Boxplot.pdf
def sharedLeaves():
    all_list = []
    
    listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]

    for e in  listes:
        weightStructure = e[0]
        weightSequence = e[1]
        sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
        
        percentageSharedLeaves = []

        for x in glob.glob("clusters_" + sub_path + "/clusters/*_completeClusters.txt"): 
            try:
                selected_transcripts = []
                tmp = x.split("/")[-1]
                tmp = tmp.split("_")[0]
                if not(splided_gene_family(tmp)):
                	continue
                intial_source2target = "initialSource/"+ tmp + "_initialsource2target.txt"
                flag = False
                file = open(intial_source2target, "r")
                lines = file.readlines()
                tmp_dict = {}
                genes = {}
                if len(lines)>3:
                    for line in lines:
                        line = line.replace("\n", "")
                        parts = line.split(" ")
                        gene_id = parts[1]
                        cds_id = parts[0]       

                        if gene_id in genes.keys():
                            genes[gene_id].append(cds_id)
                        else:
                            genes[gene_id] = [cds_id]
                                                
                        if gene_id in tmp_dict.keys():
                            tmp_dict[gene_id] += 1
                        else:
                            tmp_dict[gene_id] = 1

                    for k, v in tmp_dict.items():
                        if v > 1:
                            flag = True
                            for e in genes[k]:
                                selected_transcripts.append(e)
            
                    if flag == False:
                        pass
                    else:                           
                        fuzy = open(x, "r")
                        lines = fuzy.readlines()
                        if len(lines)>=2:
                            line = lines[1]
                            line = line.replace("\n", "")
                            seq_from_fuzzy = line.split("\t")
                            #print(seq_from_fuzzy)
                                                
                            seq_from_isosel = []
                            for record in SeqIO.parse("isoselSeq/" + tmp + "_filtered.fasta", "fasta"):
                                seq_from_isosel.append(record.id)                    

                            if len(seq_from_fuzzy) == len(seq_from_isosel):
                                shares = [e for e in seq_from_fuzzy if e in seq_from_isosel]
                                perecentage = 1.0*len(shares)/len(seq_from_fuzzy)
                                percentageSharedLeaves.append(perecentage)
                                #print(perecentage)
            except:
                pass
        all_list.append(percentageSharedLeaves)
                    

    c = "black"

    plt.boxplot(all_list,
                         vert=True,  # vertical box alignment
                         patch_artist=True,  # fill with color
                         capprops=dict(color=c),
                         whiskerprops=dict(color=c),
                         flierprops=dict(color=c, markeredgecolor=c),
                         medianprops=dict(color=c),
                         labels=[r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])  # will be used to label x-ticks
    plt.xticks(fontsize=70, rotation=50)
    plt.yticks(fontsize=70)
    plt.ylabel('Average common transcripts', fontsize=50)
    plt.legend(prop={'size': 60})    
    plt.show()

#computeStats()
sharedLeaves()
