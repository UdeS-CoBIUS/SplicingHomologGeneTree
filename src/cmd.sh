#!/bin/bash

#compute SHT gene tree for alpha = [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]. 
#This script compute the similarity matrix (structure and sequence), applied SVD
#Applied Fuzzy C-Means to group transcripts. The results are store in the directory
#clusters_0_100/fuzzyCMeans, clusters_20_80/fuzzyCMeans, clusters_40_60/fuzzyCMeans, clusters_50_50/fuzzyCMeans,
#clusters_60_40/fuzzyCMeans, clusters_80_20/fuzzyCMeans, clusters_100_0/fuzzyCMeans
#The the SHT gene for eache gene family is call locate at: clusters_xx_yy/fuzzyCMeans/[GeneFamilyId]_cluster_1_root.nhx.
#This gene is given in NHX format, because it is output by TreeBest. The python2.7 prepareTreeForConselTest.py 
#will extract the gene and write in a newick format later.
python2.7 spliceGraphMaker.py

#Recompute Ensembl gene Tree based on the loguest transcript of each. The gene build are located
#at clusters_100_0/fuzzyCMeans/[GeneFamilyID]cluster_1_root.nhx.
#This gene is given in NHX format, because it is output by TreeBest. The python2.7 prepareTreeForConselTest.py 
#will extract the gene and write in a newick format later.
python2.7 recomputeEnsemblTree.py

#This script compute transcript isoform group by using IsoSel program. And It used the 
#return group to compute isosel gene trees based on their selected isoform. The gene trees build by this program
#are located at isoselSeq/[GeneFamilyID]_root.nhx
#This gene is given in NHX format, because it is output by TreeBest. The python2.7 prepareTreeForConselTest.py 
#will extract the gene and write in a newick format later.
python2.7 isosel.py

#This script extract all gene tree (SHT, LT, LT_db, IsoSel gene trees). Create a new located
#at ensemblGeneTree/trees_xx_yy/[GeneFamilyID].nw. The first line of this file est SHT tree
#Second LT tree, third LT_db tree and fourth is IsoSel gene tree in Newick format. 
#Then, it launches consel test for all genetree and all 4 methods SHT, LT, LT_db, IsoSel
#Then, it computes the reconciliation cost of this gene trees.
python2.7 prepareTreeForConselTest.py

#Compute all species trees (LT LT_db, IsoSel, SHT [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]) 
#using astral by varying the percentage of gene used [20%, 40%, 60%, 80%, 100%]
#The results species tree are located at Astral/test_data
python2.7 prepareTreeForAstral.py

#Compute AU test and Reconiliation plot based on the produced trees
python2.7 AU_score_and_Reconciliation_score.py

#Compute Boxplot of ratios of the mean similarity between the SHT transcripts over the mean similarity between the LT transcripts.
python2.7 ratio_SHT_isosel.py

#Compute Boxplot of ratios of the mean similarity between the SHT transcripts over the mean similarity between the IsoSel transcripts.
python2.7 ratio_SHT_LT.py

#Ratio value classes Ratio percentage
python2.7 ratio_value_classes.py

#Number of CDS per gene
python2.7 transcript_per_gene.py

#Percentage of gene families per percentage of genes having a single CDS
python2.7 one_transcript_per_gene.py

#Percentage of gene families per percentage of genes having a single CDS (restricted dataset)
python2.7 one_transcript_per_gene_128.py

#Number of CDS per gene (restricted dataset)
python2.7 transcript_per_gene_restricted_128.py

#Average common transcripts SHT IsoSel (Restricted dataset)
python2.7 python2.7 shares_leaves_isosel.py

#boxplot of common transcrit when considering onnly gene having more than one transcript	
python2.7 shares_leaves_isosel_gene_having_more_cds.py

