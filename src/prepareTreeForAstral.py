#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob

def geneToSpecies(clean_original_gene_tree):
	file = open(clean_original_gene_tree, "r")
	geneToSpeciesDict = {}
	lines = file.readlines()
	for line in lines:
		line = line.replace("\n", "")
		line = line.split(" ")
		geneToSpeciesDict[line[1]] = line[0]
	return geneToSpeciesDict


def confirmedNumberOfSelectedGeneTree():
	"""
	this function replace gene if by it species name in each gene tree.
	"""
	nb_t = 0
	listes = [[1,0]]#, [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]	
	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))	
		for root, dirs, files in os.walk("ensemblGeneTree/trees_" + sub_path): 		
				for file in files:					
					trees = open("ensemblGeneTree/trees_" + sub_path + "/" + file, "r")
					lines = trees.readlines()
					listestmp = [[1,0]]#, [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		
					#listes = [[0.5,0.5]]		

					skip = False				 
					for e in  listestmp:
						weightStructure = e[0]
						weightSequence = e[1]
						sub_pathtmp = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))				
						if os.path.exists("ensemblGeneTree/trees_" + sub_pathtmp + "/" + file) == False:
							skip = True

					if len(lines) != 4 or skip:
						pass
					else:			
	
						treeFromtreeBest = PhyloTree(lines[0])
						treeFromEnsembl = PhyloTree(lines[1])
						treeFromReconstructEnsembl = PhyloTree(lines[2])
						treeFromIsosel = PhyloTree(lines[3])		
						nb = len(treeFromtreeBest.get_leaf_names())
						if(len(treeFromEnsembl.get_leaf_names()) == nb  and len(treeFromReconstructEnsembl.get_leaf_names()) == nb and len(treeFromIsosel.get_leaf_names()) == nb):
							if(len(treeFromEnsembl.get_leaf_names()) == nb  and len(treeFromReconstructEnsembl.get_leaf_names()) == nb and len(treeFromIsosel.get_leaf_names()) == nb):

								clean_original_gene_tree = "ensemblGeneTree/cleanTree/" + file.split(".")[0] + "_cleanoutput.out"

								geneToSpeciesDict = geneToSpecies(clean_original_gene_tree)					
								
								for node in treeFromtreeBest.traverse("preorder"):
									if node.is_leaf():
										node.name = geneToSpeciesDict[node.name]

								for node in treeFromEnsembl.traverse("preorder"):
									if node.is_leaf():
										node.name = geneToSpeciesDict[node.name]							

								for node in treeFromReconstructEnsembl.traverse("preorder"):
									if node.is_leaf():
										node.name = geneToSpeciesDict[node.name]						

								for node in treeFromIsosel.traverse("preorder"):
									if node.is_leaf():
										node.name = geneToSpeciesDict[node.name]								

								if len(list(set(treeFromReconstructEnsembl.get_leaf_names()))) == len(treeFromReconstructEnsembl.get_leaf_names()):
									nb_t +=1
	
	return nb_t


listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		

percentages = [0.2, 0.4, 0.6, 0.8, 1.0]

nb_tree = confirmedNumberOfSelectedGeneTree()

for p in percentages:
	prefix = str(int(p*100)) + "_"	
	for e in  listes:		
		iteration = 0
		limit = nb_tree*p
	
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))

		astral_input_tree_orthogroup = open("Astral/test_data/" + prefix + sub_path + "_inputAstralOrthoGroup.tre", "w")
		astral_input_tree_Ensembl = open("Astral/test_data/" + prefix + "inputAstralEnsembl.tre", "w")
		astral_input_tree_ReconstructEnsembl = open("Astral/test_data/" + prefix +  "inputReconstructEnsembl.tre", "w")
		astral_input_tree_IsoSel = open("Astral/test_data/" + prefix + "inputIsoSel.tre", "w")


		for root, dirs, files in os.walk("ensemblGeneTree/trees_" + sub_path): 		
				for file in files:					
					trees = open("ensemblGeneTree/trees_" + sub_path + "/" + file, "r")
					lines = trees.readlines()
					listestmp = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		
					#listes = [[0.5,0.5]]		

					skip = False				 
					for e in  listestmp:
						weightStructure = e[0]
						weightSequence = e[1]
						sub_pathtmp = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))				
						if os.path.exists("ensemblGeneTree/trees_" + sub_pathtmp + "/" + file) == False:
							skip = True

					if len(lines) != 4 or skip or iteration>limit:
						pass
					else:								
						treeFromtreeBest = PhyloTree(lines[0])
						treeFromEnsembl = PhyloTree(lines[1])
						treeFromReconstructEnsembl = PhyloTree(lines[2])
						treeFromIsosel = PhyloTree(lines[3])		
						nb = len(treeFromtreeBest.get_leaf_names())
						if(nb>2 and len(treeFromEnsembl.get_leaf_names()) == nb  and len(treeFromReconstructEnsembl.get_leaf_names()) == nb and len(treeFromIsosel.get_leaf_names()) == nb):
							iteration +=1	
							clean_original_gene_tree = "ensemblGeneTree/cleanTree/" + file.split(".")[0] + "_cleanoutput.out"

							geneToSpeciesDict = geneToSpecies(clean_original_gene_tree)					
							
							for node in treeFromtreeBest.traverse("preorder"):
								if node.is_leaf():
									node.name = geneToSpeciesDict[node.name]

							for node in treeFromEnsembl.traverse("preorder"):
								if node.is_leaf():
									node.name = geneToSpeciesDict[node.name]							

							for node in treeFromReconstructEnsembl.traverse("preorder"):
								if node.is_leaf():
									node.name = geneToSpeciesDict[node.name]						

							for node in treeFromIsosel.traverse("preorder"):
								if node.is_leaf():
									node.name = geneToSpeciesDict[node.name]								

							if len(list(set(treeFromReconstructEnsembl.get_leaf_names()))) == len(treeFromReconstructEnsembl.get_leaf_names()):
								astral_input_tree_orthogroup.write(treeFromtreeBest.write(format=9) + "\n")
								astral_input_tree_Ensembl.write(treeFromEnsembl.write(format=9) + "\n")
								astral_input_tree_ReconstructEnsembl.write(treeFromReconstructEnsembl.write(format=9) + "\n")
								astral_input_tree_IsoSel.write(treeFromIsosel.write(format=9) + "\n")

		astral_input_tree_orthogroup.close()
		astral_input_tree_Ensembl.close()
		astral_input_tree_ReconstructEnsembl.close()
		astral_input_tree_IsoSel.close()

		os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + prefix + sub_path + "_inputAstralOrthoGroup.tre " + " -o Astral/test_data/"  + prefix + sub_path + "_SpeciesTreeOrthoGroup.tre")
		os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + prefix + "inputAstralEnsembl.tre " + " -o Astral/test_data/" + prefix + "SpeciesTreeEnsembl.tre")
		os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + prefix + "inputReconstructEnsembl.tre " + " -o Astral/test_data/"  + prefix + "SpeciesTreeReconstructEnsembl.tre")
		os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + prefix + "inputIsoSel.tre " + " -o Astral/test_data/"  + prefix + "SpeciesTreeIsoSel.tre")
