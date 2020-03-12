# SplicingHomologGeneTree
SplicingHomologGeneTree: Gene Tree Construction While Accounting for Splicing Structure Similarity

### Esaie Kuitche and Aïda Ouangraoua 

### Contact Esaie.Kuitche.Kamela@USherbrooke.ca

### Requirements The program requires the following to be available

* python (2.7) and python (3.5)
* ete3 
* argparse 
* biopython 
* numpy
* pandas 
* matplotlib
* scipy
* skbio

### Usage
```
usage: spliceGraphMaker.py [-h] [-s SAVE_PATH] [-ma MACROALIGNMENT]
                           [-wst WEIGHTSTRUCTURE] [-wsq WEIGHTSEQUENCE]
                           [-i INPUTFILE]


SplicingHomologGeneTree program parameters

optional arguments:
  -h, --help            show this help message and exit
  -s SAVE_PATH, --save_path SAVE_PATH (default = spliceGraph/)
  -ma MACROALIGNMENT, --macroalignment MACROALIGNMENT (default = macroalignment/)
  -wst WEIGHTSTRUCTURE, --weightStructure WEIGHTSTRUCTURE (default = 0.5)
  -wsq WEIGHTSEQUENCE, --weightSequence WEIGHTSEQUENCE (default = 0.5)
  -i INPUTFILE, --inputFile INPUTFILE (default = None)
```

## Input files and directories

```
macroalignment/: set of macroalignment files
initialSource/: set of initialsource fast files and set of initialsource2target files
microalignment/: set of microalignment files
ressources/ensemblSpecies.tree : Ensembl species tree
GENE_TREE_ID_mapping_transcript.txt : transcriptID, GeneID, and speciesID mapping
GENE_TREE_ID_mapping_species.txt : speciesID
```
## Ouptut directory and folder
* clusters_50_50/fuzzyCMeans/
* clusters_50_50/fuzzyCMeans/GENE_TREE_ID_cluster_1_root.nhx is the tree obtained using first our approach to select transcript sequences and then gene tree is build by using treebest
* clusters_50_50/fuzzyCMeans/GENE_TREE_ID_cluster_1_root.nhx is the tree obtained using first our approach to select transcript sequences and then gene tree is build by using neighbors joining

## Running SimSpliceEvol on an example
```
python2.7 spliceGraphMaker.py
```

