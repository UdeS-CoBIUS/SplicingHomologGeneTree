import glob
import os
dict = {}
nb_gene_one_elt = 0
nb_total = 0
for x in glob.glob("clusters_0_100/fuzzyCMeans/*_cluster_1_root.nhx"): 
    tmp = x.split("/")
    tmp = tmp[-1]  
    tmp = tmp.split("_")
    id = tmp[0]
    keep = False
    initial2target = "initialSource/"+ id +"_initialsource2target.txt"
    file = open(initial2target, "r")
    lines = file.readlines()
    tmp_dict = {}
    for line in lines:
        line = line.replace("\n", "")
        parts = line.split(" ")
        gene_id = parts[1]
        cds_id = parts[0]        
        if gene_id in tmp_dict.keys():
            tmp_dict[gene_id] += 1
            keep = True
        else:
            tmp_dict[gene_id] = 1

        #used = False
        #for k, v in tmp_dict.items():
        #    if v != 1:
        #        used = True
        #if not(used):
        #    print(tmp_dict)
    nb_total += 1            
    if keep == False and len(lines) > 2:
        nb_gene_one_elt += 1

print(str(nb_total) + " gene familes \n" + str(nb_gene_one_elt) + " gene familes with genes having more than one transcript")
