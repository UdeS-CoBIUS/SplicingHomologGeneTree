import glob
import os
dict = {}
for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/initialSource/*initialsource2target.txt"): 
    tmp = x.split("/")
    tmp = tmp[-1]  
    tmp = tmp.split("_")
    id = tmp[0]
    keep = False
    if os.path.exists("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/clusters_0_100/fuzzyCMeans/" + id + "_cluster_1_root.nhx"):
        file = open(x)
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
        if keep:
            nb_gene_tmp = 0
            for k, v in tmp_dict.items():
                if v == 1:
                    nb_gene_tmp += 1
            fraction_gene_with_one_transcript = nb_gene_tmp*100.0/len(tmp_dict.keys())

            if fraction_gene_with_one_transcript in dict.keys():
                dict[fraction_gene_with_one_transcript] += 1
            else:
                dict[fraction_gene_with_one_transcript] = 1                    

for i in sorted (dict.keys()) :  
    print(i, dict[i])       
         

