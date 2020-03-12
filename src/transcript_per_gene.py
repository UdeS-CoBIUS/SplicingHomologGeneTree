import glob

dict = {}
for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/initialSource/*initialsource2target.txt"): 
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
        else:
            tmp_dict[gene_id] = 1
    for k, v in tmp_dict.items():
        if v in dict.keys():
            dict[v] += 1
        else:
            dict[v] = 1                    

nb_gene = 0
for k, v in dict.items():
    nb_gene += v           
X = []
Y = []    
for i in sorted (dict.keys()) :  
    print(i, dict[i], dict[i]*100/nb_gene)       
    X.append(i)
    Y.append(dict[i]*100/nb_gene)
print(X)    
print(Y)
         

