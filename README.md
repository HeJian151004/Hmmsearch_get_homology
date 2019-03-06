# Hmmsearch_get_homology
使用hmmsearch程序，基于预先准备好的hmmer库，搜索几个物种之间的直系同源基因。

```
import os           
import numpy as np
import re
global now_dir
now_dir = os.getcwd()        
species_list = ["Arabidopsis_thaliana","Citrus_sinensis","Malus_domestica","Vitis_vinifera","Cucumis_sativus","Glycine_max","Gossypium_raimondii","Populus_trichocarpa"]

def get_protein(IDs):
    with open(now_dir + "/" + IDs.split("|")[0] + ".fasta") as read_file:
        for each_line in read_file:
            if IDs in each_line:
                return(read_file.readline())
        

ortho_file = open(now_dir + "//myproject.proteinortho", "r")
n = 1
for each_line in ortho_file:
    if each_line:
        protein_name = "ortho" + str(n)
        n = n + 1
        with open(now_dir + "/" + protein_name + ".fasta","a") as write_file:
            for each_species in species_list:
                if each_line.count(each_species) > 4:
                    write_file.close()
                    os.remove(now_dir + "/" + protein_name + ".fasta")
                    break
                elif each_line.count(each_species) > 1:
                    pattern_protein = re.compile(each_species + "\|[a-zA-Z]+\_[0-9]+[\.]*[0-9]*")
                    pattern_protein_select1 = re.findall(pattern_protein,each_line)
                    temp_list = []
                    a = "a"
                    for each in pattern_protein_select1:
                        protein = get_protein(each)
                        if len(protein) > len(a):
                            a = protein
                            b = each
                    write_file.write(">" + protein_name + "|" + b + "\n" + a)
                else:
                    pattern_protein = re.compile(each_species + "\|[a-zA-Z]+\_[0-9]+[\.]*[0-9]*")
                    pattern_protein_select3 = re.findall(pattern_protein,each_line)[0]
                    pattern = get_protein(pattern_protein_select3)
                    write_file.write(">" + protein_name + "|" + pattern_protein_select3 + "\n" + pattern)
                                
ortho_file.close()        
```
