# Hmmsearch_get_homology
使用hmmsearch程序，基于预先准备好的hmmer库，搜索几个物种之间的直系同源基因。  

Use the hmmsearch program to search for direct homologous genes between several species based on the prepared hmmer library.  

文件夹中需要有三种类型的文件，首先是.hmm为后缀的由hmmbuild程序生成的文件，之后是.cds后缀的核酸文件，最后是.pep后缀的氨基酸文件。  

You need to have three types of files in the folder, first the hmmbuild file with the.hmm suffix, then the nucleic acid file with the.cds suffix, and finally the amino acid file with the.pep suffix.  

该pyothon程序可以找出在核酸文件以及氨基酸文件中符合".hmm"结构域模式的序列。  
The pyothon program can find sequences in nucleic acid files and amino acid files that conform to the ".hmm" domain pattern. 



```
import os

now_dir = os.getcwd()            
file_temp = os.listdir(path=now_dir)    
pep_list = []         
cds_list = []
hmm_list = []

for each in file_temp:        
    if ".pep" in each:
        pep_list.append(each)
    elif ".cds" in each:
        cds_list.append(each)
    elif ".hmm" in each:
        hmm_list.append(each)

def get_seq(each_name,phy_sequence):
    with open(now_dir + "/" + each_name, "r") as read_file:
        seq_list = []
        for each_line in read_file:           
            if phy_sequence == each_line.split(" ")[0][1:]:
                a = "a"
                while a[0] != ">":
                    a = read_file.readline()
                    if a == "":
                        break
                    elif a[0] == ">":
                        break
                    else:
                        seq_list.append(a)
        seq = ''.join(seq_list)
        return seq
                                 
def get_ortho_aln(each_hmm):
    for each_pep in pep_list:
        search_output_name = each_hmm  + "_" + each_pep
        os.system("hmmsearch " + each_hmm + " " + each_pep + " > " + search_output_name)
        with open(now_dir + "/" + search_output_name, "r") as read_file:
            for each_line in read_file:
                if "E-value  score  bias    E-value  score  bias" in each_line:
                    read_file.readline()
                    need_read_line = read_file.readline()
                    need_read_line_list = []
                    for each_element in need_read_line.split(" "):
                        if each_element != "":
                            need_read_line_list.append(each_element)
                    try:
                        if float(need_read_line_list[4]) > 500:
                            phy_sequence = need_read_line_list[8]
                            with open(now_dir + "/fa_output/" + each_hmm[:-4] + "_pep.fa", "a") as write_file1:
                                write_file1.write(">" + each_pep[:-21] + "\n" + get_seq(each_pep,phy_sequence))
                            with open(now_dir + "/fa_output/" + each_hmm[:-4] + "_cds.fa", "a") as write_file2:
                                write_file2.write(">" + each_pep[:-21] + "\n" + get_seq((each_pep[:-3] + "cds"),phy_sequence))
                        else:
                            print("low_hmmsearch_score in " + each_pep + " to search " + each_hmm)
                            with open(now_dir + "/fa_output/log.txt", "a") as log_file:
                                log_file.write("low_hmmsearch_score in " + each_pep + " to search " + each_hmm)
                    except:
                        print("low_hmmsearch_score in " + each_pep + " to search " + each_hmm)
                        with open(now_dir + "/fa_output/log.txt", "a") as log_file:
                            log_file.write("low_hmmsearch_score in " + each_pep + " to search " + each_hmm)

for each_hmm in hmm_list:
    get_ortho_aln(each_hmm)
    print(each_hmm + " have been finished search")    
```
