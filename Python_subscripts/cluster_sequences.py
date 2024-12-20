import subprocess
import os
from pathlib import Path

def open_file(filename): #Open and read in files
    with open(filename, "r") as file:
        temp = file.readlines()
    in_data = [line[:-1] for line in temp] #Remove newline from end of lines
    return in_data

def prepare_data_ssv(data):
    prepared_data = []
    for line in data:
        prepared_data.append(line.split(" "))
    return prepared_data

def get_all_possible_clusters(directory,cluster_prefix,cluster_suffix,fasta):
    d = 0
    clusters = ["empty","empty"]
    while(len(clusters) > 1):
        cluster_files = os.listdir(directory)
        cluster_file = cluster_prefix +str(d) +cluster_suffix
        outfile = directory + cluster_prefix +str(d) +cluster_suffix
        if cluster_file not in cluster_files:
            subprocess.call("swarm -t 8 %s -o %s -d %s -n -l log_dump.txt" % (fasta,outfile,d), shell=True)  #Get absolute path
        clusters = prepare_data_ssv(open_file(outfile))
        d += 1

def check_cluster_nestedness(dir,cluster_prefix,cluster_suffix):
    out_log = ""
    flag = True
    d = 1
    while flag:
        my_file1 = dir+cluster_prefix+str(d)+cluster_suffix
        my_file2 = dir+cluster_prefix+str(d+1)+cluster_suffix  
        my_path1 = Path(my_file1)
        my_path2 = Path(my_file2)
        if my_path2.is_file():
            clusters1 = prepare_data_ssv(open_file(my_file1))
            clusters2 = prepare_data_ssv(open_file(my_file2))
            for c1 in clusters1:
                nest_flag = False
                for c2 in clusters2:
                    if set(c1) <= set(c2):
                        nest_flag = True
                        break
                if not nest_flag:
                    out_log += "cluster not nested at" + str(d) +"\n"
                    out_log += " ".join(c1) + "\n" 
        else:
            out_log += "All files run to: "+ str(d) +"\n\n"
            flag=False
        d += 1
    return out_log