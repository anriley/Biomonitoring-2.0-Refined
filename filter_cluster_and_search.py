from sys import argv
import os
from Bio import SeqIO,AlignIO
import argparse
import Python_subscripts.filter_initial_data as fid
import Python_subscripts.cluster_sequences as cs
import Python_subscripts.find_taxonomy_clusters as ftc

def open_file(filename): #Open and read in files
    with open(filename, "r") as file:
        temp = file.readlines()
    in_data = [line.strip() for line in temp] #Remove newline from end of lines
    return in_data

def prepare_data_csv(data): #Make csv files list ( all lines) of lists (lines split by commas)
    prepared_data = []
    for line in data:
        prepared_data.append(line.split(",")) #For csv files split on comma
    return prepared_data

def extract_guide_information(in_file):
    for line in in_file:
        temp = line.split(":")
        if temp[0].lower() == "amplicons":
            amplicons = temp[1].split(",")
        elif temp[0].lower() == "filters":
            filters = temp[1].split(",")
            filters = [int(i) for i in filters]
        elif temp[0].lower() == "taxonomy":
            taxonomy = temp[1].split(",")
        else:
            print("Error in guide file")
            exit()
    return amplicons,filters, taxonomy

def create_directory(directory):
    if not os.path.exists(directory): #Check if directory already exist
        os.makedirs(directory)        #Create directory

def prepare_directory(directory):
    if directory[-1] == "/":
        dir_prepared = directory
        create_directory(dir_prepared)
    else:
        dir_prepared = directory + "/"
        create_directory(dir_prepared)
    return dir_prepared

def write_csv_file(outname,out_data):
    fo = open(outname,"w")
    for line in out_data:
        temp = [str(i) for i in line]
        fo.write(",".join(temp))
        fo.write("\n")
    fo.close()

#Setup run
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d","--directory", required=True,help="Use the same directory that was used as \
                    the output directory of prepare_data_for_SWARM_and_alignment_filter.py")
    ap.add_argument("-g","--guide_file", required=False,help="File that includes information on \
                    amplicons to be used, absolute filters for read count, and taxonomy assignment")
    args = ap.parse_args()
    
    in_dir = args.directory
    guide_file = open_file(args.guide_file)
    amplicon_list,filter_list,taxonomy_search = extract_guide_information(guide_file)
    taxon_level = taxonomy_search[0]
    bp_list = taxonomy_search[1:]
    #Filter data
    for amp in amplicon_list:
        esv_table_file = prepare_data_csv(open_file(in_dir + "ESV_tables/ESV_" + amp + ".csv"))
        taxonomy_file = prepare_data_csv(open_file(in_dir + "Taxonomy/taxonomy_" + amp + ".csv"))
        sequence_file = AlignIO.read(in_dir + "Sequences/nt_" + amp + ".fasta", 'fasta')
        for f in filter_list:
            esv_filter_filename = "ESV_rcf_" + str(f) + "_" + amp + ".csv"
            taxonomy_filter_filename = "taxonomy_rcf_" + str(f) + "_" + amp + ".csv"
            sequence_filter_filename = "nt_rcf_" + str(f) + "_" + amp + ".fasta"
            if esv_filter_filename not in os.listdir(in_dir + "ESV_tables/") or taxonomy_filter_filename not in os.listdir(in_dir + "Taxonomy/") \
                or sequence_filter_filename not in os.listdir(in_dir + "Sequences/"):
                filtered_esv_table_file,filtered_taxonomy_file,filtered_sequence_file = fid.filter_data(esv_table_file,taxonomy_file,sequence_file,f)
                write_csv_file(in_dir + "ESV_tables/" + esv_filter_filename,filtered_esv_table_file)
                write_csv_file(in_dir + "Taxonomy/" + taxonomy_filter_filename,filtered_taxonomy_file)
                SeqIO.write(filtered_sequence_file, in_dir + "Sequences/" + sequence_filter_filename, "fasta")

    #Cluster data and check nestedness
    nested_log = ""        
    for amp in amplicon_list:
        for f in filter_list:
            cluster_dir = in_dir + "Clusters/Clusters_rcf_" + str(f) + "_" +amp
            cluster_dir = prepare_directory(cluster_dir)
            cluster_filename_prefix = "clusters_rcf_" + str(f) + "_d_"
            cluster_filename_suffix = "_" + amp + ".txt"
            sequences_to_cluster = in_dir + "Sequences/nt_rcf_" + str(f) + "_" + amp + ".fasta"
            cs.get_all_possible_clusters(cluster_dir,cluster_filename_prefix,cluster_filename_suffix,sequences_to_cluster)
            nested_log += amp + " rcf "+ str(f) + "\n"
            nested_log += cs.check_cluster_nestedness(cluster_dir,cluster_filename_prefix,cluster_filename_suffix)
    fo = open(in_dir + "Clusters/nestedness_log.txt","w")
    fo.write(nested_log)
    fo.close()
            
    #Taxonomy search        
    taxonomy_cluster_dir = in_dir + "Taxonomy_clusters/"
    taxonomy_cluster_dir = prepare_directory(taxonomy_cluster_dir)
    for amp in amplicon_list:
        for f in filter_list:
            in_taxonomy = prepare_data_csv(open_file(in_dir + "Taxonomy/" +"taxonomy_rcf_" + str(f) + "_" + amp + ".csv"))
            cluster_dir = in_dir + "Clusters/Clusters_rcf_" + str(f) + "_" +amp + "/"
            cluster_filename_prefix = "clusters_rcf_" + str(f) + "_d_"
            cluster_filename_suffix = "_" + amp + ".txt"
            for bp in bp_list:
                taxa_isolated,taxa_split_summary,taxa_split = ftc.taxonomy_cluster_search(in_taxonomy, taxon_level, bp, cluster_dir,cluster_filename_prefix,cluster_filename_suffix)
                out_name = taxonomy_cluster_dir+"taxonomy_clusters_rcf_"+str(f)+ "_" + taxon_level + "_" + bp + "_"+ amp + "_isolated.csv"
                write_csv_file(out_name,taxa_isolated)

                out_name_split_summary = taxonomy_cluster_dir+"taxonomy_clusters_rcf_"+str(f)+ "_" + taxon_level + "_" + bp + "_"+ amp + "_split_summary.csv"
                write_csv_file(out_name_split_summary,taxa_split_summary)

                out_name_split = taxonomy_cluster_dir+"taxonomy_clusters_rcf_"+str(f)+ "_" + taxon_level + "_" + bp + "_"+ amp + "_split.csv"
                write_csv_file(out_name_split,taxa_split)

if __name__ == "__main__":
    main()