import argparse
import os
import Python_subscripts.alignment_filtering as af#COME UP WITH A BETTER NAME

#Functions for basic tasks
def open_file(filename):
    with open(filename, "r") as file: #Read in file
        temp = file.readlines()
    in_data = [line[:-1] for line in temp] #Remove newline character from lines
    return in_data

def prepare_data_csv(data): #Split csv files into lists using commas
    prepared_data = []
    for line in data:
        prepared_data.append(line.split(","))
    return prepared_data

def prepare_data_tsv(data): #Split tsv files into lists using tabs
    prepared_data = []
    for line in data:
        prepared_data.append(line.split("\t"))
    return prepared_data

def create_directory(directory):
    if not os.path.exists(directory): #Check if directory already exist
        os.makedirs(directory)

def prepare_directory(directory):
    if directory[-1] == "/":
        dir_prepared = directory
        create_directory(dir_prepared)
    else:
        dir_prepared = directory + "/"
        create_directory(dir_prepared)
    return dir_prepared

def add_dict_list(d,key,value): #Append values to lists for keys already present in dictionaries
    if key in d:
        d[key].append(value)
    else:
        d[key] = [value]
    return

def create_id_sitename_dict(metadata):
    id_sitename_dict = {}
    for line in metadata:
        add_dict_list(id_sitename_dict,line[1],line[0])
    return id_sitename_dict

def primer_indices_in_ESV_table(header_ESV,check_index,sites=0):       #Find indices of columns in the ESV table associated with specific primers
    out_index_dict = {}
    for x in range(1, len(header_ESV)):
        temp = header_ESV[x].split("_")                                #Split site name
        sitename = "_".join(temp[:-2])                                 #Change site name to match site names in metadata
        if not sites:                                                  #If no metadata is provided, do not check if a site is present within a given range
            add_dict_list(out_index_dict,temp[check_index],x)
        elif sitename in sites:
            add_dict_list(out_index_dict,temp[check_index],x)          #Use primer name in site name to separate data by primer
    return out_index_dict

def assign_zotus_to_primer(ESV_table,in_index_dict):                  
    out_zotu_dict = {}
    for line in ESV_table:                                          #For each zotu
        flag = False
        for x in range(1,len(line)):
            if int(line[x]) > 0:                                    #Find the first column where there are more than zero reads
                for key in in_index_dict:
                    if x in in_index_dict[key]:                     #Use that column index and the primer index dict to assign the zotu to a primer 
                        add_dict_list(out_zotu_dict,key,line[0])
                        flag = True                                 #Skip checking more columns once a non-zero column is found
                        break
            if flag:
                break
    return out_zotu_dict

def pull_sequences_by_primer(sequences,in_zotu_dict):                                   #Using zotus assigned to primers
    out_sequence_dict = {}
    for x in range(0,len(sequences)):       
        if sequences[x].startswith(">"):                                                #Find lines in fasta that represent headers
            for key in in_zotu_dict:
                if sequences[x][1:] in in_zotu_dict[key] and "n" not in sequences[x+1].lower(): #Identify which sequences belong to each primer and assign them to a dictionary
                    add_dict_list(out_sequence_dict,key,[sequences[x][1:],sequences[x+1]])
    return out_sequence_dict

def pull_ESV_data_by_primer(in_site_dict,in_zotu_dict,ESV_table,in_zotu_reads):
    out_ESV_dict = {}
    for key in in_site_dict:
        primer_indices = in_site_dict[key]
        primer_zotus = in_zotu_dict[key]
        ESV_outheader = [ESV_table[0][0]] + ["_".join(ESV_table[0][x].split("_")[:-2]) for x in primer_indices]                       #Setup ESV table header
        ESV_out_table = [ESV_outheader]
        for line in ESV_table:
            if line[0] in primer_zotus:                                                 #Check if current zotu belongs to the current primer
                ESV_data = [int(line[x]) for x in primer_indices]                       #Only pull site specific primers
                ESV_outline = [line[0]] + ESV_data
                ESV_out_table.append(ESV_outline)
                in_zotu_reads[line[0]] = sum(ESV_data)
        out_ESV_dict[key] = ESV_out_table
    return out_ESV_dict

def pull_ESV_data_by_primer_after_msa(in_site_dict,in_zotu_dict):
    out_ESV_dict = {}
    for key in in_site_dict:
        out_esv_table = [in_site_dict[key][0]]
        for line in in_site_dict[key][1:]:
            if line[0] in in_zotu_dict[key]:
                out_esv_table.append(line)
        out_ESV_dict[key] = out_esv_table
    return out_ESV_dict

def pull_taxonomy_data_by_primer(in_taxonomy, in_zotu_dict):
    out_taxonomy = {}
    for key in in_zotu_dict:
        primer_zotus = in_zotu_dict[key]
        current_taxonomy = [in_taxonomy[0]]
        for line in in_taxonomy[1:]:
            if line[0] in primer_zotus:
                current_taxonomy.append(line)
        out_taxonomy[key] = current_taxonomy
    return out_taxonomy

def add_read_abundance_to_zotu(in_ESV_dict,in_sequence_dict,in_taxonomy_dict,in_read_dict):
    for key in in_sequence_dict:
        for line in in_sequence_dict[key]:
            line[0] = line[0] +"_"+str(in_read_dict[line[0]])
        for line in in_ESV_dict[key][1:]:
            line[0] = line[0] +"_"+str(in_read_dict[line[0]])
        for line in in_taxonomy_dict[key][1:]:
            line[0] = line[0] +"_"+str(in_read_dict[line[0]])
    return

def write_csv_file(outname,out_data):
    fo = open(outname,"w")
    for line in out_data:
        temp = [str(i) for i in line]
        fo.write(",".join(temp))
        fo.write("\n")
    fo.close()
    return

def write_fasta_file(outname,out_data):
    fo = open(outname,"w")
    for line in out_data:
        fo.write(">"+"\n".join(line))
        fo.write("\n")
    fo.close()
    return

def output_ESV_fasta_and_taxonomy(in_ESV_dict,in_sequence_dict,in_taxonomy_dict,out_dir):
    for key in in_ESV_dict:
        out_dir_ESV = prepare_directory(out_dir + "ESV_tables")
        out_ESV_filename = out_dir_ESV + "ESV_" + key + ".csv"
        write_csv_file(out_ESV_filename,in_ESV_dict[key])
        out_dir_sequences = prepare_directory(out_dir + "Sequences")
        out_fasta_filename = out_dir_sequences + "nt_" +key +".fasta"
        write_fasta_file(out_fasta_filename,in_sequence_dict[key])
        out_dir_taxonomy = prepare_directory(out_dir + "Taxonomy")
        out_taxonomy_filename = out_dir_taxonomy + "taxonomy_" +key +".csv"
        write_csv_file(out_taxonomy_filename,in_taxonomy_dict[key])
    return

def write_summary(sum_dict,out_dir):
    out_dir_MSA = prepare_directory(out_dir + "MSA")
    outname = out_dir_MSA + "summary_test.csv"
    fo = open(outname,"w")
    fo.write("Amplicon,Before_MSA,After_MSA,Most_gap,Most_no_gap\n")
    for key in sum_dict:
        temp = [str(i) for i in sum_dict[key]]
        fo.write(key+","+",".join(temp)+"\n")
    fo.close()
    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d","--directory", required=False,help="Assumes output is directly from Metaworks: ESV.table,\
                    chimera.denoised.nonchimera.taxon, and taxonomy.csv")
    ap.add_argument("-e","--esvtable", required=False,help="Tab separated ESV table file")
    ap.add_argument("-s","--sequences",required=False,help="Sequence file in FASTA format")
    ap.add_argument("-t","--taxonomy",required=False,help="Taxonomy file in csv format")
    ap.add_argument("-o","--out_directory",required=True,help="Directory for output")
    ap.add_argument("-a","--alignment",action="store_true", help="Create a MAFFT alignment and remove sequences based on the alignment")

    args = ap.parse_args()

    out_dir = prepare_directory(args.out_directory)

    if args.directory: #Check if Metaworks output directory was provided
        in_dir = prepare_directory(args.directory)
        ESV_table_infile = prepare_data_tsv(open_file(in_dir +"ESV.table"))
        sequence_infile = open_file(in_dir +"chimera.denoised.nonchimeras.taxon")
        taxonomy_infile = prepare_data_csv(open_file(in_dir +"taxonomy.csv"))
    else: #Get manually input files
        ESV_table_infile = prepare_data_tsv(open_file(args.esvtable))
        sequence_infile = open_file(args.sequences)
        taxonomy_infile = prepare_data_csv(open_file(args.taxonomy))

    primerID_index = -1
    zotu_reads = {} #Create dictionary to hold total reads for each zotu

    #Create function  variables = args.metadata,ESV_table_in_file[0],primerID_index
    primer_site_index_dict = primer_indices_in_ESV_table(ESV_table_infile[0],primerID_index)
    primer_zotu_dict = assign_zotus_to_primer(ESV_table_infile[1:],primer_site_index_dict)
    primer_sequence_dict = pull_sequences_by_primer(sequence_infile,primer_zotu_dict)
    primer_ESV_dict = pull_ESV_data_by_primer(primer_site_index_dict,primer_zotu_dict,ESV_table_infile,zotu_reads)

    del ESV_table_infile #FREE UP MEMORY
    del sequence_infile #FREE UP MEMORY

    #Create function variables = args.alignment,primer_ESV_dict,primer_sequence_dict,zotu_reads
    if args.alignment:
        primer_zotu_dict_after_msa,primer_sequence_dict_after_msa,alignment_sum_dict = af.align_and_filter_sequences(primer_sequence_dict,out_dir)
        primer_ESV_dict_after_msa = pull_ESV_data_by_primer_after_msa(primer_ESV_dict,primer_zotu_dict_after_msa)
        primer_taxonomy_dict = pull_taxonomy_data_by_primer(taxonomy_infile, primer_zotu_dict_after_msa)
        write_summary(alignment_sum_dict,out_dir)
    else:
        print("No alignment created")
        primer_ESV_dict_after_msa = primer_ESV_dict
        primer_taxonomy_dict = pull_taxonomy_data_by_primer(taxonomy_infile, primer_zotu_dict)

    add_read_abundance_to_zotu(primer_ESV_dict_after_msa,primer_sequence_dict_after_msa,primer_taxonomy_dict,zotu_reads)
        

    output_ESV_fasta_and_taxonomy(primer_ESV_dict_after_msa,primer_sequence_dict_after_msa,primer_taxonomy_dict,out_dir)

if __name__ == "__main__":
    main()