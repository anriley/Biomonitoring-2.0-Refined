import subprocess
import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def add_dict_list(d,key,value): #Append values to lists for keys already present in dictionaries
    if key in d:
        d[key].append(value)
    else:
        d[key] = [value]
    return

def write_file(outname,out_data):
    fo = open(outname,"w")
    for line in out_data:
        fo.write(">"+"\n".join(line))
        fo.write("\n")
    fo.close()
    return

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

def align_sequences(in_sequences,primer,out_dir):
    out_dir_MSA = prepare_directory(out_dir + "MSA")
    out_fasta_for_align = out_dir_MSA + "mafft_nt_"+primer+".fasta"
    write_file("temp.fasta",in_sequences)                             #Create temporary fasta file to create the alignment
    subprocess.call("/mnt/c/Users/andre/Documents/mafft-linux64/mafft.bat --op 9999 --thread 8 temp.fasta > %s" % out_fasta_for_align, shell=True) #Run the MAFFT alignment
    aln = AlignIO.read(out_fasta_for_align, "fasta")                                  #Read in alignment file
    return aln

def find_gap_no_gap_indices(aln):
    no_gap_indices = []
    no_gap_percent = []
    gap_indices = []
    gap_percent = []
    for x in range(0,len(aln[0])):
        if aln[:,x].count("-")/len(aln[:,x]) < 0.5:                                 #Find sites in alignment with less than 50% gaps
            no_gap_indices.append(x)
            no_gap_percent.append(aln[:,x].count("-")/len(aln[:,x]))
        else:                                                                       #Find sites in alignment with more than 50% gaps
            gap_indices.append(x)
            gap_percent.append(aln[:,x].count("-")/len(aln[:,x]))
    return no_gap_indices,gap_indices,max(no_gap_percent),min(gap_percent)

def filter_sequences_using_alignment(aln,no_gap_indices,gap_indices,io_primer_sequence_dict,io_primer_zotu_dict,key):
    for record in aln:
        no_gaps = [record.seq[i] for i in no_gap_indices]
        gaps = [record.seq[i] for i in gap_indices]
        if "-" not in no_gaps and all(i=="-" for i in gaps):                        #Keep sequences that have bases where sites are less than 50% gaps and gaps where siters have more than 50% gaps
            out_seq = [record.id,"".join(no_gaps)]
            add_dict_list(io_primer_sequence_dict,key,out_seq)
            add_dict_list(io_primer_zotu_dict,key,record.id)
    return

def align_and_filter_sequences(in_sequence_dict,out_dir):
    out_primer_zotu_dict = {}
    out_primer_sequence_dict = {}
    out_sum_dict = {}
    primers = [key for key in in_sequence_dict]                                         #Get keys in list, so that keys can be cleared
    for p in primers:
        sequences_before = len(in_sequence_dict[p])##################
        aln = align_sequences(in_sequence_dict[p],p,out_dir)
        in_sequence_dict.pop(p)##########################                                                         #Remove from dictionary to free up memory
        no_gap_check,gap_check,max_no_gap,min_gap = find_gap_no_gap_indices(aln)
        filter_sequences_using_alignment(aln,no_gap_check,gap_check,out_primer_sequence_dict,out_primer_zotu_dict,p)
        sequences_after = len(out_primer_sequence_dict[p])############################################
        out_sum_dict[p] = [sequences_before,sequences_after,max_no_gap,min_gap] #Collect data on the alignment
    return out_primer_zotu_dict,out_primer_sequence_dict,out_sum_dict


