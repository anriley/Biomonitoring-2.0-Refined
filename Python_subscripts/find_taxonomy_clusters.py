def open_file(filename): #Open and read in files
    with open(filename, "r") as file:
        temp = file.readlines()
    in_data = [line[:-1] for line in temp] #Remove newline from end of lines
    return in_data

def prepare_data_ssv(data): #Make csv files list ( all lines) of lists (lines split by commas)
    prepared_data = []
    for line in data:
        prepared_data.append(line.split(" ")) #For csv files split on comma
    return prepared_data

def get_taxonomy_level_index(taxonomy_header,taxonomy): #Get index of taxonomy level in taxonomy.csv
    for x in range(0,len(taxonomy_header)):
        if taxonomy.lower() == taxonomy_header[x].lower():
            taxonomy_index = x
    return taxonomy_index

def add_dict_list(d,key,value):
    if key in d:
        d[key].append(value)
    else:
        d[key] = [value]

def get_zotu_taxonomy_dict(taxonomy_data,BP_cutoff,taxonomy_index): #Based on sBP and the taxonomy.csv determine if a haplotype/zotu is assigned to a species or not
    a_dict = {}
    a_list = []
    u_list = []
    for line in taxonomy_data[1:]:
        if len(line) > 1:
            if float(line[taxonomy_index+2]) >= float(BP_cutoff):
                add_dict_list(a_dict,line[taxonomy_index],line[0])
                a_list.append(line[0])
            else:
                u_list.append(line[0])
    return a_dict,a_list,u_list

def get_hierarchy(taxonomy_data,taxonomy_index):
    species_hier_dict = {}
    for line in taxonomy_data[1:]:
        if len(line) >1:
            species_hier_dict[line[taxonomy_index]] = [line[taxonomy_index-3],line[taxonomy_index-6],line[taxonomy_index-9],line[taxonomy_index-12]]
    return species_hier_dict

def assign_taxonomy(taxonomy_file, taxonomy_level, BP):
    taxonomy_level_index = get_taxonomy_level_index(taxonomy_file[0],taxonomy_level)
    a_dict,a_list,u_list = get_zotu_taxonomy_dict(taxonomy_file,BP,taxonomy_level_index)
    hier_dict = get_hierarchy(taxonomy_file,taxonomy_level_index)
    return a_dict,a_list,u_list,hier_dict

def find_zotus_in_one_cluster(zotus, dir, cluster_prefix, cluster_suffix):
    d = 0
    while d <= 255:
        clusters_at_d = prepare_data_ssv(open_file(dir+cluster_prefix+str(d)+cluster_suffix))
        for x in range(0,len(clusters_at_d)):
            if set(zotus).issubset(clusters_at_d[x]): #Check if all zotus are present in one cluster
                d_min = d
                cluster = clusters_at_d[x]
                return d_min,cluster #If all zotus are in one cluster, return cluster information
        d += 1
    return 256,0

def check_if_zotus_are_isolated(cluster,assigned_list,taxon_zotus):
    for zotu in cluster:
        if zotu in assigned_list and zotu not in taxon_zotus: #Check if there are zotus in the assigned list, but not assigned to the current taxon
            return False #Return false if there are zotus assigned to different taxa
    return True #Return true if all assigned zotus belong to the current taxon 

def find_d_max(dir, cluster_prefix, cluster_suffix,d_min,min_cluster):
    d_max = d_min
    d = d_min + 1
    while True:
        flag = 1 #Set flag to exit loop if d_max has been found
        clusters_at_d = prepare_data_ssv(open_file(dir+cluster_prefix+str(d)+cluster_suffix))
        for x in range(0,len(clusters_at_d)):
            if set(clusters_at_d[x]) == set(min_cluster): #check if cluster is still present in clusters at current d threshold
                d_max = d
                d += 1
                flag = 0 #Stay in while loop, d_max could be at higher d
        if flag:
            return d_max

#This needs to be fixed                
def find_d_min(cluster,taxon_zotus, dir, cluster_prefix, cluster_suffix,d):
    zotus_in_cluster = [z for z in cluster if z in taxon_zotus]
    cluster_at_d_min = cluster
    while d > 0:
        d -= 1
        flag = True        
        clusters_at_d = prepare_data_ssv(open_file(dir+cluster_prefix+str(d)+cluster_suffix))
        for x in range(0,len(clusters_at_d)):
            if set(zotus_in_cluster).issubset(clusters_at_d[x]):
                cluster_at_d_min = clusters_at_d[x]
                flag = False
        if flag:
            return d+1,cluster_at_d_min,zotus_in_cluster
    return 0,cluster_at_d_min,zotus_in_cluster

def find_zotus_across_fewest_clusters(zotus,d, dir, cluster_prefix, cluster_suffix,assigned_list):
    output = []
    taxon_zotus = zotus #get zotus of taxon
    while d > 0: #cycle through d 
        d -= 1
        clusters_at_d = prepare_data_ssv(open_file(dir+cluster_prefix+str(d)+cluster_suffix))
        clusters_w_zotus = []
        isolated_clusters = []
        mixed_clusters = []
        for x in range(0,len(clusters_at_d)):
            for z in zotus:
                if z in clusters_at_d[x] and clusters_at_d[x] not in clusters_w_zotus: #check if cluster has just taxon zotus, ignore zotus that are already isolated
                    if check_if_zotus_are_isolated(clusters_at_d[x],assigned_list,zotus):
                        isolated_clusters.append(clusters_at_d[x]) #get clusters that are isolated
                    else:
                        mixed_clusters.append(clusters_at_d[x]) # get clusters with taxon and other zotus
                    clusters_w_zotus.append(clusters_at_d[x])
        isolated_zotus = sum(isolated_clusters,[])
        zotus_not_clustered = []
        for z in zotus: #remove zotus that isolated
            if z not in isolated_zotus:
                zotus_not_clustered.append(z)
        zotus = zotus_not_clustered
        for ic in isolated_clusters:
            min_d,cluster_d_min,zotus_in_cluster = find_d_min(ic,taxon_zotus, dir, cluster_prefix, cluster_suffix,d)
            taxon_in_cluster = [i for i in taxon_zotus if i in zotus_in_cluster]
            #print(["here",taxon_in_cluster,min_d,d,zotus_in_cluster,ic])
            output.append([zotus_in_cluster,min_d,d,cluster_d_min,ic])
        if len(mixed_clusters) == 0:
            break
    return output

def taxonomy_cluster_search(in_taxonomy, taxon_level, bp, cluster_dir,cluster_filename_prefix,cluster_filename_suffix):
    assigned_dict,assigned_list,unassigned_list,hierarchy_dict = assign_taxonomy(in_taxonomy, taxon_level, bp)
    taxa_isolated = [["class","order","family","genus",taxon_level.lower(),"taxon_haplotypes","d_min","d_max","cluster_size","cluster"]]
    taxa_large_d = [[taxon_level.lower(),"zotus"]]
    taxa_split_summary = [["class","order","family","genus",taxon_level.lower(),"taxon_haplotypes","d_split","split_cluster_size","cluster"]]
    taxa_split =[["class","order","family","genus",taxon_level.lower(),"subcluster_ID","taxon_haplotypes","d_min","d_min_size",
                    "d_max","d_max_size","d_min_cluster","d_max_cluster"]]
    for taxon in assigned_dict:
        d_min,cluster_d_min = find_zotus_in_one_cluster(assigned_dict[taxon],cluster_dir,cluster_filename_prefix,cluster_filename_suffix) #For the current taxon, find the smallest cluster that contains all zotus
        if d_min > 255:
            taxa_large_d.append([taxon," ".join(assigned_dict[taxon])]) #If no cluster exist that contains all zotus with a threshold of 255 or less, but taxon in large d file
        else:
            isolated_taxon_cluster = check_if_zotus_are_isolated(cluster_d_min,assigned_list,assigned_dict[taxon]) #Evaluates to TRUE if only zotus assigned to that taxon are present in the cluster
            if isolated_taxon_cluster: 
                #For isolated clusters
                d_max = find_d_max(cluster_dir,cluster_filename_prefix,cluster_filename_suffix,d_min,cluster_d_min) #Find highest d value where cluster membership remains unchanged
                taxa_isolated.append([hierarchy_dict[taxon][3],hierarchy_dict[taxon][2],hierarchy_dict[taxon][1],hierarchy_dict[taxon][0],taxon, #Append data for the isolated cluster file
                                    len(assigned_dict[taxon]),d_min,d_max,len(cluster_d_min)," ".join(cluster_d_min)])
            else:
                #For clusters that did not isolate
                split_clusters = find_zotus_across_fewest_clusters(assigned_dict[taxon],d_min,cluster_dir,cluster_filename_prefix,cluster_filename_suffix,assigned_list)
                taxa_split_summary.append([hierarchy_dict[taxon][3],hierarchy_dict[taxon][2],hierarchy_dict[taxon][1],hierarchy_dict[taxon][0],taxon,
                                        len(assigned_dict[taxon]),d_min,len(cluster_d_min)," ".join(cluster_d_min)])
                for i in range(0,len(split_clusters)):
                    taxa_split.append([hierarchy_dict[taxon][3],hierarchy_dict[taxon][2],hierarchy_dict[taxon][1],hierarchy_dict[taxon][0],taxon,
                                    i,len(split_clusters[i][0]),split_clusters[i][1],len(split_clusters[i][3]),
                                    split_clusters[i][2],len(split_clusters[i][4])," ".join(split_clusters[i][3])," ".join(split_clusters[i][4])])
    return taxa_isolated,taxa_split_summary,taxa_split
