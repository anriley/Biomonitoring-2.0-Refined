def filter_esv_table(esv_table,filter):
    filtered_zotus=[]
    filtered_esv_table=[esv_table[0]]
    for line in esv_table[1:]:
        reads = sum([int(i) for i in line[1:]])
        if reads >= filter:
            filtered_zotus.append(line[0])
            filtered_esv_table.append(line)
    return filtered_zotus,filtered_esv_table

def filter_taxonomy(taxonomy,filtered_zotus):
    filtered_taxonomy = [taxonomy[0]]
    for line in taxonomy[1:]:
        for ft in filtered_zotus:
            if line[0] == ft:
                filtered_taxonomy.append([ft] + line[1:])
    return filtered_taxonomy

def filter_aln(sequences,filtered_zotus):
    filtered_aln = []
    for record in sequences:
        if record.id in filtered_zotus:
            filtered_aln.append(record)
    return filtered_aln

def filter_data(esv_table,taxonomy,sequences,filter):
    filtered_zotus,filtered_esv_table = filter_esv_table(esv_table,filter)
    filtered_taxonomy = filter_taxonomy(taxonomy,filtered_zotus)
    filtered_aln = filter_aln(sequences,filtered_zotus)
    return filtered_esv_table,filtered_taxonomy,filtered_aln