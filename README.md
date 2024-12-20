# Species_bound_clusters
After running MetaWorks with results option 2, you will need to take the ESV.table, taxonomy.csv, and chimera.denoised.nonchimeras.taxon files. Do not change their names.
Place these files in a directory, in this example Data/.
If you want to perform alignment filtering include the option -a, MAFFT must be installed and included in your path.
Use option -d to provide the directory of the input files. Use the option -o to provide the directory for the output files. This can be the same directory.
This code separates data by primers/amplicons, prepares data for SWARM clustering (adds read counts to sequence headers), and performs alignment filtering.
```
python3 prepare_data_for_SWARM_and_alignment_filter.py -d Data/ -o Data/ -a
```
Before generate clusters and species bound clusters, the guide file is used to input parameters:

After "Amplicons:" list marker names separated by a comma. Use the same marker names that were used in MetaWorks.

After "Filter:" provide read count filters by a comma. Filters will be less than.

After "Taxonomy:" provide the taxonomic level of assignment and the bootstrap support value to use for assignment, separated by a comma.
```
Amplicons:F230R,MLJG
Filters:100
Taxonomy:Species,0.8
```
To generate clusters and species bound clusters, SWARM must be installed and included in your path. Provide the previously used output directory in the option -d and provide the guide file in the option -g
```
python3 filter_cluster_and_search.py -d Data/ -g guide_file.txt
```
Clusters will be placed in the 
