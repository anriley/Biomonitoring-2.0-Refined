# Biomonitoring 2.0 Refined
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15546739.svg)](https://doi.org/10.5281/zenodo.15546739)

This code was created to process MetaWorks output and create clusters that are bound by taxonomic assignments.

After running MetaWorks with results option 2, you will need to take the ESV.table, taxonomy.csv, and chimera.denoised.nonchimeras.taxon files. Do not change their names.
Place these files in a directory, in this example Data/. If using Data/ included in this repository unzip ESV.zip and make sure it is named ESV.table.

The below code separates data by amplicons/markers, prepares data for SWARM clustering (adds read counts to sequence headers), and performs alignment filtering.
If you want to perform alignment filtering include the option -a, MAFFT must be installed and included in your path.
Use option -d to provide the directory of the input files. Use the option -o to provide the directory for the output files. This can be the same directory.
This code creates the directories Sequences/, ESV_table/, and Taxomony/. These directories contain the same input files split up by amplicon/marker. 
If the -a option is used, an MSA/ directory is created that contains the multiple sequence alignments that were used to filter
```
python3 prepare_data_for_SWARM_and_alignment_filter.py -d Data/ -o Data/ -a
```
SWARM must be installed and included in your path in order to run clustering.

Before generating clusters and species bound clusters, the guide file is used to input parameters:

After "Amplicons:" list amplicon/marker names separated by a comma. Use the same amplicon/marker names that were used in MetaWorks.

After "Filter:" provide read count filters by a comma. Filters will be less than.

After "Taxonomy:" provide the taxonomic level of assignment and the bootstrap support value to use for assignment, separated by a comma.

Example guide file:
```
Amplicons:F230R,MLJG
Filters:100
Taxonomy:Species,0.8
```
The below code filters data, creates clusters, and creates species bound clusters. Use the -d option to provide the input directory, use the same output directory from prepare_data_for_SWARM_and_alignment_filter.py, all output will also be placed in this directory. Use the -g option to provide the guide file. Files with the read count filter applied will be created in the Sequences/, ESV_table/, and Taxomony/ directories. Clusters will be placed in subdirectories in the Clusters/ directory, subdirectories will have a unique name based on the parameters used.
Species bound clusters will be placed in the Taxonomy_clusters/ directory. Files ending with the "isolated" suffix contain species bound clusters .
Files with the "split" and "split_summary" suffix provide information on taxonomic assignments that were not recovered by clustering.
```
python3 filter_cluster_and_search.py -d Data/ -g guide_file.txt
```
R scripts are included to show how Tables and Figures were generated. The script setup_parameters.R must be run first, followed by community_and_population_tables_dms.R

## How to cite

If you use any of the provide scripts please cite the Biomonitoring 2.0 paper:

In review

Or cite the repository:

Riley, Andrew C. (2025, May 18). Biomonitoring 2.0 Refined. Zenodo. http://doi.org/10.5281/zenodo.15546739

