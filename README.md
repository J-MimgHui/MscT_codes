The MscT_cods project contains five script files.
The codes in these scripts can help reproduce the calculation, tables, pictures, etc. of the Resource article "Single-cell transcriptomics across 2,534 microbial species reveals functional heterogeneity in the rumen microbiome".

The first file 'Reading_and_Filtering.R'
The codes in this file is used to read the microbial single-cell gene expression matrix into R to generate Seurat objects. At the same time, quality control is performed on all read-in cells, removing low-quality cells and retaining high-quality cells.
The microbial single-cell gene expression matrix (MscT_matrix) could be found in the link below.
"https://figshare.com/articles/dataset/Microbiome_single-cell_transcriptomics_reveal_functional_heterogeneity_of_metabolic_niches_covering_more_than_2_500_species_in_the_rumen/24844344"

The second file 'MscT_initial_matrix.R'
The codes in this file is used to perform benchmarking and cluster analysis.
Please note that not all of the R packages and colors shown as loaded in the script are necessarily required, and they can be selected on a case-by-case basis when using this script.

The third file 'MscT_HMACs.R'
The codes in this file is used to perform recluster analysis for HMACs (high metabolic activaty cells).

The fourth file 'MscT_basfia.R'
The codes in this file is used to perform recluster analysis for the cells of Basfia succiniciproducens species.

The fifth file 'MscT_basfia_monocle.R'
The codes in this file is used to perform pseudo-time analysis for the cells of Basfia succiniciproducens species.

Please note that raw data filtering, comparison, etc. can be done in linux system by using software according to the regular process.
Calculating the proportion of functional genes, calculating the difference between groups, etc. are routine analyses, which can be done with the R package or other statistical software.
The code will be kept up to date.
