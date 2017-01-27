# TCGA-KM-R
looking through the The Cancer Genome Atlas dataset for expression of a list of genes, and running Kaplan–Meier analysis in R, with multiple testing correction

to start, 3rd tier data including gene expression info is the easiest to work with
go to https://genome-cancer.ucsc.edu/proj/site/hgHeatmap/ or http://gdac.broadinstitute.org/ for downloading data of various cancers
without a specific type in mind, http://firebrowse.org/viewGene.html will give a summary of expression patterns

the R script will compile gene expression data with patient information, and perform Kaplan–Meier analysis (survival plot)

Step-by-step:

1)download data from UCSC

2)unzip to master folder

3)change the path in the script to the master folder

4)change the gene names to ones of your interest

5)run script, a pdf with the KM plot will be added to each sub-folder, p-values and fdr will be reported in csv file

or see variants for single gene or single cancer type
