# TCGA-KM-R
looking through the The Cancer Genome Atlas dataset, and running Kaplan–Meier analysis in R

to start, 3rd tier data including gene expression info is the easiest to work with
go to https://genome-cancer.ucsc.edu/proj/site/hgHeatmap/ or http://gdac.broadinstitute.org/ for downloading data of various cancers
without a specific type in mind, http://firebrowse.org/viewGene.html will give a summary of expression patterns

the R script will compile gene expression data with patient information, and perform Kaplan–Meier analysis (survival plot)
