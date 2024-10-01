######################################################################### 
#### Script for statistical analysis and plotting of data (figure 2) ####
####															 
#### Author: Dr. Niccolo Bassetti - niccolo.bassetti@protonmail.com
#### 
#### Principal Investigator:   Dr. Nina Fatouros - nina.fatouros@wur.nl
####
####
#### Dataset: rqtl_HR_BC1.csv
####

# Set working directories
rm(list=ls())
sessionInfo() #check packages available

setwd("write_your_path")	# set working directory @ home

## Set R library folder location
.libPaths()
Sys.getenv("R_LIBS_USER") 		#same as above

# Install R packages (update.packages(), install.packages(), library())
packageurl <- "https://cran.r-project.org/src/contrib/Archive/devtools/devtools_1.13.2.tar.gz" #write URL of package to install
install.packages(packageurl, repos=NULL, type="source") #install old version packages
install.packages("glue")

### Install and load packages
packages = c("qtl", "ggplot2", "tidyverse", "binom","dplyr","RColorBrewer")

lapply(packages, install.packages, character.only = TRUE)

lapply(packages, library, character.only = TRUE)


### 1) Load data
data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BC1.csv", crosstype = "bc", estimate.map = TRUE)


summary(data)

### 2.1) Calculate marker segregation distortion
# NOTE: NOT USED because skewness in the distribution of phenotypic data.
#       Marker distortion was instead tested with Chi-square test, see file "1_KASP-01_final_data_BC1.xlsx".
gt <- geno.table(data)
gt
  #write.table(gt, file = "marker.distortion.txt")


### 2.2) QTL analysis - single QTL analysis
# method = "hm"
data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BC1.csv", crosstype = "bc")

data = jittermap(data)
data_2 = calc.genoprob(data, step=1, error.prob=0.001)

# calculate LOD score for marker-trait association
out.em_2 <- scanone(data_2, model = "binary")
plot(out.em_2, col = "red", ylab="LOD score", alternate.chrid=TRUE)

# calculate LOD score threshold
out.em_2_perm <- scanone(data_2, n.perm = 1000)
summary(out.em_2_perm)

#  extract exact markers (remove imputation) and gives actual LOD scores
out.em_2
lod_score=mqmextractmarkers(out.em_2) 
lod_score

  #write.table(lod_score, file = "lod_score_BC1.txt") # FINAL RESULTS for Sup Table S2
