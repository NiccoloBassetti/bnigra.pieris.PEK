######################################################################### 
#### Script for statistical analysis and plotting of data (figure S2) ####
####															 
#### Author: Dr. Niccolo Bassetti - niccolo.bassetti@protonmail.com
#### 
#### Principal Investigator:   Dr. Nina Fatouros - nina.fatouros@wur.nl
####
####
#### Dataset: figure_S2_dataset.txt
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
packages=c("ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse", "binom", "dplyr", "multcompView")

lapply(packages, install.packages, character.only = TRUE)

lapply(packages, library, character.only = TRUE)


#### Prepare dataset ####
#_______________________#

mydata <- read.table(file = "figure_S2_dataset.txt",sep = "\t", header = TRUE) # !!! test sep character

str(mydata)
head(mydata)


# Reformat dataframe for plotting purpose
data_temp  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
data_temp$Treatment = factor(data_temp$Treatment, levels = unique(data_temp$Treatment)) # helps to plot x-axis variable in correct order
data_temp$Max_hr_score = as.numeric(data_temp$Max_hr_score) # helps to plot x-axis variable in correct order
data_temp


# Filter data to keep only a treatment
head(data_temp)
tail(data_temp)

data = filter(data_temp, Treatment == "Alternaria")
data = filter(data_temp, Treatment == "Eggwash")
data = filter(data_temp, Treatment == "Rhizoctonia")
data = filter(data_temp, Treatment == "Xanthomonas")
data = filter(data_temp, Treatment == "Water droplet")
data = filter(data_temp, Treatment == "Water infiltrated")

str(data)


#### Statistical analysis ####
#____________________________#

packages=c("multcomp", "multcompView", "FSA", "rcompanion")

lapply(packages, library, character.only = TRUE)


# Test for normality and homogeneity of variances
qqnorm(data$Max_hr_score)
qqline(data$Max_hr_score)
shapiro.test(data$Max_hr_score)		# test for normality (small P = NOT NORMAL)
fligner.test(data$Max_hr_score~data$Genotype)	# test for equal variances (small P = NOT EQUAL VARIANCE)



# non-parametric test
kruskal.test(data$Max_hr_score~data$Genotype)

# Data were plotted in Excel to generate figure S2.See file "figure_S2.xlsx".
