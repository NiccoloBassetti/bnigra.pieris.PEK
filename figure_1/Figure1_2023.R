##################################################################### 
#### 
#### Data analysis and graphs for Figure 1 for Bassetti et al. PEC 2023. 
####
#### Author: Dr. Lotte Caarls - lotte.caarls@wur.nl
#### 
#### Principal Investigator:   Dr. Nina Fatouros - nina.fatouros@wur.nl
####



# Data to test variation of HR response in different B.nigra accessions

require(dplyr)
getwd()
setwd("C:/Users/caarl001/OneDrive - Wageningen University & Research/PUBLICATIONS/Bassetti 2023")
file <- "RawdataFigure1.csv"
# import dataset
data <- read.csv(file,sep=";")
#rename 
data <- data  %>%  
  rename(experiment=ï..experiment)
#give me a summary of the data
summary(data)
#change to factors
data$accession <- as.factor(data$accession)
data$plant <- as.factor(data$plant)
data$treatment <- as.factor(data$treatment)

#reorder levels 
order <- levels(data$plant)
order
data$plant <- factor(data$plant,levels = c("DG1","SF47-O","SF29-O","SF19-O","DG29","DG12","SF3-O","SF25-O","SF48-O"))


#make an occurance (HR yes) column
data$HR <- as.factor(ifelse(data$score>=2,"1","0"))
data$noHR <- as.factor(ifelse(data$score<2,"1","0"))
data$score0 <- as.factor(ifelse(data$score==0,"1","0"))
data$score1 <- as.factor(ifelse(data$score==1,"1","0"))
data$score2 <- as.factor(ifelse(data$score==2,"1","0"))
data$score3 <- as.factor(ifelse(data$score==3,"1","0"))
data$score4 <- as.factor(ifelse(data$score==4,"1","0"))

#remove NAs
data <- na.omit(data)

#summary table of the data
#summarize information for accession and store in 'supdata': N, N of plants, HRy, HR ratio, mean, SD and SE
supdata <- data %>%
  group_by(data$plant) %>%
  summarize(
            class0    = length(score0[score0 == 1]),
            class1    = length(score1[score1 == 1]),
            class2    = length(score2[score2 == 1]),
            class3    = length(score3[score3 == 1]),
            class4    = length(score4[score4 == 1]),
            N    = length(score),
            HRy   = length(HR[HR == 1]),
            HRn   = length(HR[HR == 0]),
            HRratio = HRy/N,
            Mean = mean(score),
            SD = sd(score),
            SE = SD/N)
as.data.frame(supdata)
write.table(supdata,"plantresult.csv", append = FALSE, sep = ";", dec = ".",
            row.names = TRUE, col.names = TRUE)

# is there stat difference between accessions?
#### Statistical analysis ####
packages=c("multcomp", "multcompView", "FSA", "rcompanion")
lapply(packages, library, character.only = TRUE)
# Check normality of the data
qqnorm(data$score)
qqline(data$score)
shapiro.test(data$score)		# test for normality (small P = NOT NORMAL)
fligner.test(data$score~data$plant)	# test for equal variances (small P = NOT EQUAL VARIANCE)
#non-parametric test
kruskal.test(data$score~data$plant)
#sig different
# pairwise comparison: Withney-Mann U test/Wilcoxon rank sum test
PT = pairwise.wilcox.test(data$score, data$plant,
                          p.adjust.method = "BH")
PT = PT$p.value
PT1 = round(fullPTable(PT),3) # require "rcompanion"

multcompLetters(PT1,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE) # require "multcompView"
PT1

#produce graph
packagesgraph=c("ggplot2", "reshape2", "RColorBrewer")
lapply(packagesgraph, library, character.only = TRUE)

#first prepare table
#decide which variables to use in the table
rowvariable <- data$score
colvariable <- data$plant

# First make a table to count the scores per variable
table1 <- table(rowvariable,colvariable)
table1
#make a proportions table,  the second argument is if totals should come from rows or columsn 
prop <- prop.table(table1,2)
perc <- prop*100

#Make  barplot with ggpplot
#Need to convert the table into a frame
data2 <-as.data.frame.matrix(perc)
#add score as seperate variable
data2$score<-rownames(data2)
#convert to long format
data2long <- melt(data2, id.vars=c("score"), value.name = "proportion")
#paste variable as name column 2
names(data2long)[2] <-paste("treatment")
# make the plot
# the order of the factors is the "wrong" way around per default, so reverse is TRUE
display.brewer.pal(n=9,"Blues")
my.blue <- brewer.pal(n=9,"Blues")[c(2,3,5,7,9)]
bp <- ggplot() + geom_bar(aes(y = proportion,
                              x = treatment,
                              fill = score),
                          data=data2long,
                          stat="identity",
                          position=position_stack(reverse=TRUE), width = 0.6) + # change the width of the bars
  scale_fill_manual(values=my.blue,name = 'class', guide = guide_legend(reverse=TRUE)) + # change the colours and reverse the legend
  scale_y_continuous(name = "% plants in class", expand = c(0,0), limits = c(0,101)) + 
  scale_x_discrete(name = "")
bp 

require(svglite)

# save plot as image, set measurements
graphname <- "Bassetti_figure1.svg"
