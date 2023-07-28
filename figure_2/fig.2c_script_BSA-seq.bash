##################################################################### 
#### Pipeline to analyse NGS data for B. nigra BSA-01 experiment ####
####															 
#### (BSA coupled with NGS, following COSSA pipeline - https://github.com/cprodhom/CoSSA-workflows)
####
#### Author: Dr. Niccolò Bassetti - niccolo.bassetti@protonmail.com
#### 
#### Principal Investigator:   Dr. Nina Fatouros - nina.fatouros@wur.nl
####
####
#### Dataset: DNA sequencing data deposited on European Nucleotide Achieve (ENA) under accession number PRJEB64240


##########################################
#### 0 - QC and raw data manipulation ####

## Check quality raw data
# loop
nano reads_report.txt # make report reads number
for filename in *.fq; do
	wc -l $filename > cat nano reads_report.txt # "wc -l" counts lines in file
done


## Rename bulk files	
for f in *_*.fq.gz ; do mv -- "$f" "${f//_1./.1.}"; done
for f in *_*.fq.gz ; do mv -- "$f" "${f//_2./.2.}"; done

## Combine .fq reads output
for sample in $(ls *fq | rev | cut -c 6- | rev | uniq)
do 
cat $sample.1.fq >> $sample.fq # repeated for d1 and d2
cat $sample.2.fq >> $sample.fq # repeated for d1 and d2
echo "completed $sample"
done





################
#### FastQC ####
# loop through all samples
for filename in *.fq.gz
do
../Programs/FastQC/fastqc $filename # change .fq for other commands
echo "Completed"
done



#################################
#### Reads quality filtering ####
####

# NOTE. Raw reads were not trimmed because FASTQC report showed good basecall quality acroos the whole 150 bp read length

#### Trimmomatic

sample_read1="D1_L4_1"	# change name to process different sample
sample_read2="D1_L4_2"	# change name to process different sample

java -jar ../Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 -phred33 D2_L4_1.fq D2_L4_2.fq D2_L4_1.fq.trim.fq D2_L4_1.fq.trim.unpaired.fq D2_L4_2.trim.fq D2_L4_2.trim.unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75



###################################################
#### Kmer calculation and kmer sets comparison ####

### Software required
### GenomeTester4 v4.0 - https://github.com/bioinfo-ut/GenomeTester4
###
### NOTE. Don`t use v4.1, it has a bug which doesn`t allow a correct k-mer calculation


# set environment variables 

glistmaker="../GenomeTester4-Version_4_0_stable/bin/glistmaker"
glistcompare="../GenomeTester4-Version_4_0_stable/bin/glistcompare"
glistquery="../GenomeTester4-Version_4_0_stable/bin/glistquery"

## k-mer tables calculation
for sample in $(ls *fq | cut -c -7)
do
$glistmaker $sample.fq -w 31 --num_threads 60 -o $sample.raw
2> error.glistmaker.$sample.txt
echo "Completed"
done
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tcalculation of 31-mer for $sample" >> $log

	# alternative to for loop (process one sample at the time)
	sample="D6_L1_2"
	$glistmaker $sample.fq -w 31 --num_threads 60 -o $sample.raw
	2> error.glistmaker.$sample.txt
	echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tcalculation of 31-mer for $sample" >> $log

## unite k-mers set from read_1 and read_2
sample_read1="D1_L4_1"
sample_read2="D1_L4_2"
sample="D1_L4"

$glistcompare $sample_read1.raw_31.list $sample_read2.raw_31.list -u -o $sample.raw_31.list
mv $sample.raw_31.list_31_union.list $sample.raw_31_union.list
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tunion of 31-mer list for read_1 and read_2 for $sample" >> $log
	
	# filter k-mers with abundance >1
	$glistcompare $sample.raw_31.list $sample.raw_31.list -i -c 2 -o $sample.cutoff.02
	
	# filter all samples at 10 < x < 21 coverage
	$glistcompare $sample.raw_31_union.list $sample.raw_31_union.list -i -c 10 -o $sample.cutoff.10
	$glistcompare $sample.raw_31_union.list $sample.raw_31_union.list -i -c 21 -o $sample.cutoff.21 
	$glistcompare $sample.cutoff.10_31_intrsec.list $sample.cutoff.21_31_intrsec.list -d -o $sample.cutoff.10to21
	
	# check stat of k-mers list
	$glistquery $sample.cutoff.10_31_intrsec.list -stat
	$glistquery $sample.cutoff.21_31_intrsec.list -stat
	$glistquery $sample.cutoff.10to21_31_0_diff1.list -stat

	# transform binary k-mers list in text file tables
	$glistquery $sample.raw_31.list > $sample.stat.pre.filter.kmer
	
done

## R-allele (score 4)  in R_bulk from R_parent: ((R_bulk – S_bulk)-S_parent (DG1-S1)) intersect with R_parent (F1_#130) and R_donor (SF48-O1): 
# NOTE. Both S_bulk (HR=0) and (HR=1) were subtracted from R_bulk

# 1. difference between bulk k-mer sets (HR_bulk)
$glistcompare D6_L4.cutoff2_31_intrsec.list D4_L4.cutoff2_31_intrsec.list -d -o D6-D4_bulk.cutoff.02
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tR_bulk (4): D6-D4" $(ls -t1 |  head -n 1) >> $log
$glistcompare D6-D4_bulk.cutoff.2_31_0_diff1.list D5.cutoff.02_31_intrsec.list -d -o R_bulk.cutoff.02
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tR_bulk (4): D6-D4-D5" $(ls -t1 |  head -n 1) >> ../$log

	# filter kmers with abundance <10 or >21
	$glistcompare R_bulk.cutoff.02_31_0_diff1.list R_bulk.cutoff.02_31_0_diff1.list -i -c 10 -o R_bulk.cutoff.10
	$glistcompare R_bulk.cutoff.02_31_0_diff1.list R_bulk.cutoff.02_31_0_diff1.list -i -c 21 -o R_bulk.cutoff.21 
	$glistcompare R_bulk.cutoff.10_31_intrsec.list R_bulk.cutoff.21_31_intrsec.list -d -o R_bulk.cutoff.10to21
	echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tR_bulk (4): filtered 10<kmer<21" $(ls -t1 |  head -n 1) >> $log

	# transform binary k-mers list in text file tables
	$glistquery R_bulk.cutoff.02_31_0_diff1.list -stat
	$glistquery R_bulk.cutoff10to21_31_0_diff1.list -stat

	
# 2. remove DG1-S1 (S_parent)
$glistcompare R_bulk.cutoff.10to21_31_0_diff1.list D2.cutoff2_31_intrsec.list -d -o R_bulk-D2.cutoff.10to21
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tR_allele (4): remove S_parent:" $(ls -t1 |  head -n 1) >> ../$log

# 3. intersect with F1_#130 (R parent) and SF48-O1 (R donor)
$glistcompare R_bulk-D2.cutoff.10to21_31_0_diff1.list D3.cutoff.02_31_intrsec.list -i -o R_bulkiD3.cutoff.10to21
$glistcompare R_bulkiD3.cutoff.10to21_31_intrsec.list D1.cutoff.02_31_intrsec.list -i -o R_bulk.cutoff.10to21 # final result
mv $(ls -t1 |  head -n 1) RiR_p.list
output=$(ls -t1 |  head -n 1) # select "RiRp-Sp.list" as environment variable
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tR_allele (4): intersect R_parent and SF48:" $(ls -t1 |  head -n 1) >> ../$log

	## Transform the binary k-mer tables files in text files using glistquery
	$glistquery $output > $output.kmer
	echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tR_allele (4): intersect R_parent and SF48: final $output.kmer table" $(ls -t1 |  head -n 1) >> ../$log

	## Mapping the k-mers to a reference ####
	python ../../Programs/CoSSA-workflows-master/kmer_to_fastq.py $output.kmer > $output.kmer.fq # check if number k-mers is correct
	echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tR_allele (4): intersect R_parent and SF48: convert .kmer table into fastaq " $(ls -t1 |  head -n 1) >> ../$log


###########################################
#### Mapping the k-mers to a reference ####

### Softwares required 

### BWA v0.7.17 - https://bio-bwa.sourceforge.net/
### SAMTOOLS v1.16.1 - http://www.htslib.org/download/
### PICARD - v2.22.0 - https://github.com/broadinstitute/picard



# set environment variables 

bwa="../bwa-0.7.17/bwa"

fasta=Bnigra_C2.v1.genome_plastids.fasta 
kmer_set=RiRp_3
ref_genome=C2

# 1. index reference genome
$bwa index $fasta
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tindex reference genome $fasta" >> ../$log

# 2.alignment (ALN/SAMSE)
maxDiff=3 # allow a maximum of 3 mismatches to pass alignment of each kmer

$bwa aln -n $maxDiff -t 30 $fasta $kmer_set.fq > $kmer_set.$maxDiff.$ref_genome.sai
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tbwa aln -n $maxDiff $ref_genome > $(ls -t1 |  head -n 1)" >> ../$log

$bwa samse $fasta $kmer_set.$maxDiff.$ref_genome.sai $kmer_set.fq > $kmer_set.$maxDiff.$ref_genome.sam
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tbwa samse, aln -n $maxDiff $ref_genome, $(ls -t1 |  head -n 1)" >> ../$log

# 3. CHECK validity SAM file
# check https://broadinstitute.github.io/picard/explain-flags.html for SAM flags

## PICARD
# Check validity SAM/BAM file 
# for detail check ( https://gatkforums.broadinstitute.org/gatk/discussion/7571/errors-in-sam-bam-files-can-be-diagnosed-with-validatesamfile)
java -jar ../../Programs/Picard/picard.jar ValidateSamFile \
	I=$kmer_set.$maxDiff.$ref_genome.sam MODE=SUMMARY
	# I=$kmer_set.$maxDiff.$ref_genome.sam IGNORE=RECORD_MISSING_READ_GROUP MO=1000 MODE=VERBOSE

	# Add/replace read group in SAM file
	java -jar ../../Programs/Picard/picard.jar AddOrReplaceReadGroups \
	  I=$kmer_set.$maxDiff.$ref_genome.sam \
	  O=$kmer_set.$maxDiff.$ref_genome.header.sam 2>error\
	  VALIDATION_STRINGENCY=LENIENT \
	  RGID=HWNCTDSXX.4 \
	  RGLB=TruSeq3_PE \
	  RGPL=illumina \
	  RGPU=HWNCTDSXX.4 \
	  RGSM=$kmer_set
	# check if SAM/BAM file has been corrected
	samtools view -H $kmer_set.$maxDiff.$ref_genome.header.sam | grep '^@RG'		

	echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tSAM read group tag added/corrected to $kmer_set.$maxDiff.sam " >> ../$log
		

# 4. create BAM file
samtools view -Sbh -o $kmer_set.$maxDiff.$ref_genome.bam $kmer_set.$maxDiff.$ref_genome.header.sam
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tSAM to BAM conversion, aln -n $maxDiff" >> ../$log
	
	# sort BAM file
	samtools sort -o $kmer_set.$maxDiff.$ref_genome.srt.bam $kmer_set.$maxDiff.$ref_genome.bam
	echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tBAM sorted, aln -n $maxDiff" >> ../$log

# 5. Inspect quality alignment
## PICARD
# Aligment summary metrics
java -jar ../../Programs/Picard/picard.jar CollectAlignmentSummaryMetrics \
	R=$fasta \
	I=$kmer_set.$maxDiff.$ref_genome.srt.bam \
	O=$kmer_set.$maxDiff.$ref_genome.alig.metrics.txt 2>error
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tPicard AligmentSummaryMetrics for $kmer_set.$maxDiff.bam " >> ../$log



########################################
#### k-mer Alignement visualization ####


### Softwares required 

### bedtools v2.30 - https://bedtools.readthedocs.io/en/latest/ 


# set enviromental variables

fasta=Bnigra_C2.v1.genome_plastids.fasta 

interval=1Mb

kmer_set=R_extra	# CHANGE
ref_genome=C2
maxDiff=3 # set the max # of mismatches alllowed

# 1. Split genomes in 1 Mb windows and prepare a .bed file
# Convert genome .fasta file to .bed with chromosome as single entry
	# samtools faidx $fasta -o $fasta.fai
	# alternatively use this awk command:
 	# awk 'BEGIN {FS="\t"};{print $1 FS "0" FS $2}' $fasta.fai > $fasta.bed # add "BEGIN {FS="\t"};" statement for .bed tab delimited (without quotes) 

	# Separate the genome in 1Mb intervals
	../bedtools makewindows -b $fasta.bed -w 1000000 > $fasta.$interval.window.bed
	../bedtools makewindows -b $fasta.bed -w 1000000 -s 10000 > $fasta.$interval.sliding.window.bed # sliding window

# 2. Intersect alignment results with chromosome location
../bedtools intersect -a $fasta.$interval.window.bed -b $kmer_set.$maxDiff.$ref_genome.srt.bam -c -bed > $kmer_set.$maxDiff.$ref_genome.$interval.mapping.txt
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\t$kmer_set.$maxDiff alignemnt for $interval window" >> ../$log
	
	# alternatively you can use sliding window	
	../bedtools intersect -a $fasta.$interval.sliding.window.bed -b $kmer_set.$maxDiff.$ref_genome.srt.bam -c -bed > $kmer_set.$maxDiff.$ref_genome.$interval.sliding_window.mapping.txt


# 3. visualization of k-mer mapping results
########
### NOTE. The file $kmer_set.$maxDiff.$ref_genome.$interval.sliding_window.mapping.txt was imported in Excel (see file "fig.1b_BSAseq_resistant_kmers_plot.xlsx") to generate the final plots of figure 2C.
########

	
# 4. prepare .bam file for visualizationON genome browser (for example IGV)
# Index BAM file for visualization
samtools index $kmer_set.$maxDiff.$ref_genome.srt.bam
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tIndex final .srt.bam file for -n 0,1,3 " >> ../$log

samtools flagstat $kmer_set.$maxDiff.$ref_genome.srt.bam



