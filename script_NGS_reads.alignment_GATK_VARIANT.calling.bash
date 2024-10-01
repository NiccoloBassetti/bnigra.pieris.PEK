############################################################################## 
#### Pipeline for GATK variant calling for B. nigra BSA-01 experiment ########
####
####
#### Author: Dr. Niccolò Bassetti - niccolo.bassetti@protonmail.com
#### 
#### Principal Investigator:   Dr. Nina Fatouros - nina.fatouros@wur.nl
####
####
#### Dataset: DNA sequencing data deposited on European Nucleotide Achieve (ENA) under accession number PRJEB64240



##########################################
#### Mapping raw reads to a reference ####

### Softwares required 

### BWA v0.7.17 - https://bio-bwa.sourceforge.net/
### SAMTOOLS v1.16.1 - http://www.htslib.org/download/
### PICARD - v2.22.0 - https://github.com/broadinstitute/picard

export PATH=/usr/local/bin/samtools:$PATH 
export PATH=/usr/local/bin/bcftools/bin:$PATH 
export PATH=/usr/local/bin/htslib/bin:$PATH 


# Set environment variable
samtools="/home/niccolo/Programs/samtools-1.11/samtools"
bwa="/home/niccolo/Programs/bwa-0.7.17/bwa"
picard="/home/niccolo/Programs/Picard/picard.jar"

log="../B.nigra_BSA_raw_reads_mapping.log"
log_C2="../B.nigra_BSA_raw_reads_mapping_C2.log"

ref_fasta="./Bnigra_NI100.v2.nuclear+plastids.fasta"
genome="NI100"

ref_fasta="./Bnigra_C2.v1.nuclear+plastids.fasta"
genome="C2"

ref_fasta="./Bnigra_sangam.v1.nuclear+plastids.fasta"
genome="sangam"

read_1="../D6_L4+L1_1.fq.gz"
read_2="../D6_L4+L1_2.fq.gz"

sample="D6_L4+L1"
dup="dup"

# 1. index reference genome
$bwa index $ref_fasta
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tindex $genome reference genome (with plastids)" >> $log

# 2.alignment (MEM)
$bwa mem -t 20 -M $ref_fasta $read_1 $read_2 | $samtools view -Sbh - | $samtools sort - -o $sample.$genome.srt.without.header.bam

echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tbwa mem $sample $(ls -t1 | head -n 1)" >> $log

	$samtools merge $sample.$genome.srt.without.header.bam $sample.$genome.srt.without.header.bam*
	rm -r $sample.$genome.srt.without.header.bam.0*
	
# 4. CHECK validity BAM file
# check https://broadinstitute.github.io/picard/explain-flags.html for SAM flags

## PICARD
# Check validity SAM/BAM file 
# for detail check (https://gatkforums.broadinstitute.org/gatk/discussion/7571/errors-in-sam-bam-files-can-be-diagnosed-with-validatesamfile)
java -jar ../../Programs/Picard/picard.jar ValidateSamFile \
	-I $sample.$genome.srt.without.header.bam \
	-MODE SUMMARY

# Add/replace read group in SAM file
java -jar $picard AddOrReplaceReadGroups \
  I=$sample.$genome.srt.without.header.bam \
  O=$sample.$genome.srt.bam 2>error.AddOrReplaceRG.$sample.$genome.txt \
  RGID=HWNCTDSXX.4 \
  RGLB=TruSeq3_PE \
  RGPL=illumina \
  RGPU=HWNCTDSXX.4 \
  RGSM=$sample # IMPORTANT! - CHANGE based on sample name 
	
# check if SAM/BAM file has been corrected (use output file from previous command)
$samtools view -H $sample.$genome.srt.bam | grep '^@RG'		

echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tSAM read group tag added/corrected to $sample.sam " >> ../$log
		
# 5. Inspect quality alignment
## SAMTOOLS flagstat
$samtools flagstat $sample.$genome.srt.$dup.bam

picard="/home/niccolo/Programs/Picard/picard.jar"

## PICARD
# Alignment summary metrics
java -jar $picard CollectAlignmentSummaryMetrics \
	R=../$ref_fasta \
	I=$sample.$genome.srt.$dup.bam \
	O=$sample.$genome.srt.$dup.alig.metrics.txt 2>error.$sample.$genome.srt.$dup.bam.align.metrics.txt
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tPicard AligmentSummaryMetrics for $kmer_set.$maxDiff.bam " >> ../$log


########################################################
#### 			BAM alignment filtering 			####

samtools="/home/niccolo/Programs/samtools-1.11/samtools"

fasta="./Bnigra_NI100.v2.nuclear+plastids.fasta"
genome="NI100"
bedtools="../../Programs/bedtools" 
sample="D2_L4"
plastids_bed="Bnigra_plastids.fasta.bed"
nuclear_bed="Bnigra_NI100.v2.genome.fasta.bed"

locus="hr.locus"

# 0. Index FASTA/BAM for visualization
samtools faidx $ref_fasta # -o $fasta.fai # -o specifiy output file
samtools index $sample.$genome.srt.bam


# 1. Extract alignment results based on MAPQ score
samtools view -q 5 $sample.$genome.srt.bam | awk -F "\t" '{print $3"\t"$5"\t"$4"\t"$10}' > output.sam
	
	# original settings from Clorentin Clot
	samtools view -q 2 sample.srt.bam.bam -L regions_of_interest.bed | awk -F ' ' '{print $3"\t"$5"\t"$4"\t"$10}' > output.txt
	
	
# 2. Extract reads that map to plastids (method 1 - see Clorentin example)
samtools view -b -L $nuclear_bed -U $sample.plastids_reads_method1.bam $sample.$genome.srt.bam > $sample.nuclear_reads_method1.bam # method 1
	samtools index D2_nuclear_reads.bam		# index BAM for visualization
	samtools index D2_plastids_reads.bam	# index BAM for visualization

	# Alternative methods to filter BAM file (method 2 - from Biostars)
	samtools view -h $sample.$genome.srt.bam | awk '{if($3 != "NC_029182.1_mitochondrion" && $3 != "NC_030450.1_chloroplast"){print $0}}' | samtools view -Sb - > D2_nuclear_reads_method2.bam # method 2
	
	$bedfile="Bnigra_plastids.fasta.bed" # create a bedfile for intersection
	$bedtools intersect -abam $sample.$genome.srt.bam -b $bedfile > D2_plastids_reads_method3.bam # method 3


# 3. Calculate average sequencing coverage
## SAMTOOLS
# average X coverage for covered regions (from https://www.biostars.org/p/356937/)
samtools depth $sample.$genome.srt.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' # this is average for all covered regions. include parameter "-a" to output all positions

# average X coverage genome-wide
samtools view -H $sample.$genome.srt.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' # compute the total size of the genome, sostitute this to NR in above command
	# D2_L4 genome="506386272" 
samtools depth $sample.$genome.srt.bam |  awk '{sum+=$3} END { print "Average = ",sum/506386272}' # this is average for all covered regions. include parameter "-a" to output all positions

# average seq coverage per nt position (to be exported for plotting, from https://www.biostars.org/p/104063/)
samtools depth -aa $sample.$genome.srt.bam > $sample.coverage
	### select area of interest
	left_border='5905000'	# left of marker M27
	right_border='5975000'	# right of BniB03g014150.2N
	awk '{OFS="\t"; if ($1 == "B3" && $2 > 5905000 && $2 < 5975000){print}}' $sample.coverage > $sample.$locus.coverage # select coverage for a particular locus, plot
	

# 4. Extact region from genome
reformat="/home/niccolo/Programs/bbmap/reformat.sh" # @Mary
genome="C2"
fasta="Bnigra_C2.v1.genome.fasta"

bedtools getfasta -fi $fasta -bed $genome.kmer.locus.bed -name -fo result.fasta
$reformat in=result.fasta out=$genome.kmer.locus.fasta fastawrap=60 # BBmap utily to wrap FASTA in 60 letters line

# 4. Extract region from GFF file
genome="C2"
locus="kmer.locus"

bedfile="C2.kmer.locus.bed"
GFFfile="Bnigra_C2.v1.genes.gff3" # NI100

bedtools intersect -a $GFFfile -b $bedfile > $genome.$locus.genes.gff3

# 4. Extract region from BAMfile 
#bamfile="D3_L4.N100.srt.dup.bam"
bedfile_path="/home/niccolo/Genomes"

sample="D1_L4"
	sample="D2_L4"
		sample="D3_L4"
			sample="D4_L4"
				sample="D5_L4"
					sample="D6_L4+L1"
		
genome="NI100"
type_bam="srt.dup"
locus="hr.locus"

bedtools intersect -abam $sample.$genome.$type_bam.bam -b $bedfile_path/$genome.$locus.bed > $sample.$genome.$type_bam.$locus.bam
$samtools index $sample.$genome.$type_bam.$locus.bam


# 5. Subset fasta files
reformat="/home/niccolo/Programs/bbmap/reformat.sh" # @Mary

fasta="/home/niccolo/Genomes/B.rapa/Chiifu_v3.0/Brapa_genome_v3.0_pep.fasta"

region="Pbc2.pep"
gene_list="genes_list.txt"

# 5.1 seqtk
seqtk subseq $fasta $gene_list > $region.subset.fasta
$reformat in=$region.subset.fasta out=genes.fasta ignorejunk=t fastawrap=60 overwrite=t # BBmap utily to wrap FASTA in 60 letters line

# 5.2 grep 
grep -w -A 2 -f  test.txt test.fa --no-group-separator # not tested


# 6. Make BED from REF fasta
# Convert genome .fasta file to .bed with chromosome as single entry
samtools="/home/niccolo/Programs/samtools-1.11/samtools"
fasta="Bnigra_sangam.v1.nuclear+plastids.fasta"

$samtools faidx $fasta -o $fasta.fai
awk 'BEGIN {FS="\t"};{print $1 FS "0" FS $2}' $fasta.fai > $fasta.bed # add "BEGIN {FS="\t"};" statement for .bed tab delimited (without quotes) 

	# Separate the genome in 1Mb intervals
	../../Programs/bedtools makewindows -b $fasta.bed -w 1000000 > $fasta.$interval.window.bed
	../../Programs/bedtools makewindows -b $fasta.bed -w 1000000 -s 10000 > $fasta.$interval.sliding.window.bed # sliding window



								#####################
								#### SNP calling ####

cd /B.nigra_BSA_analysis/gatk_Emma/gatk_C2 # @Mary
cd /B.nigra_BSA_analysis/gatk_Emma/gatk_NI100 # @Mary

cd /media/biosystematics/data/Niccolo/B.nigra_BSA_analysis/bwa_mapping_raw_reads	# @Emma
log="../B.nigra_BSA_raw_reads_mapping.log"

genome="NI100"
ref_fasta="Bnigra_NI100.v2.nuclear+plastids.fasta"
genome_interval="../Bnigra_NI100.v2.nuclear+plastids.fasta.bed"
locus_interval="NI100.hr.locus.bed"

genome="C2"
ref_fasta="Bnigra_C2.v1.nuclear+plastids.fasta"
locus_interval="C2.hr.locus.bed"


genome="sangam"
ref_fasta="Bnigra_sangam.v1.nuclear+plastids.fasta"
locus_interval="sangam.hr.locus.bed"


locus="hr.locus" 
dup="dup"

sample="D1_L4"




#####   		GATK best practices				#####
#### 											 #### - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
### 1.  Data pre-processing for variant discovery ### - https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery


## 1.1 Mark duplicate in BAM (PICARD)
java -jar /home/niccolo/Programs/Picard/picard.jar MarkDuplicates \
      I=$sample.$genome.srt.bam \
      O=$sample.$genome.srt.$dup.bam 2>error.$sample.$genome.$dup.MarDup.txt\
      M=$sample.$genome.MarDup.metrics.txt
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tMark duplicates with Picard on $bamfile" >> ../$log

	# Check quality of BAM afterwards
	samtools flagstat $sample.$genome.srt.$dup.bam

# 1.2 BaseRecalibrator ### NOT USED ###
gatk --java-options "-Xmx4g" BaseRecalibrator \
   -I my_reads.bam \
   -R reference.fasta \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table

####                                  ####
### 2. Variant calling (joint calling) ### - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Intro - joint calling workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants
#								 - https://gatk.broadinstitute.org/hc/en-us/articles/360035890411?id=3893

# export PATH=/usr/local/bin/samtools:$PATH 
# export PATH="/media/biosystematics/data/Niccolo/Programs/tabix-0.2.6/tabix:$PATH" # @Emma
# export PATH="/media/biosystematics/data/Niccolo/Programs/gatk-4.1.9.0:$PATH" # @Emma

export PATH="/home/niccolo/Programs/gatk-4.1.9.0:$PATH" # @Mary

# 2.0 Create index and dictionary for $ref_fasta and BAM files
$samtools faidx $ref_fasta # -o $fasta.fai 	# index ref genome # -o specifiy output file
gatk CreateSequenceDictionary -R $ref_fasta  # create ref_fasta.dict file

$samtools index $sample.$genome.srt.$dup.bam	# index BAM file

genome="sangam"
ref_fasta="Bnigra_sangam.v1.nuclear+plastids.fasta"
locus_interval="sangam.hr.locus.bed"

locus="hr.locus" 
dup="dup"


# 2.1 HaplotypeCaller - (SNP calling) - https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
gatk HaplotypeCaller \
	-R ../$ref_fasta \
	-I $sample.$genome.srt.$dup.bam \
	-O $sample.$genome.$dup.g.vcf.gz 2>error.$sample.$genome.$dup.HapCall.txt \
	-bamout $sample.$genome.$dup.HapCall.bam \
	-ERC GVCF \
	-L /home/niccolo/Genomes/$locus_interval \
	--pcr-indel-model NONE # NONE=PCR free library 
	# -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation # Use for allele-specific annotation (https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622)
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tHaplotype Caller on $bamfile to get $final gVCF file" >> $log

	# Inspect error report
	less +G error.$sample.$genome.$dup.HapCall.txt # +G open less and jump to the ened of the file
	
# 2.2 GenomicsDBImport - (merge GVCFs from multiple sample) - https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport#--intervals
gatk --java-options "-Xmx10g" GenomicsDBImport \
	-V D1_L4.$genome.$dup.g.vcf.gz \
	-V D2_L4.$genome.$dup.g.vcf.gz \
	-V D3_L4.$genome.$dup.g.vcf.gz \
	-V D4_L4.$genome.$dup.g.vcf.gz \
	-V D5_L4.$genome.$dup.g.vcf.gz \
	-V D6_L4+L1.$genome.$dup.g.vcf.gz \
	--genomicsdb-workspace-path my_GVCFs_database_BSA.$genome.$dup.$locus \
	-L /home/niccolo/Genomes/$locus_interval \
	--reader-threads 20 \
	2>error.log
	# -G StandardAnnotation -G AS_StandardAnnotation # Use for allele-specific annotation (https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622)
	
# 2.3 GenotypeGVCFs (joint SNP calling on merged GVCFs) - https://gatk.broadinstitute.org/hc/en-us/articles/360035889971--How-to-Consolidate-GVCFs-for-joint-calling-with-GenotypeGVCFs
gatk  --java-options "-Xmx10g" GenotypeGVCFs \
	-R ../$ref_fasta \
	-V gendb://my_GVCFs_database_BSA.$genome.$dup.$locus \
	-L /home/niccolo/Genomes/$locus_interval \
	-O BSA.$genome.$dup.$locus.vcf 2>error.BSA.$genome.$dup.$locus.GenotypeGVCFs.log
	# -G StandardAnnotation -G AS_StandardAnnotation # Use for allele-specific annotation (https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622)
echo -e $(date +"%Y-%m-%d"$" %H:%M:%S")"\tGenotypeGVCFs to create VCF file for $final" >> $log




# Variant filtering (VQSR/variant recalibration or hard-filtering) - https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
	# tutorial variant filtering
	# -AS to activale allele-specific mode in VQSR

locus="hr.locus"

# 2.4.0 Variant to Table
gatk --java-options "-Xmx10g" VariantsToTable \
-R $ref_fasta \
-V BSA.$genome.$dup.$locus.vcf \
-F CHROM -F POS -F REF -F ALT -F QUAL \
-GF AD -GF DP -GF GQ -GF PL \
-O BSA.$genome.$dup.$locus.vcf.csv 2>error


# 2.4.1 Variant selection
# select SNPs
gatk SelectVariants \
    -V BSA.$genome.$dup.$locus.vcf \
    -select-type SNP \
    -O snps.BSA.$genome.$dup.$locus.vcf 2>error.snps # .gz for archive file type
	
# select Indels
gatk SelectVariants \
    -V BSA.$genome.$dup.$locus.vcf \
    -select-type INDEL \
    -O indels.BSA.$genome.$dup.$locus.vcf 2>error.indel ## .gz for archive file type
	
# select mixed type (multiallelic SNPs + InDels)
gatk SelectVariants \
	-V BSA.$genome.$dup.$locus.vcf \
	-select-type MIXED \
	-O mixed.BSA.$genome.$dup.$locus.vcf 2>error.mixed ## .gz for archive file type

# 2.4.1 Variant filering (based on quality scores)
# snps
gatk VariantFiltration \
	-V snps.BSA.$genome.$dup.$locus.vcf \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O snps.filtered.BSA.$genome.$dup.$locus.vcf 2>error

# indels
gatk VariantFiltration \
	-V indels.BSA.$genome.$dup.$locus.vcf \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 200.0" --filter-name "FS200" \
	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	-O indels.filtered.BSA.$genome.$dup.$locus.vcf 2>error

# mixed
gatk VariantFiltration \
	-V mixed.BSA.$genome.$dup.$locus.vcf \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 200.0" --filter-name "FS200" \
	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	-O mixed.filtered.BSA.$genome.$dup.$locus.vcf
	
	# valuate effect of filtering
	cat snps.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l
	cat snps.filtered.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l

	cat indels.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l
	cat indels.filtered.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l
	
	cat mixed.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l
	cat mixed.filtered.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l


##################################################
### 3. Variant filtering based on samples/genotype


## 1. combine all VCFs files (snps, indels, mixed) -- GATK
java -jar ../../../Programs/Picard/picard.jar MergeVcfs \
	I=snps.filtered.BSA.$genome.$dup.$locus.vcf \
	I=indels.filtered.BSA.$genome.$dup.$locus.vcf \
	I=mixed.filtered.BSA.$genome.$dup.$locus.vcf \
	O=variants.filtered.BSA.$genome.$dup.$locus.vcf 2>error
	
	# check size of final VCF (not necessary)
	cat snps.filtered.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l
	cat indels.filtered.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l
	cat mixed.filtered.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l
	cat variants.qualfilter.BSA.$genome.$dup.$locus.vcf | grep -v '^#' | wc -l

## 2. filter sites based on genetic data (Hom or Het occurence through specific samples)
vcf='variants.qualfilter.BSA.C2.dup.hr.locus.vcf' # CHANGE file name

cat $vcf | awk '{OFS="\t"; if (!(substr($11,1,3) == substr($12,1,3))) {print}}' > $vcf.f1.vcf # remove variants that are equal between parents (sample $10 and $11)
cat $vcf.f1.vcf | awk '{OFS="\t"; if (!(substr($12,1,1) == substr($12,3,1))) {print}}' > $vcf.f2.vcf # select only heterozygous in HR+ sample
cat $vcf.f2.vcf | awk '{OFS="\t"; if (substr($11,1,1) == substr($11,3,1)) {print}}' > final.$vcf # select only homozygous in HR- sample
	

	# inspect filtering results
	cat final.$vcf | grep -v '^#' | less
	# count variants (to compare effect of filtering)
	cat $vcf | grep -v '^#' | wc -l
	cat $vcf.f1.vcf | grep -v '^#' | wc -l
	cat $vcf.f2.vcf | grep -v '^#' | wc -l
	cat final.$vcf | grep -v '^#' | wc -l
	
	
	# remove temp files
	rm $vcf.f1.vcf
	rm $vcf.f2.vcf
	# get header from original VCF
	cat $vcf | grep '^#' >> genfilter.$vcf  # CHANGE to rename files
	# append filtered VCF to header
	cat final.$vcf >> genfilter.$vcf

	
## 3. filter based on position locus
chr="B3"				# NI100
start_locus="5700000"
end_locus="7100000"

chr="B3"				# C2
start_locus="6000000"
end_locus="7520000"

chr="B3"				# sangam TO FIX
start_locus="6000000"
end_locus="7520000"

locus="hr.locus" 

vcftools --vcf --chr $chr --from-bp $start_locus --to-bp $end_locus --out goi.final.variants.filtered.BSA.$genome.$locus.vcf

