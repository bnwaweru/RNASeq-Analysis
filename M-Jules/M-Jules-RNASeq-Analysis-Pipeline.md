J M RNASeq Data Analysis Pipeline
================
Bernice Waweru
Wed 31, Mar 2021

-   [Background of Experiment](#background-of-experiment)
    -   [Objective](#objective)
    -   [Design](#design)
-   [Overview of Pipeline](#overview-of-pipeline)
-   [Analysis of Data](#analysis-of-data)
    -   [1. Retrieve the data and check the
        quality](#retrieve-the-data-and-check-the-quality)
    -   [2. Mapping the reads to the
        genome](#mapping-the-reads-to-the-genome)
        -   [Summary of star log files with mapping
            statistics](#summary-of-star-log-files-with-mapping-statistics)

## Background of Experiment

#### Objective

The goal of the experiment was to find whether there are genes that are
differentially expressed between cultivars of grass species under study
in response to aluminum present in the soil.

#### Design

There were 2 cultivars (Basilisk and Marandu) each with a treatment
(aluminium /lime) and control. For each treatment, there were 4
biological samples per cultivar with 3 technical replications per
sample. For the biological sample, four root samples were collected from
one plot of each of the cultivars (Basilisk, Marandu) one treated with
lime, and one from a plot not treated with lime i.e the control. So for
each cultivar, we had 8 biological samples. From each of the eight
samples, three different nucleic acid extraction was done, i.e three lab
technical replicates, leading to 24 samples per cultivar, 12 from a plot
treated with lime and 12 from a control. Truseq Illumina kit was used to
prepare the libraries. A total of 48 samples were submitted for
sequencing, which were sequenced in two pools.

## Overview of Pipeline

RNA-seq experiments are performed with an aim to comprehend
transcriptomic changes in organisms in response to a certain treatment.

In general, analysis of RNASeq data follows several steps, mainly
depending on whether a genome and annotation is available for the
organism under study.

1.  Retrieving the raw data
2.  Check the quality of the raw data, if any cleaning is required,
    clean up the reads before proceeding
3.  Generate the transcript using a genome-guided approach or a *de
    novo* approach
4.  Estimate the abundance of features within your data i.e count the
    number of reads assigned to each feature observed in the data
5.  Differential gene expression to identify which genes are expressed
    differently between your samples depending on the treatment

## Analysis of Data

For this analysis, there is a genome and annotation available, hence we
will follow a genome-guided approach to in analyzing our data.

### 1. Retrieve the data and check the quality

The sequencing was done with the sequencing company Novagene. The data
was retrieved from an ftp site provide by the sequencing company and
stores on the [ilri cluster](https://hpc.ilri.cgiar.org/) under the
fellow’s home directory.

The first step wis to check the quality of the data. We do this using
[Fastqc](https://github.com/s-andrews/FastQC), a tool designed to parse
through high throughput sequencing data and check for potential errors.
After the report for each file are generated, we use another tool,
[multiqc](https://github.com/ewels/MultiQC), to aggregate the reports
into ine single report.

There are 48 sample files that were sequenced in paired-end, hence we
have a total of 96 fastq files. We check the quality of the files using
a script that loops over all the files. To do this we submit a [batch
script](https://hpc.ilri.cgiar.org/using-slurm) to the slurm scheduler
on the cluster.

    #!/bin/bash
    #SBATCH -p batch
    #SBATCH -w compute04
    #SBATCH -n 4
    #SBATCH -J fastqc_trim


    #============= load tools ======================================

    module load fastqc/0.11.7 multiqc/1.6 trimmomatic/0.38

    #===============================================================
    #========== define working directories =========================
    #===============================================================

    #the original data directory
    raw_data='/home/bngina/Fellows/mjules/orig_data/raw_data'

    #fastqc results directory

    fastqc_out='/home/bngina/Fellows/mjules/fastqc_dir/raw_data'

    #directory for multiqc output

    multiqc_out='/home/bngina/Fellows/mjules/multiqc_out/raw_data'

    #================================================================
    #=== use fastqc to check the quality of the raw fastq files =====
    #================================================================

    #for file in $(ls ${raw_data}/*.fq.gz) ;\
    #do echo $file ;\
    #fastqc -t 5 ${file} -o ${fastqc_out}/;\
    #done

    #================================================================
    #============== use multiqc to aggregate fastqc reports =========
    #================================================================

    multiqc ${fastqc_out}/ -o ${multiqc_out}

Below is a look at the quality bars and the adapter content;

<img src="./embedded-images/raw-mean-quality-scores.PNG" style="width:30.0%" alt="Quality scores" />
<img src="./embedded-images/raw-adapter-content.PNG" style="width:30.0%" alt="Adapter content" />

All the samples have quality score that are above 30, all of them within
the green zone of the graph, meaning that thay are pretty good. However
the adapter content within the sample files seems to be high starting at
around 55bp within the reads. The adpters need to be removed from the
reads, we use a tool called
[trimmomatic](https://github.com/timflutre/trimmomatic) to do this. It’s
not a very recent tool hence its abit slow on large files, but it works
very well in performing the trimming.

As we have 96 fastq files to clean of adapters, we use a loop command
once again in a batch script to make use of the cluster resources. After
the trimming we re-do the fastqc to check that our we have cleaned our
reads of the adapter sequences.

Continuing within the batch script we had set up initially;


    # ======== load trimmomatic tools used for trimming ======================================================

    module load trimmomatic/0.38

    # ============= directory to store output files from trimmomatic =========================================

    trim_out='/var/scratch/waweru/mjules/trimmomatic_out'

    # =========================================================================================================
    # ============= run trimmomatic ===========================================================================
    # =========================================================================================================

    for line in $(ls ${raw_data}/S21_23*1.fq.gz) ;\
    do echo ${line} ;\
    out_file=$(echo $line |cut -f 8 -d "/" | sed 's/_1.fq.gz//g') ;\
    echo $out_file ;\
    trimmomatic PE -threads 4 -trimlog ${trim_out}/mjules_trimmomatic_log_file.txt\
     -basein ${line} -baseout ${trim_out}/${out_file}_trmd.fq \
    ILLUMINACLIP:/home/bngina/Bambara/all_adapeters.fa:2:28:4 ;\
    cd /var/scratch/waweru/mjules/trimmomatic_out/ ;\
    rm *U.fq ;\
    done


    # ============= re-do the fastqc to check if the adapters have been removed and the sequence length after trimming ========
    # ==========================================================================================================================

    #=== use fastqc to check the quality of the trimmed fastq files =====
    #====================================================================


    #=========== directories to store results for the trimmed data =========
    #=======================================================================

    # ========= from fastqc ==========

    fastqc_out='/var/scratch/waweru/mjules/fastqc_out'

    # ========= from multiqc =========

    multiqc_out='/home/bngina/Fellows/mjules/multiqc_out/trmd_data'

We then take a look at the adapter content after trimming;

<img src="./embedded-images/trimmed-quality-scores.PNG" style="width:50.0%" alt="trimmed quality scores" /><img src="./embedded-images/trimmed-adater-content.PNG" style="width:50.0%" alt="trimmed adapter content" />

We see that the quality now for most sequences is closer to 35, and the
adapter content is less than 0.1% in the samples, hence our trimming
worked well.

### 2. Mapping the reads to the genome

Several mapping algorithms/tools ar available for mapping RNASeq reads
to a genome. It is however important to use an aligner that is able
handle reads spanning splice junctions, i.e [splice
aware](https://www.biostars.org/p/175454/) aligners like
[STAR](https://github.com/alexdobin/STAR) and
[HISAT](https://github.com/infphilo/hisat).

In this pipeline, we use *STAR* aligner to find where the reads in the
sample files originated from in the genome. Before beginning the
genome-guided approach to, one has to find choose a well annotated
genome from the same species or aclosely realted species of the organism
under study. An associated [genome feature
file](https://mblab.wustl.edu/GTF22.html)(`*.gff*`), indicates to the
aligner where the transcripts are located on the genome.

[This
page](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)
gives a nice overview of the star mapping algorithm. To use STAR we
first build an index of the genome, and then map the files in a loop one
after the other using the cluster resources.

    #!/bin/bash
    #SBATCH -p batch
    #SBATCH -w compute04
    #SBATCH -n 7
    #SBATCH -J star_mapping
    #SBATCH -e /home/bngina/Fellows/mjules/batch_logs/indexing.%N.%J.err
    #SBATCH -o /home/bngina/Fellows/mjules/batch_logs/indexing.%N.%J.out

    #============= load tools ======================================

    module load star/2.7.1a


    #===============================================================
    #========== define working directories =========================
    #===============================================================


    # ===== dir to store ref genome and gff file ========================

    ref_dir='/home/bngina/Fellows/mjules/ref'

    # ==== directory to store the index files =========================

    mkdir -p /home/bngina/Fellows/mjules/ref/geno_dir

    geno_dir='/home/bngina/Fellows/mjules/ref/geno_dir'


    # ===== once the data is dwonloaded we unzip the files before supplying them as input to star
    # ===========================================================================================
    cd ${ref_dir}

    gzip -d *.gz

    # ===== index the files =============================================
    # ===================================================================

    STAR --runThreadN 7 \
     --runMode genomeGenerate \
     --genomeDir ${geno_dir}\
     --limitGenomeGenerateRAM 80000000000 \
     --genomeFastaFiles Bruz_v1.min1000.fa\
     --sjdbGTFfile Bruz_v2.gff3 \
     --sjdbOverhang 100

Once the indexing is complete, we move on to map the files within our
batch script.


    # ===================================================================
    # ===== Mapping the trimmed fastq files to the indexed genome =======
    # ===================================================================


    # ===== dir with the fastq files ====================================

    fastq_dir='/var/scratch/waweru/mjules/trimmomatic_out'

    # ===== dir for the output files ====================================

    mkdir -p /var/scratch/waweru/mjules/star_out

    # ===== move into the directory where we want the output files to be stored

    cd /var/scratch/waweru/mjules/star_out/

    # ===== the mapping ==================

    for R1 in ${fastq_dir}/*trmd_1P.fq ;\
    do echo ${R1};\
    R2=$(echo ${R1} | sed 's/trmd_1P.fq/trmd_2P.fq/g');\
    echo ${R2};\
    out_name=$(echo ${R1} | cut -f 7 -d "/" | sed 's/trmd_1P.fq//g');\
    echo ${out_name};\
    STAR --runThreadN 6\
     --genomeDir ${geno_dir}\
     --readFilesIn ${R1} ${R2}\
     --outFileNamePrefix ${out_name}\
     --outSAMtype BAM SortedByCoordinate\
     --outBAMsortingThreadN 1 ;\
    done

The mapping generates alignment files in `bam` format as we requested in
the command, sorted by coordinate. Star also give nice log files for
each sample file that give a summary of the mapping. We retrieve the log
files from the cluster and store them on out computers, then process the
files within R to extract the information about how well our samples
mapped to the genome.

#### Summary of star log files with mapping statistics

``` r
setwd("C:/Users/BWaweru/OneDrive - CGIAR/Documents/Fellows/Jules_Mutabazi/RWD_Git/RNASeq-Analysis/M-Jules/")

# ===== the package we will need ===============

require(magrittr)
require(tidyverse)
require(kableExtra)

# ===== path to the star log files we retrieved from the cluster
path <- ("./data_raw/star_log_stats/")

# ===== files names of the logs

files <- list.files(path, pattern = ".final.out")

# ===== check the files

files # we have 15 files
```

    ##  [1] "S11_10_Log.final.out" "S11_11_Log.final.out" "S11_9_Log.final.out" 
    ##  [4] "S13_16_Log.final.out" "S13_18_Log.final.out" "S13_19_Log.final.out"
    ##  [7] "S14_20_Log.final.out" "S14_21_Log.final.out" "S14_22_Log.final.out"
    ## [10] "S15_2_Log.final.out"  "S15_4_Log.final.out"  "S15_5_Log.final.out" 
    ## [13] "S16_12_Log.final.out" "S16_6_Log.final.out"  "S16_7_Log.final.out"

``` r
# ===== function to extract the information we need from the log files

merge_out <- function (files) {
  df <- df <- read.delim(paste0(path, files[1]), header= F) %>% 
    filter(grepl("Number of input reads", V1) |
             grepl("Uniquely mapped reads", V1) |
             grepl("Uniquely mapped reads %", V1) |
             grepl("Average mapped length", V1) |
             grepl("Number of splices: Total", V1) |
             grepl(" Number of reads mapped to multiple loci", V1) |
             grepl("% of reads mapped to multiple loci", V1) |
             grepl("Number of reads unmapped: too many mismatches", V1) |
             grepl(" % of reads unmapped: too many mismatches", V1) |
             grepl("Number of reads unmapped: too short", V1) |
             grepl("% of reads unmapped: too short", V1) |
             grepl("Number of reads unmapped: other", V1) |
             grepl(" % of reads unmapped: other", V1) |
             grepl("Number of chimeric reads", V1) |
             grepl("% of chimeric reads", V1)) %>% 
    set_names("Var", "value")
}

#use the function to put the results in a list 

results <- lapply(files, merge_out)

# ===== next we want to generate a dataframe with the file names as a data frame, 
# ===== then we add columns of extracted information to build one large data frame

star_df <- as.data.frame(files)

star_df # the dataframe with one column of file names
```

    ##                   files
    ## 1  S11_10_Log.final.out
    ## 2  S11_11_Log.final.out
    ## 3   S11_9_Log.final.out
    ## 4  S13_16_Log.final.out
    ## 5  S13_18_Log.final.out
    ## 6  S13_19_Log.final.out
    ## 7  S14_20_Log.final.out
    ## 8  S14_21_Log.final.out
    ## 9  S14_22_Log.final.out
    ## 10  S15_2_Log.final.out
    ## 11  S15_4_Log.final.out
    ## 12  S15_5_Log.final.out
    ## 13 S16_12_Log.final.out
    ## 14  S16_6_Log.final.out
    ## 15  S16_7_Log.final.out

``` r
# add to the newly created data frame number of mapped reads
for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[1] -> star_df[i, "Number of input reads"]
}
  
# add to the newly created data frame number of uniquely mapped reads

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[2] -> star_df[i, "Uniquely mapped reads number"]
}


# add to the newly created data frame number of uniquely mapped reads

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[3] -> star_df[i, "% of Uniquely mapped reads"]
}
 
# add average mapped length

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[4] -> star_df[i, "Average mapped length"]
}

# add total number of splices

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[5] -> star_df[i, "Number of splices: Total"]
}

# add Number of reads mapped to multiple loci

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[6] -> star_df[i, "Number of reads mapped to multiple loci"]
}

# add % of reads mapped to multiple loci

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[7] -> star_df[i, "% of reads mapped to multiple loci"]
}

# add Number of reads unmapped: too many mismatches

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[8] -> star_df[i, "Number of reads unmapped: too many mismatches"]
}

# add % of reads unmapped: too many mismatches

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[9] -> star_df[i, "% of reads unmapped: too many mismatches"]
}

# add Number of reads unmapped: too short

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[10] -> star_df[i, "Number of reads unmapped: too short"]
}

# add % of reads unmapped: too short

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[11] -> star_df[i, "% of reads unmapped: too short"]
}

# add Number of reads unmapped: other

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[12] -> star_df[i, "Number of reads unmapped: other"]
}

# add % of reads unmapped: other

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[13] -> star_df[i, "% of reads unmapped: other"]
}

# add Number of chimeric reads

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[14] -> star_df[i, "Number of chimeric reads"]
}

# add % of chimeric reads

for (i in 1:15) {
  file_num <- i
  temp_df <- as.data.frame(results[file_num])
  temp_df$value[15] -> star_df[i, "% of chimeric reads"]
}


# ===== lets have a look at the generated table

kable(star_df, caption = "A summary of mapping statistics generated by STAR aligner")
```

<table>
<caption>
A summary of mapping statistics generated by STAR aligner
</caption>
<thead>
<tr>
<th style="text-align:left;">
files
</th>
<th style="text-align:left;">
Number of input reads
</th>
<th style="text-align:left;">
Uniquely mapped reads number
</th>
<th style="text-align:left;">
% of Uniquely mapped reads
</th>
<th style="text-align:left;">
Average mapped length
</th>
<th style="text-align:left;">
Number of splices: Total
</th>
<th style="text-align:left;">
Number of reads mapped to multiple loci
</th>
<th style="text-align:left;">
% of reads mapped to multiple loci
</th>
<th style="text-align:left;">
Number of reads unmapped: too many mismatches
</th>
<th style="text-align:left;">
% of reads unmapped: too many mismatches
</th>
<th style="text-align:left;">
Number of reads unmapped: too short
</th>
<th style="text-align:left;">
% of reads unmapped: too short
</th>
<th style="text-align:left;">
Number of reads unmapped: other
</th>
<th style="text-align:left;">
% of reads unmapped: other
</th>
<th style="text-align:left;">
Number of chimeric reads
</th>
<th style="text-align:left;">
% of chimeric reads
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
S11\_10\_Log.final.out
</td>
<td style="text-align:left;">
148198381
</td>
<td style="text-align:left;">
40477546
</td>
<td style="text-align:left;">
27.31%
</td>
<td style="text-align:left;">
187.98
</td>
<td style="text-align:left;">
16075157
</td>
<td style="text-align:left;">
1774133
</td>
<td style="text-align:left;">
1.20%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
105794013
</td>
<td style="text-align:left;">
71.39%
</td>
<td style="text-align:left;">
146863
</td>
<td style="text-align:left;">
0.10%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S11\_11\_Log.final.out
</td>
<td style="text-align:left;">
128738307
</td>
<td style="text-align:left;">
36381530
</td>
<td style="text-align:left;">
28.26%
</td>
<td style="text-align:left;">
187.45
</td>
<td style="text-align:left;">
14320691
</td>
<td style="text-align:left;">
1552910
</td>
<td style="text-align:left;">
1.21%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
90625497
</td>
<td style="text-align:left;">
70.40%
</td>
<td style="text-align:left;">
172855
</td>
<td style="text-align:left;">
0.13%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S11\_9\_Log.final.out
</td>
<td style="text-align:left;">
162366827
</td>
<td style="text-align:left;">
50923907
</td>
<td style="text-align:left;">
31.36%
</td>
<td style="text-align:left;">
180.85
</td>
<td style="text-align:left;">
18040026
</td>
<td style="text-align:left;">
2502676
</td>
<td style="text-align:left;">
1.54%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
108712152
</td>
<td style="text-align:left;">
66.95%
</td>
<td style="text-align:left;">
220954
</td>
<td style="text-align:left;">
0.14%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S13\_16\_Log.final.out
</td>
<td style="text-align:left;">
105013067
</td>
<td style="text-align:left;">
25131597
</td>
<td style="text-align:left;">
23.93%
</td>
<td style="text-align:left;">
200.20
</td>
<td style="text-align:left;">
9598698
</td>
<td style="text-align:left;">
1462576
</td>
<td style="text-align:left;">
1.39%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
78347699
</td>
<td style="text-align:left;">
74.61%
</td>
<td style="text-align:left;">
66404
</td>
<td style="text-align:left;">
0.06%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S13\_18\_Log.final.out
</td>
<td style="text-align:left;">
124999927
</td>
<td style="text-align:left;">
34366932
</td>
<td style="text-align:left;">
27.49%
</td>
<td style="text-align:left;">
209.84
</td>
<td style="text-align:left;">
13411811
</td>
<td style="text-align:left;">
2135647
</td>
<td style="text-align:left;">
1.71%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
88418533
</td>
<td style="text-align:left;">
70.73%
</td>
<td style="text-align:left;">
73485
</td>
<td style="text-align:left;">
0.06%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S13\_19\_Log.final.out
</td>
<td style="text-align:left;">
125748113
</td>
<td style="text-align:left;">
33838748
</td>
<td style="text-align:left;">
26.91%
</td>
<td style="text-align:left;">
193.28
</td>
<td style="text-align:left;">
12120907
</td>
<td style="text-align:left;">
2230599
</td>
<td style="text-align:left;">
1.77%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
89322723
</td>
<td style="text-align:left;">
71.03%
</td>
<td style="text-align:left;">
347855
</td>
<td style="text-align:left;">
0.28%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S14\_20\_Log.final.out
</td>
<td style="text-align:left;">
171696835
</td>
<td style="text-align:left;">
47929051
</td>
<td style="text-align:left;">
27.91%
</td>
<td style="text-align:left;">
183.03
</td>
<td style="text-align:left;">
17874994
</td>
<td style="text-align:left;">
2282105
</td>
<td style="text-align:left;">
1.33%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
121216607
</td>
<td style="text-align:left;">
70.60%
</td>
<td style="text-align:left;">
259901
</td>
<td style="text-align:left;">
0.15%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S14\_21\_Log.final.out
</td>
<td style="text-align:left;">
127785864
</td>
<td style="text-align:left;">
81582634
</td>
<td style="text-align:left;">
63.84%
</td>
<td style="text-align:left;">
189.75
</td>
<td style="text-align:left;">
26439168
</td>
<td style="text-align:left;">
5731119
</td>
<td style="text-align:left;">
4.48%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
40265626
</td>
<td style="text-align:left;">
31.51%
</td>
<td style="text-align:left;">
192403
</td>
<td style="text-align:left;">
0.15%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S14\_22\_Log.final.out
</td>
<td style="text-align:left;">
148425510
</td>
<td style="text-align:left;">
40796711
</td>
<td style="text-align:left;">
27.49%
</td>
<td style="text-align:left;">
196.43
</td>
<td style="text-align:left;">
17760347
</td>
<td style="text-align:left;">
1683542
</td>
<td style="text-align:left;">
1.13%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
105790395
</td>
<td style="text-align:left;">
71.28%
</td>
<td style="text-align:left;">
149233
</td>
<td style="text-align:left;">
0.10%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S15\_2\_Log.final.out
</td>
<td style="text-align:left;">
83594420
</td>
<td style="text-align:left;">
29312074
</td>
<td style="text-align:left;">
35.06%
</td>
<td style="text-align:left;">
205.76
</td>
<td style="text-align:left;">
7121229
</td>
<td style="text-align:left;">
4300538
</td>
<td style="text-align:left;">
5.14%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
49927333
</td>
<td style="text-align:left;">
59.73%
</td>
<td style="text-align:left;">
22219
</td>
<td style="text-align:left;">
0.03%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S15\_4\_Log.final.out
</td>
<td style="text-align:left;">
72511685
</td>
<td style="text-align:left;">
27785374
</td>
<td style="text-align:left;">
38.32%
</td>
<td style="text-align:left;">
189.03
</td>
<td style="text-align:left;">
5250200
</td>
<td style="text-align:left;">
5052474
</td>
<td style="text-align:left;">
6.97%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
39541621
</td>
<td style="text-align:left;">
54.53%
</td>
<td style="text-align:left;">
73419
</td>
<td style="text-align:left;">
0.10%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S15\_5\_Log.final.out
</td>
<td style="text-align:left;">
67732345
</td>
<td style="text-align:left;">
23022180
</td>
<td style="text-align:left;">
33.99%
</td>
<td style="text-align:left;">
183.26
</td>
<td style="text-align:left;">
7376733
</td>
<td style="text-align:left;">
1859796
</td>
<td style="text-align:left;">
2.75%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
42756211
</td>
<td style="text-align:left;">
63.13%
</td>
<td style="text-align:left;">
79168
</td>
<td style="text-align:left;">
0.12%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S16\_12\_Log.final.out
</td>
<td style="text-align:left;">
87666853
</td>
<td style="text-align:left;">
30941272
</td>
<td style="text-align:left;">
35.29%
</td>
<td style="text-align:left;">
201.48
</td>
<td style="text-align:left;">
11450824
</td>
<td style="text-align:left;">
2202614
</td>
<td style="text-align:left;">
2.51%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
54455287
</td>
<td style="text-align:left;">
62.12%
</td>
<td style="text-align:left;">
53074
</td>
<td style="text-align:left;">
0.06%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S16\_6\_Log.final.out
</td>
<td style="text-align:left;">
153276603
</td>
<td style="text-align:left;">
90848194
</td>
<td style="text-align:left;">
59.27%
</td>
<td style="text-align:left;">
187.92
</td>
<td style="text-align:left;">
32648851
</td>
<td style="text-align:left;">
4135798
</td>
<td style="text-align:left;">
2.70%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
58079802
</td>
<td style="text-align:left;">
37.89%
</td>
<td style="text-align:left;">
196952
</td>
<td style="text-align:left;">
0.13%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
<tr>
<td style="text-align:left;">
S16\_7\_Log.final.out
</td>
<td style="text-align:left;">
138475284
</td>
<td style="text-align:left;">
92828377
</td>
<td style="text-align:left;">
67.04%
</td>
<td style="text-align:left;">
189.03
</td>
<td style="text-align:left;">
34199289
</td>
<td style="text-align:left;">
4242917
</td>
<td style="text-align:left;">
3.06%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
<td style="text-align:left;">
41207987
</td>
<td style="text-align:left;">
29.76%
</td>
<td style="text-align:left;">
179052
</td>
<td style="text-align:left;">
0.13%
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.00%
</td>
</tr>
</tbody>
</table>

``` r
# ===== save dataframe into a csv file

write.csv(star_df, file = "./results/jm_star_mapping_stats.csv", quote = FALSE )
```

Looking at the results, we observe that many of the samples had very low
percentages of reads uniquely mapping to on e location on the genome.
But we also observe that a lot of these reads were not mapped because
they were too short. Maybe there is parameter we can add to instruct
star to also use the short reads?

From the [star
manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf),
the parameter that might have an effect on read length is specified
while indexing the genome;

> `--sjdbOverhang` specifies the length of the genomic sequence around
> the annotated junction to be used in constructing the splice junctions
> database. Ideally, this length should be equal to the ReadLength-1,
> where ReadLength is the length of the reads. For instance, for
> Illumina 2x100b paired-end reads, the ideal value is 100-1=99. In case
> of reads of varying length, the ideal value is max(ReadLength)-1. In
> most cases, the default value of 100 will work as well as the ideal
> value.

We specified this read length size as 100, which is the default. Maybe
we can go a bit lower to include more shorter reads? We can test this on
one file with a very low mapping rate and see if the results improve.
