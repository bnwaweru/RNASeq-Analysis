Abudance Counting
================
Bernice Waweru
Mon 12, Apr 2021

-   [Session information](#session-information)

The next thing we need to do is count how many reads are assigned to
each feature (our interest here being genes) found with the samples, i.e
find the abundance of genes within the samples.

[DESEq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
expects estimated count data as obtained, i.e data not normalized, in
the form of a matrix of integer values.

After mapping reads to trimmed reads to the genome, the count of reads
assigned to identified features/genes can be counted/estimated using
several tools including [Subread](http://bioinf.wehi.edu.au/subread/)
and [HTseq](https://htseq.readthedocs.io/en/master/). We used the
[`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/) read
summarization program within Subread.

`featureCounts` counts mapped reads for genomic features such as genes,
exons among others. The `featureCounts` program uses the gene\_id
attribute available in the [GTF](https://mblab.wustl.edu/GTF22.html)
format annotation (or the GeneID column in the SAF format annotation) to
group features into meta-features, ie. features belonging to the same
meta-feature have the same gene identifier.

To use it, we supply the annotation file and the bam files from mapping.
*featureCounts* supports multi-threading, hence we can specify the
number of cpus, to use. We also tell it which feature should be used to
group features using the `-t` flag. `-p` flag specifies that our reads
are paired, hence fragments are counted. `-g` specifies which attribute
is the gene identifier, by default this is *gene\_id*. A quick look at
the annotation file;

    cat ref_genome/pasi3.clean.gff3 | head -n 50

The result;

    LG1     AUGUSTUS        gene    2195    3661    0.01    -       .       ID=g1
    LG1     AUGUSTUS        transcript      2195    3661    0.01    -       .       ID=g1.t1;Parent=g1
    LG1     AUGUSTUS        transcription_end_site  2195    2195    .       -       .       Parent=g1.t1
    LG1     AUGUSTUS        exon    2195    2627    .       -       .       Parent=g1.t1
    LG1     AUGUSTUS        stop_codon      2292    2294    .       -       0       Parent=g1.t1
    LG1     AUGUSTUS        intron  2628    3264    0.05    -       .       Parent=g1.t1
    LG1     AUGUSTUS        CDS     2292    2627    0.04    -       0       ID=g1.t1.cds;Parent=g1.t1
    LG1     AUGUSTUS        CDS     3265    3372    0.09    -       0       ID=g1.t1.cds;Parent=g1.t1
    LG1     AUGUSTUS        exon    3265    3661    .       -       .       Parent=g1.t1
    LG1     AUGUSTUS        start_codon     3370    3372    .       -       0       Parent=g1.t1
    LG1     AUGUSTUS        transcription_start_site        3661    3661    .       -       .       Parent=g1.t1
    LG1     AUGUSTUS        gene    14193   16401   0.09    -       .       ID=g2
    LG1     AUGUSTUS        transcript      14193   16401   0.09    -       .       ID=g2.t1;Parent=g2
    LG1     AUGUSTUS        transcription_end_site  14193   14193   .       -       .       Parent=g2.t1
    LG1     AUGUSTUS        exon    14193   15024   .       -       .       Parent=g2.t1
    LG1     AUGUSTUS        stop_codon      14860   14862   .       -       0       Parent=g2.t1
    LG1     AUGUSTUS        intron  15025   16254   0.94    -       .       Parent=g2.t1
    LG1     AUGUSTUS        CDS     14860   15024   0.34    -       0       ID=g2.t1.cds;Parent=g2.t1
    LG1     AUGUSTUS        CDS     16255   16365   0.78    -       0       ID=g2.t1.cds;Parent=g2.t1
    LG1     AUGUSTUS        exon    16255   16401   .       -       .       Parent=g2.t1
    LG1     AUGUSTUS        start_codon     16363   16365   .       -       0       Parent=g2.t1
    LG1     AUGUSTUS        transcription_start_site        16401   16401   .       -       .       Parent=g2.t1

The gene identifier here is `ID=` so we use this in our script.


    # ===== count feature using featureCount tool in subread tool

    module load subread/1.6.2

    # ===== loop over the files bam files generated from mapping with star

    bam_dir='/var/scratch/waweru/linly_banda/star_out'

    for bam in ${bam_dir}/*Aligned.sortedByCoord.out.bam
    do echo ${bam}
    out_name=$(echo ${bam} | cut -f 7 -d "/" | sed 's/Aligned.sortedByCoord.out.bam//')
    echo ${out_name}
    featureCounts -T 10 \
     -p \
     -a ${ref_dir}/pasi3.clean.gff3 \
     -t gene \
     -g ID \
     -o ${bam_dir}/${out_name}counts_ID.txt ${bam}
    done

Looks like by specifying that the gene identifier as `ID` we get the
tabular result expected. Looking at the top of one of the files from
feature counts;

    cat /var/scratch/waweru/linly_banda/star_out/NASPOT11_1_counts_ID.txt | head

We see;

    # Program:featureCounts v1.6.2; Command:"featureCounts" "-T" "10" "-p" "-a" "/home/bngina/Fellows/linly_banda/ref_genome/pasi3.clean.gff3" "-t" "gene" "-g" "ID" "-o" "/var/scratch/waweru/linly_banda/star_out/NASPOT11_1_counts_ID.txt" "/var/scratch/waweru/linly_banda/star_out/NASPOT11_1_Aligned.sortedByCoord.out.bam"
    Geneid  Chr     Start   End     Strand  Length  /var/scratch/waweru/linly_banda/star_out/NASPOT11_1_Aligned.sortedByCoord.out.bam
    g1      LG1     2195    3661    -       1467    0
    g2      LG1     14193   16401   -       2209    5
    g3      LG1     16841   24171   -       7331    10
    g4      LG1     30666   32876   +       2211    0
    g5      LG1     33053   34101   -       1049    0
    g6      LG1     38933   44551   -       5619    23
    g7      LG1     45865   50081   -       4217    26
    g8      LG1     50463   53101   -       2639    0

The file has a commented header file with the details of the command we
gave, and a table with the genes within the annotation fie, the last
column being a count of how many reads were observed for each gene in
the sample file named

This is great, we can now work on how to incorporate this to work to
analyse for differentially expressed genes. We will use the R package
**DESEq2**.

#### Session information

Details of packages used for the work flow

``` r
devtools::session_info()
```

    ## - Session info ---------------------------------------------------------------
    ##  setting  value                       
    ##  version  R version 4.0.3 (2020-10-10)
    ##  os       Windows 10 x64              
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  ctype    English_United States.1252  
    ##  tz       Africa/Nairobi              
    ##  date     2021-04-12                  
    ## 
    ## - Packages -------------------------------------------------------------------
    ##  package     * version date       lib source        
    ##  assertthat    0.2.1   2019-03-21 [2] CRAN (R 4.0.3)
    ##  callr         3.5.1   2020-10-13 [2] CRAN (R 4.0.3)
    ##  cli           2.2.0   2020-11-20 [2] CRAN (R 4.0.3)
    ##  crayon        1.3.4   2017-09-16 [2] CRAN (R 4.0.3)
    ##  desc          1.2.0   2018-05-01 [2] CRAN (R 4.0.3)
    ##  devtools      2.3.2   2020-09-18 [2] CRAN (R 4.0.3)
    ##  digest        0.6.27  2020-10-24 [1] CRAN (R 4.0.3)
    ##  ellipsis      0.3.1   2020-05-15 [2] CRAN (R 4.0.3)
    ##  evaluate      0.14    2019-05-28 [2] CRAN (R 4.0.3)
    ##  fansi         0.4.2   2021-01-15 [2] CRAN (R 4.0.3)
    ##  fs            1.5.0   2020-07-31 [2] CRAN (R 4.0.3)
    ##  glue          1.4.2   2020-08-27 [2] CRAN (R 4.0.3)
    ##  htmltools     0.5.1   2021-01-12 [2] CRAN (R 4.0.3)
    ##  knitr         1.30    2020-09-22 [2] CRAN (R 4.0.3)
    ##  lifecycle     0.2.0   2020-03-06 [2] CRAN (R 4.0.3)
    ##  magrittr      2.0.1   2020-11-17 [2] CRAN (R 4.0.3)
    ##  memoise       1.1.0   2017-04-21 [2] CRAN (R 4.0.3)
    ##  pkgbuild      1.2.0   2020-12-15 [2] CRAN (R 4.0.3)
    ##  pkgload       1.1.0   2020-05-29 [2] CRAN (R 4.0.3)
    ##  prettyunits   1.1.1   2020-01-24 [2] CRAN (R 4.0.3)
    ##  processx      3.4.5   2020-11-30 [2] CRAN (R 4.0.3)
    ##  ps            1.5.0   2020-12-05 [2] CRAN (R 4.0.3)
    ##  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.0.3)
    ##  R6            2.5.0   2020-10-28 [2] CRAN (R 4.0.3)
    ##  remotes       2.2.0   2020-07-21 [2] CRAN (R 4.0.3)
    ##  rlang         0.4.10  2020-12-30 [2] CRAN (R 4.0.3)
    ##  rmarkdown     2.6     2020-12-14 [2] CRAN (R 4.0.3)
    ##  rprojroot     2.0.2   2020-11-15 [2] CRAN (R 4.0.3)
    ##  sessioninfo   1.1.1   2018-11-05 [2] CRAN (R 4.0.3)
    ##  stringi       1.5.3   2020-09-09 [2] CRAN (R 4.0.3)
    ##  stringr       1.4.0   2019-02-10 [2] CRAN (R 4.0.3)
    ##  testthat      3.0.1   2020-12-17 [2] CRAN (R 4.0.3)
    ##  usethis       2.0.0   2020-12-10 [2] CRAN (R 4.0.3)
    ##  withr         2.4.0   2021-01-16 [2] CRAN (R 4.0.3)
    ##  xfun          0.20    2021-01-06 [2] CRAN (R 4.0.3)
    ##  yaml          2.2.1   2020-02-01 [2] CRAN (R 4.0.3)
    ## 
    ## [1] C:/Users/BWaweru/OneDrive - CGIAR/Documents/R/win-library/4.0
    ## [2] C:/R/R-4.0.3/library
