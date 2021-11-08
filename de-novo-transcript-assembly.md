Assembly and quantification of *de novo* transcripts
================
Sam Buckberry
02/11/2021

For this section, you will need to have HISAT2 and Stringtie installed
and in your PATH to be accessed by the `system()` call from R.

You can build your own alignment index by following the documentation
for `hisat2-build` or download the pre-complied index to save a lot of
time from <https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz>

Load libraries

``` r
library(magrittr)
library(stringr)
library(Rsamtools)
```

Construct a manifest file for mapping all the samples

``` r
r1_list <- list.files("fastq/", pattern = "_filt_fastp_R1.fastq.gz",
                      full.names = TRUE)

r2_list <- list.files("fastq/", pattern = "_filt_fastp_R2.fastq.gz",
                      full.names = TRUE)

sam_list <- basename(r1_list) %>%
    str_replace("_filt_fastp_R1.fastq.gz", "_hisat2.sam")

dir.create("aligned_data")
sam_list <- str_c("aligned_data/", sam_list)

manifest <- data.frame(r1 = r1_list, r2 = r2_list, sam = sam_list)

write.table(x = manifest, file = "hisat_map_manifest.tsv",
            sep = "\t", quote = FALSE, col.names = FALSE,
            row.names = FALSE)
```

Align ***trimmed reads*** using HISAT2 with the option `--dta` that is
reccomended for the stringtie assembly. For *de novo* assemblies,
removing adapter sequences is much more important that if oyu are just
performing differential expression testing.

This code points to the pre-compliled indexes that will be accessible if
you are doing this workshop on one of the virtual machines with RStudio.

``` r
align_cmd <- str_c("while read i j k; do hisat2 --time --dta --threads 2 -x /home/data/hisat2/grch38_tran/genome_tran -1 $i -2 $j -S $k; done < hisat_map_manifest.tsv")

system(align_cmd)
```

Convert SAM to BAM, sort and index. Here we will use `mclapply` to speed
things up by processing in parallel. All of these samtools functions can
be performed in R using Rsamtools.

``` r
# List the hisat2 aligned SAM files
hisat_sam <- list.files(path = "aligned_data/",
                        pattern = "hisat2.sam$",
                        full.names = TRUE)

# Convert SAM to BAM
mclapply(hisat_sam, FUN = asBam, mc.cores = 4)

# List the hisat2 aligned BAM files
hisat_bam <- list.files(path = "aligned_data/",
                        pattern = "hisat2.bam$",
                        full.names = TRUE)

# Sort the hisat2 BAM files
sorted_bam <- str_replace_all(string = hisat_bam,
                                   pattern = ".bam",
                                   replacement = "_sorted")

mclapply(1:length(hisat_bam),
         FUN = function(x){sortBam(hisat_bam[x], sorted_bam[x])},
       mc.cores = 4)

# Index the sorted BAM files
mclapply(sorted_bam, indexBam, mc.cores = 4)
```

Make the reference-guided transcript assemblies for each library using
Stringtie

``` r
sorted_bam <- list.files("aligned_data/", "_hisat2_sorted.bam$",
                         full.names = TRUE)

system("pigz -d --keep reference/Homo_sapiens.GRCh38.104.chromosome.22.gtf")

run_stringtie <- function(bam){
    
    gtf <- str_replace(string = bam,
                           pattern = ".bam",
                           replacement = ".gtf")
    
    stringtie_cmd <- str_c("stringtie ", bam," -p 12 -G  reference/Homo_sapiens.GRCh38.104.chromosome.22.gtf -o ", gtf)
    
    system(stringtie_cmd)
}

mclapply(sorted_bam, run_stringtie, mc.cores = 4)
```

Merge sample assemblies into a meta-assembly of transcripts

``` r
stringtie_gtf <- list.files(path = "aligned_data/",
                            pattern = ".gtf$",
                            full.names = TRUE)

stringtie_cmd <- "stringtie aligned_data/*_hisat2_sorted.gtf --merge -G reference/Homo_sapiens.GRCh38.104.chromosome.22.gtf -p 4 -o aligned_data/merged_stringtie_transcripts.gtf"

system(stringtie_cmd)
```

Now you can use the reference guided *de novo* transcripts GTF with the
featureCounts function to get the transcripts counts as in the
`map-rna-subread` workflow in this workshop.
