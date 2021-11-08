# RNA-seq data analysis in R

With the emergence of high-throughput sequencing technologies, it has become trivial to profile the expression levels of tens of thousands of genes in a single experiment. However, the challenge in many of these experiments is often the analysis of the data rather than the generation of the data.

A whole suite of computational tools and methods are required for getting meaningful results from RNA-seq experiments. Although RNA-seq experiments can yield a lot of valuable information, challenges also exist in dealing with multiple sources of noise and bias in the data. However, most steps in RNA-seq data analysis are reasonably mature.

In this series of exercises, we will first describe how to perform quality contol on the sequence data obtained from an Illumina sequencing run, how to map these reads to a reference genome/transcriptome, derive a read count table for control-treatment differential expression analysis, how to perform functional enrichment tests, and how to visualise your the data throught your analyses.

In this workshop, example datasets have beed obtained from Calderon et al. (2019) **Landscape of stimulation-responsive chromatin across diverse human immune cells** _Nature Genetics_. [https://www.nature.com/articles/s41588-019-0505-9](https://www.nature.com/articles/s41588-019-0505-9). See below for more information.  

### Getting setup

**RStudio:** Although not essential, this workshop is designed to be run from RStudio. Therefore it is reccomended that you install RStudio to your local computer, or RStudio server on a remote system.  

Within RStudio, you can then clone this GitHub repository to your computer. 

Otherwise you can clone this repository in a bash terminal:  
`git clone https://github.com/SamBuckberry/RNAseq-workshop.git`

**Installing required packages:** This repository contains an R script that will install all the required R packages for the workshop. Install these packages by running the following command in the R console:  
`source ./install-r-packages.R`

---

### Experimental design considerations
- What RNA-seq library preperation method was used? Are the transcripts polyA selected? Were ribosomal RNA depleted? Was the library preperation protocol stranded? Was there a size selection for miRNA?
- Were spike-in controls used?
- What sequencing platform was used?
- What is the sequencing paired-end or single-end?

---

## Analysis workflow

### Pre-alignment
Follow the workflow in the `fastq-quality-control` Rmd or html file which covers:  

- FASTQ file quality assessment
- Read filtering and trimming for adapters and low-quality base calls

### Alignment 
Follow the workflow in the `map-rna-subread` Rmd or html file which covers:

- Build an alignment index Select reference genome (FASTA files) and gene models (GTF/GFF files)
- Align fastq files to index
- Inspect alignment statistics

### Post alignment
Follow the workflow in the `

- Sort and post-processing BAM/SAM files
- Perform gene and/or transcript quantificaion
- Inspect quantification metrics, and aggregate plots to indentify problematic samples or batch effects. 
- Test for differential expressed genes
- Inspect differential testing plots to identify potential issues with normalisation
- Plot differentially expressed genes 
- Perform Ontology testing to identify biological pathways and functions associated with experimental treatment(s). 

### Optional: _de novo_ Transcript asssembly
- HISAT2 alignment 
- Stringtie assembly of _de novo_ transcripts

If you do not have a reference genome or transcriptome, [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) is a highly regarded tool for de novo transcript assembly using short read illumina data. The computational requirements are extensive and its use is beyond the scope of this workshop. One decent guide (among many) to the considerations of de novo transcriptome assembly can be found [here](https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html) 

---

### Data sources
All raw fastq data analysed in this workshop is from Calderon et al. (2019) **Landscape of stimulation-responsive chromatin across diverse human immune cells** _Nature Genetics_. [https://www.nature.com/articles/s41588-019-0505-9](https://www.nature.com/articles/s41588-019-0505-9). 

Data are available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118165


**Reference genome and transcriptome (Ensembl GRCh38):**

Whole genome FASTA:  
http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Chromosome 22 FASTA:  
http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

Gene assembly:  
http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
