# HLA-typing using HLA-LA tool (Nextflow pipeline)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

The pipeline does Human Leukocyte Antigen (HLA) typing using HLA-LA.
See the [Paper](https://academic.oup.com/bioinformatics/article/35/21/4394/5426702) and the [GitHub](https://github.com/DiltheyLab/HLA-LA) repo.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker and singularity containers making installation trivial and results highly reproducible.

The pipeline takes in bam files as input. As a prerequisite, ensure your bam files are located within position 29mb -34mb on chromosome 6.
If not, you can extract the files using:
```
samtools index "path to sorted_bam"
samtools view -b "path to sorted_bam" "6:29000000-34000000" > "your output file"
```

## Installation 

1. Nextflow
```
wget -qO- https://get.nextflow.io | bash
```
2. HLA-LA
```
conda install -c bioconda hla-la
```
3. Reference graph
- Download the data package (2.3G)
```
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
```
- Index the graph, can take a few hours and might take up to 40G of memory.
```
HLA-LA.pl --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
```

## Running the pipeline
The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.

### Test data
To execute the pipeline on test dataset run:

 ```
 nextflow run nanjalaruth/HLA-typing -profile test -r main --reference_genome "path to the graph reference genome <hg19>" -resume
 ```
### Own data
Start running your own analysis either by using flags as shown below:

```
nextflow run nanjalaruth/hla_typing_using_HLA-LA -profile slurm -resume --input "*.bam" --reference_genome "path to the graph reference genome"  
```
 or run your own analysis by modifying the conf/test.config file to suit the path to your data location and then run the command as below:
 
 ```
 nextflow run nanjalaruth/hla_typing_using_HLA-LA -profile slurm -c <path to your edited config file> -resume
 ```
    
## To run the updated version of this pipeline, run:

 ```
 nextflow pull nanjalaruth/hla_typing_using_HLA-LA
 ```
    
## Arguments

### Required Arguments
| Argument  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,slurm\>                    | Configuration profile to use. Slurm is a job scheduler, you could otherwise use pbs                                       |
| --input  | \</project/\*.bam\> | Directory pattern for bam files.                                   |
| --reference_genome    | \<hg19\>              | Path to the reference genome to which the samples will be mapped |


## Support
I track open tasks using github's [issues](https://github.com/nanjalaruth/hla_typing_using_HLA-LA/issues)
