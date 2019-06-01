"# SNARE-seq-processing pipeline" 

## Dependency
[Drop-seq_tools-1.13/Picard](https://github.com/broadinstitute/Drop-seq/releases)

[Reference(mm10,hg38)](https://github.com/kundajelab/atac_dnase_pipelines#genome-data)

[.dict file](https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)

[fastx_trimmer](http://hannonlab.cshl.edu/fastx_toolkit/index.html)

## fastq file Structure 
Read1(30 bp): [12-bp barcode][9-bp UMI]TTTTTTTTT

Read2(75 bp): Nextera Read1

Read3(75 bp): Nextera Read2

## Usage
bash SNARE_processing.sh -n $name -s $species

