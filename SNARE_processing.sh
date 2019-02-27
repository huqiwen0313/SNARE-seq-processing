#!/bin/bash
# fastq file processing steps modified from SNARE_Prep.sh, https://github.com/chensong611/SNARE_prep 
# Usage: SNARE_processing.sh -n $name -s $species

dropseq_dir="/home/qiwenhu/software/Drop-seq_tools-1.13"
picard_dir="/home/qiwenhu/software/Drop-seq-1.13/lib"
atac_dir="/home/qiwenhu/software/atac_dnase_pipelines"
ref_dir="/home/qiwenhu/ref_genome/mm10_bd/mm10"

species='m'

# Barcode end base position
bEnd=12

# UMI start and end base positions
uStart=$(($bEnd+1))
uEnd=21

while getopts ":n:s:c:" options; do
    case $options in
    	n ) name=$OPTARG;;
        s ) species=$OPTARG;;
	    #c ) cells=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))

if [ $species = 'h' ]; then
    ref_fasta=hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    ref_dict=hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.dict
    species='hg38'
elif [ $species = 'm' ]; then
    ref_fasta=/mm10_no_alt_analysis_set_ENCODE.fasta
    ref_dict=/mm10_no_alt_analysis_set_ENCODE.dict
    species='mm10'
else
	echo "No reference found!"
fi

mkdir -p Reports

# reads prefiltering - filtering small cells (cutoff 3000) and output barcode summary
gunzip $name'_R1_001.fastq.gz'
gunzip $name'_R2_001.fastq.gz'
gunzip $name'_R3_001.fastq.gz'
./src/filter.barcode $name'_R1_001.fastq' $name'_R2_001.fastq' $name'_R3_001.fastq'
gzip  $name*'filter.fastq' 

#merge read1 and read2 together as a single read 
paste -d "" <(zcat $name'_R1_001.filter.fastq.gz') <(zcat $name'_R2_001.filter.fastq.gz' | sed '1~2s/.*//') > $name'_R12_001.filter.fastq'

# Make a pair-ended bam file, extract cell barcode and UMI information and store it in bam tags (XC/XM), and convert it back to pair-end fastq files
java -Xmx16g -jar $picard_dir/picard-2.18.5.jar FastqToSam FASTQ=$name'_R12_001.filter.fastq' FASTQ2=$name'_R3_001.filter.fastq.gz' SAMPLE_NAME=$name OUTPUT=/dev/stdout | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=Reports/$name'.cell_tag_report.txt' BASE_RANGE=1-$bEnd BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=Reports/$name'.molecule_tag_report.txt' BASE_RANGE=$uStart-$uEnd BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 | \
tee $name'.unaligned.tagged.bam' | \
java -Xmx16g -jar $picard_dir/picard-2.18.5.jar SamToFastq I=/dev/stdin FASTQ=$name'_R12_002.filter.fastq' SECOND_END_FASTQ=$name'_R3_001.filter.fastq'

# Remove read1 from merged single read
fastx_trimmer -f 31 -i $name'_R12_002.filter.fastq' -o $name'_R2_001.filter.fastq' -Q33

# Align reads and generate sorted bam files
bowtie2 -x $ref_dir'/bowtie2_index'$ref_fasta -1 $name'_R2_001.filter.fastq' -2 $name'_R3_001.filter.fastq' | \
samtools view -Su /dev/stdin | samtools sort - $name
samtools index $name'.bam' 

# Remove duplicated reads
samtools rmdup -s $name'.bam' $name'.unique.bam'

# Merge bam file and tagged barcode informaiton
java -Xmx16g -jar $picard_dir/picard-2.18.5.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir$ref_fasta UNMAPPED_BAM=$name'.unaligned.tagged.bam' ALIGNED_BAM=$name'.unique.bam' O=$name'.tagged.merged.bam' INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=true ATTRIBUTES_TO_RETAIN=ZI TMP_DIR=Tmp

# Call peaks and Generate count matrix
Rscript ./R/sc.count.R $name'.tagged.merged'
