#!/bin/bash

SRR="test"

# Variables that user should customize:
bowtie2_db="/mfs/gdouglas/projects/honey_bee/helpful_bee/honeybee_microbiota_pangenome/bowtie2_indices/all_species_pangenome_reference"
pseudoR_folder="/mfs/gdouglas/local/prg/pseudoR"
input_fastq_dir="/mfs/gdouglas/projects/honey_bee/insertion_seq/test_fastqs"
output_dir="/mfs/gdouglas/projects/honey_bee/insertion_seq/test_output/${SRR}"
threads=5

mkdir -p $output_dir/bowtie2_output

# Map reads, and get mapped reads in BAM output.
bowtie2 \
	-x $bowtie2_db \
	-p ${threads} \
	-1 $input_fastq_dir/${SRR}_paired.1.fastq.gz \
	-2 $input_fastq_dir/${SRR}_paired.2.fastq.gz \
	-U $input_fastq_dir/${SRR}_unpaired.1.fastq.gz \
	-U $input_fastq_dir/${SRR}_unpaired.2.fastq.gz \
	-S $output_dir/bowtie2_output/${SRR}.sam \
	2> $output_dir/bowtie2_output/${SRR}.initial_mapping_log.txt

# Convert to BAM, sort, and index.
samtools view -Sb --threads ${threads} $output_dir/bowtie2_output/${SRR}.sam > $output_dir/bowtie2_output/${SRR}.bam
samtools sort --threads ${threads} -m 4G -o $output_dir/bowtie2_output/${SRR}.sort.bam $output_dir/bowtie2_output/${SRR}.bam
rm $output_dir/bowtie2_output/${SRR}.sam $output_dir/bowtie2_output/${SRR}.bam
samtools index $output_dir/bowtie2_output/${SRR}.sort.bam

# Get unmapped reads and save to FASTA (with numbered sequence IDs).
# Note that unmapped reads are sorted by sequence ID, to ensure the numbering is consistent if this workflow is re-run.
mkdir -p $output_dir/unmapped
samtools fastq -f 4 $output_dir/bowtie2_output/${SRR}.sort.bam | \
    paste - - - - | \
    sort -k1,1 | \
    tr '\t' '\n' > $output_dir/unmapped/${SRR}.unmapped.fastq

seqkit replace -p .+ -r "seq_{nr}" $output_dir/unmapped/${SRR}.unmapped.fastq > $output_dir/unmapped/${SRR}.unmapped.numbered.fastq
rm $output_dir/unmapped/${SRR}.unmapped.fastq

seqkit fq2fa $output_dir/unmapped/${SRR}.unmapped.numbered.fastq > $output_dir/unmapped/${SRR}.unmapped.fa

# BLAST unmapped reads against IS termini database.
# Do so in split files.
mkdir -p $output_dir/unmapped_blast_results
seqkit split2 -p ${threads} --force $output_dir/unmapped/${SRR}.unmapped.fa
ls $output_dir/unmapped/${SRR}.unmapped.fa.split > $output_dir/unmapped/${SRR}.unmapped_list.txt

cat $output_dir/unmapped/${SRR}.unmapped_list.txt | xargs -P${threads} -I% blastn -db ${pseudoR_folder}/IRs.fa -query $output_dir/unmapped/${SRR}.unmapped.fa.split/% -out $output_dir/unmapped_blast_results/%_blast_results.txt -num_threads 1 -outfmt 6

cat $output_dir/unmapped_blast_results/* > $output_dir/${SRR}_unmapped_blast_results.tsv

# Remove IS termini containing section from reads and remap shortened reads to contigs.
cut -f1 $output_dir/${SRR}_unmapped_blast_results.tsv | seqkit grep -f - $output_dir/unmapped/${SRR}.unmapped.numbered.fastq | seqkit fx2tab > $output_dir/${SRR}_IS_hits.tsv
timeout 20m Rscript --vanilla ${pseudoR_folder}/custom_workflow/V9_analysis_fork.R $output_dir $SRR
seqkit tab2fx $output_dir/${SRR}_trimmed_reads.tsv > $output_dir/${SRR}_trimmed_reads.fq

bowtie2 \
	-x $bowtie2_db \
	-p ${threads} \
    --very-sensitive \
	-U $output_dir/${SRR}_trimmed_reads.fq \
	2> $output_dir/bowtie2_output/${SRR}.IS_hits_remapped_log.txt | samtools view --threads ${threads} -Sb -o $output_dir/bowtie2_output/${SRR}.IS_hits_remapped.bam

# Get sequence IDs in bowtie2 index that have at least one IS mapped to them (i.e., unique the genes or contig IDs).
samtools view -F 4 $output_dir/bowtie2_output/${SRR}.IS_hits_remapped.bam | cut -f3 | sort | uniq > $output_dir/IS_hits_remapped_refids.txt
