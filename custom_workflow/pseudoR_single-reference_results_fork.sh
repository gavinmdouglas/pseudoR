#!/bin/bash

pseudoR_folder="/mfs/gdouglas/local/prg/pseudoR"

bash ${pseudoR_folder}/post-analysis_V3.sh -s ../${sraList} -d ${pseudoR_folder} 
Rscript --vanilla ${database}/Final_Output.singleReference.R ../${sraList} final_results/contig_analysis.step1.tsv ${database} contig
Rscript --vanilla ${database}/Final_Output.singleReference.R ../${sraList} final_results/orf_analysis.step1.tsv ${database} ORF
