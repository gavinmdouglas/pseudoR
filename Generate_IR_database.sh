while getopts "i:h" flag
do
    case "${flag}" in
        i) input_seq=${OPTARG};;
        h)
		echo "This program builds an IS termini database."
;;

    esac
done

#change nucl.fa to the fasta file with IS Sequences

mkdir temp
seqkit fx2tab ${input_seq}  | sed 's/_/:/g' > temp/IS_nucl.txt

cut -f1 temp/IS_nucl.txt | awk '{print $0 "_5prime"}' > temp/name_5.txt
cut -f1 temp/IS_nucl.txt | awk '{print $0 "_3prime"}' > temp/name_3.txt
cut -f1 temp/IS_nucl.txt > temp/name_full.txt
cut -f2 temp/IS_nucl.txt > temp/seqs.txt
cut -c -150 temp/seqs.txt > temp/seqs_5_150bp.txt
paste temp/name_5.txt temp/seqs_5_150bp.txt > temp/IRs.5.txt
rev temp/seqs.txt | cut -c -150 | rev > temp/seqs_3_150bp.txt
paste temp/name_3.txt temp/seqs_3_150bp.txt > temp/IRs.3.txt
paste temp/name_full.txt temp/seqs_5_150bp.txt temp/seqs_3_150bp.txt > temp/all_IRs.txt


paste temp/name_full.txt temp/seqs_5_150bp.txt temp/seqs_3_150bp.txt > temp/combined_seqs.part1.txt
awk 'BEGIN {FS = "\t";OFS = "\t"}  { print $1, $2 $3 } ' temp/combined_seqs.part1.txt > temp/combined_seqs.part2.txt
seqkit tab2fx temp/combined_seqs.part2.txt > temp/combined_seqs.part3.fna

cd-hit-est -i temp/combined_seqs.part3.fna -o temp/combined_seqs.part4.fna -c 0.9 -n 8

seqkit fx2tab temp/combined_seqs.part4.fna > temp/combined_seqs.part5.fna
cut -f1 temp/combined_seqs.part5.fna > temp/dedupe_names.txt

awk '{print $0 "_5prime"}' temp/dedupe_names.txt > temp/dedupe_names.5.txt
awk '{print $0 "_3prime"}' temp/dedupe_names.txt > temp/dedupe_names.3.txt

grep -w -f temp/dedupe_names.5.txt temp/IRs.5.txt > temp/IRs.5.final.txt
grep -w -f temp/dedupe_names.3.txt temp/IRs.3.txt > temp/IRs.3.final.txt


cat temp/IRs.5.final.txt temp/IRs.3.final.txt > temp/IRs.txt
seqkit tab2fx temp/IRs.txt > IRs.fa
makeblastdb -in IRs.fa -dbtype 'nucl'