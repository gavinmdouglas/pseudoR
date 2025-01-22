if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("readr" %in% rownames(installed.packages()) == FALSE) {install.packages("readr")}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyr")}
library (readr, quietly = T)
library (tidyr, quietly = T)
library(dplyr, quietly = T)


print ("Begin Finding IS Circular Intermediates")

args = commandArgs(trailingOnly=TRUE)
options(dplyr.summarise.inform = FALSE)

args <- c("/mfs/gdouglas/projects/honey_bee/insertion_seq/test_SRRs.txt",
          "/mfs/gdouglas/projects/honey_bee/insertion_seq/test_output",
          "/mfs/gdouglas/projects/honey_bee/insertion_seq/test_output_combined/IS_circles.csv")

sampleList = read.table(args[1], stringsAsFactors = FALSE, header = FALSE, sep = "\t")
colnames(sampleList) = c("Sample")

intermediate_dir = args[2]

output_filepath = args[3]

combinedResults = data.frame()

for (i in seq(nrow(sampleList))){

  blast_results_filepath <- paste0(intermediate_dir, "/", sampleList$Sample[i], '/', sampleList$Sample[i], "_unmapped_blast_results.tsv")

  blast_results=readr::read_tsv(blast_results_filepath,
                               col_names = FALSE, show_col_types = FALSE)
  colnames (blast_results) = c("read_id", 'sseqid', 'pident','lngth', 'mismatch',
                               'gapopen','qstart', 'qend', 'sstart','send', 'evalue', 'bitscore')

  # this block sets up the blast_results dataframe
  blast_results = blast_results %>%
    rowwise() %>%
    mutate (IR_end = strsplit(sseqid, "_")[[1]][2])%>%
    mutate (subj_max = max(c(sstart,send)))%>%
    mutate (subj_min = min(c(sstart,send))) %>%
    mutate (IS = strsplit(sseqid, "_")[[1]][1])

  #this block finds reads that have multiple blast hits
  n_hits = blast_results %>% count (read_id, IS) %>%
    filter (n>1) %>%
    rowwise() %>%
    #the filter_key variable is a way to filter just the read_id+IS hits that had multiple blast hits, because some reads will also have hits to other
    # IS elements other than the IS that we care about.
    mutate (filter_key = paste (read_id, IS, sep='-'))
  if (nrow(n_hits) > 0) {
    blast_results = blast_results %>%
      mutate (filter_key = paste (read_id, IS, sep='-')) %>%
      filter (filter_key %in% n_hits$filter_key)
  }

  #this block finds reads that have blast hits from both the 3' and 5' termini of an IS element
  threePrime_reads = blast_results %>% filter (IR_end == "3prime") %>%
    select (read_id)
  fivePrime_reads = blast_results %>% filter (IR_end == "5prime") %>%
    select (read_id)
  blast_results = blast_results %>%
    #these first query position filters out multiple blast hits for each read_id+IS pair. This clears the way for the pivoting later and reduces redundancy.
    # It might also remove some real hits but I can't think of a way around this for the moment.
    group_by (read_id, IR_end) %>%
    filter (qstart== min (qstart)) %>%
    filter (qend == max(qend)) %>% ungroup() %>%
    # filter reads that map to both 3' and 5' termini
    filter (read_id %in% threePrime_reads$read_id) %>%
    filter (read_id %in% fivePrime_reads$read_id)

  #this block filters out blast hits that do not occur from the very termini of the IS element. I'm convinced that mapping to 3'
  # termini should map to the end (150bp) of that piece. Likewise, the 5' mapping should start at the very beginning (0bp) of the element.

  threePrime_reads = blast_results %>% filter (IR_end == "3prime") %>%
    filter (subj_max > 144) %>%
    select (read_id)
  fivePrime_reads = blast_results %>% filter (IR_end == "5prime") %>%
    filter (subj_min < 6) %>%
    select (read_id)
  #this block makes sure that the query starts at the beginning of the read, which would be weird if it didn't
  queryStart = blast_results %>%
    filter (qstart < 6) %>%
    select (read_id)
  blast_results = blast_results %>%
    filter (read_id %in% threePrime_reads$read_id) %>%
    filter (read_id %in% fivePrime_reads$read_id) %>%
    filter (read_id %in% queryStart$read_id) %>%
    #this final distinct function just clears out any duplicates or other redundancy
    distinct (read_id, sseqid, .keep_all = TRUE)

  if (nrow(blast_results) > 0) {
  blast_results = blast_results %>% select (read_id, qstart,qend, sstart,send, IR_end, subj_max, subj_min, IS)%>%
    pivot_wider(names_from = IR_end, values_from = c(qstart,qend, sstart,send, subj_max, subj_min)) %>%
    rowwise() %>%
    filter (qend_3prime != qend_5prime) %>%
    filter (qstart_3prime != qstart_5prime)

  col_order = c("read_id", "IS", "qstart_3prime","qend_3prime", "sstart_3prime", "send_3prime",
                "qstart_5prime","qend_5prime", "sstart_5prime", "send_5prime")

  blast_results <- blast_results[, col_order]

  blast_results = blast_results %>% mutate (sample = sampleList$Sample[i])
  pre_mapping_csv <- paste0(intermediate_dir, "/", sampleList$Sample[i], '/', sampleList$Sample[i], "_pre_mapping_output.csv")
  IS_annot = readr::read_csv(pre_mapping_csv ,show_col_types = FALSE) %>%
    select (read_id, read_seq, read_length)
  blast_results = blast_results %>% left_join (IS_annot) %>%
    mutate (map_length = ((qend_3prime-qstart_3prime)+(qend_5prime-qstart_5prime))) %>%
    filter (map_length >= 0.75 * read_length)

  combinedResults = rbind(combinedResults, blast_results)
  }
}

combinedResults = combinedResults %>% distinct (read_id, .keep_all = TRUE)

write_csv (combinedResults, file=output_filepath)
print ("Finished Finding IS Circular Intermediates")
