#!/usr/bin/env Rscript


## Download transcript tables
## Written by Olga Kondrashova


################## Load libraries and arguments ################################

# Load options parser
library(optparse, quietly = TRUE)

option_list = list(
  make_option(c("-g", "--gene_list"), type="character", default=NULL, 
              help="Gene list, tab-separated with two columns - Gene and RefSeq_ID", 
              metavar="FILE"),
  make_option(c("--out_refseq"), type="character", default=NULL, 
              help="Save downloaded preprocessed transcript table", 
              metavar="FILE"),
  make_option(c("--out_tx_spliceai"), type="character", default=NULL, 
              help="Save preprocessed transcript table for SpliceAI", 
              metavar="FILE"),
  make_option(c("--out_refseq_full"), type="character", default=NULL, 
              help="Optional: save downloaded ncbiRefSeq table", 
              metavar="FILE"),
  make_option(c("--refseq_full"), type="character", default=NULL, 
              help="Optional: ncbiRefSeq table pre-downloaded from ucsc", 
              metavar="FILE"),
  make_option("--ref", type="character", default="hg19",
              help="Optional: reference genome (hg19 or hg38) [default= %default]",
              metavar="STRING"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Load other required libraries
library(tidyverse, quietly=TRUE)
library(rtracklayer, quietly=TRUE)


# Load all arguments as variables
genes <- read_tsv(opt$gene_list)
output_refseq <- opt$out_refseq
output_refseq_full <- opt$out_refseq_full
output_tx_spliceai <- opt$out_tx_spliceai
refseq_full_file <- opt$refseq_full
ref_genome <- opt$ref


## Troubleshooting
#genes <- read_tsv("./example_gene_list.txt")
#output_refseq <- "./example_refseq_tx.txt"
#output_tx_spliceai <- "./example_tx_spliceai.txt"
#refseq_full_file <- "./hg19_ucsc.txt"
#output_refseq_full <- "./hg19_ucsc.txt"
#ref_genome <- "hg19"
# 
# output_refseq <- "./example_refseq_tx_hg38.txt"
# output_tx_spliceai <- "./example_tx_spliceai_hg38.txt"
# ref_genome <- "hg38"
# refseq_full_file <- "./hg38_ucsc.txt"
# output_refseq_full <- "./hg38_ucsc.txt"

################## Pre-process transcripts #####################################

cat("\nPre-processing transcript list\n")

# Using the supplied gene list with RefSeq IDs to download the full UCSC 
# ncbiRefSeq table, and then subset by transcripts of interest.
genes <- genes %>%
  separate(RefSeq_ID, into = c("RefSeq_ID_no_v"), sep = "[.]", 
           extra = "drop", remove = FALSE)

# Read or download UCSC ncbiRefSeq table
if (is.null(refseq_full_file)) {
  query <- ucscTableQuery(ref_genome, table = "ncbiRefSeq")
  refseq_table <- getTable(query)
} else {
  refseq_table <- read_tsv(refseq_full_file, col_types = cols(.default = "c"))
}

if (!is.null(output_refseq_full)){
  write_tsv(refseq_table, output_refseq_full)
}


# only use refseq ID for matching not version
refseq_table <- refseq_table %>%
  separate(name, into = c("refseq_nov","refseq_version"), 
           sep = "[.]", extra = "drop", remove = FALSE)

# subsetting refseq list by gene list
refseq_table_filtered <- refseq_table %>%
  right_join(genes, by = c("refseq_nov"="RefSeq_ID_no_v")) %>%
  filter(!str_detect(chrom, "_")) %>% 
  mutate(refseq_match = if_else(RefSeq_ID == name, TRUE, FALSE))  %>%
  mutate(chrom = str_remove(chrom, "chr"))

refseq_table_filtered %>%
  write_tsv(output_refseq)

cat("The following RefSeq ID versions have been updated")
refseq_table_filtered %>%
  filter(!refseq_match) %>%
  dplyr::rename(RefSeq_found = name, RefSeq_submitted = RefSeq_ID) %>% 
  dplyr::select(RefSeq_submitted,RefSeq_found)

refseq_table_spliceAI <- refseq_table_filtered %>%
  mutate(`#NAME` = paste0("RefSeqTx-",name),      
         CHROM = str_remove(chrom,"chr"),
         STRAND = strand,
         TX_START = txStart,
         TX_END = txEnd,
         EXON_START = exonStarts,
         EXON_END = exonEnds) %>%
  select(`#NAME`,CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END) %>%
  distinct()

refseq_table_spliceAI %>%
  write_tsv(output_tx_spliceai)
  

