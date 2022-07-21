#!/usr/bin/env Rscript


## SpliceAI output parser
## Developed by Daffodil Canson
## Implemented in R by Olga Kondrashova
## Date: 21/07/2022

###############################################################################

# Load argument parser
library(optparse)

# Arguments
option_list = list(
  make_option(c("-i", "--in_vcf"), type="character", default=NULL, 
              help="Input vcf file direct from SpliceAI", metavar="FILE"),
  make_option(c("-g", "--gene_list"), type="character", default=NULL, 
              help="Gene list, tab-separated with two columns - Gene and RefSeq_ID", 
              metavar="FILE"),
  make_option(c("-o", "--out_file"), type="character", default="out.tsv", 
              help="output file name [default= %default]", metavar="FILE"),
  make_option("--DS_AGDG_MIN_T", type="character", default=0.02, 
              help="Delta score (acceptor & donor gain) - minimum [default= %default]", 
              metavar="NUMERIC"),  
  make_option("--DS_AGDG_MAX_T", type="character", default=0.05, 
              help="Delta score (acceptor & donor gain) - maximum [default= %default]",
              metavar="NUMERIC"), 
  make_option("--GEX_size_MIN", type="character", default=25, 
              help="Gained exon size range - minimum [default= %default]",
              metavar="INT"), 
  make_option("--GEX_size_MAX", type="character", default=500, 
              help="Gained exon size range - maximum [default= %default]",
              metavar="INT"), 
  make_option("--DS_ALDL_MIN_T", type="character", default=0, 
              help="Delta score (acceptor & donor loss)  - minimum [default= %default]", 
              metavar="NUMERIC"), 
  make_option("--DS_ALDL_MAX_T", type="character", default=0.05, 
              help="Delta score (acceptor & donor loss) - maximum [default= %default]",
              metavar="NUMERIC"), 
  make_option("--AG_T", type="character", default=0.2, 
              help="Cryptic splice site - acceptor gain [default= %default]",
              metavar="NUMERIC"), 
  make_option("--DG_T", type="character", default=0.2, 
              help="Cryptic splice site - donor gain [default= %default]",
              metavar="NUMERIC"),
  make_option(c("-r", "--refseq_table"), type="character", default=NULL, 
              help="Optional: ncbiRefSeq table pre-downloaded from ucsc", 
              metavar="FILE"),
  make_option(c("--out_refseq"), type="character", default=NULL, 
              help="Optional: save downloaded ncbiRefSeq table", 
              metavar="FILE")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load required libraries
library(tidyverse, quietly=TRUE)

# Load all arguments as variables
DS_AGDG_MIN_T <- opt$DS_AGDG_MIN_T
DS_AGDG_MAX_T <- opt$DS_AGDG_MAX_T
GEX_size_MIN <- opt$GEX_size_MIN
GEX_size_MAX <- opt$GEX_size_MAX
DS_ALDL_MIN_T <-opt$DS_ALDL_MIN_T
DS_ALDL_MAX_T <- opt$DS_ALDL_MAX_T
AG_T <- opt$AG_T
DG_T <- opt$DG_T
input <- read_tsv(opt$in_vcf, comment="##", col_types = c(`#CHROM` = "c"))
genes <- read_tsv(opt$gene_list)
output_file <- opt$out_file
refseq_table_file <- opt$refseq_table
output_refseq <- opt$out_refseq

# ## Troubleshooting
# DS_AGDG_MIN_T <- 0.02
# DS_AGDG_MAX_T <- 0.05
# GEX_size_MIN <- 25
# GEX_size_MAX <- 500 
# DS_ALDL_MIN_T <- 0
# DS_ALDL_MAX_T <- 0.05
# AG_T <- 0.2
# DG_T <- 0.2
# input <- read_tsv("./example_variants.vcf",comment="##")
# genes <- read_tsv("./example_gene_list.txt")
# output_file <- "./example_variants_parsed.tsv"


###############################################################################

cat("\nPre-processing input\n")

# Splitting info field
input_splice_all <- input %>%
  separate(INFO, 
           into = c("ALLELE","SYMBOL","DS_AG","DS_AL","DS_DG",
                    "DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"), 
           sep="[|]",
           fill="right")

# Only including variants with SPLICE AI annotation & SNVs
input_splice_annot <- input_splice_all %>% 
  filter(DS_AG != "." & ALLELE != ".") %>% 
  filter(str_count(REF) == 1 & str_count(ALT) == 1)

# Keeping the rest of the variants for later to add to the final output
input_splice_other <- input_splice_all %>% 
  filter((DS_AG == "." | ALLELE == ".") | 
           (str_count(REF) != 1 | str_count(ALT) != 1))  %>%
  mutate_at(vars(c("DS_AG","DS_AL","DS_DG","DS_DL",
                   "DP_AG","DP_AL","DP_DG","DP_DL")), ~NA_real_)

###############################################################################

cat("\nPre-processing transcript list\n")


# Using the supplied gene list with RefSeq IDs to download the full UCSC 
# ncbiRefSeq table, and then subset by transcripts of interest.

genes <- genes %>%
  separate(RefSeq_ID, into = c("RefSeq_ID_no_v"), sep = "[.]", 
           extra = "drop", remove = FALSE)


if (is.null(refseq_table_file)) {
  library(rtracklayer, quietly=TRUE)
  query <- ucscTableQuery("hg19", table = "ncbiRefSeq")
  refseq_table <- getTable(query)

} else {
  refseq_table <- read_tsv(refseq_table_file, col_types = cols(.default = "c"))
  head(refseq_table)
}

if (!is.null(output_refseq)){
  write_tsv(refseq_table, output_refseq)
}


refseq_table <- refseq_table %>%
  separate(name, into = c("refseq_nov","refseq_version"), 
           sep = "[.]", extra = "drop", remove = FALSE)

refseq_table_filtered <- refseq_table %>%
  right_join(genes, by = c("refseq_nov"="RefSeq_ID_no_v")) %>%
  filter(!str_detect(chrom, "_")) %>% 
  mutate(refseq_match = if_else(RefSeq_ID == name, TRUE, FALSE))

cat("The following RefSeq ID versions have been updated")
refseq_table_filtered %>%
  filter(!refseq_match) %>%
  dplyr::rename(RefSeq_found = name, RefSeq_submitted = RefSeq_ID) %>% 
  dplyr::select(RefSeq_submitted,RefSeq_found)

refseq_table_filtered_expanded <- refseq_table_filtered %>%
  separate_rows(exonStarts, exonEnds, exonFrames, sep = ",", convert=TRUE) %>%
  filter(exonStarts != "") %>% 
  mutate(chrom_nochr = str_remove(chrom,"chr")) %>% 
  mutate(strand = if_else(strand == "-", -1, 1))

# Reformatting transcript table to have previous and next exons and 
# introns as columns (with exon and intron numbering).

refseq_table_filtered_expanded_bondaries <- refseq_table_filtered_expanded %>%
  group_by(name) %>%
  arrange(exonStarts) %>%
  mutate(exon_num_chrom = row_number(),
         eNum = if_else(strand == -1, 
                           as.integer(max(exon_num_chrom) - exon_num_chrom + 1), 
                           (exon_num_chrom))) %>%
  select(-exon_num_chrom) %>%
  dplyr::rename(eStart = exonStarts,
                eEnd = exonEnds,
                eFrame = exonFrames) %>%
  arrange(name,eNum) %>%
  mutate(prev_eStart = if_else(eNum > 1, lag(eStart), NA_integer_),
         prev_eEnd = if_else(eNum > 1, lag(eEnd), NA_integer_),
         prev_eFrame = if_else(eNum > 1, lag(eFrame), NA_integer_)) %>%
  mutate(next_eStart = if_else(eNum < exonCount, lead(eStart), NA_integer_),
         next_eEnd = if_else(eNum < exonCount, lead(eEnd), NA_integer_),
         next_eFrame = if_else(eNum < exonCount, lead(eFrame), NA_integer_)) %>%
  mutate(intronStart = case_when(eNum == 1 ~ NA_integer_, 
                                 strand == 1 ~ as.integer(prev_eEnd + 1),
                                 TRUE ~ as.integer(eEnd + 1)),
         intronEnd =  case_when(eNum == 1 ~ NA_integer_, 
                                strand == 1 ~ as.integer(eStart - 1),
                                TRUE ~ as.integer(prev_eStart - 1)),
         intron_num = if_else(eNum == 1, NA_integer_, as.integer(eNum - 1))) %>% 
  mutate(next_intronStart = case_when(eNum == exonCount ~ NA_integer_, 
                                      strand == 1 ~ as.integer(eEnd + 1),
                                      TRUE ~ as.integer(next_eEnd + 1)),
         next_intronEnd =  case_when(eNum == exonCount ~ NA_integer_, 
                                     strand == 1 ~ as.integer(next_eStart - 1),
                                     TRUE ~ as.integer(eStart - 1)),
         next_intron_num = if_else(eNum == exonCount, NA_integer_, as.integer(eNum)))

###############################################################################

cat("\nPerforming calculations\n")


# Combining variant table with transcripts table

input_splice_distance <- input_splice_annot %>%
  left_join(refseq_table_filtered_expanded, by = c("#CHROM" = "chrom_nochr",
                                                   "SYMBOL" = "Gene")) %>%
  group_by(ID) %>%
  # -1 is to account for 0-based position
  mutate(dist_exon_start = as.numeric(POS) - as.numeric(exonStarts) -1,
         dist_exon_end = as.numeric(POS) - as.numeric(exonEnds)) %>%
  mutate(exonic = if_else(dist_exon_start > 0 & dist_exon_end < 0,
                          "yes","no")) %>% 
  rowwise() %>% 
  mutate(dist_exon_closest_abs = min(abs(dist_exon_start), 
                                         abs(dist_exon_end))) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  # keeping only the exon that is the closest and not all other exons
  filter(dist_exon_closest_abs == min(dist_exon_closest_abs)) %>% 
  # adding directional distance from exon, not absolute distance
  mutate(d_from_exon = case_when(dist_exon_start <= 0 & dist_exon_end < 0 ~ abs(dist_exon_start),
                                 dist_exon_start > 0 & dist_exon_end >= 0 ~ dist_exon_end,
                                 TRUE ~ 0)) %>% 
  # add annotation for Daff's calculations
  mutate(d_250bp = if_else(d_from_exon > 250, "YES", "NO"),
         d_50bp = if_else(d_from_exon > 50, "YES", "NO")) %>%
  mutate_at(vars(c("DS_AG","DS_AL","DS_DG","DS_DL",
                   "DP_AG","DP_AL","DP_DG","DP_DL")), ~as.numeric(.)) 


# GEX prediction
output <- input_splice_distance %>%
  rowwise() %>%
  mutate(DS_AGDG_MIN = min(DS_AG, DS_DG),
         DS_AGDG_MAX = max(DS_AG, DS_DG),
         DS_AGDG_CUTOFF = if_else(DS_AGDG_MIN >= DS_AGDG_MIN_T & DS_AGDG_MAX >= DS_AGDG_MAX_T,
                                  "PASS","FAIL"),
         SS_AGDG_orientation = if_else((strand == 1 & DP_AG < DP_DG) | (strand == -1 & DP_AG > DP_DG), 
                                       TRUE, FALSE),
         GEX_size = if_else(DS_AGDG_CUTOFF == "PASS" & SS_AGDG_orientation,
                            abs(DP_AG - DP_DG) + 1, NA_real_),
         GEX_predict = if_else((GEX_size >= GEX_size_MIN  &  GEX_size <= GEX_size_MAX), 
                               "PASS", "FAIL"))  

# LEX, RET, Cryptic Acceptor and Donor prediction
output <- output %>% 
  mutate(DS_ALDL_MIN = min(DS_AL, DS_DL),
         DS_ALDL_MAX = max(DS_AL, DS_DL),
         DS_ALDL_CUTOFF = if_else(DS_ALDL_MIN >= DS_ALDL_MIN_T & DS_ALDL_MAX >= DS_ALDL_MAX_T,
                                  "PASS","FAIL"),
         SS_ALDL_orientation = if_else((strand == 1 & DP_AL < DP_DL) | (strand == -1 & DP_AL > DP_DL),
                                       TRUE, FALSE),
         LEX_predict = if_else(DS_ALDL_CUTOFF == "PASS" & SS_ALDL_orientation,
                               abs(DP_AL - DP_DL) + 1, NA_real_),
         RET_predict = if_else(DS_AL >=  DS_ALDL_MAX_T & DS_DL >= DS_ALDL_MAX_T & !SS_ALDL_orientation,
                               abs(DP_AL - DP_DL) - 1, NA_real_)) %>% 
  mutate(Cryptic_Acceptor_activ = if_else(DS_AG >= AG_T & DS_AG > DS_DG,
                                               "YES", "NO"),
         Cryptic_Donor_activ = if_else(DS_DG >= DG_T & DS_DG > DS_AG,
                                            "YES", "NO")) 

# Predicted exon sizes
output <- output %>% 
  mutate(Gained_exon_size = GEX_size,
         Lost_exon_size = LEX_predict,
         Retained_intron_size = RET_predict,
         bp_5prime = case_when(Cryptic_Acceptor_activ == "YES" & strand == 1 ~ (DP_AG - DP_AL),
                               Cryptic_Acceptor_activ == "YES" & strand == -1 ~ (DP_AL- DP_AG),
                               TRUE ~ 0),
         bp_3prime = case_when(Cryptic_Donor_activ == "YES" & strand == 1 ~ (DP_DG - DP_DL),
                               Cryptic_Donor_activ == "YES" & strand == -1 ~ (DP_DL- DP_DG),
                               TRUE ~ 0)) 

# Type of event
output <- output %>% 
  mutate(Partial_intron_retention = if_else(d_250bp == "NO" &  
                                              (Cryptic_Acceptor_activ == "YES" | 
                                                 Cryptic_Donor_activ == "YES") &
                                              ((bp_5prime < 0 & bp_5prime > -251) | 
                                                 (bp_3prime > 0 & bp_3prime < 251)), 
                                            "YES", "NO")) %>% 
  mutate(Pseudoexon_activation = if_else(d_50bp == "YES" & 
                                           GEX_predict == "PASS" & 
                                           !is.na(Gained_exon_size),
                                         "YES", "NO")) %>% 
  mutate(Partial_exon_deletion = if_else(d_50bp == "NO" & 
                                           (bp_5prime > 0 | bp_3prime < 0 ), 
                                         "YES", "NO")) %>% 
  mutate(Exon_skipping = if_else(d_50bp == "YES" | 
                                   is.na(Lost_exon_size),
                                 "NO", "YES")) %>% 
  mutate(Intron_retention = if_else(d_50bp == "YES" |
                                      is.na(Retained_intron_size),"NO","YES"))


# Summary
output <- output %>% 
  mutate(Any_splicing_abberation = if_else(Partial_intron_retention == "YES" |
                                             Pseudoexon_activation == "YES" | 
                                             Partial_exon_deletion == "YES" |
                                             Exon_skipping == "YES" | 
                                             Intron_retention == "YES", 
                                           "YES", "NO")) %>%
  mutate_at(vars(c(GEX_size, GEX_predict, LEX_predict, RET_predict)),
            ~ replace_na(as.character(.),"FAIL"))

###############################################################################


cat("\nWriting output\n")

## Save output
output_all <- output %>% 
  bind_rows(input_splice_other)# %>% 
#  arrange(`#CHROM`,POS)

output_all %>%
  write_tsv(output_file)

cat("\nSession Information:\n")

sessionInfo()
