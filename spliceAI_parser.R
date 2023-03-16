#!/usr/bin/env Rscript

## SpliceAI output parser
## Developed by Daffodil Canson
## Implemented in R by Olga Kondrashova & Aimee Davidson
## Functionality extended by Aimee Davidson
## Date: 14/12/2022

################## Load libraries and arguments ################################

# Load options parser
library(optparse, quietly = TRUE)

# Arguments
option_list = list(
  make_option(c("-i", "--in_vcf"), type="character", default=NULL, 
              help="Input vcf file direct from SpliceAI", metavar="FILE"),
  make_option(c("-r", "--refseq_table"), type="character", default=NULL, 
              help="Preprocessed UCSC ncbiRefSeq table (download_tx.R script)", 
              metavar="FILE"),
  make_option(c("-o", "--out_file"), type="character", default="out.tsv", 
              help="output file name [default= %default]", metavar="FILE"),
  make_option("--DS_AGDG_MIN_T", type="double", default=0.02, 
              help="Delta score (acceptor & donor gain) - minimum [default= %default]", 
              metavar="NUMERIC"),  
  make_option("--DS_AGDG_MAX_T", type="double", default=0.05, 
              help="Delta score (acceptor & donor gain) - maximum [default= %default]",
              metavar="NUMERIC"), 
  make_option("--GEX_size_MIN", type="double", default=25, 
              help="Gained exon size range - minimum [default= %default]",
              metavar="INT"), 
  make_option("--GEX_size_MAX", type="double", default=500, 
              help="Gained exon size range - maximum [default= %default]",
              metavar="INT"), 
  make_option("--DS_ALDL_MIN_T", type="double", default=0.02, 
              help="Delta score (acceptor & donor loss)  - minimum [default= %default]", 
              metavar="NUMERIC"), 
  make_option("--DS_ALDL_MAX_T", type="double", default=0.2, 
              help="Delta score (acceptor & donor loss) - maximum [default= %default]",
              metavar="NUMERIC"), 
  make_option("--AG_T", type="double", default=0.2, 
              help="Cryptic splice site - acceptor gain [default= %default]",
              metavar="NUMERIC"), 
  make_option("--DG_T", type="double", default=0.2, 
              help="Cryptic splice site - donor gain [default= %default]",
              metavar="NUMERIC"),
  make_option("--ref", type="character", default="hg19",
              help="Optional: reference genome (hg19 or hg38) [default= %default]",
              metavar="STRING")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load other required libraries
library(tidyverse, quietly=TRUE)
library(Biostrings, quietly = TRUE)

# Load all arguments as variables
DS_AGDG_MIN_T <- opt$DS_AGDG_MIN_T
DS_AGDG_MAX_T <- opt$DS_AGDG_MAX_T
GEX_size_MIN <- opt$GEX_size_MIN
GEX_size_MAX <- opt$GEX_size_MAX
DS_ALDL_MIN_T <-opt$DS_ALDL_MIN_T
DS_ALDL_MAX_T <- opt$DS_ALDL_MAX_T
AG_T <- opt$AG_T
DG_T <- opt$DG_T
input <- read_tsv(opt$in_vcf, comment="##",
                  col_types = c(`#CHROM` = "c"))
refseq_table <- read_tsv(opt$refseq_table,
                         col_types = cols(.default = "c"))
output_file <- opt$out_file
ref_genome <- opt$ref

# # Troubleshooting
# DS_AGDG_MIN_T <- 0.02
# DS_AGDG_MAX_T <- 0.05
# GEX_size_MIN <- 25
# GEX_size_MAX <- 500
# DS_ALDL_MIN_T <- 0.02
# DS_ALDL_MAX_T <- 0.2
# AG_T <- 0.2
# DG_T <- 0.2
# input <- read_tsv("./example_variants.vcf",comment="##",
#                   col_types = c(`#CHROM` = "c"))
# refseq_table <- read_tsv("./example_refseq_tx.txt",
#                          col_types = cols(.default = "c"))
# output_file <- "./example_variants_parsed.tsv"
# ref_genome <- "hg19"

################## Helper functions ############################################
# Helper functions amino acid predictions

# Mutate the DNA sequence to add the variant
add_variant <- function(ref,alt,varPos,extractRef,exonTable,exonSEQs) {
  # first check that the given ref and the genome reference match
  if (ref != extractRef) {
    return("reference mismatch")
  }
  # assuming that the reference is ok
  affectedExonLocale = which(exonTable$eStartAdj <= varPos & exonTable$eEnd >= varPos)
  affectedExon = exonSEQs[affectedExonLocale]
  # variant is not located within an affected sequence region
  # so don't need to make any adjustments
  if (is_empty(affectedExonLocale) == TRUE) {
    return("cannot_determine")
  } 
  # now update the exon with the variant
  adjustmentSize = (varPos - exonTable$eStartAdj[[affectedExonLocale]])+1
  adjustedExon = replaceAt(affectedExon,IRanges(adjustmentSize, adjustmentSize),alt)
  adjustedExon = as.character(adjustedExon)
  
  return(paste(affectedExonLocale,adjustedExon,sep="_"))
}


# Find the first non-identical character between two strings.
find_difference_point <- function(seqA,seqB,direction) {
  seqA = as.character(seqA)
  seqB = as.character(seqB)
  if(seqA == seqB) {
    return("identical")
  } else {
    seqAsplit = strsplit(seqA,"")[[1]]
    seqBsplit = strsplit(seqB,"")[[1]]
    if (direction == "reverse") {
      seqAsplit = rev(seqAsplit)
      seqBsplit = rev(seqBsplit)
    }
    point = suppressWarnings(which.min(seqAsplit == seqBsplit))
    if (direction == "reverse") {
      lengthAA = nchar(as.character(seqA))
      endPos = (lengthAA-point)+1
      point = endPos
    }
    return(point)
  }
}


# Format the refseq transcript, exon table
make_consensus_table <- function(transcript,refseqTab,filterNONCOD=TRUE) {
  # Extract the single transcript
  transcriptTab <- refseqTab %>%
    ungroup %>%
    filter(., `name` == transcript)
  # now correct for the zero UCSC positioning
  transcriptTab <- transcriptTab %>%
    mutate(eStartAdj = case_when(is.na(eStart) ~ NA_integer_,
                                 TRUE ~ as.integer(eStart+1)),
           intronEndAdj = case_when(is.na(intronEnd) ~ NA_integer_,
                                    TRUE ~ as.integer(intronEnd+1))) %>%
    rownames_to_column(., "rows") %>%
    mutate_at(vars(c("rows")), ~as.numeric(.))
  # fix the exon sizes
  transcriptTab <- transcriptTab %>%
    mutate(exon_size = as.integer(eEnd - eStart))
  # filter if don't want non-coding exons
  if (filterNONCOD == TRUE) {
    transcriptTab <- transcriptTab %>%
      filter(., eFrame != -1)
  }
  return(transcriptTab)
}


# Extract the intron that is retained with complete intron retention
find_retained_intron <- function(ALpos,DLpos,transcript,refseqTab) {
  intronTab <- make_consensus_table(transcript,refseqTab,filterNONCOD = FALSE)
  strand = intronTab$strand[[1]]
  if(strand == 1){
    intronnum = filter(intronTab, intronStart == DLpos+1 & intronEndAdj == ALpos-1)
  } else {
    intronnum = filter(intronTab, intronStart == ALpos+1 & intronEndAdj == DLpos-1)
  }
  # find intron and intron type
  intron = intronnum %>% pull(., intron_num)
  intron_type = intronnum %>% pull(., prev_eFrame)
  if(is_empty(intron) | is_empty(intron_type)) {
    return("NA|NA")
  } else {
    return(paste(intron,intron_type,sep="|"))
  }
}

# Extract the intron where the pseudoexon is placed, and the type of intron
find_pseudoexon_position <- function(DGpos,AGpos,strand,refseqTab,transcript) {
  # select exon info for transcript
  intronTab <- make_consensus_table(transcript,refseqTab,filterNONCOD = FALSE)
  if (strand == 1) {
    intronnum = filter(intronTab, intronStart+50 < AGpos & intronEndAdj-50 > DGpos)
  } else {
    intronnum = filter(intronTab, intronStart+50 < DGpos & intronEndAdj-50 > AGpos)
  }
  # find intron and intron type
  # exclude if branches multiple introns
  intron = intronnum %>% pull(., intron_num)
  intron_type = intronnum %>% pull(., prev_eFrame)
  if(is_empty(intron) | is_empty(intron_type) | length(intron) != 1) {
    return("NA|NA")
  } else {
    return(paste(intron,intron_type,sep="|"))
  }
}

# Extract which exon(s) have been skipped
find_exon_lost <- function(transcript,DLposition,ALposition,refseqTab) {
  exonTab = make_consensus_table(transcript,refseqTab,filterNONCOD = FALSE)
  lostExons = NA_character_
  size = NA_integer_
  # find the exons, for forward strand
  if (DLposition %in% exonTab$eStartAdj & ALposition %in% exonTab$eEnd) {
    exons = exonTab %>% filter(., DLposition == eStartAdj | ALposition == eEnd)
    exonNums = exons %>% pull(eNum)
    size = exons %>% pull(exon_size) %>% sum()
    lostExons = paste(unique(exonNums),collapse=",")
  }
  # find the exons, for reverse strand
  if (ALposition %in% exonTab$eStartAdj & DLposition %in% exonTab$eEnd) {
    exons = exonTab %>% filter(., ALposition == eStartAdj | DLposition == eEnd)
    exonNums = exons %>% pull(eNum)
    size = exons %>% pull(exon_size) %>% sum()
    lostExons = paste(unique(exonNums),collapse=",")
  }
  return(paste(lostExons,size,sep="|"))
}


# Extract the changed amino acid sequence for partial's
get_partial_SEQ <- function(transcript,consensusStart,consensusEnd,
                            partialStart,partialEnd,refseqTable,
                            frameshift,varPos,ref,alt) {
  # correct the start site, from ucsc table format
  consensusStart = consensusStart+1
  consensusTable = make_consensus_table(transcript,refseqTable,filterNONCOD = TRUE)
  # pull out the CDS start and stop
  cdsStartPos = as.integer(consensusTable$cdsStart[[1]])+1
  cdsEndPos = as.integer(consensusTable$cdsEnd[[1]])
  # check for partial end movement to before the start
  if((partialEnd-partialStart)+1 <= 0) {
    return("deletion greater than exon size")
  }
  # now manipulate the table to include the altered partial exon
  # under assumption partial has already been corrected for ucsc 1-base start
  partialTable = consensusTable
  strand = consensusTable$strand[[1]]
  if (strand == -1) {
    partialTable = partialTable %>% arrange(., desc(eStartAdj))
  }
  rowNumber = which(partialTable$eStartAdj==consensusStart | partialTable$eEnd == consensusEnd)
  # if can't find the target/affected exons then assume an entirely non-coding exon
  if (is_empty(rowNumber)) {
    return("does not affect coding region")
  }
  # now for variants altering the first or last coding exon
  # find first and last coding exon
  cdsStartRow = which(partialTable$eStart <= cdsStartPos & partialTable$eEnd >= cdsStartPos)
  cdsEndRow = which(partialTable$eStart <= cdsEndPos & partialTable$eEnd >= cdsEndPos)
  # first coding exon is altered
  if (rowNumber == cdsStartRow) {
    # for changes to start of exon (non-coding portion)
    if (partialStart != partialTable$eStartAdj[[rowNumber]] & partialEnd == partialTable$eEnd[[rowNumber]]) {
      if (partialStart > cdsStartPos) {
        # start site is altered
        return("impacts native start of stop site")
      } else if (partialStart <= cdsStartPos) {
        return ("does not affect coding region")
      } else {
        # shift exon start to coding start for aaseq prediction
        partialTable$eStartAdj[rowNumber] = cdsStartPos
      }
      # for changes to end of the exon (coding portion)
    } else if (partialStart == partialTable$eStartAdj[[rowNumber]] & partialEnd != partialTable$eEnd[[rowNumber]]) {
      # start site is deleted (any of 3 start site bases)
      if (partialEnd < cdsStartPos) {
        return("impacts native start or stop site")
      } else {
        partialTable$eStartAdj[rowNumber] = cdsStartPos
        partialTable$eEnd[rowNumber] = partialEnd
      }
    } else {
      return("cannot determine")
    }
  }
  # last coding exon is altered
  if (rowNumber == cdsEndRow) {
    # for changes to start of exon (coding portion)
    if (partialStart != partialTable$eStartAdj[[rowNumber]] & partialEnd == partialTable$eEnd[[rowNumber]]) {
      if (partialStart > cdsEndPos) {
        # stop site is altered
        return("impacts native start or stop site")
      } else {
        # shift exon end to coding end for aaseq prediction
        partialTable$eEnd[rowNumber] = cdsEndPos
        partialTable$eStartAdj[rowNumber] = partialStart
      }
      # for changes to end of the exon (non-coding portion)
    } else if (partialStart == partialTable$eStartAdj[[rowNumber]] & partialEnd != partialTable$eEnd[[rowNumber]]) {
      # stop site is deleted 
      if (partialEnd < cdsEndPos) {
        return("impacts native start of stop site")
      } else if (partialEnd >= cdsEndPos) {
        return ("does not affect coding region")
      } else {
        partialTable$eEnd[rowNumber] = cdsEndPos
      }
    } else {
      return("cannot determine")
    }
  }
  # for all other middle coding exons
  if (cdsStartRow != rowNumber & cdsEndRow != rowNumber) {
    # update the transcript table with the altered start and stop
    partialTable$eStartAdj[rowNumber] = partialStart
    partialTable$eEnd[rowNumber] = partialEnd
  }
  # correct for CDS positioning in both tables
  partialTable <- partialTable %>%
    mutate(eStartAdj = ifelse(eStartAdj < cdsStartPos, cdsStartPos, eStartAdj),
           eEnd = ifelse(eEnd > cdsEndPos, cdsEndPos, eEnd))
  consensusTable <- consensusTable %>%
    mutate(eStartAdj = ifelse(eStartAdj < cdsStartPos, cdsStartPos, eStartAdj),
           eEnd = ifelse(eEnd > cdsEndPos, cdsEndPos, eEnd))  
  # remove any but the immediate upstream exon
  partialTable <- partialTable %>% ungroup() %>%
    dplyr::slice(., ((rowNumber-1):n()))
  predictSEQ = determine_aaSEQ(partialTable,consensusTable,frameshift,varPos,ref,alt)
  return(predictSEQ)
}


# Extract the changed amino acid sequence for exon skipping
get_skip_SEQ <- function(exons,refseqTable,frameshift,transcript,varPos,ref,alt){
  if (is.na(exons)) {
    return("lost site/s do not match consensus")
  }
  consensusTable = make_consensus_table(transcript,refseqTable,filterNONCOD=FALSE)
  skipTable = consensusTable
  exonsList = purrr::flatten(str_split(exons, ","))
  exonsList = as.numeric(exonsList)
  # check for whether only non-coding exons are lost
  skippedExons = skipTable %>% filter(., eNum %in% exonsList)
  if (all(skippedExons$eFrame %in% c(-1))) {
    return("non-coding exons lost")
  }
  # check for whether start or stop is lost
  cdsstart = as.integer(skippedExons$cdsStart[[1]])+1
  cdsend = as.integer(skippedExons$cdsEnd[[1]])
  startPos = which(cdsstart >= skippedExons$eStart & cdsstart <= skippedExons$eEnd)
  endPos = which(cdsend >= skippedExons$eStart & cdsend <= skippedExons$eEnd)
  if (is_empty(startPos) == FALSE | is_empty(endPos) == FALSE) {
    return("native start or stop is lost")
  }
  # remove skipped exons
  # and also remove non-coding exons now that are accounted for
  skipTable = skipTable %>% filter(., !eNum %in% exonsList) %>%
    filter(., eFrame != -1)
  consensusTable = consensusTable %>% filter(., eFrame != -1)
  minExon = min(exonsList)-1
  strand = consensusTable$strand[[1]]
  if (strand == -1) {
    skipTable = skipTable %>% arrange(., desc(eStartAdj))
  }
  rowNumber = which(skipTable$eNum==minExon)
  # correct for CDS positioning in both tables
  skipTable <- skipTable %>%
    mutate(eStartAdj = ifelse(eStartAdj <= cdsstart, cdsstart, eStartAdj),
           eEnd = ifelse(eEnd >= cdsend, cdsend, eEnd))
  consensusTable <- consensusTable %>%
    mutate(eStartAdj = ifelse(eStartAdj <= cdsstart, cdsstart, eStartAdj),
           eEnd = ifelse(eEnd >= cdsend, cdsend, eEnd)) 
  # remove any but the immediate upstream exon
  skipTable <- skipTable %>% ungroup() %>%
    dplyr::slice(., (rowNumber:n()))
  predictSEQ = determine_aaSEQ(skipTable,consensusTable,frameshift,varPos,ref,alt)
  return(predictSEQ)
}



# Format the refseq table to reflect the pseudoexon activation impact
get_pseudo_SEQ <- function(pseudoStart,pseudoEnd,refseqTable,frameshift,
                           transcript,varPos,ref,alt) {
  if (is.na(pseudoStart) | is.na(pseudoEnd)) {
    return("gain site/s not intronic")
  }
  consensusTable = make_consensus_table(transcript,refseqTable,filterNONCOD = TRUE)
  pseudoTable = consensusTable
  strand = consensusTable$strand[[1]]
  currentChr = pseudoTable$chrom[[1]]
  # check whether pseudoexon is within coding region
  cdsstart = as.integer(min(consensusTable$cdsStart, na.rm = TRUE))+1
  cdsend = as.integer(min(consensusTable$cdsEnd, na.rm = TRUE))
  if (!(pseudoStart >= cdsstart & pseudoEnd <= cdsend)) {
    return("non-coding pseudoexon")
  }
  # add pseudoexon to the transcript table
  pseudoTable = pseudoTable %>% 
    ungroup() %>% 
    add_row(., chrom=currentChr,eStartAdj=pseudoStart,eEnd=pseudoEnd)
  # for reverse strand order the table
  if (strand == -1) {
    pseudoTable = pseudoTable %>%
      arrange(., desc(eStartAdj))
  } else {
    pseudoTable = pseudoTable %>%
      arrange(., eStartAdj)
  }
  pseudoTable = pseudoTable %>% rownames_to_column(., "new_rows")
  rowNumber = pseudoTable %>% filter(is.na(rows)) %>% 
    mutate_at(vars(c("new_rows")), ~as.numeric(.)) %>% 
    pull(., new_rows)
  # correct for CDS positioning in both tables
  pseudoTable <- pseudoTable %>%
    mutate(eStartAdj = ifelse(eStartAdj < cdsstart, cdsstart, eStartAdj),
           eEnd = ifelse(eEnd > cdsend, cdsend, eEnd))
  consensusTable <- consensusTable %>%
    mutate(eStartAdj = ifelse(eStartAdj < cdsstart, cdsstart, eStartAdj),
           eEnd = ifelse(eEnd > cdsend, cdsend, eEnd)) 
  # remove any but the immediate upstream exon
  pseudoTable <- pseudoTable %>% ungroup() %>%
    dplyr::slice(., ((rowNumber-1):n()))
  predictSEQ = determine_aaSEQ(pseudoTable,consensusTable,frameshift,varPos,ref,alt)
  return(predictSEQ)
}


# Format the refseq table to reflect the intron retention impact
get_retention_SEQ <- function(refseqTable,intron,frameshift,transcript,varPos,ref,alt) {
  if (is.na(intron)) {
    return("loss site/s do not match consensus")
  }
  consensusTable = make_consensus_table(transcript,refseqTable,filterNONCOD = TRUE)
  retentionTable = consensusTable
  strand = consensusTable$strand[[1]]
  # pull out the CDS start and stop
  cdsStartPos = as.integer(consensusTable$cdsStart[[1]])+1
  cdsEndPos = as.integer(consensusTable$cdsEnd[[1]])
  # extract appropriate positions to add the intron to transcript table
  currentChr = retentionTable$chrom[[1]]
  if (strand == 1) {
    neededStart = retentionTable %>% filter(., eNum==as.integer(intron)) %>% pull(., eEnd)
    neededEnd = retentionTable %>% filter(., eNum==as.integer(intron)+1) %>% pull(., eStartAdj)
  } else {
    neededStart = retentionTable %>% filter(., eNum==as.integer(intron)+1) %>% pull(., eEnd)
    neededEnd = retentionTable %>% filter(., eNum==as.integer(intron)) %>% pull(., eStartAdj)     
  }
  # assumes that for non-coding introns
  if (is_empty(neededStart) == TRUE | is_empty(neededEnd) == TRUE) {
    return("non-coding intron retained")
  }
  # add the intron to the transcript table
  retentionTable = retentionTable %>% 
    ungroup() %>% 
    add_row(., chrom=currentChr,eStartAdj=neededStart+1, eEnd=neededEnd-1)
  if (strand == 1) {
    retentionTable = retentionTable %>% 
      arrange(., eStartAdj) %>% 
      rownames_to_column(., "new_rows")
  } else {
    retentionTable = retentionTable %>% 
      arrange(., desc(eStartAdj)) %>% 
      rownames_to_column(., "new_rows")
  }
  rowNumber = retentionTable %>% filter(is.na(rows)) %>% 
    mutate_at(vars(c("new_rows")), ~as.numeric(.)) %>% 
    pull(., new_rows)
  retentionTable = retentionTable %>% dplyr::select(., -c("rows")) %>% 
    dplyr::rename("rows" = "new_rows") %>% 
    mutate_at(vars(c("rows")), ~as.numeric(.))
  # correct for CDS positioning in both tables
  retentionTable <- retentionTable %>%
    mutate(eStartAdj = ifelse(eStartAdj < cdsStartPos, cdsStartPos, eStartAdj),
           eEnd = ifelse(eEnd > cdsEndPos, cdsEndPos, eEnd))
  consensusTable <- consensusTable %>%
    mutate(eStartAdj = ifelse(eStartAdj < cdsStartPos, cdsStartPos, eStartAdj),
           eEnd = ifelse(eEnd > cdsEndPos, cdsEndPos, eEnd))
  # remove any but the immediate upstream exon
  retentionTable <- retentionTable %>% ungroup() %>%
    dplyr::slice(., ((rowNumber-1):n()))
  predictSEQ = determine_aaSEQ(retentionTable,consensusTable,frameshift,varPos,ref,alt)
  return(predictSEQ)
}


# Overview function to predict the amino acid sequence
determine_aaSEQ <- function(altTable,consensusTab,frameshift,varPos,ref,alt) {
  strand = consensusTab$strand[[1]]
  # remove the equivalent rows from the consensus table to match the altered table
  if (strand == -1) {
    consensusTab <- consensusTab %>% arrange(., desc(eStartAdj))
  }
  minKeepExon = min(altTable$eNum, na.rm = TRUE)
  minKeepExonRow = which(consensusTab$eNum == minKeepExon)
  currentChr = consensusTab$chrom[[1]]
  consensusTab <- consensusTab %>% ungroup() %>%
    dplyr::slice(., ((minKeepExonRow):n()))  
  # get the DNA sequences
  alteredExonDNAseqs <- getSeq(Hsapiens,paste0("chr",altTable$chrom),start=altTable$eStartAdj,end=altTable$eEnd)
  consensusExonDNAseqs <- getSeq(Hsapiens,paste0("chr",consensusTab$chrom),start=consensusTab$eStartAdj,end=consensusTab$eEnd)
  # make necessary adjustments for the variant itself
  # also check whether the reference is correct
  genomeRef = as.character(getSeq(Hsapiens,paste0("chr",currentChr),varPos,varPos))
  adjustedExonDNAseq = add_variant(ref = ref,
                                   alt = alt,
                                   varPos = varPos,
                                   extractRef = genomeRef,
                                   exonTable = altTable,
                                   exonSEQs = alteredExonDNAseqs)
  # exit as supplied reference is wrong
  if (adjustedExonDNAseq == "reference mismatch") {
    return(adjustedExonDNAseq)
  }
  # no changes needed as variant outside of altered sequence
  if (adjustedExonDNAseq != "cannot_determine") {
    tempNewExon = unlist(stringr::str_split(adjustedExonDNAseq,"_"))
    exonPos = as.numeric(tempNewExon[1])
    exonSEQ = DNAStringSet(tempNewExon[2])
    alteredExonDNAseqs[exonPos] = exonSEQ
  } 
  # reverse complement DNA sequences if on reverse strand 
  if (strand == -1) {
    alteredExonDNAseqs = reverseComplement(alteredExonDNAseqs)
    consensusExonDNAseqs = reverseComplement(consensusExonDNAseqs)
  }
  # correct the framing if the exon frame is not "0"
  if (altTable$eFrame[1] == 1) {
    alteredExonDNAseqs[1] = subseq(alteredExonDNAseqs[1],3)
    consensusExonDNAseqs[1] = subseq(consensusExonDNAseqs[1],3)
  } 
  # for frame of 2, remove leading one base
  if (altTable$eFrame[1] == 2) {
    alteredExonDNAseqs[1] = subseq(alteredExonDNAseqs[1],2)
    consensusExonDNAseqs[1] = subseq(consensusExonDNAseqs[1],2)
  }
  # now paste all the DNA sequences together
  alteredJointDNAseq <- paste(alteredExonDNAseqs,collapse="")
  consensusJointDNAseq <- paste(consensusExonDNAseqs,collapse="")
  # convert to amino acid sequence
  alteredAAseq <- suppressWarnings(translate(DNAString(alteredJointDNAseq)))
  consensusAAseq <- suppressWarnings(translate(DNAString(consensusJointDNAseq)))
  # check for whether there is a difference in protein sequence
  forwardDiff <- find_difference_point(alteredAAseq,consensusAAseq,direction="forward")
  # if no difference stop trying to predict
  if (forwardDiff == "identical") {
    return("no difference in protein sequence")
  }
  # now find the difference ins sequence
  if (forwardDiff >= 4) {
    alteredAAseq <- subseq(alteredAAseq,forwardDiff-3)
  }
  # for inframe sequences, strip off any consensus sequence from the end
  if (frameshift == "NO") {
    reverseDiff = find_difference_point(alteredAAseq,consensusAAseq,direction="reverse")
    # stop looking if identical sequence
    if (reverseDiff == "identical") {
      return("no difference in protein sequence")
    }
    if (reverseDiff < 3) {
      alteredAAseq <- subseq(alteredAAseq,1,6)
    } else {
      alteredAAseq <- subseq(alteredAAseq,1,reverseDiff+3)
    }
  }
  # now check for whether a stop has been introduced
  alteredStopPos <- regexpr(pattern="\\*",as.character(alteredAAseq))[1]
  # if frameshift and no stop, might continue reading past the consensus stop site
  if (frameshift == "YES" & alteredStopPos == -1) {
    return("protein sequence extends beyond native stop site")
  }
  # strip off all extra at end if a stop is introduced
  if (alteredStopPos != -1) {
    alteredAAseq <- subseq(alteredAAseq,1,alteredStopPos)
  }
  # now do some extra formatting for the end aa sequence
  preOutInfo = paste(as.character(alteredAAseq))
  lengthOutInfo = nchar(preOutInfo)
  if (alteredStopPos == -1) {
    # put branching square brackets around the altered sequence
    prefix = substr(x = preOutInfo, start = 1, stop = 3)
    if (lengthOutInfo >= 6) {
      suffix = substr(x = preOutInfo, start = lengthOutInfo-2, stop = lengthOutInfo)
      middle = substr(preOutInfo,4,(lengthOutInfo-3))
    } else {
      suffix = substr(x = preOutInfo, start = 4, stop = lengthOutInfo)
      middle = ""
    }
    outInfo = paste0(prefix,"[",middle,"]",suffix)
  } else {
    # put branching square brackets around the altered sequence until the stop
    outInfo = paste0(substr(preOutInfo,1,3),"[",substr(preOutInfo,4,lengthOutInfo),"]")
  }
  return(outInfo)
}

################## Genome version ##############################################

if (ref_genome == "hg19") {
  library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
} else if (ref_genome == "hg38") {
  library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
} else {
  stop("Valid genome reference was not provided (specify either hg19 or hg38)")
}

cat("\nReference genome used: ",ref_genome,"\n")

################## Load input ##################################################

cat("\nPre-processing input\n")


input_splice_all <- input %>%
  # Splitting SpliceAI predictions if multiple transcripts included
  separate_rows(INFO, sep = ",") %>%
  # Splitting info field
  separate(INFO, 
           into = c("ALLELE","SYMBOL","DS_AG","DS_AL","DS_DG",
                    "DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"), 
           sep="[|]",
           fill="right") %>%
  mutate(ALLELE = if_else((str_detect(ALLELE,"SpliceAI=") | ALLELE == "."),
                           ALLELE,
                           paste0("SpliceAI=",ALLELE))) %>%
  # Stripping "chr" in case vcf contains "chr"
  mutate(`#CHROM` = as.character(str_remove(`#CHROM`,"chr")))

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

################## Prepare transcripts #########################################

cat("\nPrepare transcripts\n")

refseq_table_expanded <- refseq_table %>%
  separate_rows(exonStarts, exonEnds, exonFrames, sep = ",", convert=TRUE) %>%
  filter(exonStarts != "") %>% 
 # mutate(chrom_nochr = str_remove(chrom,"chr")) %>% 
  mutate(strand = if_else(strand == "-", -1, 1))


# Reformatting transcript table to have previous and next exons and 
# introns as columns (with exon and intron numbering).
refseq_boundaries <- refseq_table_expanded %>%
  group_by(name) %>%
  arrange(exonStarts) %>%
  mutate(exon_num_chrom = row_number(),
         eNum = if_else(strand == -1, 
                           as.integer(max(exon_num_chrom) - exon_num_chrom + 1), 
                           (exon_num_chrom))) %>%
  dplyr::select(-exon_num_chrom) %>%
  dplyr::rename(eStart = exonStarts,
                eEnd = exonEnds,
                eFrame = exonFrames) %>%
  arrange(name,eNum) %>%
  mutate_at(vars(c("exonCount")), ~as.integer(.)) %>%
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
         next_intron_num = if_else(eNum == exonCount, NA_integer_, as.integer(eNum))) %>%
  rowwise() %>%
  mutate(exon_size = as.integer(eEnd - eStart),
         intron_size = as.integer(intronEnd - intronStart))

################## Perform calculations ########################################

cat("\nPerforming calculations\n")

# Depending on how SpliceAI was run (default or with supplied transcripts)
# the transcript table will be matched accordingly
if(str_detect(input_splice_annot$SYMBOL[1],"^RefSeqTx-")){
  refseq_boundaries <- refseq_boundaries %>%
    mutate(chrom = str_remove(chrom,"chr")) %>% # REMOVE
    mutate(SYMBOL = paste0("RefSeqTx-", name))
  
}else if(!str_detect(input_splice_annot$SYMBOL[1],"^RefSeqTx-")){
  refseq_boundaries <- refseq_boundaries %>%
    mutate(chrom = str_remove(chrom,"chr")) %>% # REMOVE
    mutate(SYMBOL = Gene)
}

# Combining variant table with transcripts table
input_splice_distance <- input_splice_annot %>%
  left_join(refseq_boundaries, by = c("#CHROM" = "chrom",
                                      "SYMBOL" = "SYMBOL")) %>%
  group_by(ID,SYMBOL) %>%
  # -1 is to account for 0-based position
  mutate(dist_exon_start = as.numeric(POS) - as.numeric(eStart) -1,
         dist_exon_end = as.numeric(POS) - as.numeric(eEnd)) %>%
  mutate(exonic = if_else(dist_exon_start > 0 & dist_exon_end < 0,
                          "yes","no")) %>% 
  rowwise() %>% 
  mutate(dist_exon_closest_abs = min(abs(dist_exon_start), 
                                         abs(dist_exon_end))) %>% 
  ungroup() %>% 
  group_by(ID,SYMBOL) %>% 
  # keeping only the exon that is the closest and not all other exons
  filter(dist_exon_closest_abs == min(dist_exon_closest_abs)) %>% 
  # adding directional distance from exon, not absolute distance
  mutate(d_from_exon = case_when(dist_exon_start <= 0 & dist_exon_end < 0 ~ abs(dist_exon_start),
                                 dist_exon_start > 0 & dist_exon_end >= 0 ~ dist_exon_end,
                                 TRUE ~ 0)) %>% 
  # add annotation for distance to closest exon
  mutate(d_250bp = if_else(d_from_exon > 250, "YES", "NO"),
         d_50bp = if_else(d_from_exon > 50, "YES", "NO")) %>%
  mutate_at(vars(c("DS_AG","DS_AL","DS_DG","DS_DL",
                   "DP_AG","DP_AL","DP_DG","DP_DL")), ~as.numeric(.)) 

# Add genomic positioning of the spliceAI prediction sites
output <- input_splice_distance %>%
  ungroup() %>%
  rowwise() %>%
  mutate(., GEO_AG = POS+DP_AG,
         GEO_AL = POS+DP_AL,
         GEO_DG = POS+DP_DG,
         GEO_DL = POS+DP_DL)

# GEX prediction
output <- output %>%
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
  rowwise() %>%
  mutate(DS_ALDL_MIN = min(DS_AL, DS_DL),
         DS_ALDL_MAX = max(DS_AL, DS_DL),
         DS_ALDL_CUTOFF = if_else(DS_ALDL_MIN >= DS_ALDL_MIN_T & DS_ALDL_MAX >= DS_ALDL_MAX_T,
                                  "PASS","FAIL"),
         SS_ALDL_orientation = if_else((strand == 1 & DP_AL < DP_DL) | (strand == -1 & DP_AL > DP_DL),
                                       TRUE, FALSE),
         LEX_predict = if_else(DS_ALDL_CUTOFF == "PASS" & SS_ALDL_orientation,
                               abs(DP_AL - DP_DL) + 1, NA_real_),
         RET_predict = if_else(DS_ALDL_CUTOFF == "PASS" & !SS_ALDL_orientation,
                               abs(DP_AL - DP_DL) - 1, NA_real_)) %>% 
  
  mutate(Cryptic_Acceptor_activation = if_else(DS_AG >= AG_T & DS_AG > DS_DG,
                                               "YES", "NO"),
         Cryptic_Donor_activation = if_else(DS_DG >= DG_T & DS_DG > DS_AG,
                                            "YES", "NO")) 

# Orientation check for partial retention and deletions
output <- output %>%
  rowwise() %>%
  mutate(Cryptic_Acceptor_orientation = if_else((strand == 1 & prev_eEnd - GEO_AG < 0) | (strand == -1 & GEO_AG - (prev_eStart+1) < 0),"PASS", "FAIL"),
         Cryptic_Donor_orientation = if_else((strand == 1 & GEO_DG - (next_eStart+1) < 0) | (strand == -1 & next_eEnd - GEO_DG < 0), "PASS", "FAIL"))


# Predicted exon sizes
output <- output %>% 
  rowwise() %>%
  mutate(Gained_exon_size = GEX_size,
         Lost_exon_size = LEX_predict,
         Retained_intron_size = RET_predict,
         bp_5prime = case_when(Cryptic_Acceptor_activation == "YES" & strand == 1 & Cryptic_Acceptor_orientation == "PASS" ~ GEO_AG-(eStart+1),
                               Cryptic_Acceptor_activation == "YES" & strand == -1 & Cryptic_Acceptor_orientation == "PASS" ~ eEnd-GEO_AG,
                               TRUE ~ 0),
         bp_3prime = case_when(Cryptic_Donor_activation == "YES" & strand == 1 & Cryptic_Donor_orientation == "PASS" ~ GEO_DG-eEnd,
                               Cryptic_Donor_activation == "YES" & strand == -1 & Cryptic_Donor_orientation == "PASS" ~ (eStart+1)-GEO_DG,
                               TRUE ~ 0))

# Placement of GEX
output <- output %>%
  rowwise() %>%
  mutate(Pseudoexon_intron_info = if_else(d_50bp == "YES" & GEX_predict == "PASS" & !is.na(Gained_exon_size),
                                          find_pseudoexon_position(DGpos = GEO_DG,
                                                                   AGpos = GEO_AG,
                                                                   strand = strand,
                                                                   refseqTab = refseq_boundaries,
                                                                   transcript = name),
                                          "NA|NA")) %>%
  separate(., col = Pseudoexon_intron_info, into = c("Pseudoexon_intron","Pseudoexon_intron_type"), sep = "[|]", remove = TRUE) %>%
  mutate(Pseudoexon_intron = na_if(Pseudoexon_intron, "NA"),
         Pseudoexon_intron_type = na_if(Pseudoexon_intron_type, "NA"))

# Type of event
output <- output %>% 
  rowwise() %>%
  mutate(Partial_intron_retention = if_else(d_250bp == "NO" &  
                                              (Cryptic_Acceptor_activation == "YES" | 
                                                 Cryptic_Donor_activation == "YES") &
                                              ((bp_5prime < 0 & bp_5prime > -251) | 
                                                 (bp_3prime > 0 & bp_3prime < 251)), 
                                            "YES", "NO")) %>% 
  mutate(Pseudoexon_activation = if_else(d_50bp == "YES" & 
                                           GEX_predict == "PASS" & 
                                           !is.na(Gained_exon_size) &
                                           !is.na(Pseudoexon_intron),
                                         "YES", "NO")) %>% 
  mutate(Partial_exon_deletion = if_else(d_50bp == "NO" & 
                                           (bp_5prime > 0 | bp_3prime < 0 ), 
                                         "YES", "NO")) %>% 
  mutate(Exon_skipping = if_else(d_50bp == "YES" | 
                                   is.na(Lost_exon_size),
                                 "NO", "YES")) %>% 
  mutate(Intron_retention = if_else(d_50bp == "YES" |
                                      is.na(Retained_intron_size),"NO","YES"))


# Summary of predictions
output <- output %>% 
  rowwise() %>%
  mutate(Any_splicing_aberration = if_else(Partial_intron_retention == "YES" |
                                             Pseudoexon_activation == "YES" | 
                                             Partial_exon_deletion == "YES" |
                                             Exon_skipping == "YES" | 
                                             Intron_retention == "YES", 
                                           "YES", "NO")) %>%
  mutate_at(vars(c(GEX_size, GEX_predict, LEX_predict, RET_predict)),
            ~ replace_na(as.character(.),"FAIL"))

################## Add additional information ##################################

# Add additional psuedoexon activation information
output <- output %>%
  rowwise() %>%
  mutate(Pseudoexon_start = case_when(Pseudoexon_activation == "YES" & 
                                        is.na(Pseudoexon_intron)==FALSE ~ as.integer(min(GEO_AG,GEO_DG)),
                                      TRUE ~ NA_integer_),
         Pseudoexon_end = case_when(Pseudoexon_activation == "YES" & 
                                      is.na(Pseudoexon_intron)==FALSE~ as.integer(max(GEO_AG,GEO_DG)),
                                    TRUE ~ NA_integer_),
         Pseudoexon_frameshift = case_when(Pseudoexon_activation == "NO" ~ as.character(NA),
                                           is.na(Pseudoexon_intron_type)==TRUE ~ as.character(NA),
                                           as.numeric(Pseudoexon_intron_type)==-1 ~ as.character(NA),
                                           (abs(Pseudoexon_start-Pseudoexon_end)+1) %% 3 != 0 ~ "YES",
                                           TRUE ~ "NO"))

# Add additional retained intron information 
output <- output %>%
  rowwise() %>%
  mutate(Retained_intron_info = ifelse(Intron_retention == "YES",
                                       find_retained_intron(ALpos = GEO_AL,
                                                            DLpos = GEO_DL,
                                                            transcript = name,
                                                            refseqTab = refseq_boundaries),"NA|NA")) %>%
  separate(., col = Retained_intron_info, into = c("Retained_intron","Retained_intron_type"), sep = "[|]", remove = TRUE) %>%
  mutate(Retained_intron = na_if(Retained_intron, "NA"),
         Retained_intron_type = na_if(Retained_intron_type, "NA")) %>%
  mutate(Intron_retention_frameshift = case_when(Intron_retention == "NO" ~ as.character(NA),
                                                 is.na(Retained_intron_type)==TRUE ~ as.character(NA),
                                                 as.numeric(Retained_intron_type)==-1 ~ as.character(NA),
                                                 (abs(GEO_AL - GEO_DL)-1) %% 3 != 0 ~ "YES",
                                                 TRUE ~ "NO"))

# Add additional exon skipping information
output <- output %>%
  rowwise() %>%
  mutate(Lost_exon_info = ifelse(Exon_skipping == "YES",
                                 find_exon_lost(transcript = name,
                                                DLposition = GEO_DL,
                                                ALposition = GEO_AL,
                                                refseqTab = refseq_boundaries),NA_character_)) %>%
  separate(., col = Lost_exon_info, into = c("Lost_exons","Lost_exons_combined_size"), sep = "[|]", remove = TRUE) %>%
  mutate(Lost_exons_combined_size = na_if(Lost_exons_combined_size, "NA"),
         Lost_exons = na_if(Lost_exons, "NA")) %>%
  mutate_at(vars(c("Lost_exons_combined_size")), ~as.numeric(.)) %>%
  mutate(Exon_skipping_frameshift = case_when(Exon_skipping == "NO" ~ as.character(NA),
                                              Exon_skipping == "YES" & is.na(Lost_exons_combined_size) ~ as.character(NA),
                                              as.integer(Lost_exons_combined_size) %% 3 != 0 ~ "YES",
                                              TRUE ~ "NO"))

# Add additional information for partial deletion and retention
# Calculate the start and end of altered exon
output <- output %>%
  rowwise() %>%
  mutate(Partial_exon_start = case_when(Partial_intron_retention != "YES" & Partial_exon_deletion != "YES" ~ NA_integer_,
                                        Partial_intron_retention == "YES" & strand == 1 & bp_5prime < 0 ~ as.integer((eStart+1)-abs(bp_5prime)),
                                        Partial_intron_retention == "YES" & strand == -1 & bp_3prime > 0 ~ as.integer((eStart+1)-bp_3prime),
                                        Partial_exon_deletion == "YES" & strand == 1 & bp_5prime > 0 ~ as.integer((eStart+1)+abs(bp_5prime)),
                                        Partial_exon_deletion == "YES" & strand == -1 & bp_3prime < 0 ~ as.integer((eStart+1)+abs(bp_3prime)),
                                        TRUE ~ as.integer(eStart+1)),
         Partial_exon_end = case_when(Partial_intron_retention != "YES" & Partial_exon_deletion != "YES" ~ NA_integer_,
                                      Partial_intron_retention == "YES" & strand == 1 & bp_3prime > 0 ~ as.integer(eEnd+bp_3prime),
                                      Partial_intron_retention == "YES" & strand == -1 & bp_5prime < 0 ~ as.integer(eEnd+abs(bp_5prime)),
                                      Partial_exon_deletion == "YES" & strand == 1 & bp_3prime < 0 ~ as.integer(eEnd-abs(bp_3prime)),
                                      Partial_exon_deletion == "YES" & strand == -1 & bp_5prime > 0 ~ as.integer(eEnd-bp_5prime),
                                      TRUE ~ as.integer(eEnd)),
         Partial_frameshift = case_when(is.na(Partial_exon_start) == TRUE & is.na(Partial_exon_end) == TRUE ~ as.character(NA),
                                        abs(bp_5prime+bp_3prime) %% 3 != 0 ~ "YES",
                                        TRUE ~ "NO"))


cat("\nExtracting amino sequence predictions\n")


# Calculate the predicted amino acid sequence changes
# for partial intron retention
output <- output %>%
  rowwise() %>%
  mutate(Partial_intron_retention_aaseq = ifelse(Partial_intron_retention == "YES",
                                                 get_partial_SEQ(transcript = name,
                                                                 consensusStart = eStart,
                                                                 consensusEnd = eEnd,
                                                                 partialStart = Partial_exon_start,
                                                                 partialEnd = Partial_exon_end,
                                                                 refseqTable = refseq_boundaries,
                                                                 frameshift = Partial_frameshift,
                                                                 varPos = POS,
                                                                 ref = REF,
                                                                 alt = ALT),"-"))
# for partial exon deletion
output <- output %>%
  rowwise() %>%
  mutate(Partial_exon_deletion_aaseq = ifelse(Partial_exon_deletion == "YES",
                                              get_partial_SEQ(transcript = name,
                                                              consensusStart = eStart,
                                                              consensusEnd = eEnd,
                                                              partialStart = Partial_exon_start,
                                                              partialEnd = Partial_exon_end,
                                                              refseqTable = refseq_boundaries,
                                                              frameshift = Partial_frameshift,
                                                              varPos = POS,
                                                              ref = REF,
                                                              alt = ALT),"-"))
# for (multi)exon skipping
output <- output %>%
  rowwise() %>%
  mutate(Exon_skipping_aaseq = ifelse(Exon_skipping == "YES",
                                      get_skip_SEQ(exons = Lost_exons,
                                                   transcript = name,
                                                   refseqTable = refseq_boundaries,
                                                   frameshift = Exon_skipping_frameshift,
                                                   varPos = POS,
                                                   ref = REF,
                                                   alt = ALT),"-"))
# for pseudoexon activation
output <- output %>%
  rowwise() %>%
  mutate(Pseudoexon_activation_aaseq = ifelse(Pseudoexon_activation == "YES",
                                              get_pseudo_SEQ(pseudoStart = Pseudoexon_start,
                                                             pseudoEnd = Pseudoexon_end,
                                                             refseqTable = refseq_boundaries,
                                                             frameshift = Pseudoexon_frameshift,
                                                             transcript = name,
                                                             varPos = POS,
                                                             ref = REF,
                                                             alt = ALT),"-"))
# for intron retention
output <- output %>%
  rowwise() %>%
  mutate(Intron_retention_aaseq = ifelse(Intron_retention == "YES",
                                         get_retention_SEQ(refseqTable = refseq_boundaries,
                                                           transcript = name,
                                                           intron = Retained_intron,
                                                           frameshift = Intron_retention_frameshift,
                                                           varPos = POS,
                                                           ref = REF,
                                                           alt = ALT),"-"))

# clean up the partial frameshift column
output <- output %>%
  rowwise() %>%
  mutate(Partial_frameshift = case_when(Partial_exon_deletion_aaseq == "deletion greater than exon size" ~ NA_character_,
                                        Partial_exon_deletion_aaseq == "does not affect coding region" ~ NA_character_,
                                        Partial_exon_deletion_aaseq == "impacts native start of stop site" ~ NA_character_,
                                        Partial_exon_deletion_aaseq == "cannot determine" ~ NA_character_,
                                        TRUE ~ Partial_frameshift))



################## Write output ################################################

cat("\nWriting output\n")

## Save output
output_all <- output %>% 
  bind_rows(input_splice_other)

# Drop columns from output
output_all <- output_all %>%
  dplyr::select(., c(`#CHROM`:DP_DL,
                     name,strand,Cryptic_Acceptor_activation,
                     Cryptic_Donor_activation,Any_splicing_aberration,
                     bp_5prime,bp_3prime,Partial_intron_retention,
                     Partial_exon_deletion,Partial_exon_start,
                     Partial_exon_end,Partial_frameshift,
                     Partial_intron_retention_aaseq,
                     Partial_exon_deletion_aaseq,Gained_exon_size,
                     Pseudoexon_activation,Pseudoexon_start,
                     Pseudoexon_end,Pseudoexon_frameshift,
                     Pseudoexon_intron,Pseudoexon_activation_aaseq,
                     Exon_skipping,Lost_exons,Exon_skipping_frameshift,
                     Exon_skipping_aaseq,Retained_intron_size,
                     Intron_retention,Retained_intron,
                     Intron_retention_frameshift,Intron_retention_aaseq)) %>%
  dplyr::rename(Used_RefSeq_Transcript = name)

output_all %>%
  write_tsv(output_file)

cat("\nSession Information:\n")

sessionInfo()
