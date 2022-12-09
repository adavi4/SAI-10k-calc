#!/usr/bin/env Rscript


## SpliceAI output parser
## Developed by Daffodil Canson
## Implemented in R by Olga Kondrashova
## Date: 21/07/2022

###############################################################################
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
    filter(., `name` == transcript)
  # Now correct table to account for CDS Start and End
  cdsStartPos <- as.integer(transcriptTab$cdsStart[[1]])
  cdsEndPos <- as.integer(transcriptTab$cdsEnd[[1]])
  cdsStartExon = filter(transcriptTab, eStart <= cdsStartPos & eEnd >= cdsStartPos) %>% 
    dplyr::select(., name, eStart, eEnd, eNum, strand, chrom)
  cdsEndExon = filter(transcriptTab, eStart <= cdsEndPos & eEnd >= cdsEndPos) %>% 
    dplyr::select(., name, eStart, eEnd, eNum, eFrame, intronStart, intronEnd, intron_num, prev_eFrame)
  # modify the existing exons
  transcriptTab <- transcriptTab %>%
    mutate(eStart = ifelse(eNum == cdsStartExon$eNum, cdsStartPos, eStart),
           eStart = ifelse(eNum == cdsEndExon$eNum, cdsEndPos, eStart),
           eFrame = ifelse(eNum == cdsEndExon$eNum, -1, eFrame),
           prev_eFrame = ifelse(eNum == cdsEndExon$eNum, NA_integer_, prev_eFrame),
           intronStart = ifelse(eNum == cdsEndExon$eNum, NA_integer_, intronStart),
           intronEnd = ifelse(eNum == cdsEndExon$eNum, NA_integer_, intronEnd),
           intron_num = ifelse(eNum == cdsEndExon$eNum, NA_integer_, intron_num))
  # add new non-coding exons
  transcriptTab <- transcriptTab %>%
    ungroup() %>%
    add_row(., eStart = cdsStartExon$eStart, eEnd = cdsStartPos, eFrame = -1, 
            eNum = cdsStartExon$eNum, strand = cdsStartExon$strand, 
            chrom = cdsStartExon$chrom) %>%
    add_row(., eStart = cdsEndExon$eStart, eEnd = cdsEndPos, 
            eFrame = cdsEndExon$eFrame, eNum = cdsEndExon$eNum, 
            prev_eFrame = cdsEndExon$prev_eFrame, 
            strand = cdsStartExon$strand, chrom = cdsStartExon$chrom,
            intronStart = cdsEndExon$intronStart, intronEnd = cdsEndExon$intronEnd,
            intron_num = cdsEndExon$intron_num) %>%
    arrange(., eStart) 
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
    intronnum = filter(intronTab, intronStart+51 <= AGpos & intronEndAdj-51 >= DGpos)
  } else {
    intronnum = filter(intronTab, intronStart+51 <= DGpos & intronEndAdj-51 >= AGpos)
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
  # pull out the CDS start and stop for some extra screens for sequence prediction
  cdsStartPos = consensusTable$cdsStart[[1]]
  cdsEndPos = consensusTable$cdsEnd[[1]]
  # perform some pre-emptive checks for partial exon deletion and partial intron retention
  if(partialEnd-partialStart <= 0) {
    return("deletion greater than exon size")
  }
  if(partialStart < cdsStartPos | partialEnd > cdsEndPos) {
    return("impacts native start or stop site")
  }
  partialTable = consensusTable
  strand = consensusTable$strand[[1]]
  if (strand == -1) {
    partialTable = partialTable %>% arrange(., desc(eStartAdj))
  }
  # for variants affecting the first or last coding exon
  if (partialStart == as.integer(cdsStartPos)) {
    consensusStart = as.integer(cdsStartPos+1)
  }
  if (partialEnd == as.integer(cdsEndPos)) {
    consensusEnd = as.integer(cdsEndPos)
  }
  rowNumber = which(partialTable$eStartAdj==consensusStart & partialTable$eEnd == consensusEnd)
  # if can't find the target/affected exons then return that the table is "empty"
  if(is_empty(rowNumber)) {
    return("does not affect coding region")
  } else {
    # update the transcript table with the altered start and stop
    partialTable$eStartAdj[rowNumber] = partialStart
    partialTable$eEnd[rowNumber] = partialEnd
  }  
  # remove any but the immediate upstream exon
  partialTable <- partialTable %>% ungroup() %>%
    dplyr::slice(., ((rowNumber-1):n()))
  predictSEQ = determine_aaSEQ(partialTable,consensusTable,frameshift,varPos,ref,alt)
  return(predictSEQ)
}


# Extract the changed amino acid sequence for exon skipping
get_skip_SEQ <- function(exons,refseqTable,frameshift,transcript,varPos,ref,alt){
  consensusTable = make_consensus_table(transcript,refseqTable,filterNONCOD=FALSE)
  skipTable = consensusTable
  exonsList = purrr::flatten(str_split(exons, ","))
  exonsList = as.numeric(exonsList)
  # check for whether coding exons are lost
  skippedExons = skipTable %>% filter(., eNum %in% exonsList)
  if (-1 %in% skippedExons$eFrame) {
    if (max(skippedExons$eFrame)>=0) {
      return("native start or stop is lost")
    } else {
      return("non-coding exons lost")
    }
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
  # remove any but the immediate upstream exon
  if (is_empty(rowNumber)) {
    return("lost site/s do not match consensus")
  }
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
  cdsstart = min(consensusTable$cdsStart, na.rm = TRUE)
  cdsend = min(consensusTable$cdsEnd, na.rm = TRUE)
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
  alteredExonDNAseqs <- getSeq(Hsapiens,altTable$chrom,start=altTable$eStartAdj,end=altTable$eEnd)
  consensusExonDNAseqs <- getSeq(Hsapiens,consensusTab$chrom,start=consensusTab$eStartAdj,end=consensusTab$eEnd)
  # make necessary adjustments for the variant itself
  # also check whether the reference is correct
  genomeRef = as.character(getSeq(Hsapiens,currentChr,varPos,varPos))
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
  make_option("--DS_ALDL_MIN_T", type="character", default=0.02, 
              help="Delta score (acceptor & donor loss)  - minimum [default= %default]", 
              metavar="NUMERIC"), 
  make_option("--DS_ALDL_MAX_T", type="character", default=0.2, 
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
              metavar="FILE"),
  make_option("--ref", type="character", default="hg19",
              help="Optional: reference genome (hg19 or hg38) [default= %default]",
              metavar="STRING")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load required libraries
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
input <- read_tsv(opt$in_vcf, comment="##", col_types = c(`#CHROM` = "c"))
genes <- read_tsv(opt$gene_list)
output_file <- opt$out_file
refseq_table_file <- opt$refseq_table
output_refseq <- opt$out_refseq
ref_genome <- opt$ref

## Troubleshooting
# DS_AGDG_MIN_T <- 0.02
# DS_AGDG_MAX_T <- 0.05
# GEX_size_MIN <- 25
# GEX_size_MAX <- 500
# DS_ALDL_MIN_T <- 0.02
# DS_ALDL_MAX_T <- 0.2
# AG_T <- 0.2
# DG_T <- 0.2
# input <- read_tsv("./example_variants.vcf",comment="##")
# genes <- read_tsv("./example_gene_list.txt")
# output_file <- "./example_variants_parsed.tsv"
# ref_genome <- "hg19"

###############################################################################

if (ref_genome == "hg19") {
  library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
} else if (ref_genome == "hg38") {
  library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
} else {
  stop("Valid genome reference was not provided (specify either hg19 or hg38)")
}

cat("\nReference genome used: ",ref_genome,"\n")

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
  query <- ucscTableQuery(ref_genome, table = "ncbiRefSeq")
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
refseq_boundaries <- refseq_table_filtered_expanded %>%
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

###############################################################################

cat("\nPerforming calculations\n")

# Combining variant table with transcripts table
input_splice_distance <- input_splice_annot %>%
  left_join(refseq_boundaries, by = c("#CHROM" = "chrom_nochr",
                                      "SYMBOL" = "Gene")) %>%
  group_by(ID) %>%
  # -1 is to account for 0-based position
  mutate(dist_exon_start = as.numeric(POS) - as.numeric(eStart) -1,
         dist_exon_end = as.numeric(POS) - as.numeric(eEnd)) %>%
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

# Predicted exon sizes
output <- output %>% 
  rowwise() %>%
  mutate(Gained_exon_size = GEX_size,
         Lost_exon_size = LEX_predict,
         Retained_intron_size = RET_predict,
         bp_5prime = case_when(Cryptic_Acceptor_activation == "YES" & strand == 1 ~ GEO_AG-(eStart+1),
                               Cryptic_Acceptor_activation == "YES" & strand == -1 ~ eEnd-GEO_AG,
                               TRUE ~ 0),
         bp_3prime = case_when(Cryptic_Donor_activation == "YES" & strand == 1 ~ GEO_DG-eEnd,
                               Cryptic_Donor_activation == "YES" & strand == -1 ~ (eStart+1)-GEO_DG,
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

###############################################################################

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
                                      TRUE ~ as.integer(eEnd)))

# Now correct for partial start, ends for the CDS start, end
output <- output %>%
  rowwise() %>%
  mutate_at(vars(c("cdsStart","cdsEnd")), ~as.integer(.)) %>%
  mutate(Partial_exon_start = case_when(is.na(Partial_exon_start) ~ Partial_exon_start,
                                        Partial_exon_start <= (cdsStart+1) ~ as.integer(cdsStart+1),
                                        TRUE ~ Partial_exon_start),
         Partial_exon_end = case_when(is.na(Partial_exon_end) ~ Partial_exon_end,
                                      Partial_exon_end >= cdsEnd ~ cdsEnd,
                                      TRUE ~ Partial_exon_end)) %>%
  mutate(Partial_frameshift = case_when(is.na(Partial_exon_start) == TRUE & is.na(Partial_exon_end) == TRUE ~ as.character(NA),
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


###############################################################################

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
