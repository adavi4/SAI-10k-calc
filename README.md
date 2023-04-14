# SAI-10k-calc

============

## 1. About
SpliceAI post-processing calculator. Requires two files: 
* VCF file - preferred (or SpliceAI output - VCF)
* List of genes and transcripts of interest


If you use the parser, please cite our manuscript: Canson DM, Davidson AL et al. Bioinformatics 2023 https://doi.org/10.1093/bioinformatics/btad179

**Authors:**  
Daffodil Canson - developed the parser   
Aimee Davidson - implemented the parser in R and added extra functionality   
Olga Kondrashova - implemented the parser in R  
*QIMR Berghofer Medical Research Insititute*


## 2. Requirements
**Software Requirements:**  

* R version 4.1.2
* R packages:  
	- tidyverse
	- optparse
	- rtracklayer
	- Biostrings
	- BSgenome.Hsapiens.UCSC.hg19 (or BSgenome.Hsapiens.UCSC.hg38)
* htslib 1.9 (for running SpliceAI)
* SpliceAI 1.3.1 (for running SpliceAI)


## 3. Usage and Program Options

**Step 1.** Run `download_tx.R` script to download and pre-process transcript tables:

`Rscript download_tx.R -g example_gene_list.txt --out_refseq example_refseq_tx_hg19.txt --out_tx_spliceai example_spliceai_tx_hg19.txt`

Full options are available by running:

`Rscript download_tx.R -h`

Options:
	-g FILE, --gene_list=FILE
		Gene list, tab-separated with two columns - Gene and RefSeq_ID

	--out_refseq=FILE
		Save downloaded preprocessed transcript table

	--out_tx_spliceai=FILE
		Save preprocessed transcript table for SpliceAI

	--out_refseq_full=FILE
		Optional: save downloaded ncbiRefSeq table

	--refseq_full=FILE
		Optional: ncbiRefSeq table pre-downloaded from ucsc

	--ref=STRING
		Optional: reference genome (hg19 or hg38) [default= hg19]

	-h, --help
		Show this help message and exit


**Step 2.** Run SpliceAI using the pre-processed transcript table (example_spliceai_tx_hg19.txt) from the first script. Check the example script `run_splice_ai.sh`.   
  
It is recommended to run [SpliceAI](https://github.com/Illumina/SpliceAI) with a custom annotation file (pre-processed transcript table) using `-A` flag, because predictions for specific transcripts may differ from GENCODE V24 canonical transcript predictions (SpliceAI default). 

**Step 3.** Run `spliceAI_parser.R` script using the pre-processed transcript table from Step 1 (example_refseq_tx_hg19.txt) and the SpliceAI output vcf file.

`Rscript spliceAI_parser.R -i example_variants.vcf -r example_refseq_tx_hg19.txt -o example_variants_parsed.tsv`

Full options are available by running:

`Rscript spliceAI_parser.R -h`


 Options:
	-i FILE, --in_vcf=FILE
		Input vcf file direct from SpliceAI

	-r FILE, --refseq_table=FILE
		Preprocessed UCSC ncbiRefSeq table (download_tx.R script)

	-o FILE, --out_file=FILE
		output file name [default= out.tsv]

	--DS_AGDG_MIN_T=NUMERIC
		Delta score (acceptor & donor gain) - minimum [default= 0.02]

	--DS_AGDG_MAX_T=NUMERIC
		Delta score (acceptor & donor gain) - maximum [default= 0.05]

	--GEX_size_MIN=INT
		Gained exon size range - minimum [default= 25]

	--GEX_size_MAX=INT
		Gained exon size range - maximum [default= 500]

	--DS_ALDL_MIN_T=NUMERIC
		Delta score (acceptor & donor loss)  - minimum [default= 0.02]

	--DS_ALDL_MAX_T=NUMERIC
		Delta score (acceptor & donor loss) - maximum [default= 0.2]

	--AG_T=NUMERIC
		Cryptic splice site - acceptor gain [default= 0.2]

	--DG_T=NUMERIC
		Cryptic splice site - donor gain [default= 0.2]

	--ref=STRING
		Optional: reference genome (hg19 or hg38) [default= hg19]

	-h, --help
		Show this help message and exit
