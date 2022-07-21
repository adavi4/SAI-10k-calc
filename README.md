# SpliceAI Parser

============

## 1. About
SpliceAI post-processing calculator. Requires two files: 
* SpliceAI output
* List of genes and transcripts of interest


If you use the parser, please cite our pre-print: [TBD]

**Authors:**  
Daffodil Canson - developed the parser   
Olga Kondrashova - implemented the parser in R  
*QIMR Berghofer Medical Research Insititute*


## 2. Requirements
**Software Requirements:**  

* R version 4.1.2
* R packages:  
	- tidyverse
	- optparse
	- rtracklayer


## 3. Usage and Program Options

Run spliceAI_parser.R script:

`Rscript spliceAI_parser.R -i example_variants.vcf -g example_gene_list.txt -o example_variants_parsed.tsv`



Full options are available by running:

`Rscript spliceAI_parser.R -h`


 Options:
 

	
	-i FILE, --in_vcf=FILE
		Input vcf file direct from SpliceAI

	-g FILE, --gene_list=FILE
		Gene list, tab-separated with two columns - Gene and RefSeq_ID

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
		Delta score (acceptor & donor loss)  - minimum [default= 0]

	--DS_ALDL_MAX_T=NUMERIC
		Delta score (acceptor & donor loss) - maximum [default= 0.05]

	--AG_T=NUMERIC
		Cryptic splice site - acceptor gain [default= 0.2]

	--DG_T=NUMERIC
		Cryptic splice site - donor gain [default= 0.2]

	-r FILE, --refseq_table=FILE
		Optional: ncbiRefSeq table pre-downloaded from ucsc

	--out_refseq=FILE
		Optional: save downloaded ncbiRefSeq table

	-h, --help
		Show this help message and exit
