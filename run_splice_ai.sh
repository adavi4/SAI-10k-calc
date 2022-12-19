#!/bin/sh


## Requirements:
## htslib (v1.9) and SpliceAI (v1.3.1)


# tabix index vcf file
bgzip -c example_variants.vcf > example_variants.vcf.gz
tabix -p vcf example_variants.vcf.gz

# run spliceAI
spliceai -I example_variants.vcf.gz \
-O ./example_variants.spliceAI.vcf \
-R ./hg19.fa \
-A ./example_spliceai_tx.txt \
-D 4999
