#!/bin/bash

tblastn -query GCF_000182965.3_ASM18296v3_orthologs.fasta -db GCA_000149445.2_ASM14944v2_genomic.fna -evalue .00001 -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue frames'
