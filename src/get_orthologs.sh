#!/bin/bash

blastp -query fungus_conserved.fasta -db ./GCF_000182965.3_ASM18296v3_protein.faa -max_target_seqs 1 -outfmt 6 | cut  -f 2 | xargs  -I {} sh -c "esearch -db protein -query "{}" | efetch -format fasta" > GCF_000182965.3_ASM18296v3_orthologs.fasta
