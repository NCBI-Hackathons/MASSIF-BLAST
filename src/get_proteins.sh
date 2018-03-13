#!/bin/bash

cut -d ',' -f 1 conserved_genes.csv | esearch -db protein -query "{0}" | efetch -format fasta > conserved_proteins.fasta
