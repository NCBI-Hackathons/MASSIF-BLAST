#!/bin/bash

wget `esearch -db assembly -query "Candida Albicans[Organism] AND latest[filter]" | \
efetch -format docsum | \
xtract -pattern DocumentSummary -element FtpPath_GenBank | \
awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`
