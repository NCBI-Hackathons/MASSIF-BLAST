# src Bash script run commands

## get_genome_assemblies.sh

* Input strings must be in double quotes.

Options:
  * `-db | --database` for database to search
  * `-q | --query` for query to use
  * `-f | --format`, optional, `efetch` format to use. Default: `docsum`
  * `-fs | --file-string`, optional, end of file string to download. Default: `_genomic.fna.gz`

Example:
```
./get_genome_assemblies.sh -db "assembly" \
-q "Candida Albicans[Organism] AND latest[filter]" \
-f "docsum" \
-fs "_genomic.fna.gz"
```
