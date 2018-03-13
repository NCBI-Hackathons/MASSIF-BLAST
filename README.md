# Small Genome Assembly Repair Toolkit - SGART
A toolkit of different modules to quickly assess and repair badly assembled genomes.

## Introduction
Over the past 10 years, the quality of genome sequencing and assembly has improved with the advancement of both experimental and analysis techniques. However, badly assembled genomes are still used in research, particularly if the assembly is old or if the genome is complex. We offer a suite of modules that quickly assesses the quality of a genome assembly and improves it all in one pipeline. The first module searches for frameshifts in highly conserved gene of that species within the user's genome assembly of interest. The second module annotates potential exogenous virus contamination in the genome assembly of interest. The third module uses other sequencing data (DNA and/or RNA) to expand areas of the genome that did not assemble well.

## Module Workflow

## Setup
### System
### Software Packages
### Data files

## Module 1 - Assess assembly quality
### Case Usages
*

### Workflow

### Operation
` `

## Module 2 - Annotate exogenous viral DNA
### Case Usages
* Genome of interest

### Workflow

### Usage
` `

## Module 3 - Revamp genome with extra sequencing data
### Case Usages
* Genome of interest is old or has poor coverage.
* Use RNA-seq to improve assembly. 
* No existing reference, user performs de-novo assembly, and then uses dna-seq or rna-seq to improve assembly.
* Genome of interest has large number of gaps or unassigned contigs.

### Workflow

#### Inputs
* Genome assembly that needs to be improved
* Dna-seq or Rna-seq data

#### Improving assembly with dna-seq

* Align raw data to existing assembly
  - Create a BLAST db
```
sudo ./ncbi-magicblast-1.3.0/bin/makeblastdb -in /data/Candina_Albicans/assemblies/GCF_000182965.3_ASM18296v3_genomic.fna -dbtype nucl -parse_seqids
```
	- Use Magicblast for alignment
```
./ncbi-magicblast-1.3.0/bin/magicblast -db /data/Candina_Albicans/assemblies/GCF_000182965.3_ASM18296v3_genomic.fna -sra SRR3593469 -splice F -no_unaligned -num_threads 4 -out SRR3593469_into_GCF_000182965
```

* Convert SAM->BAM, Sort, and index SAM file
```
samtools sort /data/sra_BAM/SRR3593469_into_GCF_000182965 -O BAM -o /data/sra_BAM/test_example/SRR3593469_into_GCF_000182965.sorted.bam

samtools index /data/sra_BAM/test_example/SRR3593469_into_GCF_000182965.sorted.bam
```

* Use Pilon to improve assembly
```
java -jar /software/pilon-1.22.jar \
--genome /data/Candina_Albicans/assemblies/GCF_000182965.3_ASM18296v3_genomic.fna \
--frags /data/sra_BAM/test_example/SRR3593469_into_GCF_000182965.sorted.bam \
--outdir /data/sra_BAM/test_example/pilon_output/
```

` `
