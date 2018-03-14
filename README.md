# MASSIF BLAST
Modular ASSembly Improvement Framework using BLAST

## Introduction
Over the past 10 years, the quality of genome sequencing and assembly has improved with the advancement of both experimental and analysis techniques. However, badly assembled genomes are still used in research, particularly if the assembly is old or if the genome is complex. We offer a suite of modules that quickly assesses the quality of a genome assembly and improves it all in one pipeline. The first module searches for frameshifts in highly conserved gene of that species within the user's genome assembly of interest. The second module uses other sequencing data (DNA and/or RNA) to expand areas of the genome that did not assemble well. A user can run the full pipeline including all modules, some modules, or just one module depending on their needs. Other modules, such as an annotator of exogenous virus contamination or a polyploid detector, can be easily incorporated later on.

## Module Workflow

## Setup
### System
Depending on the module, the size of the genome of interest, and the amount of supporting data, the pipeline can require a significant amount of computational resources. We recommend using a computer that has at least 4 processor cores, 16 Gb of RAM and approximately 3 Gb hard drive space for each analyzed run.

### Software Packages
### Data files

## Module 1 - Assess assembly quality
### Use Cases
* 

### Workflow

### Operation
` `

## Module 2 - Assembly Improving Pipeline using DNA-sequencing data
### Use Cases
* Genome of interest is old or has poor coverage.
* No existing reference, user performs de-novo assembly, and then uses dna-seq or rna-seq to improve assembly.
* Genome of interest has large number of gaps or unassigned contigs.

### Workflow

![Workflow](https://github.com/NCBI-Hackathons/assemblyrepair/mod-2_workflow.png.jpg)

### Output file

## Module 3 - Assembly Improving Pipeline using RNA-sequencing data
### Use Cases
* Similar to Module 2, but with RNA-seq data.

### Workflow

![Workflow](https://github.com/NCBI-Hackathons/assemblyrepair/mod-3_workflow.png.jpg)

### Output file

## Future Directions

## Module 4 - Annotate exogenous viral DNA
### Use Case
Viral DNA -- use VirusFriends NCBI Hackathon LINK

## People
* Chris Ball <>
* Thomas Dyar <>
* Jessica Maia <>
* Kyle Roell <>
* Kimiko Suzuki <sksuzuki@ad.unc.edu>

## Reference
