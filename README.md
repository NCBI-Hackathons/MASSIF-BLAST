# MASSIF BLAST
Modular ASSembly Improvement Framework using BLAST

## Introduction
Over the past 10 years, the quality of genome sequencing and assembly has improved with the advancement of both experimental and analysis techniques. However, badly assembled genomes are still used in research, particularly if the assembly is old or if the genome is complex. We offer a suite of modules that quickly assesses the quality of a genome assembly and improves it all in one pipeline. The first module searches for frameshifts in highly conserved gene of that species within the user's genome assembly of interest. The second module uses other DNA-sequencing data to repair areas of the genome that did not assemble well. The third module does the same as the second, but with RNA-sequencing data. A user can run the full pipeline including all modules, some modules, or just one module depending on their needs. Other modules, such as an annotator of exogenous virus contamination or a polyploid detector, can be easily incorporated later on.

## Modules


## Setup
### System
Depending on the module, the size of the genome of interest, and the amount of supporting data, the pipeline can require a significant amount of computational resources. We recommend using a computer that has at least 4 processor cores, 16 Gb of RAM, and approximately 3 Gb hard drive space.

### Software Packages
* Please see our [docker README.md](https://github.com/NCBI-Hackathons/MASSIF-BLAST/tree/master/docker/base/README.md)

## Module 1 - Quickly assess assembly quality
This module takes a list of conserved genes which can be at a phylum taxonomic level (or lower) and compares the protein sequence of a reference genome to the assembly of interest. Poor assemblies will likely exhibit frameshifts or incomplete protein sequences. The output file is a tab-deliminated file comparing its protein quality how many of the conserved genes had poor protein sequence in the assembly of interest.

### Use Cases
* Comparison of user's assembly to a reference genome.
* Assessment of old genome assembly that interests a user.

### Workflow
![Mod 1 Workflow](https://github.com/NCBI-Hackathons/assemblyrepair/blob/master/mod-1_workflow.png)

## Module 2 - Assembly Improving Pipeline using DNA-sequencing data
### Use Cases
* Genome of interest is old or has poor coverage.
* No existing reference, user performs de-novo assembly, and then uses DNA-seq to improve assembly.
* Genome of interest has large number of gaps or unassigned contigs.

### Workflow
![Mod 2 Workflow](https://github.com/NCBI-Hackathons/assemblyrepair/blob/master/mod-2_workflow.png)

## Module 3 - Assembly Improving Pipeline using RNA-sequencing data
### Use Cases
* Similar to Module 2, but with RNA-seq data.

### Workflow
![Mod 3 Workflow](https://github.com/NCBI-Hackathons/assemblyrepair/blob/master/mod-3_workflow.png)

## Future Directions

## Module 4 - Annotate exogenous viral DNA
### Use Case
* Viral DNA -- use VirusFriends NCBI Hackathon LINK

## People
* Chris Ball - <christopherball@rti.org>
* Thomas Dyar - <thomas.dyar@biogen.com>
* Jessica Maia - <jessica.maia@bd.com>
* Kyle Roell - <krroell@ncsu.edu>
* Kimiko Suzuki - <sksuzuki@ad.unc.edu>

## Reference
