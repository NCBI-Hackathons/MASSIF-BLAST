class: CommandLineTool
cwlVersion: v1.0

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerImageId: assemblyscripts:latest

baseCommand: ["bash", "/workdir/get_genome_assemblies.sh"]

inputs:
  - id: database
    type: string
    inputBinding:
      position: 1
      prefix: -db
  - id: query
    type: string
    inputBinding:
      position: 2
      prefix: -q
  - id: format
    type: string?
    inputBinding:
      position: 3
      prefix: -f
  - id: file_string
    type: string
    inputBinding:
      position: 4
      prefix: -fs

outputs:
  output:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*_genomic.fna.gz"
