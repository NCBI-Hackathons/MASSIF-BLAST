class: CommandLineTool
cwlVersion: v1.0

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerImageId: assemblyscripts:latest

baseCommand: ./get_genome_assemblies.sh

inputs:
  - id: database
    type: String
    inputBinding:
      position: 1
      prefix: -db
  - id: query
    type: String
    inputBinding:
      position: 2
      prefix: -q
  - id: format
    type: string[]?
    inputBinding:
      position: 3
      prefix: -f
  - id: file-string
    type: string[]?
    inputBinding:
      position: 4
      prefix: -fs

outputs:
  - id: outs
    type: File
    outputBinding:
      glob: ["*${file-string}"]
