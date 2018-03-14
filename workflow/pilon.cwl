class: CommandLineTool
cwlVersion: v1.0

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerImageId: assemblyrepair:latest

baseCommand: ["java", "-Xmx8G", "-jar", "/lib/pilon-1.22.jar"]

inputs:
  - id: genome
    type: File
    inputBinding:
      position: 1
      prefix: --genome
  - id: frags
    type: File
    inputBinding:
      position: 2
      prefix: --frags
    secondaryFiles: [.bai]
  - id: output
    type: string
    inputBinding:
      position: 3
      prefix: --output
  - id: options
    type: string
    inputBinding:
      position: 4

outputs:
  - id: outs
    type: File
    outputBinding:
      glob: ["${output}*"]

