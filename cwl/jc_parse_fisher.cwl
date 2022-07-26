#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/snv_germline_granite:VERSION

baseCommand: [python3, /usr/local/bin/higlass_joint_parser.py]

inputs:
  - id: input
    type: File
    inputBinding:
      prefix: -i
    doc: expect the path to the sample vcf gz file

  - id: gnomAD
    type: string[]
    inputBinding:
      prefix: -g
    doc: list of gnomAD version(s) to use as control (v2 and/or v3)

  - id: probands
    type: string[]
    inputBinding:
      prefix: -p
    doc: list of proband sample IDs (cases in case vs control)

  - id: outputfile
    default: 'joint_called_higlass.vcf'
    type: string
    inputBinding:
      prefix: -o
    doc: base name of output vcf gz file

outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.outputfile + ".gz")
    secondaryFiles:
      - .tbi

doc: |
  run higlass_joint_parser.py to clean annotations and generate Fisher's exact scores for HiGlass viewing of a jointly-called vcf file
