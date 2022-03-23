#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/snv_germline_granite:VERSION

baseCommand: [reformat_vcf.sh]

inputs:
  - id: input
    type: File
    inputBinding:
      position: 1
    doc: expect the path to the sample vcf gz file

  - id: outputfile
    default: 'reformatted.vcf'
    type: string
    inputBinding:
      position: 2
    doc: name of output intermediate vcf file

outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.outputfile + ".gz")
    secondaryFiles:
      - .tbi


doc: |
  run portal_reformat_vcf.py |
  run integrated vcf integrity check
