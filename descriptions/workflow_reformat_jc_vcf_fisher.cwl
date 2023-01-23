cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: input_vcf
    type: File
    secondaryFiles:
      - .tbi
    doc: expect the path to the sample vcf gz file

  - id: proband_list
    type: File
    doc: expect the path to a list of sample IDs for probands (cases in case vs control)

  - id: gnomAD
    type: string[]
    doc: list of gnomAD version(s) to use as control (v2 and/or v3)

outputs:
  vcf:
    type: File
    outputSource: jc_parse_fisher/output

steps:
  reformat_vcf:
    run: reformat_vcf.cwl
    in:
      input:
        source: input_vcf
    out:  [output]

  jc_parse_fisher:
    run: jc_parse_fisher.cwl
    in:
      input:
        source: reformat_vcf/output
      proband_list:
        source: proband_list
      gnomAD:
        source: gnomAD
    out: [output]

doc: |
  run reformat_vcf.sh to identify worst consequence for each transcript and run internal integrity check |
  run higlass_joint_parser.py to clean annotations and calculate Fisher's exact scores for each variant for HiGlass
