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

outputs:
  vcf:
    type: File
    outputSource: jc_parse_fisher/output

  vcf-check:
    type: File
    outputSource: integrity-check/output

steps:
  reformat_vcf:
    run: reformat_vcf.cwl
    in:
      input:
        source: input_vcf
    out:  [output]

  integrity-check:
    run: vcf-integrity-check.cwl
    in:
      input:
        source: reformat_vcf/output
    out: [output]

  jc_parse_fisher:
    run: jc_parse_fisher.cwl
    in:
      input:
        source: reformat_vcf/output
      probands:
        source: proband_list
    out: [output]

doc: |
  run reformat_vcf.py to identify worst consequence for each transcript |
  run an integrity check on the output vcf gz |
  run higlass_joint_parser.py to clean annotations and calculate Fisher's exact scores for each variant for HiGlass
