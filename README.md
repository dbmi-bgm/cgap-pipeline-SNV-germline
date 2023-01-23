<img src="https://github.com/dbmi-bgm/cgap-pipeline/blob/master/docs/images/cgap_logo.png" width="200" align="right">

# CGAP Pipeline for Germline Single Nucleotide Variants and small INDELs

This repository contains components for the CGAP pipeline for Single Nucleotide Variants (SNVs) and small insertions/delitions (INDELs) in germline data:

  * CWL workflow descriptions
  * CGAP Portal *Workflow* and *MetaWorkflow* objects
  * CGAP Portal *Software*, *FileFormat*, and *FileReference* objects
  * ECR (Docker) source files, which allow for creation of public Docker images (using `docker build`) or private dynamically-generated ECR images (using [*cgap pipeline utils*](https://github.com/dbmi-bgm/cgap-pipeline-utils/) `pipeline_deploy`)

The pipeline starts from analysis-ready `bam` files and produces `g.vcf` and `vcf` files containing calls for SNVs and small INDELs as output.
For more details check [*documentation*](https://cgap-pipeline-main.readthedocs.io/en/latest/Pipelines/Downstream/SNV_germline/index-SNV_germline.html "SNV germline").
