<img src="https://github.com/dbmi-bgm/cgap-pipeline/blob/master/docs/images/cgap_logo.png" width="200" align="right">

# CGAP Pipeline for Germline Single-Nucleotide Variants and small INDELs

This repository contains components for the CGAP pipeline for germline single-nucleotide variants (SNVs) and small INDELs:

  * CWL workflows
  * CGAP Portal Workflows and MetaWorkflows objects
  * ECR (Docker) source files, which allow for creation of public Docker images (using `docker build`) or private dynamically-generated ECR images (using [*cgap pipeline utils*](https://github.com/dbmi-bgm/cgap-pipeline-utils/) `deploy_pipeline`)

The pipeline starts from analysis ready `bam` files and produces `g.vcf` and `vcf` files containing calls for SNVs and small INDELs as output.
For more details check [*documentation*](https://cgap-pipeline-main.readthedocs.io/en/latest/Pipelines/Downstream/SNV_germline/index-SNV_germline.html "SNV germline").

### Version Updates

#### v1.0.0
* v27 -> v1.0.0, we are starting a new more comprehensive versioning system
* Added some change in metaworkflows to accomodate the changes in foursight
* Re-organized and updated docker components

#### v27
* This repo starts from the v26 release of [*cgap-pipeline*](https://github.com/dbmi-bgm/cgap-pipeline) and contains the second half of the original pipeline (starting from `bam` files and producing `g.vcf` and `vcf` files)
* Changes in repo structure to allow for compatibility with new pipeline organization
* Changes to the hg19 liftover/HGVSG workflow now allow for HGVSg characterizations to be generated from hg19 liftover coordinates
* New cwls and workflows added for annotation of jointly called `vcf` files
