<img src="https://github.com/dbmi-bgm/cgap-pipeline/blob/master/docs/images/cgap_logo.png" width="200" align="right">

# CGAP Single-Nucleotide Variants Pipeline, Germline

* This repo contains components for CGAP single-nucleotide variants pipeline for germline mutations
  * CWL
  * CGAP Portal Workflows and Metaworkflow
  * ECR (Docker) source files, which allow for creation of public Docker images (using `docker build`) or private dynamically-generated ECR images (using `cgap pipeline utils` (https://github.com/dbmi-bgm/cgap-pipeline-utils/) `deploy_pipeline`)

For more details check [*documentation*](https://cgap-pipeline-master.readthedocs.io/en/latest/Pipelines/Downstream/SNV_germline/index-SNV_germline.html "SNV germline documentation")

### Version updates

#### v27

* Initial release based on v26 pipeline
* Changes in repo structure to allow for compatibility with new pipeline organization
