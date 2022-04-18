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
