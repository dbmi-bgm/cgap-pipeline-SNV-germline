#######################################################################
#     Basic image
#######################################################################
FROM cgap/cgap-ubuntu2004-py-38:0.0.1
MAINTAINER Phil Grayson (phil.d.grayson@gmail.com), Michele Berselli (berselli.michele@gmail.com)

#######################################################################
#     General updates & installing necessary Linux components
#       - libncurses-dev, libbz2-dev, liblzma-dev used by samtools
#       - parallel is still used by mpileupCounts at the moment
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    parallel

#######################################################################
#     Setting working env
#######################################################################
WORKDIR /usr/local/bin

## granite VERSION
ARG GRANITE_VERSION=0.2.0

#######################################################################
#     Software
#######################################################################
## conda install
RUN conda install -y scipy==1.5.3 pyliftover==0.4 && \
    conda clean -a -y -f

## samtools (used by granite mpileupCounts)
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xjf samtools-1.9.tar.bz2 && \
    rm -rf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure && \
    make && \
    cd .. && \
    ln -s samtools-1.9 samtools

## granite
RUN pip install granite-suite==${GRANITE_VERSION}

#######################################################################
#     Scripts
#######################################################################
## sort-bgzip-vcf
COPY sort-bgzip-vcf.sh .
RUN chmod +x sort-bgzip-vcf.sh

## cgap-scripts
##  - depth_filter, samplegeno and duplines_vcf, liftover_hg19
RUN git clone https://github.com/dbmi-bgm/cgap-scripts/ && \
    cd cgap-scripts && \
    git checkout 4550b2b && \
    chmod +x depth_filter.py && \
    chmod +x samplegeno.py && \
    chmod +x duplines_vcf.py && \
    chmod +x liftover_hg19.py && \
    chmod +x portal_reformat_vcf.py && \
    cd ..

## dbSNP_ID_fixer
COPY scripts/dbSNP_ID_fixer.py .
RUN chmod +x dbSNP_ID_fixer.py

## hgvsg_creator
COPY scripts/hgvsg_creator.py .
RUN chmod +x hgvsg_creator.py

## higlass_joint_parser
COPY scripts/higlass_joint_parser.py .
RUN chmod +x higlass_joint_parser.py

## reformat_vcf
COPY reformat_vcf.sh .
RUN chmod +x reformat_vcf.sh

## mpileupCounts parallel (from granite repo)
RUN wget https://raw.githubusercontent.com/dbmi-bgm/granite/${GRANITE_VERSION}/docker/mpileupCounts-parallel.sh && \
    chmod +x mpileupCounts-parallel.sh

## tar-novoCaller
COPY tar-novoCaller.sh .
RUN chmod +x tar-novoCaller.sh

## merge-sort-vcf
COPY merge-sort-vcf.sh .
RUN chmod +x merge-sort-vcf.sh

## samplegeno
COPY samplegeno.sh .
RUN chmod +x samplegeno.sh

## parallel_dbSNP_ID_fixer
COPY parallel_dbSNP_ID_fixer.sh .
RUN chmod +x parallel_dbSNP_ID_fixer.sh

## uniq_variants
COPY uniq_variants.sh .
RUN chmod +x uniq_variants.sh

#######################################################################
#     Setting env variables
#######################################################################
ENV PATH=/usr/local/bin/cgap-scripts/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH

## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]
