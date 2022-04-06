#######################################################################
#     Basic image
#######################################################################
FROM cgap/cgap-ubuntu2004-py-38:0.0.1
MAINTAINER Phil Grayson (phil.d.grayson@gmail.com), Michele Berselli (berselli.michele@gmail.com)

#######################################################################
#     General updates & installing necessary Linux components
#######################################################################
RUN apt-get update -y && apt-get install -y \
    locales \
    tabix

#######################################################################
#     Setting working env
#######################################################################
WORKDIR /usr/local/bin

#######################################################################
#     Software
#######################################################################
## conda install
RUN conda install -c bioconda -y gatk4==4.2.5.0 && \
    conda clean -a -y -f

#######################################################################
#     Scripts
#######################################################################
## HaplotypeCaller parallel
COPY HaplotypeCaller-parallel.sh .
RUN chmod +x HaplotypeCaller-parallel.sh

## HaplotypeCaller
COPY HaplotypeCaller.sh .
RUN chmod +x HaplotypeCaller.sh

## CombineGVCFs parallel
COPY CombineGVCFs-parallel.sh .
RUN chmod +x CombineGVCFs-parallel.sh

## GenotypeGVCFs parallel
COPY GenotypeGVCFs-parallel.sh .
RUN chmod +x GenotypeGVCFs-parallel.sh

#######################################################################
#     Setting env variables
#######################################################################
## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]
