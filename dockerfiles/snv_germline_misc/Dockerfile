#######################################################################
#     Basic image
#######################################################################
FROM 4dndcic/ubuntu16.04-miniconda-python36:v1
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

## granite VERSION
ARG GRANITE_VERSION=0.2.0

#######################################################################
#     Software
#######################################################################
## conda install
RUN conda install -c bioconda -y peddy==0.4.7 \
    python=3.6.8 pysam==0.15.3 pytabix==0.0.2 pyfaidx==0.5.8 tqdm==4.48.2 && \
    conda clean -a -y -f

## bamsnap
RUN pip install bamsnap-cgap==0.3.0

## granite
RUN pip install granite-suite==${GRANITE_VERSION}

#######################################################################
#     Scripts
#######################################################################
## peddy
COPY peddy.sh .
RUN chmod +x peddy.sh

#######################################################################
#     Setting env variables
#######################################################################
## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]
