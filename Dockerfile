FROM ubuntu:16.04

#####################################################################
# python 2.7
#####################################################################

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
    build-essential \
    ca-certificates \
    gcc \
    git \
    libpq-dev \
    make \
    python-pip \
    python2.7 \
    python2.7-dev \
    ssh \
    && apt-get autoremove \
    && apt-get clean

RUN pip install -U "setuptools==3.4.1"
RUN pip install -U "pip==1.5.4"
RUN pip install -U "Mercurial==2.9.1"
RUN pip install -U "virtualenv==1.11.4"

RUN apt-get update \
  && apt-get install -y build-essential \
  wget \
  unzip \
  bzip2 \
  software-properties-common

#####################################################################
# bedtools
#####################################################################

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz && \
tar -zxvf bedtools-2.25.0.tar.gz && \
cd bedtools2 && \
make

RUN apt-get install -y autoconf

ENV PATH $PATH:/bedtools2/bin 


#####################################################################
# Install Anaconda (python 2.7 version)
#####################################################################


RUN wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh
RUN bash Anaconda2-4.2.0-Linux-x86_64.sh -b -p anaconda
RUN rm Anaconda2-4.2.0-Linux-x86_64.sh
ENV PATH /anaconda/bin:$PATH

#####################################################################
# BLAST
#####################################################################

RUN conda install -c bioconda blast

#####################################################################
# BWA
#####################################################################

RUN conda install -c bioconda bwa=0.7.15

#####################################################################
# JELLYFISH
#####################################################################

ENV JELLYFISH_VERSION "2.2.6"


RUN apt-get update && apt-get install -y \
	build-essential \
	tar


ADD https://github.com/gmarcais/Jellyfish/releases/download/v$JELLYFISH_VERSION/jellyfish-$JELLYFISH_VERSION.tar.gz \
	jellyfish.tar.gz
RUN mkdir /src/ && \
	tar zxvf jellyfish.tar.gz -C /src/ --strip-components=1 && \
	cd /src/ && \
	./configure && \
	make && \
	make install


ENV PATH="/src/bin/:$PATH" \
	LD_LIBRARY_PATH="/src/lib/:$LD_LIBRARY_PATH" \
	MANPATH="/src/share/man:$MANPATH" \
	PKG_CONFIG_PATH="/src/lib/pkgconfig:$PKG_CONFIG_PATH"

#####################################################################
# SALMON
#####################################################################

#RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz
#RUN tar xzvf /salmon-1.3.0_linux_x86_64.tar.gz
#ENV PATH /salmon-latest_linux_x86_64/bin/:$PATH

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/salmon-0.14.1_linux_x86_64.tar.gz
RUN tar xzvf salmon-0.14.1_linux_x86_64.tar.gz
ENV PATH /salmon-latest_linux_x86_64/bin/:$PATH

#####################################################################
# Installazione trinity
#####################################################################

RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.7.0-PRERELEASE.tar.gz
RUN tar -zxvf Trinity-v2.7.0-PRERELEASE.tar.gz
RUN cd trinityrnaseq-Trinity-v2.7.0-PRERELEASE/
RUN cd /trinityrnaseq-Trinity-v2.7.0-PRERELEASE
WORKDIR /trinityrnaseq-Trinity-v2.7.0-PRERELEASE
RUN apt-get install -y rsync
RUN make
RUN make plugins
RUN make install

ENV PATH /trinityrnaseq-Trinity-v2.7.0-PRERELEASE:$PATH


#####################################################################
# Bowtie2
#####################################################################

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip/download -O bowtie2-2.3.4.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.1-linux-x86_64.zip && \
    mv bowtie2-2.3.4.1-linux-x86_64/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-2.3.4.1-linux-x86_64


#####################################################################
# JAVA 8
#####################################################################

RUN apt update
RUN add-apt-repository ppa:openjdk-r/ppa
RUN apt-get update
RUN apt-get install -y openjdk-8-jre

#####################################################################
# SAMTOOLS
#####################################################################

RUN apt-get install -y libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev
RUN \
  wget -c https://github.com/samtools/htslib/archive/1.4.tar.gz && \
  tar -zxvf 1.4.tar.gz && \
  rm 1.4.tar.gz && \
  mv htslib-1.4 htslib && \
  cd htslib && \
  autoreconf && \
  ./configure && \
  make && \
  make install
# Installazione samtools
RUN \
  wget -c https://github.com/samtools/samtools/archive/1.4.tar.gz && \
  tar -zxvf 1.4.tar.gz && \
  cd samtools-1.4 && \
  make && \
  make install

ENV PATH $PATH:samtools-1.4/samtools 


#####################################################################
# VIR
#####################################################################


WORKDIR /home


RUN git clone https://github.com/epischedda/ViR.git

WORKDIR /home/ViR


