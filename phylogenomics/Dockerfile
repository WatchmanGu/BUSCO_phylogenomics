from ubuntu
RUN apt-get update --fix-missing 
RUN yes | apt-get install build-essential vim python3 git wget hmmer ncbi-blast+ libboost-iostreams-dev bamtools libbamtools-dev lib32z1-dev mafft autoconf flex bison libtool
RUN git clone https://github.com/scapella/trimal.git
RUN git clone https://gitlab.com/ezlab/busco.git
WORKDIR /trimal/source
RUN make
ENV PATH=$PATH:/trimal/source
RUN yes | apt-get install raxml
WORKDIR /busco
RUN python3 setup.py install --user
ADD config.ini /busco/config/config.ini
WORKDIR /
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.tar.gz
RUN tar -zxvf augustus-3.3.tar.gz
WORKDIR /augustus
RUN make
ENV PATH=$PATH:/augustus/bin:/augustus/scripts
ENV AUGUSTUS_CONFIG_PATH=/augustus/config
WORKDIR /
RUN wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
RUN tar zxvf newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
RUN cp /newick-utils-1.6/src/nw_* /usr/local/bin
ADD run.sh /
ADD *.py /
