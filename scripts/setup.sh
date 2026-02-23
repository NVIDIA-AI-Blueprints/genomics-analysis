#!/bin/bash 

# Install bcftools 
wget https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2
tar xjf bcftools-1.23.tar.bz2
cd bcftools-1.23 && \
    ./configure && \
    make && \
    make install

# Install CodonFM Python requirements
pip install -r CodonFM/requirements.txt