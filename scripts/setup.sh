#!/bin/bash 

apt install y bcftools

# Install CodonFM Python requirements
pip install -r CodonFM/requirements.txt

# Pull Docker containers
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
