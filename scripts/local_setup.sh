#!/bin/bash 

# This setup script should be used for local development

# Install apt dependencies 
sudo apt install -y bcftools pigz tabix samtools

# Import submodules 
git submodule update --init --recursive

# Install Python dependencies
pip install jupyterlab
pip install -r CodonFM/requirements.txt

# Pull Docker containers to avoid cluttering the notebook output 
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1