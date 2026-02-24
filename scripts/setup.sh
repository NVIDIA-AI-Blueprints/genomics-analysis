#!/bin/bash 

# Install apt dependencies 
sudo apt install -y bcftools tree pigz tabix

REPO_PATH=/home/ubuntu/genomics-analysis

# Import submodules (not properly imported by brev)
cd ${REPO_PATH} && git submodule update --init --recursive

# Install Python dependencies
curl -sSL https://bootstrap.pypa.io/get-pip.py | python3
python3 -m pip install -r ${REPO_PATH}/CodonFM/requirements.txt

# Pull Docker containers to avoid cluttering the notebook output 
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1