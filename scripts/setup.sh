#!/bin/bash 

# Install apt dependencies 
sudo apt install -y bcftools pigz tabix samtools

# Each cloud provider uses a different home dir 
REPO_PATH=/home/$(whoami)/genomics-analysis
PYTHON_PATH=/home/$(whoami)/.venv/bin/python3

# Import submodules 
cd ${REPO_PATH} && git submodule update --init --recursive

# Install Python dependencies
curl -sSL https://bootstrap.pypa.io/get-pip.py | ${PYTHON_PATH}
${PYTHON_PATH} -m pip install -r ${REPO_PATH}/CodonFM/requirements.txt

# Pull Docker containers to avoid cluttering the notebook output 
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1