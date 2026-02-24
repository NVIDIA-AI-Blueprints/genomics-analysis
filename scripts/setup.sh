#!/bin/bash 

# Install additional dependencies 
sudo apt install -y bcftools tree 

# For the CodonFM Submodule (not properly loaded into brev)
# 3. Switch to the main branch 
cd CodonFM && git switch main 

# Pull Docker containers to avoid cluttering the notebook output 
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
