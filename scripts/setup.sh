#!/bin/bash 

# Install dependencies 
sudo apt install -y bcftools tree pigz

# Import submodules (not properly imported by brev)
git submodule init
git submodule update --recursive --remote

# Pull Docker containers to avoid cluttering the notebook output 
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1