#!/bin/bash 

# Install additional dependencies 
sudo apt install -y bcftools tree 

# Pull Docker containers to avoid cluttering the notebook output 
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
