#!/bin/bash
set -euo pipefail

# SPDX-FileCopyrightText: Copyright (c) 2024-2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

# Make directories for data and references files
DATA_DIR=$1
mkdir -p $DATA_DIR/ref 

# Download the data 
# Note: The tarball downloads unnecessary files and can be streamlined in the future 
cd $DATA_DIR && \
    wget -nc --progress=bar:force 2>&1 ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz && \
    wget -nc --progress=bar:force 2>&1 ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz && \
    wget -nc --progress=bar:force 2>&1 https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz 

# Pre-process the data 
cd $DATA_DIR && \
    pigz -dc parabricks_sample.tar.gz | tar xvf - && \
    mv parabricks_sample/Ref/* ref