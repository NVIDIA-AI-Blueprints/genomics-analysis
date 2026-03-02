# SPDX-FileCopyrightText: Copyright (c) 2024-2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

#!/bin/bash 

# Make directories for data and references files 
DATA_DIR=$(pwd)/data
mkdir -p $DATA_DIR/ref

# Download the exome files 
cd $DATA_DIR && \
    wget -nc --progress=bar:force 2>&1 ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz && \
    wget -nc --progress=bar:force 2>&1 ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz && \
    cd - 

# Download the reference files 
# Note: This downloads unnecessary files and can be streamlined in the future 
wget -nc --progress=bar:force 2>&1 -O parabricks_sample.tar.gz https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz
pigz -dc parabricks_sample.tar.gz | tar xf -
mv parabricks_sample/Ref $DATA_DIR/ref
